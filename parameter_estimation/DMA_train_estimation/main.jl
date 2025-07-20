## load the packages used in the estimation
# plotting
using PyPlot
rc("text", usetex=true) # use LaTeX in figures
using myPlot

# data manipulation (loading, writing, etc)
using CSV
using Dates
using DataFrames # for loading and saving the results, e.g. CSV.write("measured_psd.csv",DataFrame(PSD_meas'); writeheader=false)
using Printf

# scientific package from the official Julia repositories
using LinearAlgebra
using Statistics
using DSP
# using SpecialMatrices
using Polynomials
using StatsBase

# implemented scientific packages
using StochProc # for computing the covariance matrix
using AeroMec2  # GDE discretization and error estimation
using utilsFun  # for the softMax functions
using EKF       # Kalman filter and smoother

println("loading packages: done!")


##
## define flags to select the data
##
FLAG_0000_00 = true
FLAG_0000_01 = false
FLAG_0000_02 = false
FLAG_0000_03 = false
FLAG_0000_04 = false

FLAG_1952_02 = false
FLAG_1802_01 = false

FLAG_1906_03 = false

# all cases need 0-padding of the data to emulate a better time resolution so that the size resolution is not too bad
need_padding = true

# plotting flags: flip them to true only for extensive plotting
PLOT_MORE      = false
PLOT_EVEN_MORE = false

path_to_data = "/path/to/data/folder/" # "../../data/" # 
path_to_data = "../../data/" # 

##################################################
###    read data and init some variables       ###
##################################################

if FLAG_1952_02
    GR_EST_1952_02 = true
    input_folder = string(path_to_data,"1952_02/")
else
    GR_EST_1952_02 = false
end
if FLAG_0000_00
    input_folder = string(path_to_data,"0000_00/")
end
if FLAG_0000_01
    input_folder = string(path_to_data,"0000_01/")
end
if FLAG_0000_02
    input_folder = string(path_to_data,"0000_02/")
end
if FLAG_0000_03
    input_folder = string(path_to_data,"0000_03/")
end
if FLAG_0000_04
    input_folder = string(path_to_data,"0000_04/")
end

if FLAG_1802_01
    organic_run = true
    input_folder = string(path_to_data,"1802_01/")
    GR_EST_1802_01 = true # false
else
    organic_run = false
    GR_EST_1802_01 = false
end

if FLAG_1906_03
    input_folder = string(path_to_data,"1906_03/")
    GR_EST_1906_03 = true
else
    GR_EST_1906_03 = false
end

GT_loaded = (!FLAG_1952_02 & !FLAG_1802_01 & !FLAG_1906_03) #WARNING need more memory for this!!!


# output folder: the results will be written at this location
folder = "results/all_channels/test_that_will_be_deleted/"
mkpath(folder)


##
## load measurement operator
##
deltaT = @elapsed include("load_model.jl")
# model_meas_dim = length(diameter_model);
println(@sprintf "time to load the measurement model: %fs" deltaT)

##
## load data
##
deltaT = @elapsed include("load_measurements.jl")
println(@sprintf "time to load the data: %fs" deltaT)

# compute bin widths
cst_r_raw = mean(diameter_model[2:end]./diameter_model[1:end-1]);
# convert to the right units for the evolution model (the GDE is implemented in terms of concentration, not density)
H_DMATRAIN = H_DMATRAIN./log10(cst_r_raw);
nb_ch = size(H_DMATRAIN,1);


##
## shift error
##
dH1 = similar(H_DMATRAIN);
dH2 = similar(H_DMATRAIN);
dH  = similar(H_DMATRAIN);
d_idx_h = 1
dH1[:,1:d_idx_h]       .=0.0
dH2[:,end-d_idx_h+1:end] .=0.0
dH1[:,d_idx_h+1:end]   = abs.(H_DMATRAIN[:,1:end-d_idx_h]-H_DMATRAIN[:,d_idx_h+1:end]) # log10(cst_r_raw)*
dH2[:,1:end-d_idx_h] = abs.(H_DMATRAIN[:,d_idx_h+1:end]-H_DMATRAIN[:,1:end-d_idx_h]) # log10(cst_r_raw)*
dH = ((dH1.>dH2).*dH1+(dH1.<=dH2).*dH2);

d_mod1 = similar(diameter_model);
d_mod2 = similar(diameter_model);
d_mod1[1:d_idx_h]      .=0.0
d_mod2[end-d_idx_h+1:end] .=0.0
d_mod1[d_idx_h+1:end]  = abs.(diameter_model[1:end-d_idx_h]-diameter_model[d_idx_h+1:end])
d_mod2[1:end-d_idx_h] = abs.(diameter_model[d_idx_h+1:end]-diameter_model[1:end-d_idx_h])
d_mod = ((d_mod1.>d_mod2).*d_mod1+(d_mod1.<=d_mod2).*d_mod2);
diameter_model_full = copy(diameter_model)


##
## down sample the measurement model
##
# cropping boundaries and down-sampling factor
idx_s_0 = 231 # only consider the evolution model where there is infomartion in the measurement (from the lower end of the smallest channel)
idx_s_z = 998
d_idx   = 48 # you don't want to try a finner resolution (d_idx<48) because the evolution model because unstable de to the growth rate (ds>2.0GR*dt)
if FLAG_1906_03
    d_idx   = 24
end
d_idx_half = floor(Int64,d_idx/2) # correct the centroid of the bin
diameter_model = diameter_model[(idx_s_0+d_idx_half):d_idx:(idx_s_z+d_idx_half)];
cst_r_model = mean(diameter_model[2:end]./diameter_model[1:end-1]);

# find the region of influence of each channel (where a channel has more weight than all other)
idx_bin_low = zeros(Int64,nb_ch);
idx_bin_hig = zeros(Int64,nb_ch);
idx_bin_low[1] = 1;
idx_bin_hig[1] = findlast(diameter_model.<sqrt(diameter_data[1]*diameter_data[2]));
for i in 2:nb_ch-1
    global idx_bin_low, idx_bin_hig
    idx_bin_low[i] = findlast(diameter_model.<sqrt(diameter_data[i-1]*diameter_data[i]))
    idx_bin_hig[i] = findlast(diameter_model.<sqrt(diameter_data[i]*diameter_data[i+1]))
end
idx_bin_low[end] = findlast(diameter_model.<sqrt(diameter_data[end-1]*diameter_data[end]))
idx_bin_hig[end] = length(diameter_model);



# because of truncation error in the file, some values are interpreted as negative while they are just very small value, which are 0s
H_DMATRAIN = H_DMATRAIN.*(H_DMATRAIN.>0.0);
H_DMATRAIN_FULL = copy(H_DMATRAIN)
H_DMATRAIN = H_DMATRAIN[:,idx_s_0:idx_s_z];
dH_FULL    = copy(dH)
dH = dH[:,idx_s_0:idx_s_z];
(meas_dim,model_meas_dim) = size(H_DMATRAIN)
model_meas_dim = length(collect(idx_s_0:idx_s_z))
H_DMATRAIN_SMALL = zeros(Cdouble,meas_dim,floor(Int64,model_meas_dim/d_idx));
dH_SMALL = zeros(Cdouble,meas_dim,floor(Int64,model_meas_dim/d_idx));
for i in 1:meas_dim
    for j in 1:convert(Int64,model_meas_dim/d_idx)
        range_col = d_idx*(j-1)+1:d_idx*j
        H_DMATRAIN_SMALL[i,j] = mean(H_DMATRAIN[i,range_col])
        dH_SMALL[i,j] = mean(dH[i,range_col])
    end
end
H_DMATRAIN = 1.0e0*H_DMATRAIN_SMALL;
(meas_dim,model_meas_dim) = size(H_DMATRAIN)
if GT_loaded
    idx_p_simu = zeros(Int64,model_meas_dim)
    for i in 1:model_meas_dim
        idx_p_simu[i] = findfirst(dp.>diameter_model[i])
    end
end

# set the dimensions: try to set a fixed state size
nbin  = length(diameter_model)                              # Number of particle size bins
model_dim = nbin+2nbin+nbin+2  # +4                      # [],        dimension of the model
meas_dim  = length(diameter_data)             # [],        dimension of the observation space
d0 = diameter_model[1] # diameter_data[1] # 1.5e-9 #
v0 = (pi/6.0)*(d0^3)

# compute the diameter, volume and constant ratio for the model state
cst_v = ((pi/6.0)*(diameter_model[end]^3)/v0)^(1.0/(nbin-1.0))
cst_r = cst_v^(1.0/3.0)
volume = v0*(cst_v.^(collect(0:nbin-1)));
diameter = ((6.0/pi)*volume).^(1.0/3.0);
delta = diameter*(sqrt(cst_r)-1.0/sqrt(cst_r));



##################################################
####   load the Kalman Filter implementation   ###
##################################################
myWSKF = EKF.wsKF(model_dim,meas_dim,n_samp,false,true);
myWSKF.t_samp = t_samp;
if need_padding
    myWSKF.time_dep_noise_model = true # evolution covariance matrix change with time
    myWSKF.time_dep_meas_op = true     # measurement covariance matrix change with time
end






##################################################
###             define the model               ###
##################################################

# definition of the characteristic constants (set parameter for nondimensionalization)
gamma0 = 1.0e-4                    # the maximum wall loss rate observed in the data or a fraction of the maximum
t0     = 600.0                     # the characteristic time t0=1.0/alpha0
x0     = 1.0                       # set to 1 for no normalization of the density
GR0    = 2.8e-13                   # the maximum condensation rate observed in the data or a fraction of the maximum
J0     = 0.5                       # normalizing constant for the nucleation rate
if FLAG_1906_03
    GR0    = 5.0*2.7778e-13
    J0     = 10.0
end

##
## the evolution model
##
include("model_mod.jl")

# allocation
x00 = Array{Cdouble,1}(undef,model_dim);
gam00_vec = Array{Cdouble,1}(undef,model_dim);
var_model = Array{Cdouble,2}(undef,model_dim,model_dim);
var_noise_diag = Array{Cdouble,1}(undef,model_dim);
var_noise = Array{Cdouble,2}(undef,model_dim,model_dim);

x_fil_all = Array{Cdouble,2}(undef,model_dim,n_samp);
x_pre_all = Array{Cdouble,2}(undef,model_dim,n_samp);
o_fil_all = Array{Cdouble,3}(undef,model_dim,model_dim,n_samp);
o_pre_all = Array{Cdouble,3}(undef,model_dim,model_dim,n_samp);
x_smo_all = similar(myWSKF.x_fil_all);
o_smo_all = similar(myWSKF.o_fil_all);

percentiles_cond_smo       = Array{Cdouble}(undef,nbin,n_samp,2);
percentiles_wall_smo       = Array{Cdouble,3}(undef,nbin,n_samp,2);
percentiles_nuc_smo        = Array{Cdouble,2}(undef,n_samp,2);
percentiles = [15;85]
Ns_per = 100

##
## initialization of the Kalman Filter
##
# how much do we trust the model (init)
fill!(var_model,0.0)
var_model[R_psd,R_psd] = gamma_psd

# parameter model
var_model[R_cond_all,R_cond_all] = gamma_cond
var_model[R_loss,R_loss] = gamma_loss
var_model[R_nuc,R_nuc]   = gamma_nuc


#how much do we trust the data (init)
var_noise_diag = data_var(y_all[:,1])
var_noise = diagm(var_noise_diag)

# initialize the measurement model for time independent measurement model
if !myWSKF.time_dep_meas_op
    myWSKF.H_me = set_measurement_jacobian!(myWSKF.H_me)
end

# initial state
x00[R_psd]      .= 0.0
x00[R_cond_all] .= -1.0
# Dilution is for sure known quite well, CLOUD is well-mixed (there are several studies on that), the dilution rate should be chamber flow/chamber volume: chamber flow 270 lpm, chamber volume 26.1 m3, so that makes Q_dil = 1.72e-4 s-1. Also wall losses are experimentally verified and certainly within 20 %. You can approximate them by: k_wall = 0.00163/ dp(nm). If that helps, we should use this information. Ultimately we are interested in GR and J and not if we reproduce the loss numbers which we have measured several times.
loss_rate = 1.72e-4 .+ 0.00163./(1.0e9diameter)
if organic_run
    loss_rate = 1.58e-4 .+ 0.00163./(1.0e9diameter) # 1802 01
end
if GT_loaded
    loss_rate = wall_rate_expected[idx_p_simu]
end
x00[R_loss] .= wall_rate_inv(loss_rate./gamma0)
x00[R_nuc]  .= 0.0

# diagonal of the covariance matrix of the initial state
gam00_vec[R_psd]       = (10.0*(1.0e9delta)).^2                    # concentration for a uniform density of 10 cm^{-3} nm^{-1}
gam00_vec[R_cond_all] .= 1.0^2                                     # dimensionless growth
gam00_vec[R_loss]     .= (0.1wall_rate_inv(loss_rate./gamma0)).^2  # dimensionless losses
gam00_vec[R_nuc]      .= 1.0^2                                     # dimensionless nucleaion rate

############################################################
###             fix interval smoother                   ####
############################################################

deltaT = @elapsed myWSKF,x_smo_all,o_smo_all = EKF.KF_FIS(x00,gam00_vec,var_model,var_noise,t0,myWSKF, y_all,n_samp)
println(@sprintf "time to run the FIKS estimation: %fs" deltaT)

# alias (just to make the code more readable, especially in my_display.jl)
x_fil_all = myWSKF.x_fil_all
x_pre_all = myWSKF.x_pre_all
o_fil_all = myWSKF.o_fil_all
o_pre_all = myWSKF.o_pre_all
gam_v_all = myWSKF.gam_v_all

deltaT = @elapsed for i in 1:n_samp
    # condensation rate
    percentiles_cond_smo[:,i,:] = [CGR(x_smo_all[R_cond,i]-sqrt.(diag(o_smo_all[R_cond,R_cond,i]))); CGR(x_smo_all[R_cond,i]+sqrt.(diag(o_smo_all[R_cond,R_cond,i])))]
    # wall loss
    for j in 1:nbin
        percentiles_wall_smo[j,i,:] = [wall_rate(x_smo_all[R_loss_init-1+j,i]-sqrt(o_smo_all[R_loss_init-1+j,R_loss_init-1+j,i])); wall_rate(x_smo_all[R_loss_init-1+j,i]+sqrt(o_smo_all[R_loss_init-1+j,R_loss_init-1+j,i]))]
    end
    # nucleation rate
    percentiles_nuc_smo[i,:] = StochProc.percentile_estimation(Nucleation_rate,x_smo_all[R_nuc_init,i],sqrt(o_smo_all[R_nuc_init,R_nuc_init,i]),percentiles,Ns_per) # not quite sure this method is still relevant
end
println(@sprintf "time to compute the physically relevant errors from the estimated covariances: %fs" deltaT)

# print some constants
println(" ")
println("Normalization constants used in the GDE")
println(" ")
println("condensation growth rate,  GR0: ", GR0,    " m.s^-1")
println("condensation growth rate,   J0: ", J0,     " #.cm^-3.s^-1")
println("wall losses rate,       gamma0: ", gamma0, " s^-1")
println("number concentration,       x0: ", x0,     " #.cm^-3")
println("characteristic time,        t0: ", t0,     " s")
println("time step,                  dt: ", dt,     " s")


# start plotting
if FLAG_1952_02
    deltaT = @elapsed include("plot_1952_02.jl")
    println(@sprintf "time to plot the results: %fs" deltaT)
end
if FLAG_1802_01
    deltaT = @elapsed include("plot_1802_01.jl")
    println(@sprintf "time to plot the results: %fs" deltaT)
end
if GT_loaded
    deltaT = @elapsed include("plot_0000.jl")
    println(@sprintf "time to plot the results: %fs" deltaT)
end
if FLAG_1906_03
    deltaT = @elapsed include("plot_1906_03.jl")
    println(@sprintf "time to plot the results: %fs" deltaT)
end

# plot more figures if required
if PLOT_EVEN_MORE
    deltaT = @elapsed include("plot_operator.jl")
    println(@sprintf "time to plot measurement operator: %fs" deltaT)
end

