# TODO: implement same estimation as for the DMA train

# This file is licensed under the MIT "Expat" License:

# Copyright (c) 2020: Matthew Ozon.

# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:

# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.



using PyPlot
# set matplotlib to allow fior the use of LaTeX formula in the graphs (side effect: the font is set to LaTeX's default)
rc("text", usetex=true)
using myPlot
using CSV
using DataFrames
using Printf
using LinearAlgebra
using Statistics
using XLSX

using AeroMeas        # measurement simulation


# data with a different measurement model
input_folder = "../../data_vine/"

# load data 
xf_exp_data = XLSX.readxlsx(joinpath(input_folder,"concentration_data_vine.xlsx")); 


df_data = DataFrame(XLSX.gettable(xf_exp_data["concentration"]))
raw_data = Matrix{Cdouble}(Matrix{Cdouble}(df_data)');

df_diameter = DataFrame(XLSX.gettable(xf_exp_data["centroid_diamter"]))
diameter_data = Vector{Cdouble}(df_diameter.values)


df_time = DataFrame(XLSX.gettable(xf_exp_data["measurement_time"]))
t_samp = Vector{Cdouble}(df_time.values)
n_samp = length(t_samp)         # number of time sample
dt = t_samp[2]-t_samp[1]


# data of interest:
raw_data_interest = raw_data[143:end,:] # take only the data after the first 143 samples

# diameter data of interest:
diameter_data_interest = diameter_data[143:end]
cst_r_data = mean(diameter_data_interest[2:end]./diameter_data_interest[1:end-1])
volume_data_interest = (pi/6.0)*(diameter_data_interest.^3)
cst_v_data = mean(volume_data_interest[2:end]./volume_data_interest[1:end-1])
min_pc,max_pc = extrema(raw_data_interest)
s = @sprintf "number concentration (%1.2e,%1.2e)" min_pc max_pc
displayLogData2D(1,t_samp/3600.0,diameter_data_interest,raw_data_interest,max(min_pc,0.001max_pc),max_pc,_title=s,_colorbar_label="concentration [cm\$^{-3}\$]")
tight_layout()



# compute the measurement operator


###############################################################
#   simulation of the evolution of the number concentration   #
###############################################################
# diameter
nbin = 10000 # 2500                                                 # number of discretization point (size wise)
# d0 = 13.85e-9 #14.1e-9                                      # smallest size in simulation
d0 = 12.88e-9 #14.1e-9                                      # smallest size in simulation
d0 = diameter_data_interest[1]/cst_r_data
dmax = 1000.0e-9 #  1.0e-6                                  # biggest size in simulation
cst_r = (dmax/d0)^(1.0/(nbin-1.0))                          # constant ratio of the log-scale (actually exp-scale, but that's what people are used to)
size_scale = cst_r.^(collect(0:nbin-1))  
d = d0*size_scale                                           # discretization points (centroid of the virtual bins)
delta = (cst_r-1.0).*d0.*cst_r.^(collect(0:nbin-1).-0.5)    # width of each virtual bin

# set the channel center and width for which the initial state had been measured
# nbin_meas   = 111                                       # number of channels
nbin_meas = length(diameter_data_interest) # 230
d0_meas     = diameter_data_interest[1] # 13.0e-9 # 14.1e-9                                   # smallest bin centroid
dmax_meas   = diameter_data_interest[end] # 0.79863e-6                                 # biggest bin centroid
cst_r_meas  = (dmax_meas/d0_meas)^(1.0/(nbin_meas-1))   # constant ratio of the measured bin
d_meas = Array{Cdouble}(undef,nbin_meas)                # centroid of the measured channels
d_meas[1]=d0_meas
for i in 2:nbin_meas
    d_meas[i]=d_meas[i-1]*cst_r_meas
end
delta_meas = (cst_r_meas-1.0)*d0_meas*cst_r_meas.^(collect(0:nbin_meas-1).-0.5); # width of the measured channels

###########################################################
#                  convert to density                     #
###########################################################
# PSD is the particle size distribution (dN/dv [# cm^{-3} nm^{-1}])



###########################################################
#         simulation of the measurement process           #
###########################################################

# first create the measurement model
q_a = 0.3;             # [L min^{-1}], aerosol sample flow in the DMA
q_sh = 4.8;            # [L min^{-1}], sheath flow

# a few constant for the SMPS
s50imp=1.0e-6;         # [m],  cut off size of the impactor
delta50imp=0.1e-6;     # [m],  selectivity of the impactor WARNING: this parameter is not given in the data, it needs to be checked
s50cpc=4.0e-9          # [m],  cut off size of the CPC WARNING: this parameter is not given in the data, it needs to be checked
delta50cpc=2.0e-9      # [m],  selectivity of the CPC... more like the spread of the CPC's attenuation at the lower end of the spectrum WARNING: this parameter is not given in the data, it needs to be checked
T0 = 293.0             # [K],  temperature of the carrier gas
Pr0=1.0e5              # [Pa], pressure of the carrier gas
Nq0 = -10 # -1 #       # [],   signed number of charges
# compute the transfer function (the efficiency of the channels)
H = SMPS3936_transfer_function(d,d_meas,cst_r_meas;s50imp=s50imp,delta50imp=delta50imp,s50cpc=s50cpc, delta50cpc=delta50cpc,T=T0,Pr=Pr0,Nq=Nq0,q_a=q_a,q_sh=q_sh);



H_ds = H.*delta'; # middle Riemann integration
nbin = 40 
cst_r_disc  = (dmax_meas/d0_meas)^(1.0/(nbin-1))
dp_disc = d0_meas*cst_r_disc.^(collect(0:nbin-1))
# idx_min = zeros(Int64,length(d_meas));
# idx_max = zeros(Int64,length(d_meas));
idx_min = zeros(Int64,nbin);
idx_max = zeros(Int64,nbin);
# the simulation size discretization start at a bigger size than the smallest size measured
idx_min[1] = 1
idx_max[1] = findfirst(d.>=(d_meas[1]*sqrt(cst_r_meas)))
# for i in 2:nbin_meas
#     idx_min[i] = findlast(d.<(d_meas[i]/sqrt(cst_r_meas)))
#     idx_max[i] = findfirst(d.>=(d_meas[i]*sqrt(cst_r_meas)))
# end
idx_min[1] = 1
idx_max[1] = findfirst(d.>=(dp_disc[1]*sqrt(cst_r_disc)))
for i in 2:nbin
    idx_min[i] = findlast(d.<(dp_disc[i]/sqrt(cst_r_disc)))
    idx_max[i] = findfirst(d.>=(dp_disc[i]*sqrt(cst_r_disc)))
end
# for each channel, compute the average efficiency
H_avg = zeros(Float64,nbin_meas,nbin);
for i in 1:nbin_meas
    for j in 1:nbin
        H_avg[i,j] = sum(H_ds[i,idx_min[j]:idx_max[j]])/delta_meas[j]
    end
end

min_H_avg,max_H_avg = extrema(H_avg); # extrema(H.*delta_s'); #
# s = @sprintf "SMPS transfer function avg (%1.2f,%1.2f)" min_H_avg max_H_avg
s = @sprintf "transfer function avg (%1.2f,%1.2f)" min_H_avg max_H_avg
# displayLogData2D(3000,1.0e9.*d_meas,1.0e9.*d_meas,H_avg,max(0.001max_H_avg,min_H_avg),max_H_avg,_title=s,_colorbar_label="channel efficiency []")
displayLogData2D(3000,1.0e9.*dp_disc,1.0e9.*d_meas,H_avg,max(0.001max_H_avg,min_H_avg),max_H_avg,_title=s,_colorbar_label="channel efficiency []")
xscale("log")
xlabel("diameter [nm]")
ylabel("channel center [nm]")
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
# if SAVE_FIG
#     savefig("kernel_function_contour_plot_avg.png")
#     savefig("kernel_function_contour_plot_avg.pdf")
# end

figure();
semilogx(1.0e9.*dp_disc,H_avg')
xlabel("diameter [nm]")
ylabel("channel efficiency []")
title("SMPS transfer functions: average")
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
# if SAVE_FIG
#     savefig("kernel_function_avg.png")
#     savefig("kernel_function_avg.pdf")
# end



# estimate the concentration from the raw data using the channel efficiency
# D = diagm(0=>-ones(Float64,nbin_meas-2), 1=>2ones(Float64,nbin_meas-2),2=>-ones(Float64,nbin_meas-2))[1:end-2,:];
D = diagm(0=>-ones(Float64,nbin-2), 1=>2ones(Float64,nbin-2),2=>-ones(Float64,nbin-2))[1:end-2,:];
y_data = inv(H_avg'*H_avg + (1.0e1)^2*D'*D)*H_avg'*raw_data_interest[:,100]

# at this stage, the reconstructed concentration is in [# cm^{-3}] (the quantity of interest is more likely the log density)

figure(); 
plot(diameter_data_interest,raw_data_interest[:,1],label="concentration [cm\$^{-3}\$]"); 
# plot(diameter_data_interest,y_data,label="estimated concentration [cm\$^{-3}\$]"); 
plot(diameter_data_interest,H_avg*y_data,label="estimated concentration [cm\$^{-3}\$]"); 
xscale("log")


# compute for all instant in time the reconstructed concentration
# y_data_all = zeros(Float64,nbin_meas,size(raw_data_interest,2));
y_data_all = zeros(Float64,nbin,size(raw_data_interest,2));
for t in 1:size(raw_data_interest,2)
    y_data_all[:,t] = inv(H_avg'*H_avg + (1.0e1)^2*D'*D)*H_avg'*raw_data_interest[:,t]
end
y_data_all[y_data_all.<0.0] .= 0.0


min_pc,max_pc = extrema(y_data_all)
s = @sprintf "number concentration (%1.2e,%1.2e)" min_pc max_pc
# displayLogData2D(1000,t_samp/3600.0,diameter_data_interest,y_data_all,max(min_pc,0.001max_pc),max_pc,_title=s,_colorbar_label="concentration [cm\$^{-3}\$]")
displayLogData2D(1000,t_samp/3600.0,dp_disc,y_data_all,max(min_pc,0.001max_pc),max_pc,_title=s,_colorbar_label="concentration [cm\$^{-3}\$]")
tight_layout()




# set the dimensions: try to set a fixed state size
COAGULATION = false
# nbin  = 40 #  length(diameter_data_interest)               # Number of particle size bins
model_dim = nbin+2+nbin+2  # +4                      # [],        dimension of the model
meas_dim  = length(diameter_data_interest)             # [],        dimension of the observation space
d0 = diameter_data_interest[1] # 1.5e-9 #
dmax = diameter_data_interest[end]
# v0 = (pi/6.0)*(d0^3)

# compute the diameter, volume and constant ratio for the model state
# cst_v = (volume_data[end]/v0)^(1.0/(nbin-1.0))
# cst_r = cst_v^(1.0/3.0)
cst_r = (dmax/d0)^(1.0/(nbin-1.0))
r_vec = cst_r.^(collect(0:nbin-1));
diameter = d0*(cst_r.^(collect(0:nbin-1)));
delta = diameter*(sqrt(cst_r)-1.0/sqrt(cst_r));

# the identity matrix
Id=Matrix{Cdouble}(I,model_dim,model_dim);


##################################################
###             define the model               ###
##################################################

# definition of the characteristic constants (set parameter for nondimensionalization)
GR0    = 5.0*2.8e-13              # the maximum condensation rate observed in the data or a fraction of the maximum
gamma0 = 3.0e-5                   # the maximum wall loss rate observed in the data or a fraction of the maximum
t0     = 2.0*3600.0               # the characteristic time t0=1.0/alpha0
x0     = 1.0e2                    # a fixed number concentration 10^4 #cm^-3
J0     = 0.2                      # normalizing constant for the nucleation rate

# if false
# TODO: define y_all and COAGULATION_GAIN
COAGULATION_GAIN = false
y_all = raw_data_interest
need_padding = true
if need_padding
    dt_data = dt
    n_samp_data = n_samp
    pad_factor = floor(Int64,0.5dt_data/30.0)
    # allocate a new time array
    n_samp = (n_samp_data-1)*(pad_factor+1)+1 # (n_samp_data)*(pad_factor+1)+1 # 
    t_samp_start = t_samp[1];
    t_samp = t_samp_start .+ (dt_data/(pad_factor+1))*collect(0:1:((n_samp_data-1)*(pad_factor+1)));
    # t_samp = t_samp_start .+ (dt_data/(pad_factor+1))*collect(0:1:((n_samp_data)*(pad_factor+1)));
    dt = dt_data/(pad_factor+1)
    # allocate a new data array
    y_all_data = copy(y_all)
    y_all = zeros(Cdouble,meas_dim,n_samp)
    y_all[:,1:pad_factor+1:end] = y_all_data
end

##################################################
####   load the Kalman Filter implementation   ###
##################################################
using EKF
myWSKF = wsKF(model_dim,meas_dim,n_samp,false,true);
myWSKF.t_samp = t_samp;

# the evolution and measurement models
V_cm3_sample = 90.0
x0_dmps = 1.0 # 100.0
include("model_mod.jl") #TODO: go to StochProc and create Toeplitz matrix generator to remove SpecialMatrices package


# TODO: check the following code


##################################################
###    initialization of the Kalman Filter     ###
##################################################
# initial state variables
x00            = Array{Cdouble,1}(undef,model_dim);
gam00_vec      = Array{Cdouble,1}(undef,model_dim);
var_model      = Array{Cdouble,2}(undef,model_dim,model_dim);
var_noise_diag = Array{Cdouble,1}(undef,model_dim);
var_noise      = Array{Cdouble,2}(undef,model_dim,model_dim);

# how much do we trust the model
fill!(var_model,0.0)
var_model[R_psd,R_psd] = gamma_psd;

# wall deposition rate model (for the first bin, the rate varies in [0,1])
var_model[R_cond,R_cond] = gamma_cond
var_model[R_loss,R_loss] = gamma_loss
var_model[R_nuc,R_nuc]   = gamma_nuc

#how much do we trust the data
var_noise_diag = data_var(y_all[:,1]) # approximation of the actual covariance of the size distribution: poisson by gaussian
var_noise = diagm(var_noise_diag) 

# it is important to initialize the measurement model beforehand so that we can compute the initial guess
myWSKF.H_me = set_measurement_jacobian!(myWSKF.H_me) # time independent, must be initialized 

# draw the initial state
x00[R_psd] = (myWSKF.H_me'*(y_all[:,1]/x0))[R_psd] # WARNING: this might not be the best initial guess, if one knows a better initial concentration state, one should use it.
x00[R_cond] .= 0.0 
x00[R_loss] .= 2.0 
x00[R_nuc]  .= 0.0 

# initial model uncertainties
sig_model  = 1.0
# gam00_vec_model = var_noise_diag.*(x0*var_noise_diag.>sig_model)+(sig_model/x0)*(x0*var_noise_diag.<=sig_model) # just to make sure that the initial covariance is larger than the covariance of the model and strictly positive
gam00_vec_model = 200.0*ones(Float64,nbin)
gam00_vec[R_psd]   = 4.0gam00_vec_model
gam00_vec[R_cond] .= 1.0^2
gam00_vec[R_loss] .= (40.0gamma0*t0)^2
gam00_vec[R_nuc]  .= 1.0^2

############################################################
###             fix interval smoother                   ####
############################################################
# variables used for the KF (NOTE: there's no more need to allocate the following variables, it could be commented out)
x_fil_all = Array{Cdouble,2}(undef,model_dim,n_samp);
x_pre_all = Array{Cdouble,2}(undef,model_dim,n_samp);
o_fil_all = Array{Cdouble,3}(undef,model_dim,model_dim,n_samp);
o_pre_all = Array{Cdouble,3}(undef,model_dim,model_dim,n_samp);
x_smo_all = similar(myWSKF.x_fil_all);
o_smo_all = similar(myWSKF.o_fil_all);

# run the FIKS
@elapsed myWSKF,x_smo_all,o_smo_all = KF_FIS(x00,gam00_vec,var_model,var_noise,t0,myWSKF, y_all/x0,n_samp)

# alias
x_fil_all = myWSKF.x_fil_all;
x_pre_all = myWSKF.x_pre_all;
o_fil_all = myWSKF.o_fil_all;
o_pre_all = myWSKF.o_pre_all;



# percentiles results
percentiles = [15;85]
percentiles_cond_pre       = Array{Cdouble,2}(undef,n_samp,2);
percentiles_cond_fil       = Array{Cdouble,2}(undef,n_samp,2);
percentiles_cond_smo       = Array{Cdouble,2}(undef,n_samp,2);
percentiles_wall_pre       = Array{Cdouble,3}(undef,nbin,n_samp,2);
percentiles_wall_fil       = Array{Cdouble,3}(undef,nbin,n_samp,2);
percentiles_wall_smo       = Array{Cdouble,3}(undef,nbin,n_samp,2);
percentiles_nuc_pre        = Array{Cdouble,2}(undef,n_samp,2);
percentiles_nuc_fil        = Array{Cdouble,2}(undef,n_samp,2);
percentiles_nuc_smo        = Array{Cdouble,2}(undef,n_samp,2);
percentile_smo             = Array{Cdouble,3}(undef,n_samp,2,model_dim);

# compute percentiles or +/- 1 sigma
@elapsed for i in 1:n_samp
    # condensation rate
    percentiles_cond_pre[i,:] = [CGR(x_pre_all[R_cond_init,i]-sqrt(o_pre_all[R_cond_init,R_cond_init,i])); CGR(x_pre_all[R_cond_init,i]+sqrt(o_pre_all[R_cond_init,R_cond_init,i]))]
    percentiles_cond_fil[i,:] = [CGR(x_fil_all[R_cond_init,i]-sqrt(o_fil_all[R_cond_init,R_cond_init,i])); CGR(x_fil_all[R_cond_init,i]+sqrt(o_fil_all[R_cond_init,R_cond_init,i]))]
    percentiles_cond_smo[i,:] = [CGR(x_smo_all[R_cond_init,i]-sqrt(o_smo_all[R_cond_init,R_cond_init,i])); CGR(x_smo_all[R_cond_init,i]+sqrt(o_smo_all[R_cond_init,R_cond_init,i]))]
    # wall loss
    for j in 1:nbin
        percentiles_wall_pre[j,i,:] = [wall_rate(x_pre_all[R_loss_init-1+j,i]-sqrt(o_pre_all[R_loss_init-1+j,R_loss_init-1+j,i])); wall_rate(x_pre_all[R_loss_init-1+j,i]+sqrt(o_pre_all[R_loss_init-1+j,R_loss_init-1+j,i]))]
        percentiles_wall_fil[j,i,:] = [wall_rate(x_fil_all[R_loss_init-1+j,i]-sqrt(o_fil_all[R_loss_init-1+j,R_loss_init-1+j,i])); wall_rate(x_fil_all[R_loss_init-1+j,i]+sqrt(o_fil_all[R_loss_init-1+j,R_loss_init-1+j,i]))]
        percentiles_wall_smo[j,i,:] = [wall_rate(x_smo_all[R_loss_init-1+j,i]-sqrt(o_smo_all[R_loss_init-1+j,R_loss_init-1+j,i])); wall_rate(x_smo_all[R_loss_init-1+j,i]+sqrt(o_smo_all[R_loss_init-1+j,R_loss_init-1+j,i]))]
    end
    # nucleation rate
    percentiles_nuc_pre[i,:] = [Nucleation_rate(x_pre_all[R_nuc_init,i]-sqrt(o_pre_all[R_nuc_init,R_nuc_init,i])); Nucleation_rate(x_pre_all[R_nuc_init,i]+sqrt(o_pre_all[R_nuc_init,R_nuc_init,i]))]
    percentiles_nuc_fil[i,:] = [Nucleation_rate(x_fil_all[R_nuc_init,i]-sqrt(o_fil_all[R_nuc_init,R_nuc_init,i])); Nucleation_rate(x_fil_all[R_nuc_init,i]+sqrt(o_fil_all[R_nuc_init,R_nuc_init,i]))]
    percentiles_nuc_smo[i,:] = [Nucleation_rate(x_smo_all[R_nuc_init,i]-sqrt(o_smo_all[R_nuc_init,R_nuc_init,i])); Nucleation_rate(x_smo_all[R_nuc_init,i]+sqrt(o_smo_all[R_nuc_init,R_nuc_init,i]))]

end


# compute the evolution of the size distribution using the estimated parameters
x_sim_all = Array{Cdouble,2}(undef,nbin,n_samp);
x_sim_all[:,1] = x_smo_all[R_psd,1]
for u in 2:n_samp
    x_sim_all[:,u] = iter!(dx_coag,dx_cond,dx_nuc,dx_wall,ws,x_sim_all[:,u-1],t_samp[u]-t_samp[u-1],CGR(x_smo_all[R_cond_init,u-1]),Nucleation_rate(x_smo_all[R_nuc_init,u-1]),wall_rate(x_smo_all[R_loss,u-1]))
end



# display and save figures
try
    include("my_display.jl")
catch msgError
    println(msgError)
    println("unable to print figure for folder: ", folder)
    println(" ")
    println(" ")
end
