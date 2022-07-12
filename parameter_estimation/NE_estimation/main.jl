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

# NOTE: you need to run the data_simulation/nucleation_event/main.jl script to simulate the Ground Truth

##################################################
###    read data and init some variabeles      ###
##################################################
COAGULATION = true;
COAGULATION_GAIN = (true & COAGULATION)

# loading the ground truth or not
GT_loaded = false #WARNING need more memory for this!!!


# data with a different measurement model
input_folder = "../../data_simulation/nucleation_event/data/CPC_5_L_min/"
folder = "coagulation_loss_and_gain/multicharge/CPC_5_L_min/bis_with_pgf/x0_smps_100/"
V_cm3_sample = 90.0
# input_folder = "../../data_simulation/nucleation_event/data/CPC_005_L_min/"
# folder = "coagulation_loss_and_gain/multicharge/CPC_005_L_min/bis_with_pgf/x0_smps_100/"
# V_cm3_sample = 0.9

##WARNING: defining a threshold number concentration to account for the measurement model errors, the stochasticity not steming from the counting process
x0_dmps = 100.0 # 100.0 # 1000.0 # 10000.0 #  for the case of bad error modelling

SAVE_PGF = false

println("creating path: ", folder)
mkpath(folder);

# time
df_t_samp = CSV.File(string(input_folder,"time_hour.csv"); delim=",", header=false) |> DataFrame
t_samp = dropdims(3600.0*Matrix{Cdouble}(df_t_samp),dims=1);
n_samp = length(t_samp)         # number of time sample
dt = t_samp[2]-t_samp[1]



# diameters
df_diameter_data = CSV.File(string(input_folder,"diameter.csv"); delim=",", header=false) |> DataFrame
diameter_data = dropdims(Matrix{Cdouble}(df_diameter_data),dims=1);
volume_data = (pi/6.0)*(diameter_data.^3)
cst_v_data = mean(volume_data[2:end]./volume_data[1:end-1])

# particle concentration
df_y_all = CSV.File(string(input_folder,"particle_conc_with_noise.csv"); delim=",", header=false) |> DataFrame
y_all = Matrix{Cdouble}(Matrix{Cdouble}(df_y_all)');
min_pc,max_pc = extrema(y_all)
s = @sprintf "number concentration (%1.2e,%1.2e)" min_pc max_pc
displayLogData2D(1,t_samp/3600.0,diameter_data,y_all,max(min_pc,0.001max_pc),max_pc,_title=s,_colorbar_label="concentration [cm\$^{-3}\$]")
savefig(string(folder,"raw_data.png"))
savefig(string(folder,"raw_data.pdf"))
if SAVE_PGF
    savefig(string(folder,"raw_data.pgf"))
end
close("all")


# load the measurement model
df_H_avg = CSV.File(string(input_folder,"H_avg.csv"); header=false) |> DataFrame
H_avg = Matrix{Cdouble}(df_H_avg); # convert(Matrix{Cdouble},df_H_avg);
# H_avg = H_avg[49:92,49:92];
min_H,max_H = extrema(H_avg)
s = @sprintf "channel efficiency (%1.2e,%1.2e)" min_H max_H
displayLogData2D(1111,1.0e9diameter_data,1.0e9diameter_data,H_avg,max(min_H,0.001max_H),max_H,_title=s,_colorbar_label="channel efficiency []")
xscale("log")
xlabel("diameter [nm]")
ylabel("channel center [nm]")
savefig(string(folder,"H_avg.png"))
savefig(string(folder,"H_avg.pdf"))
if SAVE_PGF
    savefig(string(folder,"H_avg.pgf"))
end
close("all")


# set the dimensions: try to set a fixed state size
nbin  = length(diameter_data)                              # Number of particle size bins
model_dim = nbin+2+nbin+2  # +4                      # [],        dimension of the model
meas_dim  = length(diameter_data)             # [],        dimension of the observation space
d0 = diameter_data[1] # 1.5e-9 #
v0 = (pi/6.0)*(d0^3)

# compute the diameter, volume and constant ratio for the model state
cst_v = (volume_data[end]/v0)^(1.0/(nbin-1.0))
cst_r = cst_v^(1.0/3.0)
r_vec = cst_r.^(collect(0:nbin-1));
volume = Array{Cdouble,1}(undef,nbin);
volume = v0*(cst_v.^(collect(0:nbin-1)));
diameter = ((6.0/pi)*volume).^(1.0/3.0);
delta = diameter*(sqrt(cst_r)-1.0/sqrt(cst_r));


# load the expected values of the parameters
if GT_loaded
    df_tp                 = CSV.File(string(input_folder,"time_hour_parameter.csv"); delim=",", header=false) |> DataFrame
    tp                    = dropdims(Matrix{Cdouble}(df_tp),dims=1);
    df_dp                 = CSV.File(string(input_folder,"diameter_parameter.csv"); delim=",", header=false) |> DataFrame
    dp                    = dropdims(Matrix{Cdouble}(df_dp),dims=1);
    df_condensation_rate  = CSV.File(string(input_folder,"condensation_rate.csv"); delim=",", header=false) |> DataFrame
    condensation_rate_all = Matrix{Cdouble}(Matrix{Cdouble}(df_condensation_rate)');
    df_nucleation_rate    = CSV.File(string(input_folder,"nucleation_rate.csv"); delim=",", header=false) |> DataFrame
    nucleation_rate_all   = dropdims(Matrix{Cdouble}(df_nucleation_rate),dims=1);
    df_wall_rate          = CSV.File(string(input_folder,"wall_rate_constant.csv"); delim=",", header=false) |> DataFrame
    wall_rate_expected    = dropdims(Matrix{Cdouble}(df_wall_rate),dims=1);
    df_PSD_simu           = CSV.File(string(input_folder,"simulated_psd.csv"); delim=",", header=false) |> DataFrame
    PSD_simu              = Matrix{Cdouble}(Matrix{Cdouble}(df_PSD_simu)');
    cst_r_p = mean(dp[2:end]./dp[1:end-1])
    delta_p = dp.*(sqrt(cst_r_p)-1.0/sqrt(cst_r_p))
    idx_t = Array{Int64}(undef,n_samp)
    for i in 1:n_samp
        idx_t[i] = findfirst((tp*3600.0).>=t_samp[i])
    end
    idx_d = Array{Int64}(undef,length(diameter_data))
    for i in 1:length(diameter_data)
        idx_d[i] = findfirst(dp.>=diameter_data[i])
    end
end



# the identity matrix
Id=Matrix{Cdouble}(I,model_dim,model_dim);






##################################################
####   load the Kalman Filter implementation   ###
##################################################
using EKF
myWSKF = wsKF(model_dim,meas_dim,n_samp,false,true);
myWSKF.t_samp = t_samp;








##################################################
###             define the model               ###
##################################################

# definition of the characteristic constants (set parameter for nondimensionalization)
GR0    = 5.0*2.8e-13              # the maximum condensation rate observed in the data or a fraction of the maximum
gamma0 = 3.0e-5                   # the maximum wall loss rate observed in the data or a fraction of the maximum
t0     = 2.0*3600.0               # the characteristic time t0=1.0/alpha0
x0     = 1.0e2                    # a fixed number concentration 10^4 #cm^-3
J0     = 0.2                      # normalizing constant for the nucleation rate

# the evolution and measurement models
include("model_mod.jl") #TODO: go to StochProc and create Toeplitz matrix generator to remove SpecialMatrices package





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
gam00_vec_model = var_noise_diag.*(x0*var_noise_diag.>sig_model)+(sig_model/x0)*(x0*var_noise_diag.<=sig_model) # just to make sure that the initial covariance is larger than the covariance of the model and strictly positive
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

