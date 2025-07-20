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
# set the font in the plots
# rc("font",family="serif",serif="Computer Modern Roman")
rc("text", usetex=true)
using myPlot
using CSV
using Printf
using LinearAlgebra
using Statistics
using DataFrames

##################################################
###    read data and init some variables       ###
##################################################

# NOTE: you need to run the data_simulation/steady_state/main.jl script to simulate the Ground Truth

# used for the evolution model
COAGULATION = true;
COAGULATION_GAIN = (true & COAGULATION)

# loading the ground truth or not
GT_loaded = false #WARNING need more memory for this!!!

# data with a different measurement model
# input_folder = "/home/user_name/Data/steady_state_SMPS3936/no_coagulation/singlecharge/cloud_50_bins_J_5/steady_state_CPC_2000_cm3/boundary/wall/120s/"
# folder = "no_coagulation/singlecharge/cloud_50_bins_J_5/steady_state_CPC_2000_cm3/boundary/wall/120s/"
# V_cm3_sample = 2000.0
# input_folder = "/home/user_name/Data/steady_state_SMPS3936/no_coagulation/singlecharge/cloud_50_bins_J_5/steady_state_CPC_200_cm3/boundary/wall/120s/"
# folder = "no_coagulation/singlecharge/cloud_50_bins_J_5/steady_state_CPC_200_cm3/boundary/wall/120s/"
# V_cm3_sample = 200.0
# input_folder = "/home/user_name/Data/steady_state_SMPS3936/no_coagulation/singlecharge/cloud_50_bins_J_5/steady_state_CPC_20_cm3/boundary/wall/120s/"
# folder = "no_coagulation/singlecharge/cloud_50_bins_J_5/steady_state_CPC_20_cm3/boundary/wall/120s/"
# V_cm3_sample = 20.0
# input_folder = "/home/user_name/Data/steady_state_SMPS3936/no_coagulation/singlecharge/cloud_50_bins_J_5/steady_state_CPC_2_cm3/boundary/wall/120s/"
# folder = "no_coagulation/singlecharge/cloud_50_bins_J_5/steady_state_CPC_2_cm3/boundary/wall/120s/"
# V_cm3_sample = 2.0

# input_folder = "/home/user_name/Data/steady_state_SMPS3936/no_coagulation/monocharge/cloud_50_bins_J_5/steady_state_CPC_2000_cm3/boundary/wall/120s/"
# folder = "no_coagulation/monocharge/cloud_50_bins_J_5/steady_state_CPC_2000_cm3/boundary/wall/120s/"
# V_cm3_sample = 2000.0
# input_folder = "/home/user_name/Data/steady_state_SMPS3936/no_coagulation/monocharge/cloud_50_bins_J_5/steady_state_CPC_200_cm3/boundary/wall/120s/"
# folder = "no_coagulation/monocharge/cloud_50_bins_J_5/steady_state_CPC_200_cm3/boundary/wall/120s/"
# V_cm3_sample = 200.0
# input_folder = "/home/user_name/Data/steady_state_SMPS3936/no_coagulation/monocharge/cloud_50_bins_J_5/steady_state_CPC_20_cm3/boundary/wall/120s/"
# folder = "no_coagulation/monocharge/cloud_50_bins_J_5/steady_state_CPC_20_cm3/boundary/wall/120s/"
# V_cm3_sample = 20.0
# input_folder = "/home/user_name/Data/steady_state_SMPS3936/no_coagulation/monocharge/cloud_50_bins_J_5/steady_state_CPC_2_cm3/boundary/wall/120s/"
# folder = "no_coagulation/monocharge/cloud_50_bins_J_5/steady_state_CPC_2_cm3/boundary/wall/120s/"
# V_cm3_sample = 2.0

## no coagulation and multicharge
# input_folder = "/home/user_name/Data/steady_state_SMPS3936/no_coagulation/multicharge/cloud_50_bins_J_5/steady_state_CPC_2000_cm3/boundary/wall/120s/"
# folder = "no_coagulation/multicharge/cloud_50_bins_J_5/steady_state_CPC_2000_cm3/boundary/wall/120s/"
# V_cm3_sample = 2000.0
# input_folder = "/home/user_name/Data/steady_state_SMPS3936/no_coagulation/multicharge/cloud_50_bins_J_5/steady_state_CPC_200_cm3/boundary/wall/120s/"
# folder = "no_coagulation/multicharge/cloud_50_bins_J_5/steady_state_CPC_200_cm3/boundary/wall/120s/"
# V_cm3_sample = 200.0
# input_folder = "/home/user_name/Data/steady_state_SMPS3936/no_coagulation/multicharge/cloud_50_bins_J_5/steady_state_CPC_20_cm3/boundary/wall/120s/"
# folder = "no_coagulation/multicharge/cloud_50_bins_J_5/steady_state_CPC_20_cm3/boundary/wall/120s/"
# V_cm3_sample = 20.0
# input_folder = "/home/user_name/Data/steady_state_SMPS3936/no_coagulation/multicharge/cloud_50_bins_J_5/steady_state_CPC_2_cm3/boundary/wall/120s/"
# folder = "no_coagulation/multicharge/cloud_50_bins_J_5/steady_state_CPC_2_cm3/boundary/wall/120s/"
# V_cm3_sample = 2.0
## coagulation loss and gain, and multicharge
# input_folder = "/home/user_name/Data/steady_state_SMPS3936/coagulation_loss_and_gain/multicharge/cloud_50_bins_J_5/steady_state_CPC_2000_cm3/boundary/wall/120s/"
# folder = "coagulation_loss_and_gain/multicharge/cloud_50_bins_J_5/steady_state_CPC_2000_cm3/boundary/wall/120s/bis/"
# V_cm3_sample = 2000.0
# input_folder = "/home/user_name/Data/steady_state_SMPS3936/coagulation_loss_and_gain/multicharge/cloud_50_bins_J_5/steady_state_CPC_200_cm3/boundary/wall/120s/"
# folder = "coagulation_loss_and_gain/multicharge/cloud_50_bins_J_5/steady_state_CPC_200_cm3/boundary/wall/120s/bis_with_pgf/x0_smps_100/"
# V_cm3_sample = 200.0
# input_folder = "/home/user_name/Data/steady_state_SMPS3936/coagulation_loss_and_gain/multicharge/cloud_50_bins_J_5/steady_state_CPC_20_cm3/boundary/wall/120s/"
# folder = "coagulation_loss_and_gain/multicharge/cloud_50_bins_J_5/steady_state_CPC_20_cm3/boundary/wall/120s/"
# V_cm3_sample = 20.0
# input_folder = "/home/user_name/Data/steady_state_SMPS3936/coagulation_loss_and_gain/multicharge/cloud_50_bins_J_5/steady_state_CPC_2_cm3/boundary/wall/120s/"
# folder = "coagulation_loss_and_gain/multicharge/cloud_50_bins_J_5/steady_state_CPC_2_cm3/boundary/wall/120s/bis_with_pgf/x0_smps_100/"
# V_cm3_sample = 2.0


# ##WARNING: defining a threshold number concentration to account for the measurement model errors, the stochasticity not steming from the counting process
x0_dmps = 100.0 # 100.0 # 1000.0 # 10000.0 # for the case of bad error modelling

input_folder = "../../data_simulation/steady_state/cloud_50_bins_J_5/steady_state_CPC_200_cm3/boundary/wall/120s/"
folder = "coagulation_loss_and_gain/multicharge/cloud_50_bins_J_5/steady_state_CPC_200_cm3/boundary/wall/120s/x0_smps_100/test_that_will_be_deleted/"
V_cm3_sample = 200.0

SAVE_PGF = false

println("creating path: ", folder)
mkpath(folder);

# time
df_t_samp = CSV.File(string(input_folder,"time_hour.csv"); header=false,ntasks=1)|> DataFrame
t_samp = 3600.0dropdims(Matrix{Cdouble}(df_t_samp),dims=1);
Tmax = length(t_samp)
DT = 1
idx_T = 1:DT:Tmax # 1:length(t_samp) # 600:length(t_samp)
t_samp = t_samp[idx_T];
n_samp = length(t_samp)         # number of time sample
dt = t_samp[2]-t_samp[1]


# diameters
df_diameter_data = CSV.File(string(input_folder,"diameter.csv"); header=false,ntasks=1)|> DataFrame
diameter_data = dropdims(Matrix{Cdouble}(df_diameter_data),dims=1);
idx_s = findfirst(diameter_data.>=8.5*1.0e-9)
if idx_s!=nothing
    diameter_data = diameter_data[1:idx_s];
end
volume_data = (pi/6.0)*(diameter_data.^3)
cst_v_data = mean(volume_data[2:end]./volume_data[1:end-1])

# particle concentration
df_y_all = CSV.File(string(input_folder,"particle_conc_with_noise.csv"); header=false)|> DataFrame
y_all = Matrix{Cdouble}(Matrix{Cdouble}(df_y_all)');
if idx_s!=nothing
    y_all = y_all[1:idx_s,idx_T];
else
    y_all = y_all[:,idx_T];
end
min_pc,max_pc = extrema(y_all)
s = @sprintf "number concentration (%1.2e,%1.2e)" min_pc max_pc
displayLogData2D(1,t_samp/3600.0,diameter_data,y_all,max(0.01max_pc,min_pc),max_pc,_title=s,_colorbar_label="concentration [cm\$^{-3}\$]")





# the measurement model
df_H_avg = CSV.File(string(input_folder,"H_avg.csv"); header=false)|> DataFrame
H_avg = Matrix{Cdouble}(df_H_avg);
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

# set the dimensions: try to set a fixed state size
nbin  = length(diameter_data)                              # Number of particle size bins
model_dim = nbin+2nbin+nbin+2                              # dimension of the model
meas_dim  = length(diameter_data)                          # dimension of the observation space
d0 = diameter_data[1]
v0 = (pi/6.0)*(d0^3)

# compute the diameter, volume and constant ratio for the model state
cst_v = (volume_data[end]/v0)^(1.0/(nbin-1.0))
cst_r = cst_v^(1.0/3.0)
r_vec = cst_r.^(collect(0:nbin-1));
volume = Array{Cdouble}(undef,nbin);
volume = v0*(cst_v.^(collect(0:nbin-1)));
diameter = ((6.0/pi)*volume).^(1.0/3.0);
delta = diameter*(sqrt(cst_r)-1.0/sqrt(cst_r));



# NOTE: you need to run the data_simulation/steady_state/main.jl script to simulate the Ground Truth
if GT_loaded
    df_tp                 = CSV.File(string(input_folder,"time_hour_parameter.csv"); header=false,ntasks=1)|> DataFrame
    tp                    = dropdims(Matrix{Cdouble}(df_tp),dims=1);
    df_dp                 = CSV.File(string(input_folder,"diameter_parameter.csv"); header=false,ntasks=1)|> DataFrame
    dp                    = dropdims(Matrix{Cdouble}(df_dp),dims=1);
    df_G                  = CSV.File(string(input_folder,"condensation_rate.csv"); header=false,ntasks=1)|> DataFrame
    condensation_rate_all = Matrix{Cdouble}(Matrix{Cdouble}(df_G)');
    df_J                  = CSV.File(string(input_folder,"nucleation_rate.csv"); header=false,ntasks=1)|> DataFrame
    nucleation_rate_all   = dropdims(Matrix{Cdouble}(df_J),dims=1);
    df_λ                  = CSV.File(string(input_folder,"wall_rate_constant.csv"); header=false,ntasks=1)|> DataFrame
    wall_rate_expected    = dropdims(Matrix{Cdouble}(df_λ),dims=1);
    df_psd                = CSV.File(string(input_folder,"simulated_psd.csv"); header=false,ntasks=1)|> DataFrame
    PSD_simu              = Matrix{Cdouble}(Matrix{Cdouble}(df_psd)');

    # find sampling indices
    idx_t = Array{Int64}(undef,n_samp)
    for i in 1:n_samp
        idx_t[i] = findfirst((tp*3600.0).>=t_samp[i])
    end
    idx_d = Array{Int64}(undef,length(diameter_data))
    for i in 1:length(diameter_data)
        idx_d[i] = findfirst(dp.>=diameter_data[i]/sqrt(cst_r)) # lower end of the simulation intervals
    end
    # brutally sample the parameters
    tp = tp[idx_t]
    dp = dp[idx_d]
    cst_r_p = mean(dp[2:end]./dp[1:end-1])
    delta_p = dp.*(sqrt(cst_r_p)-1.0/sqrt(cst_r_p))
    PSD_simu = PSD_simu[idx_d,idx_t];
    condensation_rate_all = condensation_rate_all[idx_d,idx_t]
    nucleation_rate_all = nucleation_rate_all[idx_t]
    wall_rate_expected = wall_rate_expected[idx_d]

    # plot the parameters used for the simulation with a good resolution
    titleCond = "condensation rate"
    cblabelCond = "condensation rate [m s\$^{-1}\$]"
    minCond = max(5.0e-14,minimum(condensation_rate_all))
    maxCond = maximum(condensation_rate_all)
    # displayLogData2D(790,tp[1:idx_t[end]],dp,condensation_rate_all[:,1:idx_t[end]],minCond,maxCond,_title=titleCond,_colorbar_label=cblabelCond)
    displayLogData2D(790,tp,1.0e9dp,condensation_rate_all,minCond,maxCond,_title=titleCond,_colorbar_label=cblabelCond)
    ylabel("diameter [nm]")
    savefig(string(folder,"cond_simu.png"))
    savefig(string(folder,"cond_simu.pdf"))
    if SAVE_PGF
        savefig(string(folder,"cond_simu.pgf"))
    end

    figure(43)
    # plot(tp[1:idx_t[end]],nucleation_rate_all[1:idx_t[end]])
    plot(tp,nucleation_rate_all)
    s = @sprintf "nucleation rate"
    title(s)
    xlabel("time [h]")
    ylabel("nucleation rate [cm\$^{-3}\$ s\$^{-1}\$]")
    grid(true,which="major",ls="-")
    grid(true,which="minor",ls="-",alpha=0.5)
    # xlim(tp[1],tp[idx_t[end]])
    xlim(tp[1],tp[end])
    savefig(string(folder,"nuc_simu.png"))
    savefig(string(folder,"nuc_simu.pdf"))
    if SAVE_PGF
        savefig(string(folder,"nuc_simu.pgf"))
    end



    figure(44)
    loglog(dp*10.0^9,wall_rate_expected)
    s = @sprintf "wall loss rate"
    title(s)
    xlabel("diameter [nm]")
    ylabel("loss rate [s\$^{-1}\$]")
    grid(true,which="major",ls="-")
    grid(true,which="minor",ls="-",alpha=0.5)
    xlim(dp[1]*10.0^9,dp[end]*10.0^9)
    ylim(2e-5,2e-3)
    savefig(string(folder,"loss_simu.png"))
    savefig(string(folder,"loss_simu.pdf"))
    if SAVE_PGF
        savefig(string(folder,"loss_simu.pgf"))
    end
end
close("all")




# the identity matrix
Id=Matrix{Cdouble}(I,model_dim,model_dim);

# compute the wall loos rate (assumed to be a known parameter of the chamber model)
if (!GT_loaded)
    wall_rate_expected = (1.31e-12*ones(Cdouble,nbin)./(d0*cst_r.^(collect(0:nbin-1))))
end



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
GR0    = 4.0*2.0*2.8e-13    # the maximum condensation rate observed in the data or a fraction of the maximum
gamma0 = 4.0e-4             # the maximum wall loss rate observed in the data or a fraction of the maximum
t0     = 20.0               # the characteristic time t0=1.0/alpha0
x0     = 4e2                # a fixed number concentration 10^4 #cm^-3
J0     = 5.0                # normalizing constant for the nucleation rate

# the evolution model
include("model_mod.jl")

x00 = Array{Cdouble}(undef,model_dim);
gam00_vec = Array{Cdouble}(undef,model_dim);
var_model = Array{Cdouble}(undef,model_dim,model_dim);
var_noise_diag = Array{Cdouble}(undef,model_dim);
var_noise = Array{Cdouble}(undef,model_dim,model_dim);


############################################################
###          initialization of the Kalman Filter         ###
############################################################
# how much do we trust the model
fill!(var_model,0.0)
var_model[R_psd,R_psd] = natasha.gamma_psd

# evolution model covariance matrix
var_model[R_cond_all,R_cond_all] = natasha.gamma_cond
var_model[R_loss,R_loss]         = natasha.gamma_loss
var_model[R_nuc,R_nuc]           = natasha.gamma_nuc

#how much do we trust the data
var_noise_diag = data_var(y_all[:,1]/x0) # approximation of the actual covariance of the size distribution: poisson by gaussian
var_noise = diagm(var_noise_diag)

# it is important to initialize the measurement model before starting so that we can compute the initial guess
myWSKF.H_me = set_measurement_jacobian!(myWSKF.H_me) # time independent, must be initialized

# draw the initial state
x00[R_psd]  .= 0.0
x00[R_cond_all] .= -2.0
x00[R_loss]  =  wall_rate_inv(wall_rate_expected/gamma0)
x00[R_nuc]  .= -2.0

# I don't know how to explain those values
gam00_vec_model       = 4.0*data_var(y_all[:,1]/x0)
gam00_vec[R_psd]      = gam00_vec_model
gam00_vec[R_cond_all] = gamma_c[floor.(Int64,collect(1:0.5:nbin+0.5))].^2
gam00_vec[R_loss]     = (0.1wall_rate_inv(wall_rate_expected/gamma0)).^2
gam00_vec[R_nuc]     .= gamma_j.^2

############################################################
###             fix interval smoother                   ####
############################################################
@elapsed myWSKF,x_smo_all,o_smo_all = KF_FIS(x00,gam00_vec,var_model,var_noise,t0,myWSKF, y_all/x0,n_samp);

# alias
x_fil_all = myWSKF.x_fil_all;
x_pre_all = myWSKF.x_pre_all;
o_fil_all = myWSKF.o_fil_all;
o_pre_all = myWSKF.o_pre_all;




############################################################
###       compute percentiles or +/- 1 sigma            ####
############################################################
N_per = 2
Ns_per = 1000
percentiles = [15;85] # 15th and 85th percentiles if independent random variables and diagonal transformation
p07 = 0.7             # the region of highest probability covering a probability of 0.7
percentiles_cond_pre       = Array{Cdouble}(undef,nbin,n_samp,N_per);
percentiles_cond_fil       = Array{Cdouble}(undef,nbin,n_samp,N_per);
percentiles_cond_smo       = Array{Cdouble}(undef,nbin,n_samp,N_per);
percentiles_wall_pre       = Array{Cdouble}(undef,nbin,n_samp,N_per);
percentiles_wall_fil       = Array{Cdouble}(undef,nbin,n_samp,N_per);
percentiles_wall_smo       = Array{Cdouble}(undef,nbin,n_samp,N_per);
percentiles_nuc_pre        = Array{Cdouble}(undef,n_samp,N_per);
percentiles_nuc_fil        = Array{Cdouble}(undef,n_samp,N_per);
percentiles_nuc_smo        = Array{Cdouble}(undef,n_samp,N_per);
for i in 1:n_samp
    # condensation rate
    percentiles_cond_pre[:,i,:] = [CGR(x_pre_all[R_cond,i]-sqrt.(diag(o_pre_all[R_cond,R_cond,i]))); CGR(x_pre_all[R_cond,i]+sqrt.(diag(o_pre_all[R_cond,R_cond,i])))]
    percentiles_cond_fil[:,i,:] = [CGR(x_fil_all[R_cond,i]-sqrt.(diag(o_fil_all[R_cond,R_cond,i]))); CGR(x_fil_all[R_cond,i]+sqrt.(diag(o_fil_all[R_cond,R_cond,i])))]
    percentiles_cond_smo[:,i,:] = [CGR(x_smo_all[R_cond,i]-sqrt.(diag(o_smo_all[R_cond,R_cond,i]))); CGR(x_smo_all[R_cond,i]+sqrt.(diag(o_smo_all[R_cond,R_cond,i])))]
    # wall loss
    for j in 1:nbin
        percentiles_wall_pre[j,i,:] = [wall_rate(x_pre_all[R_loss_init-1+j,i]-sqrt(o_pre_all[R_loss_init-1+j,R_loss_init-1+j,i])); wall_rate(x_pre_all[R_loss_init-1+j,i]+sqrt(o_pre_all[R_loss_init-1+j,R_loss_init-1+j,i]))]
        percentiles_wall_fil[j,i,:] = [wall_rate(x_fil_all[R_loss_init-1+j,i]-sqrt(o_fil_all[R_loss_init-1+j,R_loss_init-1+j,i])); wall_rate(x_fil_all[R_loss_init-1+j,i]+sqrt(o_fil_all[R_loss_init-1+j,R_loss_init-1+j,i]))]
        percentiles_wall_smo[j,i,:] = [wall_rate(x_smo_all[R_loss_init-1+j,i]-sqrt(o_smo_all[R_loss_init-1+j,R_loss_init-1+j,i])); wall_rate(x_smo_all[R_loss_init-1+j,i]+sqrt(o_smo_all[R_loss_init-1+j,R_loss_init-1+j,i]))]
    end
    # nucleation rate
    percentiles_nuc_pre[i,:] = percentile_estimation(Nucleation_rate,x_pre_all[R_nuc_init,i],sqrt(o_pre_all[R_nuc_init,R_nuc_init,i]),percentiles,Ns_per)
    percentiles_nuc_fil[i,:] = percentile_estimation(Nucleation_rate,x_fil_all[R_nuc_init,i],sqrt(o_fil_all[R_nuc_init,R_nuc_init,i]),percentiles,Ns_per)
    percentiles_nuc_smo[i,:] = percentile_estimation(Nucleation_rate,x_smo_all[R_nuc_init,i],sqrt(o_smo_all[R_nuc_init,R_nuc_init,i]),percentiles,Ns_per)
end

############################################################
### simulation of the psd using the estimated parameters ###
############################################################
x_sim_all           = Array{Cdouble}(undef,nbin,n_samp);
x_sim_all[:,1] = x_smo_all[R_psd,1]
for u in 2:n_samp
    x_sim_all[:,u] = iter!(dx_coag,dx_cond,dx_nuc,dx_wall,ws,x_sim_all[:,u-1],t_samp[u]-t_samp[u-1],CGR(x_smo_all[R_cond,u-1]),Nucleation_rate(x_smo_all[R_nuc_init,u-1]),wall_rate(x_smo_all[R_loss,u-1]))
end

############################################################
###           display and save figures                  ####
############################################################
#
try
    include("my_display.jl")
catch msgError
    println(msgError)
    println("unable to print figure for folder: ", folder)
    println(" ")
    println(" ")
end
