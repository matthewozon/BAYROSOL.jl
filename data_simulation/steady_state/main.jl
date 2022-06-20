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


# This piece of code simulates the evolution of an aerosol system and the measurement process. The resulting simulation and parameters may be saved (cf. end of the code)

using PyPlot          # plotting package
# set matplotlib to allow fior the use of LaTeX formula in the graphs (side effect: the font is set to LaTeX's default)
rc("text", usetex=true)
using myPlot          # my plotting (makes it easy to plot size distributions)
using AeroMec2        # aerosol system simulation
using AeroMeas        # measurement simulation
using Statistics      # some basic stat function
using Distributions   # noise distributions
using Printf          # duh



using CSV             # CSV file interaction package (read and write csv files)
using DataFrames      # framework to organize any data (used for saving simulated data, etc)

SAVE_FIG = true
SAVE_DATA = true
FULL_SAVE = true

if false
###############################################################
#   simulation of the evolution of the number concentration   #
###############################################################
# diameter
nbin = 1731                                             # number of discretization point (size wise)
d0 = 0.87e-9                                            # smallest size in simulation
dmax = 1.0003467874623244e-8;                           # biggest size in simulation
cst_r = (dmax/d0)^(1.0/(nbin-1.0))                      # constant ratio of the log-scale (actually exp-scale, but that's what people are used to)
size_scale = cst_r.^(collect(0:nbin-1))
d = d0*size_scale                                       # discretization points (centroid of the virtual bins)
delta = (cst_r-1.0)*d0*cst_r.^(collect(0:nbin-1).-0.5)  # width of each virtual bin


########################################################
##         setup the aerosol system here              ##
########################################################
# Aerosol system (what mechanisms must be simulated)
gamma_c = 1.0e-5                 # characteristic value of the wall loss rate
GR_c = 5.0*2.777777777777778e-13 # characteristic value of the condensational growth rate (not a rate, but that's what they call it)
J_c = 5.0                        # max value of the nucleation rate # realistic values are in the rage 0.01 to 0.1, the others are rookie numbers
t_c = 300.0                      # characteristic time of change in the distribution
x_c = 1.0e4                      # characteristic number concentration of the distribution
ws = AeroSys(d, gamma_c = gamma_c, GR_c = GR_c, J_c = J_c, t_c = t_c, x_c = x_c);
ws.is_los = true # false         # linear loss mechanism
ws.is_coa = true # true         # coagulation mechanism# #WARNING: it can be long
ws.is_con = true                 # condensation growth mechanism
ws.is_nuc = true                 # nucleation mechanism
ws.is_coa_gain = true # true    # the default is false when the caogulation mechanism is active, you must explicitly activate the coagulation gain if you want both the loss and gain mechanisms to be computed #WARNING: it can be very long in the computation and the initialization of the indices
if ws.is_coa                     # if the coagulation mechanism is activated, it some more inits are needed (computing the caogulation coefficients and some computation indices to speed up the gain mechanism)
    coagulation_coefficient!(ws)
    ws.beta0 = 1.0e6*ws.beta0 # conversion of units: from m^3 s^{-1} to cm^3 s^{-1}
    if ws.is_coa_gain
        init_coagulation_loop_indices!(ws) #WARNING: this step may take time too
    end
end



# define the parameters
# parameterization of the parameters
#   - condensation
GR_array = tanh.(0.17*(ws.d*1.0e9.+0.8))
function CGR(time_dep::Cdouble)
    time_dep*GR_array
end


#   - linear losses
scale_factor_loss = 0.01
wall_rate = (1.31e-12*ones(nbin)./(ws.d0*size_scale.^1.0))/ws.gamma0


########################################################
##    simulate the evolution of the aerosolsystem     ##
########################################################
dt = 3.0 # [s]
n_samp = convert(Int64,round(2.25*3600/dt)) # convert(Int64,round(6.0*3600/dt)) # 1000
X = zeros(Cdouble,ws.nbin,n_samp)
dx_coa = zeros(Cdouble,ws.nbin)
dx_con = zeros(Cdouble,ws.nbin)
dx_nuc = zeros(Cdouble,ws.nbin)
dx_los = zeros(Cdouble,ws.nbin)
Jt = zeros(n_samp)
GR_t = zeros(Cdouble,ws.nbin,n_samp)
time_dep_cond = zeros(Cdouble,n_samp)
X[:,1] = 4.5e0*ones(ws.nbin)/ws.x0 
X[:,1] = 1.0e9.*delta/ws.x0
t_init = 0
t_trans = convert(Int64,round(.5*3600/dt)) # 600
t_end = convert(Int64,round(6.0*3600/dt))

xi = zeros(Cdouble,nbin)
zeta = zeros(Cdouble,nbin)

for t in 2:n_samp
    # nucleation parameter
    if ws.is_nuc
        if ((t>=t_init) & (t<=t_end))
            if ((t<t_trans) & (t>=t_init))
                j = 0.5*(1.0-cos(pi*t/t_trans))
            else
                j = 1.0
            end
        else
            j = 0.0
        end
    else
        j = 0.0
    end
    Jt[t] = ws.J0*j


    # condensation parameter
    if ws.is_con
        if ((t>=t_init) & (t<=t_end))
            if ((t<t_trans) & (t>=t_init))
                time_dep_cond[t] = 0.5*(1.0-cos(pi*t/t_trans))
                global zeta = CGR(time_dep_cond[t])
            else
                time_dep_cond[t] = 1.0
                global zeta = CGR(time_dep_cond[t])
            end
        else
            time_dep_cond[t] = 0.0
            global zeta[:] = 0.0
        end
    else
        time_dep_cond[t] = 0.0
        global zeta[:] = 0.0
    end
    GR_t[:,t] = ws.GR0*zeta

    # wall loss parameter
    if ws.is_los
        global xi = wall_rate
    else
        global xi = zeros(nbin)
    end

    # time evolution
    X[:,t] = iter!(dx_coa,dx_con,dx_nuc,dx_los,ws,X[:,t-1],dt,zeta,j,xi)
end


# downsample for display
X_display = X[1:40:end,1:32:end] # should be smoothed out before downsampling, but...
displayLogData2D(9,dt*collect(0:32:n_samp-1)/3600.0,d[1:40:end],ws.x0*X_display,max(0.001ws.x0*maximum(X),ws.x0*minimum(X)),ws.x0*maximum(X),_colorbar_label="concentration [\$\\log \\#\$ cm\$^{-3}\$]")
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
if SAVE_FIG
    savefig("simulated_number_concentration.pdf")
    savefig("simulated_number_concentration.png")
end

# plot theparameters
figure(31); loglog(d*1.0e9,xi); title("wall rate: size")
figure(32); loglog(d*1.0e9,GR_t[:,end]); title("condensation rate: size") # 50.0
figure(34); plot(dt*collect(0:n_samp-1)/3600.0,Jt); title("nucleation rate: time")



###########################################################
#                  convert to density                     #
###########################################################
PSD = 1.0e-9ws.x0*X./delta;
PSD_display = PSD[1:40:end,1:32:end] # should be smoothed out before downsampling, but...
# displayPSD(10,dt*collect(0:32:n_samp-1)/3600.0,d[1:40:end],PSD_display)
# displayPSD(11,dt*collect(0:16:n_samp-1)/3600.0,d[1:4:end],PSD[1:4:end,1:16:end])
s = @sprintf "size distribution (%1.2e,%1.2e)" minimum(1.0e-9ws.x0*X./delta) maximum(1.0e-9ws.x0*X./delta)
displayLogData2D(10,dt*collect(0:32:n_samp-1)/3600.0,d[1:40:end],PSD_display,max(0.001maximum(PSD_display),minimum(PSD_display)),maximum(PSD_display),_title=s,_colorbar_label="concentration [\$\\log \\#\$ cm\$^{-3}\$]")
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
displayLogData2D(11,dt*collect(0:16:n_samp-1)/3600.0,d[1:4:end],PSD[1:4:end,1:16:end],max(0.001maximum(PSD_display),minimum(PSD_display)),maximum(PSD_display),_title=s,_colorbar_label="concentration [\$\\log \\#\$ cm\$^{-3}\$]")
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)

PSD_log = ws.x0*X./log10(cst_r);
PSD_log_display = PSD_log[1:4:end,1:16:end]; # should be smoothed out before downsampling, but...
min_psd_log,max_psd_log = extrema(PSD_log);
s = @sprintf "dma simulation (%1.2f,%1.2f)" min_psd_log max_psd_log
# displayLogData2D(21,dt*collect(0:16:n_samp-1)/3600.0,d[1:4:end],PSD_log_display,max(0.001max_psd_log,min_psd_log),max_psd_log,_title=s,_colorbar_label="log density [# cm\$^{-3}\$ log\$_{10}\$(r)]")
displayLogData2D(21,dt*collect(0:16:n_samp-1)/3600.0,d[1:4:end],PSD_log_display,max(0.001max_psd_log,min_psd_log),max_psd_log,_title=s,_colorbar_label="log density [\$\\#\$ cm\$^{-3}\\log_{10}\$(r)]")
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)


########################################################
##                measured sizes                      ##
########################################################
nbin_meas = 50 # 21 # 12 #  # number of measured number concentration
d0_meas = 1.1e-9 # 9.0e-10
dmax_meas = 8.5e-9
cst_r_meas = (dmax_meas/d0_meas)^(1.0/(nbin_meas-1.0))
d_meas = d0_meas*cst_r_meas.^(collect(0:nbin_meas-1));
delta_meas = (cst_r_meas-1.0)*d0_meas*cst_r_meas.^(collect(0:nbin_meas-1).-0.5);  # width of each bin

#WARNING: an SMPS3936 is not supposed to operate at such small size, so I artificially improve it

########################################################
##           create measurement model                 ##
########################################################
q_a = 0.3;             # [L min^{-1}], aerosol sample flow in the DMA
q_sh = 5.0q_a          # [L min^{-1}], sheath flow 

# a few constant for the SMPS
s50imp=1.0e-6;         # [m],  cut off size of the impactor
delta50imp=0.1e-6;     # [m],  selectivity of the impactor
s50cpc=0.8e-9          # [m],  cut off size of the CPC
delta50cpc=0.1e-9      # [m],  selectivity of the CPC... more like the spread of the CPC's attenuation at the lower end of the spectrum
T0 = 293.0             # [K],  temperature of the carrier gas
Pr0=1.0e5              # [Pa], pressure of the carrier gas
Nq0 = -10 # -1 #       # [],   signed number of charges
#WARNING: here I use a multiplication factor because the SMPS that I modeled is not supposed to be used at so small sizes and so the charging probabilities are too small
H = 5.0*SMPS3936_transfer_function(d,d_meas,cst_r_meas;s50imp=s50imp,delta50imp=delta50imp,s50cpc=s50cpc, delta50cpc=delta50cpc,T=T0,Pr=Pr0,Nq=Nq0,q_a=q_a,q_sh=q_sh);

# get a beter resolution for computing the average measurement operator
nbin_more = 1000000; 
d0_more = 0.7e-9; 
dmax_more = 1.1e-8; 
cst_r_more = (dmax_more/d0_more)^(1.0/(nbin_more-1.0));
size_scale_more = cst_r_more.^(collect(0:nbin_more-1));
d_more = d0_more*size_scale_more;
delta_more = (cst_r_more-1.0)*d0_more*cst_r_more.^(collect(0:nbin_more-1).-0.5) ; # width of each bin
H_more = 5.0*SMPS3936_transfer_function(d_more,d_meas,cst_r_meas;s50imp=s50imp,delta50imp=delta50imp,s50cpc=s50cpc, delta50cpc=delta50cpc,T=T0,Pr=Pr0,Nq=Nq0,q_a=q_a,q_sh=q_sh);




# plot the transfer function
min_H,max_H = extrema(H); 
s = @sprintf "SMPS transfer function (%1.2f,%1.2f)" min_H max_H
displayLogData2D(2002,1.0e9.*d,1.0e9.*d_meas,H,max(0.001max_H,min_H),max_H,_title=s,_colorbar_label="channel efficiency []")
xscale("log")
xlabel("diameter [nm]")
ylabel("channel center [nm]")
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
if SAVE_FIG
    savefig("kernel_function_contour_plot.png")
    savefig("kernel_function_contour_plot.pdf")
end

figure();
semilogx(1.0e9.*d,H')
xlabel("diameter [nm]")
ylabel("channel efficiency []")
title("SMPS transfer functions")
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
if SAVE_FIG
    savefig("kernel_function.png")
    savefig("kernel_function.pdf")
end

# compute an average operator
# get index tables... it's not exactly good because if there is not enough discretization sizes, it bugs... fix later
H_ds_more = H_more.*delta_more';            # middle Riemann integration
idx_min_more = zeros(Int64,length(d_meas));
idx_max_more = zeros(Int64,length(d_meas));
# the simulation size discretization start at a bigger size than the smallest size measured
# idx_min[1] = 1
# idx_max[1] = findfirst(d.>=(d_meas[1]*sqrt(cst_r_meas)))
for i in 1:length(d_meas) #WARNING: maybe the first min index will not exists, uncomment the above to fix it and start the loop at the second index
    idx_min_more[i] = findlast(d_more.<(d_meas[i]/sqrt(cst_r_meas)))
    idx_max_more[i] = findfirst(d_more.>=(d_meas[i]*sqrt(cst_r_meas)))
end
# for each channel, compute the average efficiency
H_avg_more = zeros(Cdouble,length(d_meas),length(d_meas));
for i in 1:length(d_meas)
    for j in 1:length(d_meas)
        H_avg_more[i,j] = sum(H_ds_more[i,idx_min_more[j]:idx_max_more[j]])/delta_meas[j]
    end
end

min_H_avg_more,max_H_avg_more = extrema(H_avg_more); # extrema(H.*delta_s'); #
s = @sprintf "SMPS transfer function avg more (%1.2f,%1.2f)" min_H_avg_more max_H_avg_more
displayLogData2D(3005,1.0e9.*d_meas,1.0e9.*d_meas,H_avg_more,max(0.001max_H_avg_more,min_H_avg_more),max_H_avg_more,_title=s,_colorbar_label="channel efficiency []")
# displayLogData2D(2000,1.0e9.*s_test,1.0e9.*s_meas,H.*delta_s',max(0.01max_H,min_H),max_H,_title=s,_colorbar_label="channel efficiency []")
xscale("log")
xlabel("diameter [nm]")
ylabel("channel center [nm]")
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
if SAVE_FIG
    savefig("kernel_function_contour_plot_avg_more.png")
    savefig("kernel_function_contour_plot_avg_more.pdf")
end

figure();
semilogx(1.0e9.*d_meas,H_avg_more')
xlabel("diameter [nm]")
ylabel("channel efficiency []")
title("SMPS transfer functions: average")
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
if SAVE_FIG
    savefig("kernel_function_avg_more.png")
    savefig("kernel_function_avg_more.pdf")
end
end
###########################################################
#         simulation of the measurement process           #
###########################################################

phi_a = (0.05*1000.0/60.0)*[1.0; 10.0; 100.0; 1000.0] # sample flow rate passing in front of the optical detector
dt_meass = [120.0]
volumes_dmps = (phi_a/nbin_meas)*dt_meass[1] # volume used for counting in each channel

for i_t in 1:length(dt_meass)
    for i_v in 1:length(volumes_dmps)
        close("all")

        dt_meas   = dt_meass[i_t]
        n_meas    = convert(Int64,floor(n_samp*dt/dt_meas)) # length of the measurement time series
        X_dmps    = zeros(Cdouble,nbin_meas,n_meas);
        dma_count = zeros(Cdouble,nbin_meas,n_meas); # number of particles passing in front of the CPC sensor
        Jt_meas   = zeros(n_meas)
        volume_count = volumes_dmps[i_v] # phi_a*dt_meas_channel

        for i in collect(1:n_meas)
            # time interval indices
            istart = 1+convert(Int64,round((i-1)*dt_meas/dt));
            iend = convert(Int64,round(i*dt_meas/dt));
            PSD_time_mean = 1.0e9*dropdims(mean(PSD[:,istart:iend],dims=2),dims=2)
            dma_count[:,i] = volume_count*H*(PSD_time_mean.*delta)
            for j in 1:nbin_meas
                X_dmps[j,i] = rand(Poisson(dma_count[j,i]))/volume_count
            end
            Nt = length(collect(istart:iend))
            Dt = dt*ones(Nt)
            Jt_meas[i] = mean(Jt[istart:iend])
        end


        if ws.is_los
            folder =  string("cloud_",nbin_meas,"_bins_J_5/steady_state_CPC_",floor(Int64,volume_count),"_cm3/boundary/wall/",floor(Int64,dt_meas),"s/"); # string("here/") ; #
        else
            folder =  string("cloud_",nbin_meas,"_bins_J_5/steady_state_CPC_",floor(Int64,volume_count),"_cm3/boundary/no_wall/",floor(Int64,dt_meas),"s/"); # string("here/") ; #
        end
        println("creating path: ", folder)
        mkpath(folder);

        displayLogData2D(150,dt_meas*collect(0:n_meas-1)/3600.0,d_meas,X_dmps,max(10.0,minimum(X_dmps)),maximum(X_dmps),_colorbar_label="concentration [\$\\log \\#\$ cm\$^{-3}\$]")
        tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
        if SAVE_FIG
            savefig(string(folder,"measured_number_concentrations.png"))
            savefig(string(folder,"measured_number_concentrations.pdf"))
        end


        min_dmps = minimum(X_dmps)#
        max_dmps = maximum(X_dmps)
        local s = @sprintf "number concentration (%1.2e,%1.2e)" min_dmps max_dmps
        displayLogData2D(15,dt_meas*collect(0:n_meas-1)/3600.0,d_meas,X_dmps,max(0.001max_dmps,min_dmps),max_dmps;_title=s,_colorbar_label="concentration [\$\\#\$ cm\$^{-3}\$]")
        tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
        if SAVE_FIG
            savefig(string(folder,"data_number_concentration.png"))
            savefig(string(folder,"data_number_concentration.pdf"))
        end


        # delta_meas = (cst_r_meas-1.0)*d0_meas*cst_r_meas.^(collect(0:nbin_meas-1).-0.5);
        PSD_meas = 1.0e-9X_dmps./delta_meas;
        min_dmps = minimum(PSD_meas)#
        max_dmps = maximum(PSD_meas)
        s = @sprintf "particle size distribution (%1.2e,%1.2e)" min_dmps max_dmps
        displayLogData2D(16,dt_meas*collect(0:n_meas-1)/3600.0,d_meas,PSD_meas,max(0.001max_dmps,min_dmps),max_dmps;_title=s,_colorbar_label="density [\$\\#\$ cm\$^{-3}\$ nm\$^{-1}\$]")
        tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
        if SAVE_FIG
            savefig(string(folder,"data_density.png"))
            savefig(string(folder,"data_density.pdf"))
        end


        PSD_meas_log = X_dmps./log10(cst_r_meas);
        min_dmps_log = minimum(PSD_meas_log)#
        max_dmps_log = maximum(PSD_meas_log)
        s = @sprintf "log size distribution (%1.2e,%1.2e)" min_dmps_log max_dmps_log
        displayLogData2D(17,dt_meas*collect(0:n_meas-1)/3600.0,d_meas,PSD_meas_log,max(0.001max_dmps_log,min_dmps_log),max_dmps_log;_title=s,_colorbar_label="density [\$\\#\$ cm\$^{-3}\$ log\$_{10}\$(r)\$^{-1}\$]")
        tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
        if SAVE_FIG
            savefig(string(folder,"data_dNdlog10Dp.png"))
            savefig(string(folder,"data_dNdlog10Dp.pdf"))
        end

        ###########################################################
        #                    save the data                        #
        ###########################################################
        if SAVE_DATA
            # the output of the DMPS
            CSV.write(string(folder,"time_hour.csv"),DataFrame(dt_meas*collect(0:n_meas-1)'/3600.0,:auto); writeheader=false)
            # writecsv("time_hour.csv",dt_meas*collect(0:n_meas-1)/3600.0)
            CSV.write(string(folder,"particle_conc_with_noise.csv"),DataFrame(X_dmps',:auto); writeheader=false)
            # writecsv("particle_conc_with_noise.csv",X_dmps')
            # writecsv("diameter.csv",d_meas)
            CSV.write(string(folder,"diameter.csv"),DataFrame(d_meas',:auto); writeheader=false)
            # measurement operator
            CSV.write(string(folder,"H.csv"),DataFrame(H,:auto); writeheader=false)
            # average measurement operator
            CSV.write(string(folder,"H_avg_more.csv"),DataFrame(H_avg_more,:auto); writeheader=false) # H_avg_more is a better approximation of the averaged measurement operator

            # the parameters used for the simulation
            if FULL_SAVE
                # writecsv("diameter_parameter.csv",d)
                CSV.write(string(folder,"diameter_parameter.csv"),DataFrame(d',:auto); writeheader=false)
                # writecsv("time_hour_parameter.csv",dt*collect(0:n_samp-1)/3600.0)
                t_samp = dt*collect(0:n_samp-1)/3600.0;
                CSV.write(string(folder,"time_hour_parameter.csv"),DataFrame(t_samp',:auto); writeheader=false)
                # writecsv("wall_rate_constant.csv",ws.gamma0*wall_rate)
                CSV.write(string(folder,"wall_rate_constant.csv"),DataFrame(ws.gamma0*wall_rate',:auto); writeheader=false)
                # writecsv("condensation_rate.csv",GR_t')
                CSV.write(string(folder,"condensation_rate.csv"),DataFrame(GR_t',:auto); writeheader=false)
                # writecsv("nucleation_rate.csv",Jt)
                CSV.write(string(folder,"nucleation_rate.csv"),DataFrame(Jt',:auto); writeheader=false)

                # the simulated PSD
                CSV.write(string(folder,"simulated_psd.csv"),DataFrame(PSD',:auto); writeheader=false)
            else # down sample so that it is possible to save in a limited amount of time
                d_idx_s = 10
                d_idx_t = 10
                CSV.write(string(folder,"diameter_parameter.csv"),DataFrame(d[1:d_idx_s:end]',:auto); writeheader=false)
                # writecsv("time_hour_parameter.csv",dt*collect(0:n_samp-1)/3600.0)
                t_samp = dt*collect(0:n_samp-1)/3600.0;
                CSV.write(string(folder,"time_hour_parameter.csv"),DataFrame(t_samp[1:d_idx_t:end]',:auto); writeheader=false)
                # writecsv("wall_rate_constant.csv",ws.gamma0*wall_rate)
                CSV.write(string(folder,"wall_rate_constant.csv"),DataFrame(ws.gamma0*wall_rate[1:d_idx_s:end]',:auto); writeheader=false)
                # writecsv("condensation_rate.csv",GR_t')
                CSV.write(string(folder,"condensation_rate.csv"),DataFrame(GR_t[1:d_idx_s:end,1:d_idx_t:end]',:auto); writeheader=false)
                # writecsv("nucleation_rate.csv",Jt)
                CSV.write(string(folder,"nucleation_rate.csv"),DataFrame(Jt[1:d_idx_t:end]',:auto); writeheader=false)

                # the simulated PSD
                CSV.write(string(folder,"simulated_psd.csv"),DataFrame(PSD[1:d_idx_s:end,1:d_idx_t:end]',:auto); writeheader=false)
            end
        end
    end
end
