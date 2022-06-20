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
using Interpolations  # duh
using DSP # for the conv function
using Statistics      # some basic stat function
using Distributions   # noise distributions
using Printf          # duh

using CSV             # CSV file interaction package (read and write csv files)
using DataFrames      # framework to organize any data (used for saving simulated data, etc)

SAVE_FIG = true
SAVE_DATA = true
FULL_SAVE = true #WARNING: this may take a while the first time DataFrame is called on a large array (array bigger than some threshold)

###############################################################
#   simulation of the evolution of the number concentration   #
###############################################################
# diameter
nbin = 2500                                                 # number of discretization point (size wise)
d0 = 13.85e-9 #14.1e-9                                      # smallest size in simulation
dmax = 1000.0e-9 #  1.0e-6                                  # biggest size in simulation
cst_r = (dmax/d0)^(1.0/(nbin-1.0))                          # constant ratio of the log-scale (actually exp-scale, but that's what people are used to)
size_scale = cst_r.^(collect(0:nbin-1))  
d = d0*size_scale                                           # discretization points (centroid of the virtual bins)
delta = (cst_r-1.0).*d0.*cst_r.^(collect(0:nbin-1).-0.5)    # width of each virtual bin

# initial number concentrations for the measurement bins (it comes from some data, but ti could be taken as any value, I kept these value because it make the simulation look kinda real)
X0 = [0; 0; 59.0045; 0; 54.8117; 52.82; 50.9679; 98.354; 285.046; 188.869; 792.826; 582.369; 626.335; 1019.61; 1479.65; 2277.28; 2799.62; 2866.33; 3509.41; 4372.08; 5839.68; 5983.4; 7068.27; 7076.28; 7920.85; 8457.28; 9317.7; 10807.6; 10916.9; 10935.2; 10919.1; 11919.5; 13788.8; 14062.6; 15027.2; 15468.6; 15899.3; 16438.1; 16576.3; 16804.7; 15644.5; 16952; 16713.5; 16863.2; 16962; 17533; 15734.6; 16757.6; 16618.3; 15461.9; 16196.9; 14204.3; 12828.1; 13315.6; 12644.8; 12885.9; 12625.8; 12293.4; 10548.6; 10060.7; 10333.1; 9435.94; 8894.24; 8271.11; 7890.6; 8487.54; 7988.06; 6541.18; 7251.68; 5729.07; 6501.54; 5653.86; 5054.88; 5227.32; 4593.64; 4316.54; 4387.3; 3499.14; 3210.47; 2858.49; 2694.92; 2664.87; 2284.26; 2322.8; 1957.48; 1667.56; 1724.54; 1299.65; 1266.1; 1564.71; 968.998; 737.446; 450.526; 614.288; 522.534; 375.052; 521.08; 345.823; 248.171; 145.536; 131.3; 141.238; 0; 28.485; 0; 58.7499; 0; 0; 0; 0; 0]
X0 = conv(X0,[1.; 4.0; 6.0; 4.0; 1.0])[3:end-2]/16.0;   # smoothing of these initial value (mostly just for beauty)
X0 = X0.*(X0.>0.0)./100.0                               # exclude potential negative values
# set the channel center and width for which the initial state had been measured
nbin_meas   = 111                                       # number of channels
d0_meas     = 14.1e-9                                   # smallest bin centroid
dmax_meas   = 0.7365e-6                                 # biggest bin centroid
cst_r_meas  = (dmax_meas/d0_meas)^(1.0/(nbin_meas-1))   # constant ratio of the measured bin
d_meas = Array{Cdouble}(undef,nbin_meas)                # centroid of the measured channels
d_meas[1]=d0_meas
for i in 2:nbin_meas
    d_meas[i]=d_meas[i-1]*cst_r_meas
end
delta_meas = (cst_r_meas-1.0)*d0_meas*cst_r_meas.^(collect(0:nbin_meas-1).-0.5); # width of the measured channels
u0 = X0./delta_meas                                                              # compute the density (well... that's a horrible approximation, but why not?)
# interpolate the initial value to match the simulation's discretization points (cf Interpolations package)
knots = (d_meas,)
u0_itp = interpolate(knots, u0, Gridded(Linear()))
u0_sim = zeros(Cdouble,nbin);
i1 = findfirst(d.>d_meas[1]);
i2 = findfirst(d.>d_meas[end])-1;
u0_sim[i1:i2] = u0_itp(d[i1:i2])
u0_sim = u0_sim.*(u0_sim.>0.0)
X0_sim = u0_sim.*delta # compute the number concentrations of the initial state for the simulation

# plot the initial state (before and after interpolation
figure(1234567)
scatter(1.0e9d_meas,1.0e-9u0)
semilogx(1.0e9d[i1:i2],1.0e-9u0_itp(d[i1:i2]))
semilogx(1.0e9d,1.0e-9X0_sim./delta)
title("initial distribution")
xlabel("diameter [nm]")
ylabel("size density [\$\\#\$ cm\$^{-3}\$ nm\$^{-1}\$]")

########################################################
## change the measured sizes (or not, it's up to you) ##
########################################################
nbin_meas   = 111 # 40 #
d0_meas     = 14.1e-9
dmax_meas   = 0.7365e-6
cst_r_meas  = (dmax_meas/d0_meas)^(1.0/(nbin_meas-1))
d_meas = Array{Cdouble}(undef,nbin_meas)
d_meas[1]=d0_meas
for i in 2:nbin_meas
    d_meas[i]=d_meas[i-1]*cst_r_meas
end
delta_meas = (cst_r_meas-1.0)*d0_meas*cst_r_meas.^(collect(0:nbin_meas-1).-0.5);


########################################################
##         setup the aerosol system here              ##
########################################################
# Aerosol system (what mechanisms must be simulated)
gamma_c = 1.0e-5                  # characteristic value of the wall loss rate
GR_c = 5.0*2.777777777777778e-13  # characteristic value of the condensational growth rate (not a rate, but that's what they call it)
J_c = 0.2                         # max value of the nucleation rate
t_c = 180.0                       # characteristic time of change in the distribution
x_c = 1.0e4                       # characteristic number concentration of the distribution
ws = AeroSys(d, gamma_c = gamma_c, GR_c = GR_c, J_c = J_c, t_c = t_c, x_c = x_c);  # create an AeroSys object (a structure that encompasses all the necessary values for the simulation)
ws.is_los = true  # false         # linear loss mechanism
ws.is_coa = true # false         # coagulation mechanism# #WARNING: it can be long
ws.is_con = true  # false         # condensation growth mechanism
ws.is_nuc = true  # false         # nucleation mechanism
ws.is_coa_gain = true # true     # the default is false when the caogulation mechanism is active, you must explicitly activate the coagulation gain if you want both the loss and gain mechanisms to be computed #WARNING: it can be very long in the computation and the initialization of the indices
if ws.is_coa                      # if the coagulation mechanism is activated, it some more inits are needed (computing the caogulation coefficients and some computation indices to speed up the gain mechanism)
    coagulation_coefficient!(ws)
    ws.beta0 = 1.0e6*ws.beta0 # conversion of units: from m^3 s^{-1} to cm^3 s^{-1}
    if ws.is_coa_gain
        init_coagulation_loop_indices!(ws) #WARNING: this step may take time too
    end
end

# size dependence of the growth rate
GR_array = ones(Cdouble,nbin) #
zeta_dep = GR_array
zeta = zeros(Cdouble,nbin) # just a variable for comupation purposes

# size dependence of the wall losses
wall_rate = WallDepositionRate(ws)/ws.gamma0
xi = wall_rate # just a variable for comupation purposes


########################################################
##    simulate the evolution of the aerosolsystem     ##
########################################################
# evolution
dt = 3.0 # [s]
n_samp = convert(Int64,round(15.0*3600/dt)) 
X = zeros(Cdouble,ws.nbin,n_samp)
dx_coa = zeros(Cdouble,ws.nbin)
dx_con = zeros(Cdouble,ws.nbin)
dx_nuc = zeros(Cdouble,ws.nbin)
dx_los = zeros(Cdouble,ws.nbin)
Jt = zeros(n_samp)
GR_t = zeros(Cdouble,ws.nbin,n_samp)
X[:,1] = X0_sim/ws.x0 # 3.5e2*ones(ws.nbin)/ws.x0
t_init = 0
t_cond_init = convert(Int64,round(5.0*3600/dt))                 # time of the vapour injection
t_cond_end  = convert(Int64,round(10.0*3600/dt))                # time when the injection has stopped
t_end       = convert(Int64,round(15.0*3600/dt))
for t in 2:n_samp
    # nucleation parameter
    if ws.is_nuc
        if ((t>=t_cond_init) & (t<=t_cond_end))
            j = (0.5*(1.0-cos(2.0pi*(t-t_cond_init)/(t_cond_end-t_cond_init))))^4
        else
            j = 0.0
        end
    else
        j = 0.0
    end
    Jt[t] = ws.J0*j


    # condensation parameter
    if ws.is_con
        if ((t>=t_cond_init) & (t<=t_cond_end))
            global zeta = 3.0*(1.0-cos(2.0pi*(t-t_cond_init)/(t_cond_end-t_cond_init)))*GR_array
        else
            fill!(zeta,0.0)
        end
    else
        fill!(zeta,0.0)
    end
    GR_t[:,t] = ws.GR0*zeta

    # time evolution
    X[:,t] = iter!(dx_coa,dx_con,dx_nuc,dx_los,ws,X[:,t-1],dt,zeta,j,xi)
end


# downsample for display
X_display = X[1:40:end,1:32:end]; # should be smoothed out before downsampling, but...
# displayPSD(9,dt*collect(0:32:n_samp-1)/3600.0,d[1:40:end],ws.x0*X_display)
displayLogData2D(9,dt*collect(0:32:n_samp-1)/3600.0,d[1:40:end],ws.x0*X_display,max(0.001ws.x0*maximum(X),ws.x0*minimum(X)),ws.x0*maximum(X),_colorbar_label="concentration [\$\\log \\#\$ cm\$^{-3}\$]")
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)

# plot theparameters
figure(31); loglog(d*1.0e9,ws.gamma0*xi); title("wall rate: size"); tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
figure(32); loglog(d*1.0e9,ws.GR0*zeta_dep); title("condensation rate: size"); tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
figure(34); plot(dt*collect(0:n_samp-1)/3600.0,Jt); title("nucleation rate: time"); tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)




###########################################################
#                  convert to density                     #
###########################################################
PSD = 1.0e-9ws.x0*X./delta; 
PSD_display = PSD[1:40:end,1:32:end]; # should be smoothed out before downsampling, but...
displayLogData2D(10,dt*collect(0:32:n_samp-1)/3600.0,d[1:40:end],PSD_display,max(0.001maximum(PSD_display),minimum(PSD_display)),maximum(PSD_display),_colorbar_label="concentration [\$\\log \\#\$ cm\$^{-3}\$]")
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)



###########################################################
#         simulation of the measurement process           #
###########################################################

# first create the measurement model
q_a = 0.3;             # [L min^{-1}], aerosol sample flow in the DMA
q_sh = 3.0;            # [L min^{-1}], sheath flow

# a few constant for the SMPS
s50imp=1.0e-6;         # [m],  cut off size of the impactor
delta50imp=0.1e-6;     # [m],  selectivity of the impactor
s50cpc=4.0e-9          # [m],  cut off size of the CPC
delta50cpc=2.0e-9      # [m],  selectivity of the CPC... more like the spread of the CPC's attenuation at the lower end of the spectrum
T0 = 293.0             # [K],  temperature of the carrier gas
Pr0=1.0e5              # [Pa], pressure of the carrier gas
Nq0 = -10 # -1 #       # [],   signed number of charges
# compute the transfer function 
H = SMPS3936_transfer_function(d,d_meas,cst_r_meas;s50imp=s50imp,delta50imp=delta50imp,s50cpc=s50cpc, delta50cpc=delta50cpc,T=T0,Pr=Pr0,Nq=Nq0,q_a=q_a,q_sh=q_sh);

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

# get index tables... it's not exactly good because if there is not enough discretization sizes, it bugs... fix later
H_ds = H.*delta'; # middle Riemann integration
idx_min = zeros(Int64,length(d_meas));
idx_max = zeros(Int64,length(d_meas));
# the simulation size discretization start at a bigger size than the smallest size measured
idx_min[1] = 1
idx_max[1] = findfirst(d.>=(d_meas[1]*sqrt(cst_r_meas)))
for i in 2:length(d_meas)
    idx_min[i] = findlast(d.<(d_meas[i]/sqrt(cst_r_meas)))
    idx_max[i] = findfirst(d.>=(d_meas[i]*sqrt(cst_r_meas)))
end
# for each channel, compute the average efficiency
H_avg = zeros(Cdouble,length(d_meas),length(d_meas));
for i in 1:length(d_meas)
    for j in 1:length(d_meas)
        H_avg[i,j] = sum(H_ds[i,idx_min[j]:idx_max[j]])/delta_meas[j]
    end
end

min_H_avg,max_H_avg = extrema(H_avg); # extrema(H.*delta_s'); #
s = @sprintf "SMPS transfer function avg (%1.2f,%1.2f)" min_H_avg max_H_avg
displayLogData2D(3000,1.0e9.*d_meas,1.0e9.*d_meas,H_avg,max(0.001max_H_avg,min_H_avg),max_H_avg,_title=s,_colorbar_label="channel efficiency []")
xscale("log")
xlabel("diameter [nm]")
ylabel("channel center [nm]")
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
if SAVE_FIG
    savefig("kernel_function_contour_plot_avg.png")
    savefig("kernel_function_contour_plot_avg.pdf")
end

figure();
semilogx(1.0e9.*d_meas,H_avg')
xlabel("diameter [nm]")
ylabel("channel efficiency []")
title("SMPS transfer functions: average")
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
if SAVE_FIG
    savefig("kernel_function_avg.png")
    savefig("kernel_function_avg.pdf")
end


# simulated measurement as time and space integral of PSD*chanel_efficiency
dt_meas   = 120.0 # 600.0                           # time interval between measurements
n_meas    = convert(Int64,floor(n_samp*dt/dt_meas)) # length of the measurement time series
X_dmps    = zeros(Cdouble,nbin_meas,n_meas);        # measured concentrations
dma_count = zeros(Cdouble,nbin_meas,n_meas);        # number of particles passing in front of the CPC sensor
phi_a     = 0.05*1000.0/60.0                        # cm3/s: flux of aerosol sample
dt_meas_channel = dt_meas/nbin_meas                 # s: time the CPC counts per channel
volume_count = phi_a*dt_meas_channel                # effective volume of the sample used for counting
# for each measurement time, compute the measurement as if the whole array would be acquired at once
for i in collect(1:n_meas)
    # time interval indices
    istart = 1+convert(Int64,round((i-1)*dt_meas/dt));
    iend = convert(Int64,round(i*dt_meas/dt));
    PSD_time_mean = 1.0e9*dropdims(mean(PSD[:,istart:iend],dims=2),dims=2)  # time average of the particle size density
    dma_count[:,i] = volume_count*H*(PSD_time_mean.*delta)                  # number of particles passing in front of the CPC counting sensor
    # simulate the counting noise
    for j in 1:nbin_meas
        X_dmps[j,i] = rand(Poisson(dma_count[j,i]))/volume_count
    end
    # X_dmps[:,i] = DMPS_gaus(1.0e9PSD[:,istart:iend],d,delta,d_meas,cst_r_meas,volume_dmps) # an older version with a not so good measurement model
end

# plot the evolution of the first 12 measured bins
figure(300)
legend_text = []
for i in 1:12
    semilogy(dt_meas*collect(0:n_meas-1)/3600.0,X_dmps[i,:])
    toto = @sprintf "%1.2f nm" 1.0e9d_meas[i]
    global legend_text = [legend_text; toto]
end
title("bin evolution")
legend(legend_text)

# contour plot: concentration
min_dmps = minimum(X_dmps)#
max_dmps = maximum(X_dmps)
s = @sprintf "number concentration (%1.2e,%1.2e)" min_dmps max_dmps
displayLogData2D(15,dt_meas*collect(0:n_meas-1)/3600.0,d_meas,X_dmps,max(0.001max_dmps,min_dmps),max_dmps;_title=s,_colorbar_label="concentration [\$\\#\$ cm\$^{-3}\$]")
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
if SAVE_FIG
    savefig("data_number_concentration.png")
    savefig("data_number_concentration.pdf")
end

# contour plot: size density
PSD_meas = 1.0e-9X_dmps./delta_meas;
min_dmps = minimum(PSD_meas)#
max_dmps = maximum(PSD_meas)
s = @sprintf "particle size distribution (%1.2e,%1.2e)" min_dmps max_dmps
displayLogData2D(16,dt_meas*collect(0:n_meas-1)/3600.0,d_meas,PSD_meas,max(0.001max_dmps,min_dmps),max_dmps;_title=s,_colorbar_label="density [\$\\#\$ cm\$^{-3}\$ nm\$^{-1}\$]")
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
if SAVE_FIG
    savefig("data_density.png")
    savefig("data_density.pdf")
end

# countour plot: dNdlog10Dp
PSD_meas_log = X_dmps./log10(cst_r_meas);
min_dmps_log = minimum(PSD_meas_log)#
max_dmps_log = maximum(PSD_meas_log)
s = @sprintf "log size distribution (%1.2e,%1.2e)" min_dmps_log max_dmps_log
displayLogData2D(17,dt_meas*collect(0:n_meas-1)/3600.0,d_meas,PSD_meas_log,max(0.001max_dmps_log,min_dmps_log),max_dmps_log;_title=s,_colorbar_label="density [\$\\#\$ cm\$^{-3}\$ log\$_{10}\$(r)\$^{-1}\$]")
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
if SAVE_FIG
    savefig("data_dNdlog10Dp.png")
    savefig("data_dNdlog10Dp.pdf")
end

###########################################################
#                    save the data                        #
###########################################################

if SAVE_DATA
    # the output of the DMPS
    CSV.write("time_hour.csv",DataFrame(dt_meas*collect(0:n_meas-1)'/3600.0,:auto); writeheader=false)
    # writecsv("time_hour.csv",dt_meas*collect(0:n_meas-1)/3600.0)
    CSV.write("particle_conc_with_noise.csv",DataFrame(X_dmps',:auto); writeheader=false)
    # writecsv("particle_conc_with_noise.csv",X_dmps')
    # writecsv("diameter.csv",d_meas)
    CSV.write("diameter.csv",DataFrame(d_meas',:auto); writeheader=false)
    # measurement operator
    CSV.write("H.csv",DataFrame(H,:auto); writeheader=false)
    # average measurement operator
    CSV.write("H_avg.csv",DataFrame(H_avg,:auto); writeheader=false)

    # the parameters used for the simulation
    if FULL_SAVE
        # writecsv("diameter_parameter.csv",d)
        CSV.write("diameter_parameter.csv",DataFrame(d',:auto); writeheader=false)
        # writecsv("time_hour_parameter.csv",dt*collect(0:n_samp-1)/3600.0)
        t_samp = dt*collect(0:n_samp-1)/3600.0;
        CSV.write("time_hour_parameter.csv",DataFrame(t_samp',:auto); writeheader=false)
        # writecsv("wall_rate_constant.csv",ws.gamma0*wall_rate)
        CSV.write("wall_rate_constant.csv",DataFrame(ws.gamma0*wall_rate',:auto); writeheader=false)
        # writecsv("condensation_rate.csv",GR_t')
        CSV.write("condensation_rate.csv",DataFrame(GR_t',:auto); writeheader=false)
        # writecsv("nucleation_rate.csv",Jt)
        CSV.write("nucleation_rate.csv",DataFrame(Jt',:auto); writeheader=false)

        # the simulated PSD
        CSV.write("simulated_psd.csv",DataFrame(PSD',:auto); writeheader=false) 
    else # down sample so that it is possible to save in a limited amount of time
        d_idx_s = 10
        d_idx_t = 10
        CSV.write("diameter_parameter.csv",DataFrame(d[1:d_idx_s:end]',:auto); writeheader=false)
        # writecsv("time_hour_parameter.csv",dt*collect(0:n_samp-1)/3600.0)
        t_samp = dt*collect(0:n_samp-1)/3600.0;
        CSV.write("time_hour_parameter.csv",DataFrame(t_samp[1:d_idx_t:end]',:auto); writeheader=false)
        # writecsv("wall_rate_constant.csv",ws.gamma0*wall_rate)
        CSV.write("wall_rate_constant.csv",DataFrame(ws.gamma0*wall_rate[1:d_idx_s:end]',:auto); writeheader=false)
        # writecsv("condensation_rate.csv",GR_t')
        CSV.write("condensation_rate.csv",DataFrame(GR_t[1:d_idx_s:end,1:d_idx_t:end]',:auto); writeheader=false)
        # writecsv("nucleation_rate.csv",Jt)
        CSV.write("nucleation_rate.csv",DataFrame(Jt[1:d_idx_t:end]',:auto); writeheader=false)

        # the simulated PSD
        CSV.write("simulated_psd.csv",DataFrame(PSD[1:d_idx_s:end,1:d_idx_t:end]',:auto); writeheader=false)
    end
end
