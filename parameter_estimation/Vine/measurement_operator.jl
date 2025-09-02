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
nbin = 10000 # 2500                                                 # number of discretization point (size wise)
# d0 = 13.85e-9 #14.1e-9                                      # smallest size in simulation
d0 = 12.88e-9 #14.1e-9                                      # smallest size in simulation
dmax = 1000.0e-9 #  1.0e-6                                  # biggest size in simulation
cst_r = (dmax/d0)^(1.0/(nbin-1.0))                          # constant ratio of the log-scale (actually exp-scale, but that's what people are used to)
size_scale = cst_r.^(collect(0:nbin-1))  
d = d0*size_scale                                           # discretization points (centroid of the virtual bins)
delta = (cst_r-1.0).*d0.*cst_r.^(collect(0:nbin-1).-0.5)    # width of each virtual bin

# set the channel center and width for which the initial state had been measured
# nbin_meas   = 111                                       # number of channels
nbin_meas = 230
d0_meas     = 13.0e-9 # 14.1e-9                                   # smallest bin centroid
dmax_meas   = 0.79863e-6                                 # biggest bin centroid
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
# if SAVE_FIG
#     savefig("kernel_function_contour_plot.png")
#     savefig("kernel_function_contour_plot.pdf")
# end

figure();
semilogx(1.0e9.*d,H')
xlabel("diameter [nm]")
ylabel("channel efficiency []")
title("SMPS transfer functions")
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
# if SAVE_FIG
#     savefig("kernel_function.png")
#     savefig("kernel_function.pdf")
# end

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
# s = @sprintf "SMPS transfer function avg (%1.2f,%1.2f)" min_H_avg max_H_avg
s = @sprintf "transfer function avg (%1.2f,%1.2f)" min_H_avg max_H_avg
displayLogData2D(3000,1.0e9.*d_meas,1.0e9.*d_meas,H_avg,max(0.001max_H_avg,min_H_avg),max_H_avg,_title=s,_colorbar_label="channel efficiency []")
xscale("log")
xlabel("diameter [nm]")
ylabel("channel center [nm]")
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
# if SAVE_FIG
#     savefig("kernel_function_contour_plot_avg.png")
#     savefig("kernel_function_contour_plot_avg.pdf")
# end

figure();
semilogx(1.0e9.*d_meas,H_avg')
xlabel("diameter [nm]")
ylabel("channel efficiency []")
title("SMPS transfer functions: average")
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
# if SAVE_FIG
#     savefig("kernel_function_avg.png")
#     savefig("kernel_function_avg.pdf")
# end


# simulated measurement as time and space integral of PSD*chanel_efficiency
dt = 3.0 # [s]
n_samp = convert(Int64,round(15.0*3600/dt)) 
dt_meas   = 120.0 # 600.0                           # time interval between measurements
n_meas    = convert(Int64,floor(n_samp*dt/dt_meas)) # length of the measurement time series
phi_a     = 0.05*1000.0/60.0                        # cm3/s: flux of aerosol sample
dt_meas_channel = dt_meas/nbin_meas                 # s: time the CPC counts per channel
volume_count = phi_a*dt_meas_channel                # effective volume of the sample used for counting
# number of particles passing in front of the CPC counting sensor
# dma_count[:,i] = volume_count*H*(PSD_time_mean.*delta)                  

# save measurement operator
CSV.write("diameter.csv",DataFrame(d_meas',:auto); writeheader=false)
# measurement operator
CSV.write("diameter_discretization.csv",DataFrame(d',:auto); writeheader=false)
CSV.write("H.csv",DataFrame(H,:auto); writeheader=false)
# average measurement operator
CSV.write("diameter_discretization_avg.csv",DataFrame(d_meas',:auto); writeheader=false)
CSV.write("H_avg.csv",DataFrame(H_avg,:auto); writeheader=false)

