# This file is licensed under the MIT "Expat" License:

# Copyright (c) 2022: Matthew Ozon.

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

SAVE_FIG  = false
SAVE_DATA = false


##
## load parameters
##

data_folder = "simulation/"

load_dia = dropdims(Array{Cdouble}(CSV.File(string(data_folder,"mean_d.csv"); header=false) |> DataFrame),dims=1);
load_tim = 3600.0*dropdims(Array{Cdouble}(CSV.File(string(data_folder,"time.csv"); header=false) |> DataFrame),dims=1);
load_X0  = dropdims(Array{Cdouble}(CSV.File(string(data_folder,"aerodist_initial.csv"); header=false) |> DataFrame),dims=1);

# growth rate #TO CHECK
# load_CGR = Matrix{Cdouble}(CSV.File(string(data_folder,"condensation_rate.csv"); header=false) |> DataFrame)
load_cond_length = dropdims(Matrix{Cdouble}(CSV.File(string(data_folder,"cond_area.csv"); header=false) |> DataFrame),dims=1);
load_diff_vapor  =  CSV.File(string(data_folder,"diff_vapor.csv"); header=false)[1];

GR = (1.0/16.0)*4.0*load_diff_vapor.Column1[1]*load_cond_length./(load_dia.^2);

# coagulation coefficient
df_coa = CSV.File(string(data_folder,"coag_kernel.csv"); header=false) |> DataFrame # I have no clue why the last column is interpreted as Strings instead of Float64
for i in 1:length(load_dia)
    if (typeof(df_coa[!,i][1])==String15)
        df_coa[!,i] = parse.(Cdouble,df_coa[!,i])
    end
end
load_coa = Matrix{Cdouble}(df_coa);

# wall loss rate
load_wal  = dropdims(Array{Cdouble}(CSV.File(string(data_folder,"wallrate_beta.csv"); header=false) |> DataFrame),dims=1);


###############################################################
#   simulation of the evolution of the number concentration   #
###############################################################
# diameter
nbin = length(load_dia)                                 # number of discretization point (size wise)
d0 = load_dia[1];                                       # smallest size in simulation
dmax = load_dia[end];                                   # biggest size in simulation
cst_r = (dmax/d0)^(1.0/(nbin-1.0))                      # constant ratio of the log-scale (actually exp-scale, but that's what people are used to)
size_scale = cst_r.^(collect(0:nbin-1))
d = load_dia;                                           # discretization points (centroid of the virtual bins)
delta = (cst_r-1.0)*d0*cst_r.^(collect(0:nbin-1).-0.5)  # width of each virtual bin


########################################################
##         setup the aerosol system here              ##
########################################################
# Aerosol system (what mechanisms must be simulated)
gamma_c = 1.0                          # unit conversion of the wall loss rate
GR_c = 2.777777777777778e-13           # unit conversion of the condensation growth rate 
ws = AeroSys(d, gamma_c = gamma_c, GR_c = GR_c);
ws.is_los = true # false         # linear loss mechanism
ws.is_coa = true # true         # coagulation mechanism# #WARNING: it can be long
ws.is_con = true                 # condensation growth mechanism
ws.is_nuc = false                 # nucleation mechanism
ws.is_coa_gain = false # true    # the default is false when the caogulation mechanism is active, you must explicitly activate the coagulation gain if you want both the loss and gain mechanisms to be computed #WARNING: it can be very long in the computation and the initialization of the indices
if ws.is_coa                     # if the coagulation mechanism is activated, it some more inits are needed (computing the caogulation coefficients and some computation indices to speed up the gain mechanism)
    if true # not computing the coefficient, just plug them in the AeroSys object
        ws.beta[:,:] = 1.0e0load_coa';
        ws.beta0 = 1.0e0 # conversion of units: from m^3 s^{-1} to cm^3 s^{-1}
    else
        coagulation_coefficient!(ws)
        ws.beta0 = 1.0e6*ws.beta0 # conversion of units: from m^3 s^{-1} to cm^3 s^{-1}
    end
    if ws.is_coa_gain
        init_coagulation_loop_indices!(ws) #WARNING: this step may take time too
    end
end


########################################################
##    simulate the evolution of the aerosolsystem     ##
########################################################
n_samp = length(load_tim)
X = zeros(Cdouble,ws.nbin,n_samp)
dx_coa = zeros(Cdouble,ws.nbin)
dx_con = zeros(Cdouble,ws.nbin)
dx_nuc = zeros(Cdouble,ws.nbin)
dx_los = zeros(Cdouble,ws.nbin)
GR_t = zeros(Cdouble,ws.nbin,n_samp)
X[:,1] = load_X0.*log10(cst_r); # ./(1.0e9delta)
load_cond_length

for t in 2:n_samp
    # maximum time step alowed so that the scheme converges without oscillations
    local dt = load_tim[t]-load_tim[t-1];
    local max_dt = 0.5maximum(abs.(1.0 ./(load_wal + 1.0e-9GR./delta)));
    local n_sub_iter = floor(Int64,dt/max_dt)+1; # if dt>dt_max, it is necessary to break the iteration into several sub-iterations with smaller time step
    # time evolution
    X[:,t] = iter!(dx_coa,dx_con,dx_nuc,dx_los,ws,X[:,t-1],dt/n_sub_iter,GR,0.0,load_wal)
    # if dt>max_dt, run the sub iterations
    for _ in 2:n_sub_iter
        X[:,t] = iter!(dx_coa,dx_con,dx_nuc,dx_los,ws,X[:,t],dt/n_sub_iter,GR,0.0,load_wal)
    end
end





##
## plot results and parameters
##

# simulated concentrations
PSDlog = X./log10(cst_r);
displayLogData2D(91,load_tim/3600.0,d,PSDlog,max(0.001maximum(PSDlog),minimum(PSDlog)),maximum(PSDlog),_colorbar_label="log distribution \$\\frac{dN}{d\\log_{10}D_p}\$ [\$\\log \\#\$ cm\$^{-3}\$]")
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
if SAVE_FIG
    savefig(string(data_folder,"simulated_number_concentration.pdf"))
    savefig(string(data_folder,"simulated_number_concentration.png"))
end
if SAVE_DATA
    CSV.write(string(data_folder,"simulated_data_psd_log.csv"),DataFrame(PSDlog',:auto); writeheader=false)
end


# coagualtion coefficient
minCoa,maxCoa = extrema(load_coa)
fig,ax,pcm = imshowData(11,1.0e9d,1.0e9d,load_coa;_norm=:NoNorm,_vmin=0.0,_vmax=maxCoa,_edgecolors="face",_shading="None")
ax.ticklabel_format(axis="both",style="sci",scilimits=(-1,2),useOffset=true)
xlabel("diameter [nm]",fontsize=14)
ylabel("diameter [nm]",fontsize=14)
xscale("log")
yscale("log")
xticks(fontsize=14); 
yticks(fontsize=14);
rc("ytick",color="white")
cax   = fig.add_axes([0.83, .49, 0.03, 0.4])
cbar  = fig.colorbar(pcm, orientation="vertical", cax=cax, shrink=0.6)
cbar.formatter.set_powerlimits((-1,2))
cbar.update_ticks()
cbar.set_label("coagulation coefficient [cm\$^{3}\$ s\$^{-1}\$]",fontsize=14, color="white")
cbar.ax.tick_params(labelsize=14)
cbar.formatter.set_powerlimits((-1,2))
cbar.ax.yaxis.offsetText.set_size(14)
cbar.ax.yaxis.set_tick_params(color="white")
cbar.outline.set_edgecolor("white")
rc("ytick",color="black")
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
if SAVE_FIG
    savefig(string(data_folder,"coagulation_coefficient.pdf"))
    savefig(string(data_folder,"coagulation_coefficient.png"))
end

# growth rate
# cblabelCond = "condensation rate [nm h\$^{-1}\$]"
# minCond,maxCond = extrema(load_CGR)
# fig,ax,pcm = imshowData(10,load_tim[1:end-1]/3600.0,1.0e9d,load_CGR;_norm=:NoNorm,_vmin=minCond,_vmax=maxCond,_edgecolors="face",_shading="None")
# ax.ticklabel_format(axis="both",style="sci",scilimits=(-1,2),useOffset=true)
# xlabel("time [h]",fontsize=14)
# ylabel("diameter [nm]",fontsize=14)
# yscale("log")
# xticks(fontsize=14); 
# yticks(fontsize=14);
# rc("ytick",color="white")
# cax   = fig.add_axes([0.84, .5, 0.03, 0.4])
# cbar  = fig.colorbar(pcm, orientation="vertical", cax=cax, shrink=0.6)
# cbar.formatter.set_powerlimits((-1,2))
# cbar.update_ticks()
# cbar.set_label("growth rate [nm h\$^{-1}\$]",fontsize=14, color="white")
# cbar.ax.tick_params(labelsize=14)
# cbar.formatter.set_powerlimits((-1,2))
# cbar.ax.yaxis.offsetText.set_size(14)
# cbar.ax.yaxis.set_tick_params(color="white")
# cbar.outline.set_edgecolor("white")
# rc("ytick",color="black")
# tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)
# if SAVE_FIG
#     savefig(string(data_folder,"growth_rate.pdf"))
#     savefig(string(data_folder,"growth_rate.png"))
# end