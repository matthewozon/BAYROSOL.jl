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

SAVE_FIG = false
SAVE_DATA = false
FULL_SAVE = false


###############################################################
#   simulation of the evolution of the number concentration   #
###############################################################
# diameter
nbin = 100 # 1731                                             # number of discretization point (size wise)
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
dt = 30.0 # [s]
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
                if (t>180)
                    j = 0.0
                else
                    j = 1.0
                end
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
                if (t>180)
                    time_dep_cond[t] = 1.0
                    global zeta = CGR(time_dep_cond[t])
                    # zeta[50:end] = -zeta[50:end];
                    # zeta[1:49] = -zeta[1:49];
                    zeta[:] = (1.0 .- 2.0 ./(1.0 .+ exp.((0.5.+collect(-(nbin/2.0-1.0):nbin/2.0))/0.5))).*zeta;
                else
                   time_dep_cond[t] = 1.0
                   global zeta = CGR(time_dep_cond[t])
                end
                
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
if nbin>500
    X_display = X[1:40:end,1:32:end] # should be smoothed out before downsampling, but...
    displayLogData2D(9,dt*collect(0:32:n_samp-1)/3600.0,d[1:40:end],ws.x0*X_display,max(0.001ws.x0*maximum(X),ws.x0*minimum(X)),ws.x0*maximum(X),_colorbar_label="concentration [\$\\log \\#\$ cm\$^{-3}\$]")
else
    # X_display = X[:,1:32:end] # should be smoothed out before downsampling, but...
    # displayLogData2D(9,dt*collect(0:32:n_samp-1)/3600.0,d,ws.x0*X_display,max(0.001ws.x0*maximum(X),ws.x0*minimum(X)),ws.x0*maximum(X),_colorbar_label="concentration [\$\\log \\#\$ cm\$^{-3}\$]")
    X_display = X
    displayLogData2D(9,dt*collect(0:n_samp-1)/3600.0,d,ws.x0*X_display,max(0.001ws.x0*maximum(X),ws.x0*minimum(X)),ws.x0*maximum(X),_colorbar_label="concentration [\$\\log \\#\$ cm\$^{-3}\$]")
end
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
if SAVE_FIG
    savefig("simulated_number_concentration.pdf")
    savefig("simulated_number_concentration.png")
end

# plot theparameters
figure(31); loglog(d*1.0e9,xi); title("wall rate: size")
figure(32); loglog(d*1.0e9,GR_t[:,end]); title("condensation rate: size") # 50.0
figure(34); plot(dt*collect(0:n_samp-1)/3600.0,Jt); title("nucleation rate: time")
