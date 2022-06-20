using PyPlot
using myPlot
using AeroMec2
using AeroMeas
using CSV
using DataFrames
using Statistics



###############################################################
#   simulation of the evolution of the number concentration   #
###############################################################
# diameter
nbin = 3000
d0 = 0.87e-9 #1.0e-9
dmax = 60.0e-9 #  1.0e-6
cst_r = (dmax/d0)^(1.0/(nbin-1.0))
size_scale = cst_r.^(collect(0:nbin-1))
d = d0*size_scale
delta = (cst_r-1.0)*d0*cst_r.^(collect(0:nbin-1).-0.5)  # width of each bin


# Aerosol system (what mechanisms must be simulated)
gamma_c = 1.0e-5                 # characteristic value of the wall loss rate
GR_c = 5.0*2.777777777777778e-13 # characteristic value of the condensational growth rate (not a rate, but that's what they call it)
J_c = 1.0e2 # 50.0               # max value of the nucleation rate
t_c = 300.0                      # characteristic time of change in the distribution
x_c = 1.0e4                      # characteristic number concentration of the distribution
ws = AeroSys(d, gamma_c = gamma_c, GR_c = GR_c, J_c = J_c, t_c = t_c, x_c = x_c);
ws.is_los = true
ws.is_coa = false
ws.is_con = true
ws.is_nuc = true
ws.is_coa_gain = false # true # #WARNING: it can be very long in the computation and the initialization of the indices

if ws.is_coa
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

#   - nucleation rate
sig_mod_j = 0.1
scale_factor_nuc = 2.0/sig_mod_j # 2.0/sqrt(sig_nuc)
function Nucleation_rate(j::Union{Cdouble,Array{Cdouble,1}})
    # softMaxA(j,scale_factor_nuc)
    # logistic(j,-0.01,2.05,scale_factor_nuc)
    softMaxA(j,scale_factor_nuc)
end

#   - linear losses
scale_factor_loss = 0.01
wall_rate = (1.31e-12*ones(nbin)./(ws.d0*size_scale.^1.0))/ws.gamma0


# evolution
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
X[:,1] = 4.5e0*ones(ws.nbin)/ws.x0 # 3.5e2*ones(ws.nbin)/ws.x0
t_init = 0
t_trans = 600
t_end = convert(Int64,round(6.0*3600/dt))

xi = zeros(Cdouble,nbin)
zeta = zeros(Cdouble,nbin)

for t in 2:n_samp
    # nucleation parameter
    if ws.is_nuc
        if ((t>=t_init) & (t<=t_end))
            if ((t<t_trans) & (t>=t_init))
                j = 0.5*(1.0-cos(pi*t/t_trans)) # Nucleation_rate(softMaxInvA(0.5*(1.0-cos(pi*t/t_trans))+1.0e-14,scale_factor_nuc))
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
                zeta = CGR(time_dep_cond[t])
            else
                time_dep_cond[t] = 1.0
                zeta = CGR(time_dep_cond[t])
            end
        else
            time_dep_cond[t] = 0.0
            zeta[:] .= 0.0
        end
    else
        time_dep_cond[t] = 0.0
        zeta[:] .= 0.0
    end
    GR_t[:,t] = ws.GR0*zeta

    # wall loss parameter
    if ws.is_los
        xi = wall_rate
    else
        xi = zeros(nbin)
    end

    # time evolution

    X[:,t] = iter!(dx_coa,dx_con,dx_nuc,dx_los,ws,X[:,t-1],dt,zeta,j,xi)
end


# downsample for display
X_display = X[1:40:end,1:32:end] # should be smoothed out before downsampling, but...

# displayPSD(9,dt*collect(0:32:n_samp-1)/3600.0,d[1:40:end],ws.x0*X_display)
displayLogData2D(9,dt*collect(0:32:n_samp-1)/3600.0,d[1:40:end],ws.x0*X_display,max(10.0,ws.x0*minimum(X)),ws.x0*maximum(X))
savefig("simulated_number_concentration.pdf")
savefig("simulated_number_concentration.png")
displayLogData2D(8,dt*collect(0:32:n_samp-1)/3600.0,d[1:80],ws.x0*X[1:80,1:32:end],max(10.0,ws.x0*minimum(X)),ws.x0*maximum(X))

figure(31); loglog(d*1.0e9,ws.gamma0*wall_rate); title("wall rate: size")
figure(32); loglog(d*1.0e9,GR_t[:,end]); title("condensation rate: size") # 50.0
figure(34); plot(dt*collect(0:n_samp-1)/3600.0,Jt); title("nucleation rate: time") # Nucleation_rate(Jt) # 10.0


# idx_paper = [1; 11; 20; 28; 35; 42; 55; 66; 77] # for nbin = 400
idx_paper = [1; 50; 95; 136; 173; 209; 273; 329; 381] # for nbin = 2000
idx_paper = [1; 71; 135; 194; 248; 298; 390; 472; 545] # for nbin = 3000
idx_paper_str = ["1 nm", "1.1 nm", "1.2 nm", "1.3 nm", "1.4 nm", "1.5 nm", "1.7 nm", "1.9 nm", "2.1 nm"]
figure(100)
for i in idx_paper # 1:11:77 # 1:12
    semilogy(dt*collect(0:n_samp-1)/3600.0,ws.x0*X[i,:])
end
title("bin evolution: number concentration")
legend(idx_paper_str)




###########################################################
#                  convert to density                     #
###########################################################
PSD = ws.x0*X./repeat(delta*1.0e9,1,n_samp);
PSD_display = PSD[1:40:end,1:32:end] # should be smoothed out before downsampling, but...
displayPSD(10,dt*collect(0:32:n_samp-1)/3600.0,d[1:40:end],PSD_display)






#TODO: reconsider the case with 100 bins because the simulation diverges for that size density and time step (use also a hundred bins with a time step of one minute)

###########################################################
#         simulation of the measurement process           #
###########################################################

# measurement
# nbin_meas = 100 #43   # number of measured number concentration
# d0_meas = 1.1e-9
# dmax_meas = 12.0e-9 # 52.0e-9 #  1.0e-6
# cst_r_meas = (dmax_meas/d0_meas)^(1.0/(nbin_meas-1.0))
# d_meas = d0_meas*cst_r_meas.^(collect(0:nbin_meas-1))
# delta_meas = (cst_r_meas-1.0)*d0_meas*cst_r_meas.^(collect(0:nbin_meas-1)-0.5)  # width of each bin
nbin_meas = 50   # number of measured number concentration
d0_meas = 9.0e-10 # 1.1e-9 # 9.2e-10 # 1.1e-9
dmax_meas = 8.5e-9
cst_r_meas = (dmax_meas/d0_meas)^(1.0/(nbin_meas-1.0))
d_meas = d0_meas*cst_r_meas.^(collect(0:nbin_meas-1))
delta_meas = (cst_r_meas-1.0)*d0_meas*cst_r_meas.^(collect(0:nbin_meas-1).-0.5)  # width of each bin



figure(200)
d_test = collect(d0:d0/100:dmax);
for k in 1:nbin_meas
    semilogx(d_test,AeroMeas.chanel_efficiency_gaus(d_meas[k],cst_r_meas,d_test))
end


figure(201)
d_test = collect(d0:d0/100:dmax);
for k in 1:nbin_meas
    scatter(d_meas[k]*1.0e9,sum(AeroMeas.chanel_efficiency_gaus(d_meas[k],cst_r_meas,d)))
end

dt_meass = [6.0; 30.0; 60.0; 300.0]
N_meass = 1
volumes_dmps = [0.25; 1.0; 2.0; 1.0e3; 1.0e9] # [cm^3]
N_vols = 1

SAVE_CSV = false


for i_t in 1:N_meass
    for i_v in 1:N_vols
        close("all")
        # simulated measurement as time and space integral of PSD*chanel_efficiency
        dt_meas = dt_meass[i_t] # 60.0 # 300.0 # 30.0 # 6.0 # 180.0 # 20.0  # time interval between measurements
        n_meas = convert(Int64,floor(n_samp*dt/dt_meas)) # length of the measurement time series
        X_dmps = zeros(Cdouble,nbin_meas,n_meas);
        Jt_meas = zeros(n_meas)
        volume_dmps = volumes_dmps[i_v]

        for i in collect(1:n_meas)
            # time interval indices
            istart = 1+convert(Int64,round((i-1)*dt_meas/dt));
            iend = convert(Int64,round(i*dt_meas/dt));
            Nt = length(collect(istart:iend))
            Dt = dt*ones(Nt)
            Jt_meas[i] = mean(Jt[istart:iend])
            X_dmps[:,i] = AeroMeas.DMPS_gaus(1.0e9PSD[:,istart:iend],d,delta,d_meas,cst_r_meas,volume_dmps)
        end

        figure(300)
        for i in 1:min(10,nbin_meas)
            semilogy(dt_meas*collect(0:n_meas-1)/3600.0,X_dmps[i,:]) #semilogy
        end
        title("bin evolution")

        if ws.is_los
            folder =  string("cloud_",nbin_meas,"_bins/steady_state_CPC_",floor(Int64,volume_dmps),"_cm3/boundary/wall/",floor(Int64,dt_meas),"s/"); # string("here/") ; #
        else
            folder =  string("cloud_",nbin_meas,"_bins/steady_state_CPC_",floor(Int64,volume_dmps),"_cm3/boundary/no_wall/",floor(Int64,dt_meas),"s/"); # string("here/") ; #
        end
        println("creating path: ", folder)
        mkpath(folder);

        displayLogData2D(150,dt_meas*collect(0:n_meas-1)/3600.0,d_meas,X_dmps,max(100.0,minimum(X_dmps)),maximum(X_dmps))
        savefig(string(folder,"measured_number_concentrations.png"))
        savefig(string(folder,"measured_number_concentrations.pdf"))

        displayLogData2D(151,dt_meas*collect(0:n_meas-1)/3600.0,d_meas,1.0e-9X_dmps./delta_meas,max(1000.0,1.0e-9minimum(X_dmps./delta_meas)),1.0e-9maximum(X_dmps./delta_meas))
        savefig(string(folder,"measured_densities.png"))
        savefig(string(folder,"measured_densities.pdf"))


        ###########################################################
        #                    save the data                        #
        ###########################################################
        if SAVE_CSV

            # the output of the DMPS
            CSV.write(string(folder,"time_hour.csv"),DataFrame((dt_meas*collect(0:n_meas-1)/3600.0)'),writeheader=false)
            CSV.write(string(folder,"particle_conc_with_noise.csv"),DataFrame(X_dmps'),writeheader=false)
            CSV.write(string(folder,"diameter.csv"),DataFrame(d_meas'),writeheader=false)

            # the parameters used for the simulation
            CSV.write(string(folder,"diameter_parameter.csv"),DataFrame(d'),writeheader=false)
            CSV.write(string(folder,"time_hour_parameter.csv"),DataFrame((dt*collect(0:n_samp-1)/3600.0)'),writeheader=false)
            CSV.write(string(folder,"wall_rate_constant.csv"),DataFrame(ws.gamma0*wall_rate'),writeheader=false)
            CSV.write(string(folder,"condensation_rate.csv"),DataFrame(GR_t'),writeheader=false)
            CSV.write(string(folder,"nucleation_rate.csv"),DataFrame(Jt'),writeheader=false)
        end
    end
end

