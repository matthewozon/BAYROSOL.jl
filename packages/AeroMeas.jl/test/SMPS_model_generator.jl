using AeroMeas      # load the measurement model (you need to add the path of the source code to the juliarc or equivalent (some file where the Julia packages are stored))
using PyPlot        # plots
# set the font in the plots
fm = PyPlot.matplotlib.font_manager.json_load("/home/mattoz/.cache/matplotlib/fontlist-v310.json")
fm.findfont("serif", rebuild_if_missing=false)
fm.findfont("serif", fontext="afm", rebuild_if_missing=false)
rc("font",family="serif",serif="Computer Modern Roman")
rc("text", usetex=true)
using Distributions # for random number generator such as Poisson
using Printf
using myPlot

# define the desired measurement size: the bin centers and the bin widths
# N_meas = 20;                                                   # number of bins
# s_meas_0 = 7.0e-9                                              # [m], minimum size
# s_meas_inf = 0.5e-6                                            # [m], maximum size
# cst_r_meas = (s_meas_inf/s_meas_0).^(1/(N_meas-1));            # constant ratio
cst_r_meas = 1.0366399037168847
s_meas_0 = 10.9e-9
s_meas_inf = 685.4e-9
N_meas = 116 # 10 #
# cst_r_meas = (s_meas_inf/s_meas_0).^(1/(N_meas-1));
s_meas = s_meas_0*cst_r_meas.^(collect(0:N_meas-1));           # [m], bin centers
delta_s_meas = s_meas*(sqrt(cst_r_meas)-1.0/sqrt(cst_r_meas)); # [m], bin widths of the measured sizes # TODO: compute the size bandwidths based on the mobility bandwidths
# sig_s_meas   = delta_s_meas/sqrt(2.0pi);                       # [m], standard deviation of the channel efficiency associated to the bin widths

# define the mobilities associated to the bin centers
k_meas = AeroMeas.mobility_from_size_and_charge(s_meas,1;T=293.0,Pr=1.0e5);               # [m s^{-1} (V m^{-1})^{-1}], mobilities
# compute the actual bandwidth using the equation B-6 of the user manual of the SMPS3936 (TSI 3081 DMA + 3775 CPC): dZ = (q_a/q_sh)*Z
# delta_k_meas = abs.(AeroMeas.dk_from_ds_and_charge(sig_s_meas,s_meas,1;T=293.0,Pr=1.0e5)); # mobility width associated to the size bin width (I set it positive because we know the density is positive)
q_a = 0.3;
q_sh = 3.0; # 3.0; #
delta_k_meas = 0.5*(q_a/q_sh)*k_meas;

# a few constant for the SMPS
s50imp=1.0e-6;         # [m], cut off size of the impactor
delta50imp=0.1e-6;     # [m], selectivity of the impactor
phi0 = 1.0 # 1000.0    # [cm^{3} s^{-1}], flow rate at the inlet of the particle counter (CPC). It is the sample flow (q_a), not the sheath flow (q_{sh})
delta_t0 = 1.0 # 30.0  # [s], the integration time of the CPC
s50cpc=4.0e-9          # [m], cut off size of the CPC
delta50cpc=2.0e-9      # [m], selectivity of the CPC... more like the spread of the CPC's attenuation at the lower end of the spectrum
T0 = 293.0             # [K], temperature of the carrier gas
Pr0=1.0e5              # [Pa], pressure of the carrier gas

# define the discretization points: for plotting purposes and accuracy
# N_mod = 5000; # 1000;                              # number of discretization points
# s_0 = 1.0e-9; # 1.0e-9                           # [m], minimum simulated size
# s_inf = 1.0e-6; # 1.0e-6                         # [m], maximum simulated size
N_mod = 5000 # 5000 # 1000;
s_0 = 7.0e-9 # 1.0e-9
s_inf = 1.5e-6 # 1.0e-6
# or match the measured channel centers
# s_0 = 10.9e-9
# s_inf = 685.4e-9
# N_mod = 116
cst_r = (s_inf/s_0).^(1/(N_mod-1));              # constant ratio of simulated size
s_test = s_0*cst_r.^(collect(0:N_mod-1));        # [m], the sizes
# s_test[:] = s_meas[:] # TODO: get rid of that horrible thing!
delta_s = s_test*(sqrt(cst_r)-1.0/sqrt(cst_r));  # [m], the width of the simulated bins

# define the size density at the inlet of the SMPS (hereafter referred to as the ground truth)
idx_test = ((s_test.<148.0e-9) .| (s_test.>210.0e-9))
# idx_test = ((s_test.<115.0e-9) .| (s_test.>184.0e-9))
# idx_test = ((s_test.<100.0e-9) .| (s_test.>180.0e-9))
# idx_test = ((s_test.<40.0e-9) .| (s_test.>100.0e-9))
# idx_test = ((s_test.<94.0e-9) .| (s_test.>160.0e-9))
# idx_test = []
u_test = 1.0e13ones(Cdouble,length(s_test)); # [# cm^{-3} nm^{-1}], it is a number concentration density (concentration per particle size)
u_test[idx_test] .= 0.0;
u_meas = 1.0e13ones(Cdouble,length(s_meas)); # [# cm^{-3} nm^{-1}], what would go to the counting section of the SMPS if it had a 100% efficiency

# simulate the impactor: first module of the SMPS: remove particles larger than a given size
u_imp = AeroMeas.impactor(u_test,s_test;s50=s50imp,delta50=delta50imp);

# simulate the charging process: only need one polarity (negative is better in this case because the bias is toward having more negatively charged particles) and only the Nq0 first charges per particles because in the size range considered, it does not make sens to consider more than Nq0~10 elementary charges per particles (according to Wiedensohler's technical not J. Aero. Sci. V19, N3, p387, particles smaller than 20 nm only carry at most 1 charge and particles smaller than 70 nm carry at most 2 charges (and particles smaller than 5 nm are mostl neutral, less than 5% of charged particles)). Note that for the size range [1,1000] nm, Nq0=10 and more does not give any difference at the outlet of the SMPS.
Nq0 = -10 # it seems like a reasonable approximation that I seem to remember from a few papers (and also it works in practice)
R_chargep,Kp,u_q_kp = AeroMeas.neutralizer_Kr_85(u_imp,s_test;Nq=Nq0,T=T0,Pr=Pr0);
figure(); semilogx(abs.(Kp'),u_q_kp'); legend(["-1","-2","-3","-4","-5","-6","-7","-8","-9","-10"]); title("mobility density"); xlabel("electrical mobility [m\$^{2}\$ V\$^{-1}\$ s\$^{-1}\$]"); ylabel("density [cm\$^{-3}\$ (m\$^{2}\$ V\$^{-1}\$ s\$^{-1}\$)\$^{-1}\$]")

# simulate the third module of the SMPS: Differential Mobility Analyzer (DMA)
u_dma = AeroMeas.DMA_3010(u_test,s_test,sign(Nq0)*k_meas,delta_k_meas;s50=s50imp,delta50=delta50imp,Nq=Nq0,T=T0,Pr=Pr0);
figure(); semilogx(s_test,u_test); semilogx(s_test,u_imp); semilogx(s_test,u_dma'); title("size density"); ylim(0.0); legend(["sample at the inlet","after impactor","after DMA"]); xlabel("diameter [m]"); ylabel("density [cm\$^{-3}\$ m\$^{-1}\$]")

# simulate the density that actually reach the CPC counting chamber
u_cpc = AeroMeas.CPC_density(u_dma,s_test;s50=s50cpc,delta50=delta50cpc);

# simulate the counting process: count for a give amount of time
Ncount = similar(k_meas);
[Ncount[i] = AeroMeas.nb_count(u_cpc[i,:],delta_s,phi0,delta_t0) for i in 1:size(u_cpc,1)] # number of particle counted in each bin

# simulate the errors happening while counting: Poisson noise
Y = similar(k_meas);
[Y[i] = (1.0/(phi0*delta_t0))*rand(Poisson(Ncount[i])) for i in 1:length(Ncount)]

# plot
figure(); semilogx(s_test,u_test); semilogx(s_test,u_imp); semilogx(s_test,u_cpc'); title("size density"); ylim(0.0); legend(["sample at the inlet","after impactor","after CPC"]); xlabel("diameter [m]"); ylabel("density [cm\$^{-3}\$ m\$^{-1}\$]")
figure(); semilogx(s_meas,(1.0/(phi0*delta_t0))*Ncount./delta_s_meas); semilogx(s_meas,Y./delta_s_meas); semilogx(s_test,u_test); title("number concentration density"); xlabel("size [m]"); ylabel("density [cm\$^{-3}\$ m\$^{-1}\$]"); ylim(0.0); legend(["counted concentration density","counted concentration density with noise","input density"]); ylim(0.0)
figure(); semilogx(s_meas,(1.0/(phi0*delta_t0))*Ncount); semilogx(s_meas,Y); semilogx(s_meas,u_meas.*delta_s_meas); title("number concentration"); xlabel("size [m]"); ylabel("density [cm\$^{-3}\$]"); legend(["counted concentration","counted concentration with noise","ground truth"]); ylim(0.0)



# do the same as before, but everything at once (it is possible to add the data inversion from Hoppel 1978
inv_data = false; # select wether or not you want to inverse the data using Hoppel 1978
Nqinv = -3;       # if data inversion is selected, one needs to specify the maximum number of charges per particle to be considered (it should not exceed 3 in absolute value otherwise it is converging for a very small size range)
Yw,u_dmaw,u_cpcw = AeroMeas.SMPS3936(u_test,s_test,delta_s,s_meas,cst_r_meas,phi0,delta_t0;s50imp=s50imp,delta50imp=delta50imp,s50cpc=s50cpc, delta50cpc=delta50cpc,T=T0,Pr=Pr0,Nq=Nq0,inv_data=false,Nq_inv=-6,q_a=q_a,q_sh=q_sh)
H = AeroMeas.SMPS3936_transfer_function(s_test,s_meas,cst_r_meas;s50imp=s50imp,delta50imp=delta50imp,s50cpc=s50cpc, delta50cpc=delta50cpc,T=T0,Pr=Pr0,Nq=Nq0,q_a=q_a,q_sh=q_sh)
# plot
figure(); semilogx(s_test,u_test); semilogx(s_test,u_imp); semilogx(s_test,u_dmaw'); title("size density (SMPS 3936)"); xlabel("diameter [m]"); ylabel("density [cm\$^{-3}\$ m\$^{-1}\$]"); ylim(0.0); legend(["sample at the inlet","after impactor","after DMA"]);
figure(); semilogx(s_test,u_test); semilogx(s_test,u_imp); semilogx(s_test,u_cpcw'); title("size density (SMPS 3936)"); xlabel("diameter [m]"); ylabel("density [cm\$^{-3}\$ m\$^{-1}\$]"); ylim(0.0); legend(["sample at the inlet","after impactor","after CPC"]);
figure(); semilogx(s_meas,Yw); semilogx(s_meas,u_meas.*delta_s_meas); title("number concentration (SMPS 3936)"); xlabel("diameter [m]"); ylabel("number concentration [cm\$^{-3}\$]"); legend(["measured by SMPS 3936","desired measurements"])


min_H,max_H = extrema(H); # extrema(H.*delta_s'); #
s = @sprintf "SMPS transfer function (%1.2f,%1.2f)" min_H max_H
displayLogData2D(2000,1.0e9.*s_test,1.0e9.*s_meas,H,max(0.001max_H,min_H),max_H,_title=s,_colorbar_label="channel efficiency []")
# displayLogData2D(2000,1.0e9.*s_test,1.0e9.*s_meas,H.*delta_s',max(0.01max_H,min_H),max_H,_title=s,_colorbar_label="channel efficiency []")
xscale("log")
xlabel("diameter [nm]")
ylabel("channel center [nm]")


# get index tables... it's not exactly good because if there is not enough discretization sizes, it bugs... fix later
H_ds = H.*delta_s'; # middle Riemann integration
idx_min = zeros(Int64,length(s_meas));
idx_max = zeros(Int64,length(s_meas));
for i in 1:length(s_meas)
    idx_min[i] = findlast(s_test.<(s_meas[i]/sqrt(cst_r_meas)))
    idx_max[i] = findfirst(s_test.>=(s_meas[i]*sqrt(cst_r_meas)))
end
# for each channel, compute the average efficiency
H_avg = zeros(Cdouble,length(s_meas),length(s_meas));
for i in 1:length(s_meas)
    for j in 1:length(s_meas)
        H_avg[i,j] = sum(H_ds[i,idx_min[j]:idx_max[j]])/delta_s_meas[j]
    end
end

min_H_avg,max_H_avg = extrema(H_avg); # extrema(H.*delta_s'); #
s = @sprintf "SMPS transfer function avg (%1.2f,%1.2f)" min_H_avg max_H_avg
displayLogData2D(2001,1.0e9.*s_meas,1.0e9.*s_meas,H_avg,max(0.001max_H_avg,min_H_avg),max_H_avg,_title=s,_colorbar_label="channel efficiency []")
# displayLogData2D(2000,1.0e9.*s_test,1.0e9.*s_meas,H.*delta_s',max(0.01max_H,min_H),max_H,_title=s,_colorbar_label="channel efficiency []")
xscale("log")
xlabel("diameter [nm]")
ylabel("channel center [nm]")



H_meas = AeroMeas.SMPS3936_transfer_function(s_meas,s_meas,cst_r_meas;s50imp=s50imp,delta50imp=delta50imp,s50cpc=s50cpc, delta50cpc=delta50cpc,T=T0,Pr=Pr0,Nq=Nq0,q_a=q_a,q_sh=q_sh)


min_H_meas,max_H_meas = extrema(H_meas);
s = @sprintf "SMPS transfer function meas (%1.2f,%1.2f)" min_H_meas max_H_meas
displayLogData2D(2002,1.0e9.*s_meas,1.0e9.*s_meas,H_meas,max(0.001max_H_meas,min_H_meas),max_H_meas,_title=s,_colorbar_label="channel efficiency []")
xscale("log")
xlabel("diameter [nm]")
ylabel("channel center [nm]")
