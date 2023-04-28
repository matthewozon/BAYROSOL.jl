#TODO: create tests

using Test
using AeroMeas
using Distributions


# define the desired measurement size: the bin centers and the bin widths
cst_r_meas   = 1.0366399037168847;                             # constant ratio
s_meas_0     = 10.9e-9                                         # [m], minimum size
s_meas_inf   = 685.4e-9                                        # [m], maximum size
N_meas       = 116                                             # number of channels
s_meas       = s_meas_0*cst_r_meas.^(collect(0:N_meas-1));     # [m], bin centers
delta_s_meas = s_meas*(sqrt(cst_r_meas)-1.0/sqrt(cst_r_meas)); # [m], bin widths of the measured sizes 

# define the mobilities associated to the bin centers
k_meas = AeroMeas.mobility_from_size_and_charge(s_meas,1;T=293.0,Pr=1.0e5); # [m s^{-1} (V m^{-1})^{-1}], mobilities

# compute the bandwidth using the equation B-6 of the user manual of the SMPS3936 (TSI 3081 DMA + 3775 CPC): dZ = (q_a/q_sh)*Z
q_a          = 0.3;                   # [L/min], aerosol sample flow rate
q_sh         = 3.0;                   # [L/min], sheath flow rate
delta_k_meas = 0.5*(q_a/q_sh)*k_meas; # [m s^{-1} (V m^{-1})^{-1}], 

# a few constant for the SMPS
s50imp=1.0e-6;         # [m], cut off size of the impactor
delta50imp=0.1e-6;     # [m], selectivity of the impactor
phi0 = 1.0 # 1000.0    # [cm^{3} s^{-1}], flow rate at the inlet of the particle counter (CPC). 
delta_t0 = 1.0 # 30.0  # [s], the integration time of the CPC
s50cpc=4.0e-9          # [m], cut off size of the CPC
delta50cpc=2.0e-9      # [m], selectivity of the CPC... more like the spread of the CPC's attenuation at the lower end of the spectrum
T0 = 293.0             # [K], temperature of the carrier gas
Pr0=1.0e5              # [Pa], pressure of the carrier gas

# define the discretization points: for plotting purposes and accuracy
N_mod = 500#0 # 5000 # 1000;                     # number of discretization points
s_0 = 7.0e-9 # 1.0e-9                            # [m], minimum simulated size
s_inf = 1.5e-6 # 1.0e-6                          # [m], maximum simulated size
cst_r = (s_inf/s_0).^(1/(N_mod-1));              # constant ratio of simulated size
s_test = s_0*cst_r.^(collect(0:N_mod-1));        # [m], the sizes
delta_s = s_test*(sqrt(cst_r)-1.0/sqrt(cst_r));  # [m], the width of the simulated bins

# define the size density at the inlet of the SMPS 
u_test = 1.0e13ones(Cdouble,length(s_test));     # [cm^{-3} nm^{-1}], it is a number concentration density (concentration per particle size)

function test_impactor()
  # impactor outlet distribution
  u_imp = AeroMeas.impactor(u_test,s_test;s50=s50imp,delta50=delta50imp)
  # check that the computed density is positive
  all(u_imp.>=0.0)
end

function test_charge_conditioner()
  # impactor outlet distribution
  u_imp = AeroMeas.impactor(u_test,s_test;s50=s50imp,delta50=delta50imp);
  # simulate the charging process
    # negative
  Nq0 = -10 # maximum number of charge per particles
  R_chargep,Kp,u_q_kp = AeroMeas.neutralizer_Kr_85(u_imp,s_test;Nq=Nq0,T=T0,Pr=Pr0);
  cond1 = ((R_chargep==-1:-1:Nq0) & all(Kp.<=0.0) & all(u_q_kp.>=0.0))
    # positive
  Nq0 = 10 # maximum number of charge per particles
  R_chargep,Kp,u_q_kp = AeroMeas.neutralizer_Kr_85(u_imp,s_test;Nq=Nq0,T=T0,Pr=Pr0);
  cond2 = ((R_chargep==1:Nq0) & all(Kp.>=0.0) & all(u_q_kp.>=0.0))

  # test results
  cond1 & cond2
end

function test_DMA()
  # simulate the third module of the SMPS: Differential Mobility Analyzer (DMA)
  Nq0 = -10 # maximum number of charge per particles
  u_dma = AeroMeas.DMA_3010(u_test,s_test,sign(Nq0)*k_meas,delta_k_meas;s50=s50imp,delta50=delta50imp,Nq=Nq0,T=T0,Pr=Pr0);
  # check that the computed density is positive
  all(u_dma.>=0.0)
end

function test_CPC()
  # simulate the third module of the SMPS: Differential Mobility Analyzer (DMA)
  Nq0 = -10 # maximum number of charge per particles
  u_dma = AeroMeas.DMA_3010(u_test,s_test,sign(Nq0)*k_meas,delta_k_meas;s50=s50imp,delta50=delta50imp,Nq=Nq0,T=T0,Pr=Pr0);
  # simulate the density that actually reaches the CPC counting chamber
  u_cpc = AeroMeas.CPC_density(u_dma,s_test;s50=s50cpc,delta50=delta50cpc);
  # check that the computed density is positive
  all(u_cpc.>=0.0)
end

@testset "AeroMeas test bench" begin
  @test test_impactor()
  @test test_charge_conditioner()
  @test test_DMA()
  @test test_CPC()
end
