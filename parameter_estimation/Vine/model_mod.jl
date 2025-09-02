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


using SpecialMatrices
using Polynomials
using StatsBase
using StochProc # for computing the covariance matrix 
using AeroMec2
using utilsFun # for the softMax functions

# create an AeroSys workspace
ws = AeroSys(diameter); # diameter is already defined in the main at this point
# set up the aerosol system for this problem
ws.t0     = t0
ws.x0     = x0
ws.gamma0 = gamma0 # characteristic loss rate
ws.GR0    = GR0    # characteristic condensational growth rate
ws.J0     = J0     # characteristic nucleation rate
ws.is_coa = COAGULATION  # no coagulation
ws.is_con = true   # condensation
ws.is_los = true   # linear losses
ws.is_nuc = true   # nucleation source term
ws.is_coa_gain = COAGULATION_GAIN # true #WARNING: it can be very long in the computation and the initialization of the indices

if ws.is_coa
    coagulation_coefficient!(ws)
    ws.beta0 = 1.0e6*ws.beta0 # conversion of units: from m^3 s^{-1} to cm^3 s^{-1}
    if ws.is_coa_gain
        init_coagulation_loop_indices!(ws) #WARNING: this step may take time too
    end
end

# condensation growth: constant computed many times
scale_GR = Array{Cdouble,1}(undef,nbin);
scale_GR = (cst_r.^(2 .-collect(1:nbin)))/(cst_r-1.0)

# parameter rate evolution


######################################################
#             condensational growth rate             #
######################################################
T_cond      = 3.0*2.0pi*10.0*60.0; # mind the unit # I consider that the characteristic time is the same for each bins... which is kinda true
xij         = 0.95;  # damping coeficient
r1pr2_cond  = 2.0*(1.0-2.0pi*xij*dt/T_cond);
r1mr2_cond  = 1.0-(4.0pi*xij*dt/T_cond)+(4.0*pi^2*(dt/T_cond)^2);
B_cond_time = zeros(Cdouble,2,2)
for ii in 1:1
    B_cond_time[2ii-1,2ii-1] = r1pr2_cond 
    B_cond_time[2ii-1,2ii]   = -r1mr2_cond 
    B_cond_time[2ii,2ii-1]   = 1.0
    B_cond_time[2ii,2ii]     = 0.0
end
sig_mod_c   = 1.0 
sig_cond    = (1.0-(r1pr2_cond-r1mr2_cond)^2)*(sig_mod_c^2)
gamma_cond  = [sig_cond 0.0; 0.0 0.0]

######################################################
#                 wall loss rate                     #
######################################################
#   - deterministic model
r_loss_time        = 1.0 # 0.9999 # 1.0 #0.99 # 0.99     # root of the  process
B_loss_time        = r_loss_time*Matrix{Cdouble}(I,nbin,nbin)
#   - stochastic model
sig_loss           = 0.001^2 #0.01^2
D_loss             = sig_loss*ones(nbin) # the expected variances of every wall loss rate
f_array_loss       = exp.(-collect(0.0:nbin-1.0)/50.0).^2 # 100.0
gamma_loss         = covariance_process_array(f_array_loss,D_loss)


######################################################
#                  nucleation rate                   #
######################################################
# model of the nucleation rate evolution
T_nuc      = 3.0*2.0pi*10.0*60.0; # mind the unit... TODO!!! may be defined in terms of dt
xij        = 0.95;  # damping coeficient
r1pr2_nuc  = 2.0*(1.0-2.0pi*xij*dt/T_nuc);
r1mr2_nuc  = 1.0-(4.0pi*xij*dt/T_nuc)+(4.0*pi^2*(dt/T_nuc)^2);
B_nuc_time = [r1pr2_nuc  -r1mr2_nuc; 1.0 0.0] # 2 roots
sig_mod_j  = 1.0 
sig_nuc    = (1.0-(r1pr2_nuc-r1mr2_nuc)^2)*(sig_mod_j^2)
gamma_nuc  = [sig_nuc 0.0; 0.0 0.0]


######################################################
#                GDE uncertainty                     #
######################################################
sig_mod_n = 2.0*0.005 .+ 0.0*dropdims(0.01*((dt/t0)*(maximum(y_all,dims=2)/x0)/(5.0*60.0/t0)) .+ 1.0e-5,dims=2)
gamma_psd = diagm(sig_mod_n.^2);




# parameter ranges
R_psd  = 1:nbin

R_cond_init = R_psd[end] + 1
R_cond_end = R_cond_init + 2 - 1
R_cond = R_cond_init:R_cond_end

R_loss_init = R_cond[end] + 1  
R_loss_end = R_loss_init + nbin-1
R_loss = R_loss_init:R_loss_end

R_nuc_init = R_loss[end] + 1 
R_nuc_end = R_nuc_init + 1
R_nuc  = R_nuc_init:R_nuc_end


# temporary buffer for model computation
dx_coag = Array{Cdouble,1}(undef,nbin);
dx_cond = Array{Cdouble,1}(undef,nbin);
dx_wall = Array{Cdouble,1}(undef,nbin);
dx_nuc  = Array{Cdouble,1}(undef,nbin);

F_coa = Array{Cdouble,2}(undef,nbin,nbin);
F_con = Array{Cdouble,2}(undef,nbin,nbin);
F_los = Array{Cdouble,2}(undef,nbin,nbin);

# MODEL

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# condensation growth rate
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
scale_factor_cgr = 5.0/sig_mod_c
function CGR(zeta::Union{Cdouble,Array{Cdouble,1}}) #,epsilon::Cdouble=0.0) # zeta is the condensational growth rate parameterization
    softMaxA(zeta,scale_factor_cgr)
end

function CGR2D(zeta::Array{Cdouble,1})
    Z = zeros(Cdouble,nbin,n_samp)
    for i in 1:size(Z,2)
        Z[:,i] = CGR(zeta[i])
    end
    Z
end

function CGRderiv(zeta::Cdouble) # zeta is the condensational growth rate parameterization
    softMaxDerivA(zeta,scale_factor_cgr)
end



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# wall deposition rate
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
scale_factor_loss = 1.0
function wall_rate(xi::Union{Cdouble,Array{Cdouble,1}})
    softMaxA(xi,scale_factor_loss)
end

function wall_rate_Deriv(xi::Union{Cdouble,Array{Cdouble,1}})
    softMaxDerivA(xi,scale_factor_loss)
end

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# wall deposition rate
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
scale_factor_nuc = 5.0/sig_mod_j
function Nucleation_rate(j::Union{Cdouble,Array{Cdouble,1}})
    softMaxA(j,scale_factor_nuc)
end

function NucleationDeriv(j::Union{Cdouble,Array{Cdouble,1}})
    softMaxDerivA(j,scale_factor_nuc)
end




# the known variables: dx_wall::Array{Cdouble,1} r_los::Cdouble,mu_w::Array{Cdouble,1},sig_w::Array{Cdouble,2},nbin::Int64,model_dim::Int64
function EKF.apply_model!(x_fil_::Array{Cdouble,1},dt_::Cdouble,x_pre_::Array{Cdouble,1}) #,sig_w_::Array{Cdouble,2})
    # current state: x_fil
    # next state:    x_pre
    x_pre_[R_psd] = iter!(dx_coag,dx_cond,dx_nuc,dx_wall,ws,x_fil_[R_psd],t0*dt_,CGR(x_fil_[R_cond_init]),Nucleation_rate(x_fil_[R_nuc_init]),wall_rate(x_fil_[R_loss]))
    # growth and loss rate evolution
    x_pre_[R_cond] = B_cond_time*x_fil_[R_cond]
    x_pre_[R_loss] = B_loss_time*x_fil_[R_loss]
    x_pre_[R_nuc]  = B_nuc_time*x_fil_[R_nuc]

    # return
    x_pre_
end


# Jacobian
function EKF.set_jacobian!(x_fil_::Array{Cdouble,1},dt_::Cdouble,F_ev_::Array{Cdouble,2})
    # reset
    fill!(F_ev_,0.0)
    fill!(F_coa,0.0)
    fill!(F_con,0.0)
    fill!(F_los,0.0)
    F_ev_[R_psd,R_psd] = jacobian_GDE(F_coa,F_con,F_los,ws,x_fil_[R_psd],dt_*ws.t0,CGR(x_fil_[R_cond_init]),wall_rate(x_fil_[R_loss]))

    # condensation rate: derivative w.r.t. the condensation growth rate # the components used in the evolution of the number concentrations are... every second element of the condensation vector
    cgr_deriv = CGRderiv(x_fil_[R_cond_init]);
    x_factor = ws.scale_GR.*x_fil_[R_psd]
    x_diff = x_factor[1:nbin-1]-x_factor[2:nbin]
    F_ev_[1,R_cond_init] = F_ev_[1,R_cond_init] -  dt_*(t0*GR0/d0)*x_factor[1]*cgr_deriv
    F_ev_[2:nbin,R_cond_init] =  F_ev_[2:nbin,R_cond_init] + dt_*(t0*GR0/d0)*x_diff*cgr_deriv

    # wall loss rate: derivative w.r.t. the wall deposition rate
    for idx_m in 1:nbin
        F_ev_[idx_m,R_loss_init-1+idx_m] = F_ev_[idx_m,R_loss_init-1+idx_m] - dt_*gamma0*t0*x_fil_[idx_m]*wall_rate_Deriv(x_fil_[R_loss_init-1+idx_m])
    end

    # nucleation rate: derivative w.r.t. the nucleation rate parameter
    F_ev_[1,R_nuc_init] = F_ev_[1,R_nuc_init] + dt_*(ws.J0*ws.t0/ws.x0)*NucleationDeriv(x_fil_[R_nuc_init])

    # condensation and wall deposition rate evolution
    F_ev_[R_cond,R_cond] = B_cond_time
    F_ev_[R_loss,R_loss] = B_loss_time
    F_ev_[R_nuc,R_nuc]   = B_nuc_time
    F_ev_
end

function EKF.update_jacobian!(x_fil_::Array{Cdouble,1},dt_::Cdouble,F_ev_::Array{Cdouble,2})
    # for safety... but this is useless
    F_ev_ = set_jacobian!(x_fil_,dt_,F_ev_)
    F_ev_
end

function EKF.set_measurement_jacobian!(H_me_::Array{Cdouble,2})
    fill!(H_me_,0.0)
    H_me_[1:meas_dim,1:meas_dim] = H_avg
    H_me_
end


function EKF.apply_perfect_measure!(x_pre_::Array{Cdouble,1},H_me_::Array{Cdouble,2})
    H_me_*x_pre_
end


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# error models
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function EKF.update_model_covariance()
    # update the model covariance matrix
end

function data_var(data::Array{Cdouble,1})
    # poisson counting error + modeling error (in terms of density is probably easier to understand)
    (data/(x0*V_cm3_sample)) .+ (x0_dmps/x0)/(x0*V_cm3_sample) # 0th order approximation error
end

function EKF.update_data_covariance(diagPoisson::Array{Cdouble,1})
    # update the data covariance matrix
    diagm(data_var(diagPoisson))
end
