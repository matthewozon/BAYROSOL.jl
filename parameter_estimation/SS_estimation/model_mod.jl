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
using utilsFun

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
ws.is_coa_gain = COAGULATION_GAIN

if ws.is_coa
    coagulation_coefficient!(ws)
    ws.beta0 = 1.0e6*ws.beta0 # conversion of units: from m^3 s^{-1} to cm^3 s^{-1}
    if ws.is_coa_gain
        init_coagulation_loop_indices!(ws) #WARNING: this step may take time too
    end
end

# parameter rate evolution

include("type_nasha.jl")
T_cond    = 5.0*60.0;
T_nuc     = 5.0*60.0;
sig_mod_n = 0.01*ones(Cdouble,nbin);
sig_mod_c = 1.0*tanh.(0.17*(diameter*1.0e9.+0.8));
sig_mod_l  = 0.001 # 0.05 # 0.0001;
sig_mod_j = 1.0;
gamma_c = 2.0sig_mod_c       # condensation rate initial std
gamma_j = 2.0sig_mod_j
natasha = NASHA(nbin,dt;
    n_cond=nbin,n_loss=nbin,n_nuc=1,
    K_cond=2,K_loss=0,K_nuc=2,
    T_cond=T_cond,T_nuc=T_nuc,
    cor_len_psd=50.0,cor_len_cond=50.0,cor_len_loss=10.0,
    sig_psd=sig_mod_n,sig_cond=sig_mod_c,sig_loss=sig_mod_l,sig_nuc=sig_mod_j)
R_psd = natasha.R_psd
R_cond_all = natasha.R_cond
R_cond = natasha.R_cond[1:2:end]
R_cond_init = natasha.R_cond[1]
R_loss = natasha.R_loss
R_loss_init = natasha.R_loss[1]
R_nuc = natasha.R_nuc
R_nuc_init = natasha.R_nuc[1]


# parameterization of the parameters
#   - condensation
scale_factor_cgr = 2.0./sig_mod_c
function CGR(zeta::Array{Cdouble,1})
    softMaxA(scale_factor_cgr.*zeta,1.0)./scale_factor_cgr
end
function CGR2D(zeta::Array{Cdouble,2})
    Z = similar(zeta)
    for i in 1:size(Z,2)
        Z[:,i] = softMaxA(scale_factor_cgr.*zeta[:,i],1.0)./scale_factor_cgr
    end
    Z
end
function CGRderiv(zeta::Array{Cdouble,1})
    softMaxDerivA(scale_factor_cgr.*zeta,1.0) # ./scale_factor_cgr
end
#   - linear losses
scale_factor_loss = 0.01 # keep this number small enough for stability reasons
function wall_rate(xi::Union{Cdouble,Array{Cdouble,1}})
    softMaxA(xi,scale_factor_loss)
end
function wall_rate_Deriv(xi::Union{Cdouble,Array{Cdouble,1}})
    softMaxDerivA(xi,scale_factor_loss)
end
function wall_rate_inv(xi::Union{Cdouble,Array{Cdouble,1}})
    softMaxInvA(xi,scale_factor_loss)
end
#   - nucleation rate
scale_factor_nuc = 2.0/sig_mod_j # 2.0/sqrt(sig_nuc)
function Nucleation_rate(j::Union{Cdouble,Array{Cdouble,1}})
    # softMaxA(j,scale_factor_nuc)
    # logistic(j,-0.01,2.05,scale_factor_nuc)
    softMaxA(j,scale_factor_nuc)
end
function NucleationDeriv(j::Union{Cdouble,Array{Cdouble,1}})
    # softMaxDerivA(j,scale_factor_nuc)
    # logisticDeriv(j,-0.01,2.05,scale_factor_nuc)
    softMaxDerivA(j,scale_factor_nuc)
end

# temporary buffer for model computation
dx_coag = Array{Cdouble}(undef,nbin);
dx_cond = Array{Cdouble}(undef,nbin);
dx_wall = Array{Cdouble}(undef,nbin);
dx_nuc  = Array{Cdouble}(undef,nbin);

F_coa = Array{Cdouble}(undef,nbin,nbin);
F_con = Array{Cdouble}(undef,nbin,nbin);
F_los = Array{Cdouble}(undef,nbin,nbin);

# the known variables: dx_wall::Array{Cdouble,1} r_los::Cdouble,mu_w::Array{Cdouble,1},sig_w::Array{Cdouble,2},nbin::Int64,model_dim::Int64
function EKF.apply_model!(x_fil_::Array{Cdouble,1},dt_::Cdouble,x_pre_::Array{Cdouble,1})
    # current state: x_fil
    # next state:    x_pre
    x_pre_[R_psd] = iter!(dx_coag,dx_cond,dx_nuc,dx_wall,ws,x_fil_[R_psd],t0*dt_,CGR(x_fil_[R_cond]),Nucleation_rate(x_fil_[R_nuc_init]),wall_rate(x_fil_[R_loss]))
    # growth and loss rate evolution
    x_pre_[R_cond_all] = natasha.B_cond_time*x_fil_[R_cond_all]
    x_pre_[R_loss]     = natasha.B_loss_time*x_fil_[R_loss]
    x_pre_[R_nuc]      = natasha.B_nuc_time*x_fil_[R_nuc]

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
    F_ev_[R_psd,R_psd] = jacobian_GDE(F_coa,F_con,F_los,ws,x_fil_[R_psd],dt_*ws.t0,CGR(x_fil_[R_cond]),wall_rate(x_fil_[R_loss]))

    # condensation rate: derivative w.r.t. the condensation growth rate # the components used in the evolution of the number concentrations are... every second element of the condensation vector
    cgr_deriv = CGRderiv(x_fil_[R_cond]);
    x_factor = ws.scale_GR.*x_fil_[R_psd]
    x_diff = x_factor[1:nbin-1]-x_factor[2:nbin]
    F_ev_[1,R_cond_init] = F_ev_[1,R_cond_init] -  dt_*(t0*GR0/d0)*x_factor[1]*cgr_deriv[1]
    for idx_m in 2:nbin
        F_ev_[idx_m,R_cond_init-1+(2idx_m-1)] = F_ev_[idx_m,R_cond_init-1+(2idx_m-1)] -  dt_*(t0*GR0/d0)*x_factor[idx_m]*cgr_deriv[idx_m]
        F_ev_[idx_m,R_cond_init-1+(2idx_m-1-2)] = F_ev_[idx_m,R_cond_init-1+(2idx_m-1-2)] +  dt_*(t0*GR0/d0)*x_factor[idx_m-1]*cgr_deriv[idx_m-1]
    end

    # wall loss rate: derivative w.r.t. the wall deposition rate
    for idx_m in 1:nbin
        F_ev_[idx_m,R_loss_init-1+idx_m] = F_ev_[idx_m,R_loss_init-1+idx_m] - dt_*gamma0*t0*x_fil_[idx_m]*wall_rate_Deriv(x_fil_[R_loss_init-1+idx_m])
    end

    # nucleation rate: derivative w.r.t. the nucleation rate parameter
    F_ev_[1,R_nuc_init] = F_ev_[1,R_nuc_init] + dt_*(ws.J0*ws.t0/ws.x0)*NucleationDeriv(x_fil_[R_nuc_init])

    # condensation and wall deposition rate evolution
    F_ev_[R_cond_all,R_cond_all] = natasha.B_cond_time
    F_ev_[R_loss,R_loss]         = natasha.B_loss_time
    F_ev_[R_nuc,R_nuc]           = natasha.B_nuc_time
    F_ev_
end

function EKF.update_jacobian!(x_fil_::Array{Cdouble,1},dt_::Cdouble,F_ev_::Array{Cdouble,2})
    # for safety... but this is useless
    F_ev_ = set_jacobian!(x_fil_,dt_,F_ev_)
    F_ev_
end

function EKF.set_measurement_jacobian!(H_me_::Array{Cdouble,2})
    fill!(H_me_,0.0)
    H_me_[1:meas_dim,1:meas_dim] = H_avg # diagm(ones(meas_dim)) #
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
    # poisson counting error + modeling error
    (data/(x0*V_cm3_sample)) .+ (x0_dmps/x0)/(x0*V_cm3_sample) # 0th order approximation error
end

function EKF.update_data_covariance(diagPoisson::Array{Cdouble,1})
    # update the data covariance matrix
    diagm(data_var(diagPoisson))
end
