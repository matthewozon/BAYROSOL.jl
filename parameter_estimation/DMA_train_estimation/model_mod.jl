# using SpecialMatrices
# using Polynomials
# using StatsBase
# using StochProc # for computing the covariance matrix
# using AeroMec2
# using utilsFun # for the softMax functions

# create an AeroSys workspace
ws = AeroMec2.AeroSys(diameter); # diameter is already defined in the main at this point
# set up the aerosol system for this problem
ws.t0     = t0
ws.x0     = x0
ws.gamma0 = gamma0 # characteristic loss rate
ws.GR0    = GR0    # characteristic condensational growth rate
ws.J0     = J0     # characteristic nucleation rate
ws.is_coa = false  # no coagulation
ws.is_con = true   # condensation
ws.is_los = true   # linear losses
ws.is_nuc = true   # nucleation source term

# condensation growth: constant computed many times
scale_GR = Array{Cdouble,1}(undef,nbin);
scale_GR = (cst_r.^(2 .-collect(1:nbin)))/(cst_r-1.0)

# parameter rate evolution


######################################################
#             condensational growth rate             #
######################################################
T_cond      = 0.5*2.0pi*10.0*60.0;  # time constant
if FLAG_1906_03
    T_cond      = 0.1*2.0pi*10.0*60.0;
end
xij         = 0.95;  # damping coeficient
r1pr2_cond  = 2.0*(1.0-2.0pi*xij*dt/T_cond);
r1mr2_cond  = 1.0-(4.0pi*xij*dt/T_cond)+(4.0*pi^2*(dt/T_cond)^2);
B_cond_time = zeros(Cdouble,2nbin,2nbin)
for ii in 1:nbin
    B_cond_time[2ii-1,2ii-1] = r1pr2_cond
    B_cond_time[2ii-1,2ii]   = -r1mr2_cond
    B_cond_time[2ii,2ii-1]   = 1.0
    B_cond_time[2ii,2ii]     = 0.0
end
sig_mod_c                       = ones(Cdouble,nbin) # 1.0 .+ 0.0*2.0*tanh.(0.17*(diameter*1.0e9.+0.8));
if FLAG_1906_03
    sig_mod_c                       = 0.25*ones(Cdouble,nbin)
end
gamma_c                         = 2.0sig_mod_c
## weak correlation
# cor_len                         = nbin/2.0; # short
# cor_len                         = nbin/1.0; # long
# f_array_cond                    = exp.(-collect(0.0:nbin-1.0)/cor_len).^2
## strong correlation
cor_len                         = nbin/1.0
## short
# f_array_cond                    = exp.(-(collect(0.0:nbin-1.0)/(0.333cor_len)).^2).^2
## long
f_array_cond                    = exp.(-(collect(0.0:nbin-1.0)/(0.4cor_len)).^2).^2

gamma_proc                      = StochProc.covariance_process_array(f_array_cond,sig_mod_c.^2)
gamma_proc_XXL                  = zeros(Cdouble,2*nbin,2*nbin);
gamma_proc_XXL[1:2:end,1:2:end] = gamma_proc
gamma_proc_XXL[1:2:end,2:2:end] = gamma_proc
gamma_proc_XXL[2:2:end,1:2:end] = gamma_proc
gamma_proc_XXL[2:2:end,2:2:end] = gamma_proc
gamma_cond                      = gamma_proc_XXL - B_cond_time*gamma_proc_XXL*B_cond_time'
# force the uncorrelated variables to 0
gamma_cond[1:2:end,2:2:end]    .= 0.0
gamma_cond[2:2:end,2:2:end]    .= 0.0
gamma_cond[2:2:end,1:2:end]    .= 0.0



######################################################
#                 wall loss rate                     #
######################################################
#   - deterministic model
r_loss_time        = 1.0  # root of the  process
B_loss_time        = r_loss_time*Matrix{Cdouble}(I,nbin,nbin)
#   - stochastic model
sig_loss           = 0.001^2
D_loss             = sig_loss*ones(nbin) # the expected variances of every wall loss rate
correlarion_len    = nbin/10.0
f_array_loss       = exp.(-collect(0.0:nbin-1.0)/correlarion_len).^2 # 100.0
gamma_loss         = StochProc.covariance_process_array(f_array_loss,D_loss)


# model of the nucleation rate evolution
T_nuc      = 0.5*2.0pi*10.0*60.0; # time constant
if FLAG_1906_03
    T_nuc      = 0.1*2.0pi*10.0*60.0;
end
xij        = 0.95;  # damping coeficient
r1pr2_nuc  = 2.0*(1.0-2.0pi*xij*dt/T_nuc);
r1mr2_nuc  = 1.0-(4.0pi*xij*dt/T_nuc)+(4.0*pi^2*(dt/T_nuc)^2);
B_nuc_time = [r1pr2_nuc  -r1mr2_nuc; 1.0 0.0] # 2 roots
sig_mod_j  = 1.0
sig_nuc    = (1.0-(r1pr2_nuc-r1mr2_nuc)^2)*(sig_mod_j^2)
gamma_nuc  = [sig_nuc 0.0; 0.0 0.0]

# initial GDE covariance matrix (is update later on at each iteration)
n0                = 5.0*ones(Cdouble,length(diameter_model));
sig_mod_n         = n0.*(1.0e9*delta);
D_psd             = sig_mod_n.^2 # the expected variances of the number concentrations
correlarion_psd   = nbin/50.0
f_array_psd       = exp.(-collect(0.0:nbin-1.0)/correlarion_psd).^2
gamma_psd         = StochProc.covariance_process_array(f_array_psd,D_psd)

# parameter ranges
R_psd  = 1:nbin

R_cond_init = R_psd[end] + 1
R_cond_end  = R_cond_init + 2nbin - 1
R_cond_all  = R_cond_init:R_cond_end
R_cond      = R_cond_init:2:R_cond_end

R_loss_init = R_cond_end + 1
R_loss_end = R_loss_init + nbin -1
R_loss = R_loss_init:R_loss_end

R_nuc_init = R_loss_end + 1
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

scale_factor_cgr = 5.0./sig_mod_c
function CGR(zeta::Array{Cdouble,1})
    utilsFun.softMaxA(scale_factor_cgr.*zeta,1.0)./scale_factor_cgr
end
function CGR2D(zeta::Array{Cdouble,2})
    Z = similar(zeta)
    for i in 1:size(Z,2)
        Z[:,i] = utilsFun.softMaxA(scale_factor_cgr.*zeta[:,i],1.0)./scale_factor_cgr
    end
    Z
end
function CGRderiv(zeta::Array{Cdouble,1})
    utilsFun.softMaxDerivA(scale_factor_cgr.*zeta,1.0) 
end



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# wall deposition rate
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
scale_factor_loss = 1.0
function wall_rate(xi::Union{Cdouble,Array{Cdouble,1}})
    utilsFun.softMaxA(xi,scale_factor_loss)
end

function wall_rate_Deriv(xi::Union{Cdouble,Array{Cdouble,1}})
    utilsFun.softMaxDerivA(xi,scale_factor_loss)
end
function wall_rate_inv(xi::Union{Cdouble,Array{Cdouble,1}})
    utilsFun.softMaxInvA(xi,scale_factor_loss)
end

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# wall deposition rate
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
scale_factor_nuc = 5.0/sig_mod_j
function Nucleation_rate(j::Union{Cdouble,Array{Cdouble,1}})
    utilsFun.softMaxA(j,scale_factor_nuc)
end

function NucleationDeriv(j::Union{Cdouble,Array{Cdouble,1}})
    utilsFun.softMaxDerivA(j,scale_factor_nuc)
end


# the known variables: dx_wall::Array{Cdouble,1} r_los::Cdouble,mu_w::Array{Cdouble,1},sig_w::Array{Cdouble,2},nbin::Int64,model_dim::Int64
function EKF.apply_model!(x_fil_::Array{Cdouble,1},dt_::Cdouble,x_pre_::Array{Cdouble,1}) #,sig_w_::Array{Cdouble,2})
    # current state: x_fil
    # next state:    x_pre
    x_pre_[R_psd] = AeroMec2.iter!(dx_coag,dx_cond,dx_nuc,dx_wall,ws,x_fil_[R_psd],t0*dt_,CGR(x_fil_[R_cond]),Nucleation_rate(x_fil_[R_nuc_init]),wall_rate(x_fil_[R_loss]))
    # growth and loss rate evolution
    x_pre_[R_cond_all] .= B_cond_time*x_fil_[R_cond_all]
    x_pre_[R_loss]     = B_loss_time*x_fil_[R_loss]
    x_pre_[R_nuc]      .= B_nuc_time*x_fil_[R_nuc]

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
    F_ev_[R_psd,R_psd] =  AeroMec2.jacobian_GDE(F_coa,F_con,F_los,ws,x_fil_[R_psd],dt_*ws.t0,CGR(x_fil_[R_cond]),wall_rate(x_fil_[R_loss]))

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
    F_ev_[R_cond_all,R_cond_all] = B_cond_time
    F_ev_[R_loss,R_loss]         = B_loss_time
    F_ev_[R_nuc,R_nuc]           = B_nuc_time
    F_ev_
end

function EKF.update_jacobian!(x_fil_::Array{Cdouble,1},dt_::Cdouble,F_ev_::Array{Cdouble,2})
    # for safety... but this is useless
    F_ev_ = EKF.set_jacobian!(x_fil_,dt_,F_ev_)
    F_ev_
end

function EKF.set_measurement_jacobian!(H_me_::Array{Cdouble,2})
    fill!(H_me_,0.0)
    if need_padding
        if (mod(myWSKF.ikf,pad_factor+1)==1)
            H_me_[1:meas_dim,1:nbin] = H_DMATRAIN
        end
    else
        H_me_[1:meas_dim,1:nbin] = H_DMATRAIN
    end
    H_me_
end

function EKF.update_measurement_jacobian!(H_me_::Array{Float64,2})
    EKF.set_measurement_jacobian!(H_me_)
end

function EKF.apply_perfect_measure!(x_pre_::Array{Cdouble,1},H_me_::Array{Cdouble,2})
    H_me_*x_pre_
end


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# error models
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
err_model     = zeros(Cdouble,nb_ch);
meas_disc_err = zeros(Cdouble,nb_ch);
meas_mode_err = zeros(Cdouble,nb_ch);
e_gr = 0.5 # this can be evaluated from the KF estimate, and it work pretty much the same
reg_background = 0.001
reg_backcground_count = 1.0
n_est   = zeros(Cdouble,nbin);
n_deriv = zeros(Cdouble,nbin);
e_g_kf  = zeros(Cdouble,nbin);
e_l_kf  = zeros(Cdouble,nbin);
part_flux_err = zeros(Cdouble,nbin);
loss_err_kf   = zeros(Cdouble,nbin);
discretization_error = zeros(Cdouble,nbin);

# some storage
ERR_DISC = zeros(Cdouble,n_samp_data,nbin);
ERR_GROW = zeros(Cdouble,n_samp_data,nbin);
ERR_LOSS = zeros(Cdouble,n_samp_data,nbin);


function EKF.update_model_covariance()
    # update the model covariance matrix
    if (mod(myWSKF.ikf,pad_factor+1)==1)
        # time index in the non-padded space
        idx_idx = floor(Int64,myWSKF.ikf/(pad_factor+1))+1 # min(n_samp_data,)

        # get some estimate of the uncertainty in the growth rate and nucleation (in the right units)
        e_g_kf[:] = (1.0e9GR0*CGR(sqrt.(diag(myWSKF.o_fil[R_cond,R_cond]))))                      # [nm s^{-1}] growth rate error
        e_j_kf    = J0*Nucleation_rate(sqrt(myWSKF.o_fil[R_nuc_init,R_nuc_init]))                  # [cm^{-3} s^{-1}] nucleation rate error
        part_flux_err[:] = 0.5*AeroMec2.cond_err(dt, e_g_kf,e_j_kf,x0*abs.(myWSKF.x_fil[R_psd]),1.0e9delta) # [cm^{-3}] error due to errors in the nucleation and growth rate

        # get some estimation of the uncertainty in the loss rate
        e_l_kf[:]      = gamma0.*wall_rate(sqrt.(diag(myWSKF.o_fil[R_loss,R_loss])))
        loss_err_kf[:] = 0.5*AeroMec2.loss_err(dt, e_l_kf,x0*myWSKF.x_fil[R_psd],1.0e9delta)

        # compute an estimate of the density and its first derivative
        n_est[:]   = x0*myWSKF.x_fil[R_psd]./(1.0e9delta)
        n_deriv[:] = abs.([(n_est[1:end-1]-n_est[2:end])./delta[1:end-1]; (n_est[end-1]-n_est[end])./delta[end]])
        n_deriv[:] = conv(n_deriv,[1.; 4.0; 6.0; 4.0; 1.0])[3:end-2]/16.0;
        n_deriv[:] = n_deriv.*(n_deriv.>0.0)

        # discretization errors
        discretization_error[:] = AeroMec2.discretization_err(dt,n_deriv,GR0*CGR(myWSKF.x_fil[R_cond]),1.0e9delta)

        # store each contribution of the error
        ERR_GROW[idx_idx,:] = part_flux_err
        ERR_LOSS[idx_idx,:] = loss_err_kf
        ERR_DISC[idx_idx,:] = discretization_error

        # put all the contribution together and compute the covariance matrix
        myWSKF.gam_w[R_psd,R_psd] = StochProc.covariance_process_array(f_array_psd,(part_flux_err.^2) .+ (loss_err_kf.^2) .+ (discretization_error.^2) .+ (reg_background)^2)
    end
    myWSKF.gam_w
end

function data_var(data::Array{Cdouble,1})
    if (mod(myWSKF.ikf,pad_factor+1)==1)
        # time index in the non-padded space
        idx_idx = floor(Int64,myWSKF.ikf/(pad_factor+1))+1 # min(n_samp_data,)

        # compute an estimate of the density and its first derivative
        n_est[:]   = 1.0e-9x0*myWSKF.x_fil[R_psd]./delta
        n_deriv[:] = abs.([(n_est[1:end-1]-n_est[2:end])./delta[1:end-1]; (n_est[end-1]-n_est[end])./delta[end]])
        n_deriv[:] = conv(n_deriv,[1.; 4.0; 6.0; 4.0; 1.0])[3:end-2]/16.0;
        n_deriv[:] = n_deriv.*(n_deriv.>0.0) # the convolution may introduce some small negative values

        # discretization error
        meas_disc_err[:] = 0.5*H_DMATRAIN*(1.0e9*delta.*(n_deriv.*delta)) # discretization error

        # error due to errors in the measurement model
        if GT_loaded
            meas_mode_err[:] = 0.5*(0.2*H_DMATRAIN+dH_SMALL)*(1.0e9*delta.*n_est) # assuming 20% error in the operator values
        else
            meas_mode_err[:] = 0.5*(0.5*H_DMATRAIN+dH_SMALL)*(1.0e9*delta.*n_est) # assuming 50% error in the operator values
        end

        # overall modelization error
        err_model[:] = (meas_disc_err+meas_mode_err).^2 .+ (reg_backcground_count).^2
    end
    # overall uncertainty
    data[:] + err_model[:];
end

function EKF.update_data_covariance(diagPoisson::Array{Cdouble,1})
    # update the data covariance matrix
    diagm(data_var(diagPoisson))
end
