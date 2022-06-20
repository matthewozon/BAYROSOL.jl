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

# building block of the time evolution model for the parameters
function time_evolution_matrix(K::Int64,ndim::Int64,dt::Cdouble,T_char::Cdouble;zeta_damp::Cdouble=0.95)
    # K:         order of the process, 0 or 1 or 2 (more than 2 is risky so I don't use it) (here I use 0 for the random walk)
    # ndim:      dimension of the process (1 for a scalar process like the nucleation and more than one for, e.g., the condensation rates)
    # T_char:    characteristic time of the process
    # zeta_damp: damping coeficient in [0,1)
    # dt:        discretization time step

    B_time = zeros(Cdouble,max(1,K)*ndim,max(1,K)*ndim)
    if (K==1) # first order process
        r1 = min(1.0,max(0.0,1.0 - abs(dt/T_char))) # just making sure that 1) the system does not go backward and 2) that it is stable
        for ii in 1:ndim
            B_time[ii,ii] = r1
        end
    elseif (K==2) # second order process
        T_char_star = 2.0pi*T_char; # just an artifact due to the way I wrote the second order differential equation
        # compute the sum and product of the root of the discrete second order sytem
        # WARNING: if T_char<=dt, it's not gonna make sense, I should safeguard this computation
        r1pr2 = 2.0*(1.0-2.0pi*zeta_damp*dt/T_char_star);
        r1mr2 = 1.0-(4.0pi*zeta_damp*dt/T_char_star)+(4.0*pi^2*(dt/T_char_star)^2);
        for ii in 1:ndim
            B_time[2ii-1,2ii-1] = r1pr2
            B_time[2ii-1,2ii]   = -r1mr2
            B_time[2ii,2ii-1]   = 1.0
            B_time[2ii,2ii]     = 0.0
        end
    else # assume random walk K=0
        for ii in 1:ndim
            B_time[ii,ii] = 1.0
        end
    end
    B_time
end

function covariance_matrix_time(K::Int64,ndim::Int64,sig_proc::Union{Cdouble,Array{Cdouble,1}},B_time::Union{Cdouble,Array{Cdouble,2}};cor_len::Cdouble=1.0)
    # compute the (co)variance of the process
    gamma_proc = zeros(Cdouble,ndim,ndim);
    if (ndim==1)
        gamma_proc[1,1] = sig_proc^2
    else
        if cor_len>0.0
            if (length(sig_proc)==1)
                f_array_cond    = exp.(-collect(0.0:ndim-1.0)/cor_len).^2
                gamma_proc[:,:] = covariance_process_array(f_array_cond,sig_proc^2*ones(Cdouble,ndim))
            else
                f_array_cond    = exp.(-collect(0.0:ndim-1.0)/cor_len).^2
                gamma_proc[:,:] = covariance_process_array(f_array_cond,sig_proc.^2)
            end
        else
            if (length(sig_proc)==1)
                gamma_proc[:,:] = (sig_proc^2)*Matrix{Cdouble}(I,ndim,ndim)
            else
                gamma_proc[:,:] = diagm(sig_proc.^2)
            end
        end
    end

    # compute the (co)variance of the noise term in the stochastic process model
    if (K==0) # if one does not want to account for the time evolution
        gamma_proc_time = gamma_proc
    elseif (K==1)
        gamma_proc_time = gamma_proc - B_time*gamma_proc*B_time'
    else # (K==2)
        gamma_proc_XXL = zeros(Cdouble,2*ndim,2*ndim);
        gamma_proc_XXL[1:2:end,1:2:end] = gamma_proc
        gamma_proc_XXL[1:2:end,2:2:end] = gamma_proc
        gamma_proc_XXL[2:2:end,1:2:end] = gamma_proc
        gamma_proc_XXL[2:2:end,2:2:end] = gamma_proc
        gamma_proc_time                 = gamma_proc_XXL - B_time*gamma_proc_XXL*B_time'
        # MUST DO: force the uncorrelated variables to 0
        gamma_proc_time[1:2:end,2:2:end] .= 0.0
        gamma_proc_time[2:2:end,2:2:end] .= 0.0
        gamma_proc_time[2:2:end,1:2:end] .= 0.0
    end
    #
    gamma_proc_time
end


mutable struct NASHA
    # the main obejcts (cf AeroMec2.jl and EKF.jl)
    # aero::AeroSys
    # estk::EKF

    # some ranges that will locate the different physical quantities in the array
    R_psd::Union{UnitRange{Int64},StepRange{Int64,Int64}}  # size distribution range (WARNING: size distribution refers to number concentration here)
    R_cond::Union{UnitRange{Int64},StepRange{Int64,Int64}} # condensation/evaporation
    R_loss::Union{UnitRange{Int64},StepRange{Int64,Int64}} # linear loss
    R_nuc::Union{UnitRange{Int64},StepRange{Int64,Int64}}  # nucleation

    # order of the evolution processes
    K_cond::Int64 # condensation stochastic process order
    K_loss::Int64 # loss
    K_nuc::Int64  # nucleation

    # time evolution matrices
    B_cond_time::Array{Cdouble,2}
    B_loss_time::Array{Cdouble,2}
    B_nuc_time::Array{Cdouble,2}

    # covariance matrices #WARNING: I assume that each mechanism is independent from one another (maybe condensation and nucleation are not)
    gamma_psd::Array{Cdouble,2}
    gamma_cond::Array{Cdouble,2}
    gamma_loss::Array{Cdouble,2}
    gamma_nuc::Array{Cdouble,2}

    function NASHA(nbin_::Int64,dt_::Cdouble;
        n_cond::Int64=0,n_loss::Int64=0,n_nuc::Int64=0,
        K_cond::Int64=2,K_loss::Int64=0,K_nuc::Int64=2,
        T_cond::Cdouble=Inf,T_loss::Cdouble=Inf,T_nuc::Cdouble=Inf,
        cor_len_psd::Cdouble=-1.0,cor_len_cond::Cdouble=-1.0,cor_len_loss::Cdouble=-1.0,cor_len_nuc::Cdouble=-1.0,
        sig_psd::Union{Cdouble,Array{Cdouble,1}}=0.0,sig_cond::Union{Cdouble,Array{Cdouble,1}}=0.0,sig_loss::Union{Cdouble,Array{Cdouble,1}}=0.0,sig_nuc::Union{Cdouble,Array{Cdouble,1}}=0.0)
        # psd range
        R_psd = 1:nbin_
        # condensation range
        R_cond_init = R_psd[end] + 1
        R_cond_end  = R_cond_init + max(1,K_cond)*n_cond - 1
        R_cond = R_cond_init:R_cond_end
        # loss range
        R_loss_init = R_cond_end + 1
        R_loss_end = R_loss_init + max(1,K_loss)*n_loss - 1
        R_loss = R_loss_init:R_loss_end
        # nucleation range
        R_nuc_init = R_loss_end + 1
        R_nuc_end = R_nuc_init + max(1,K_nuc)*n_nuc - 1
        R_nuc  = R_nuc_init:R_nuc_end
        # time evolution matrices
        B_cond_time = time_evolution_matrix(K_cond,n_cond,dt_,T_cond;zeta_damp=0.95)
        B_loss_time = time_evolution_matrix(K_loss,n_loss,dt_,T_loss;zeta_damp=0.95)
        B_nuc_time  = time_evolution_matrix(K_nuc,n_nuc,dt_,T_nuc;zeta_damp=0.95)
        # covariance matrices
        if (length(sig_psd)==1)
            gamma_psd = covariance_matrix_time(0,nbin_,sig_psd*ones(Cdouble,nbin_),0.0;cor_len=cor_len_psd)
        else
            gamma_psd = covariance_matrix_time(0,nbin_,sig_psd,0.0;cor_len=cor_len_psd)
        end
        if (length(sig_cond)==n_cond)
            gamma_cond = covariance_matrix_time(K_cond,n_cond,sig_cond,B_cond_time;cor_len=cor_len_cond)
        else
            if (length(sig_cond)==1)
                gamma_cond = covariance_matrix_time(K_cond,n_cond,sig_cond*ones(Cdouble,n_cond),B_cond_time;cor_len=cor_len_cond)
            else
                # I don't know what to do
                throw("NASHA: don't know how to build the condensation rate covariance matrix")
            end
        end
        if (length(sig_loss)==n_loss)
            gamma_loss = covariance_matrix_time(K_loss,n_loss,sig_loss,B_loss_time;cor_len=cor_len_loss)
        else
            if (length(sig_loss)==1)
                gamma_loss = covariance_matrix_time(K_loss,n_loss,sig_loss*ones(Cdouble,n_loss),B_loss_time;cor_len=cor_len_loss)
            else
                # I don't know what to do
                throw("NASHA: don't know how to build the loss rate covariance matrix")
            end
        end
        if (length(sig_nuc)==n_nuc)
            gamma_nuc = covariance_matrix_time(K_nuc,n_nuc,sig_nuc,B_nuc_time;cor_len=cor_len_nuc)
        else
            if (length(sig_nuc)==1)
                gamma_nuc = covariance_matrix_time(K_nuc,n_nuc,sig_nuc*ones(Cdouble,n_nuc),B_nuc_time;cor_len=cor_len_nuc)
            else
                # I don't know what to do
                throw("NASHA: don't know how to build the nucleation rate covariance matrix")
            end
        end
        # now you can instanciate the structure
        new(R_psd,R_cond,R_loss,R_nuc,
        K_cond,K_loss,K_nuc,
        B_cond_time,B_loss_time,B_nuc_time,
        gamma_psd,gamma_cond,gamma_loss,gamma_nuc)
    end

    # function AeroSys(_d::Array{Cdouble,1};beta_c::Cdouble=0.0, GR_c::Cdouble=0.0, J_c::Cdouble=0.0, gamma_c::Cdouble=0.0, t_c::Cdouble=1.0, x_c::Cdouble=1.0, scale::String="log")
    # function wsKF(_mod::Int64, _meas::Int64,_time::Int64, _noise_model::Bool, _noise_data::Bool)
end
