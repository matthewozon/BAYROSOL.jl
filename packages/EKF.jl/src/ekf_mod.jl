#
# ekf_mod.jl --
#
# ekf_mod.jl is part of the Julia implementation of the Kalman Filter in its Extended version.
#
#------------------------------------------------------------------------------
#
# This file is part of the EKF module which is licensed under the MIT "Expat" License:
#
# Copyright (C) 2018-2019,  Matthew Ozon.
#
#------------------------------------------------------------------------------

#TODO: one could add the relative tolerance to pinv argument list rtol = sqrt(eps(real(float(one(eltype(M)))))), and why not an absolute tolerance atol (but which one?)

# one iteration of the extended Kalman Filter
function one_step(myWSKF_::wsKF,dt_::Cdouble)
    #++++++++++++++++++++
    # model updates
    #++++++++++++++++++++

    # the jacobian of the evolution equations
    myWSKF_.F_ev = update_jacobian!(myWSKF_.x_fil,dt_,myWSKF_.F_ev) # first order expansion in x_fil

    # the model covariance matrix: how precise is the model
    if myWSKF_.time_dep_noise_model
        myWSKF_.gam_w = update_model_covariance()
    end

    # the data covariance matrix: how much can we trust the data
    if myWSKF_.time_dep_noise_data
       myWSKF_.gam_v = update_data_covariance(myWSKF_.y_t)
    end

    # the measurement operator
    if myWSKF_.time_dep_meas_op
        myWSKF_.H_me = update_measurement_jacobian!(myWSKF_.H_me)
    end

    #++++++++++++++++++++
    # PREDICTION
    #++++++++++++++++++++
    # predicted state and covariance
    myWSKF_ = prediction!(myWSKF_,dt_)

    #++++++++++++++++++++
    # KALMAN GAIN
    #++++++++++++++++++++
    myWSKF_ = kalman_gain(myWSKF_)

    #++++++++++++++++++++
    # FILTERING
    #++++++++++++++++++++
    # filtered state and variance
    myWSKF_ = filtering!(myWSKF_,dt_)
    myWSKF_
end

# prediction step: state, covriance and measure
function prediction!(myWSKF_::wsKF,dt_::Cdouble)
    myWSKF_.x_pre = apply_model!(myWSKF_.x_fil,dt_,myWSKF_.x_pre)
    myWSKF_.o_pre = myWSKF_.F_ev*myWSKF_.o_fil*(myWSKF_.F_ev') + myWSKF_.gam_w
    myWSKF_.y_pre = apply_perfect_measure!(myWSKF_.x_pre,myWSKF_.H_me)
    myWSKF_
end

# compute the Kalman gain
function kalman_gain(myWSKF_::wsKF)
    myWSKF_.K_m = myWSKF_.o_pre*(myWSKF_.H_me')*pinv(myWSKF_.H_me*myWSKF_.o_pre*(myWSKF_.H_me') + myWSKF_.gam_v) #TODO why not try with / instead of pinv or inv?
    myWSKF_
end

# filtering step
function filtering!(myWSKF_::wsKF,dt_::Cdouble)
    myWSKF_.x_fil = myWSKF_.x_pre + myWSKF_.K_m*(myWSKF_.y_t-myWSKF_.y_pre)
    myWSKF_.o_fil = (Matrix{Cdouble}(I,myWSKF_.model_dim,myWSKF_.model_dim)-myWSKF_.K_m*myWSKF_.H_me)*myWSKF_.o_pre
    myWSKF_
end


function one_step_fixed_interval(myWSKF_::wsKF, i::Int64,x_smo_::Array{Cdouble,1},o_smo_::Array{Cdouble,2}, tmpMat::Array{Cdouble,2}) #TODO update the jacobian of the evolution model at each step instead of the stored ones
    # compute smoothing gain
    tmpMat = myWSKF_.o_fil_all[:,:,i]*myWSKF_.F_ev_all[:,:,i+1]'*pinv(myWSKF_.o_pre_all[:,:,i+1]) # it is not efficient to allocate a new matrix at each iteration

    # smoothing
    x_smo_next = myWSKF_.x_fil_all[:,i]   + tmpMat*(x_smo_-myWSKF_.x_pre_all[:,i+1])
    o_smo_next = myWSKF_.o_fil_all[:,:,i] + tmpMat*(o_smo_-myWSKF_.o_pre_all[:,:,i+1])*tmpMat'

    # return
    x_smo_next,o_smo_next
end











# the usual Kalman Filter

function KF(x00::Array{Cdouble,1},gam00_vec::Array{Cdouble,1},var_model::Array{Cdouble,2},var_noise::Array{Cdouble,2},t0_::Cdouble,myWSKF_::wsKF, y_all_::Array{Cdouble,2},n_samp_::Int64)
    # set the initial state of the Kalman Filter
    # initial state (x_fil) and covariance matrix (o_fil)
    myWSKF_.x_fil = x00
    myWSKF_.o_fil = diagm(gam00_vec)

    # this operation shouldn't be useful
    # myWSKF_.x_pre[:] .= 0.0
    # myWSKF_.o_pre[:,:] .= 0.0

    # noise and model covariances (if time dependent, updated at each step)
    if !myWSKF_.time_dep_noise_data
        myWSKF_.gam_v = var_noise
    end
    # if !myWSKF_.time_dep_noise_model
        myWSKF_.gam_w = var_model
    # end

    # before setting the jacobian of the model, x_fil must be initialized (first order expansion of the model arround x_fil)
    myWSKF_.F_ev = set_jacobian!(myWSKF_.x_fil,(myWSKF_.t_samp[2]-myWSKF_.t_samp[1])/t0_,myWSKF_.F_ev)
    # measurement operator
    if !myWSKF_.time_dep_meas_op
        myWSKF_.H_me = set_measurement_jacobian!(myWSKF_.H_me)
    end

    # start recursion
    for myWSKF_.ikf in 1:n_samp_
        # write the iteration number
        # println("iteration: ", i, " over ", n_samp_)

        # update time step
        if myWSKF_.ikf==n_samp_
            dt = (myWSKF_.t_samp[myWSKF_.ikf]-myWSKF_.t_samp[myWSKF_.ikf-1])/t0_
        else
            dt = (myWSKF_.t_samp[myWSKF_.ikf+1]-myWSKF_.t_samp[myWSKF_.ikf])/t0_
        end

        # update the observation
        myWSKF_.y_t = y_all_[:,myWSKF_.ikf]

        # run one step of the Kalman filter
        myWSKF_ = one_step(myWSKF_,dt)
        myWSKF_.x_fil_all[:,myWSKF_.ikf] = myWSKF_.x_fil
        myWSKF_.x_pre_all[:,myWSKF_.ikf] = myWSKF_.x_pre
        myWSKF_.o_fil_all[:,:,myWSKF_.ikf] = myWSKF_.o_fil
        myWSKF_.o_pre_all[:,:,myWSKF_.ikf] = myWSKF_.o_pre
        myWSKF_.F_ev_all[:,:,myWSKF_.ikf]  = myWSKF_.F_ev
        myWSKF_.H_me_all[:,:,myWSKF_.ikf]  = myWSKF_.H_me
        myWSKF_.gam_w_all[:,:,myWSKF_.ikf]  = myWSKF_.gam_w
        myWSKF_.gam_v_all[:,:,myWSKF_.ikf]  = myWSKF_.gam_v
    end

    # return the workspace
    myWSKF_
end



function KF_FIS(x00::Array{Cdouble,1},gam00_vec::Array{Cdouble,1},var_model::Array{Cdouble,2},var_noise::Array{Cdouble,2},t0_::Cdouble,myWSKF_::wsKF, y_all_::Array{Cdouble,2},n_samp_::Int64)
    # usual KF
    myWSKF_ = KF(x00,gam00_vec,var_model,var_noise,t0_,myWSKF_, y_all_,n_samp_)

    # Fixed Interval Smoother
    # init
    x_smo_all = similar(myWSKF_.x_fil_all)
    o_smo_all = similar(myWSKF_.o_fil_all)
    x_smo_all[:,n_samp_]   = myWSKF_.x_fil_all[:,n_samp_]
    o_smo_all[:,:,n_samp_] = myWSKF_.o_fil_all[:,:,n_samp_]
    tmpMat = similar(myWSKF_.o_fil)
    # run
    for myWSKF_.ikf in collect((n_samp_-1):-1:1)
        # write the iteration number
        # println("smoother iteration: ", i, " over ", n_samp_)

        # smooth
        x_smo_all[:,myWSKF_.ikf], o_smo_all[:,:,myWSKF_.ikf] = one_step_fixed_interval(myWSKF_,myWSKF_.ikf,x_smo_all[:,myWSKF_.ikf+1],o_smo_all[:,:,myWSKF_.ikf+1],tmpMat)
    end

    # return
    myWSKF_, x_smo_all, o_smo_all
end
