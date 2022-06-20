#
# type.jl --
#
# type.jl is part of the Julia implementation of the Kalman Filter in its Extended version.
#
#------------------------------------------------------------------------------
#
# This file is part of the EKF module which is licensed under the MIT "Expat" License:
#
# Copyright (C) 2018-2019,  Matthew Ozon.
#
#-----------------------------------------------------------------------------

# this file is part of the EKF module. It encodes the type which contains the workspace of the Kalman Filter

mutable struct wsKF

    # dimension of the problem
    model_dim::Int64            # evolution model dimension
    meas_dim::Int64             # measurement model dimension
    n_samp::Int64               # number of time samples

    # variables of the problem
    x_fil::Array{Cdouble,1}     # current state vector: filter step
    x_pre::Array{Cdouble,1}     # current state vector: prediction step
    o_fil::Array{Cdouble,2}     # current state covariance matrix: filter step
    o_pre::Array{Cdouble,2}     # current state covariance matrix: prediction step
    y_pre::Array{Cdouble,1}     # predicted measurement of the current state


    # variables for KF
    y_t::Array{Cdouble,1}       #  measurement observation

    # jacobian of the evolution and measurement model
    F_ev::Array{Cdouble,2}      # linear model of the pendulum system evolution
    H_me::Array{Cdouble,2}      # measurement device model

    # calculation intermediates
    K_m::Array{Cdouble,2}       # Kalman gain

    # covariances of the models
    gam_w::Array{Cdouble,2}     # covariance of the evolution model
    sig_w::Array{Cdouble,2}     # a matrix that satisfies gam_w=sig_w*sig_w or gam_w=sig_w'*sig_w  # TODO: this variable is not used, it should be deleted (though it could be usefull for the EnKF)
    gam_v::Array{Cdouble,2}     # covariance matrix of the measures


    # memory of the time series
    t_samp::Array{Cdouble,1}  # time index
    x_fil_all::Array{Cdouble,2} # should be n_samp+1 because it should include the initial guess
    x_pre_all::Array{Cdouble,2}
    o_fil_all::Array{Cdouble,3} # idem: should be n_samp+1 because it should include the initial guess
    o_pre_all::Array{Cdouble,3}
    F_ev_all::Array{Cdouble,3}  # time series of the Jacobian matrices
    H_me_all::Array{Cdouble,3}  # time series of the measurement operator matrices

    # error model time memory
    gam_w_all::Array{Cdouble,3}     # covariance of the evolution model
    sig_w_all::Array{Cdouble,3}     # a matrix that satisfies gam_w=sig_w*sig_w or gam_w=sig_w'*sig_w # TODO: this variable is not used, it should be deleted (though it could be usefull for the EnKF)
    gam_v_all::Array{Cdouble,3}     # covariance matrix of the measures

    # Kalman filter flags
    time_dep_noise_model::Bool  # if gam_w depends on time (or data)
    time_dep_noise_data::Bool   # if gam_v depends on time (or data)
    time_dep_meas_op::Bool      # if the measurement operator depends on time

    # iteration index
    ikf::Int64


    # default ctor (it is not really meaningful)
    function wsKF() #
        new(0,0,0,
            Array{Cdouble,1}(undef,0), Array{Cdouble,1}(undef,0), Array{Cdouble,2}(undef,0,0), Array{Cdouble,2}(undef,0,0), Array{Cdouble,1}(undef,0),
            Array{Cdouble,1}(undef,0), Array{Cdouble,2}(undef,0,0), Array{Cdouble,2}(undef,0,0), Array{Cdouble,2}(undef,0,0),
            Array{Cdouble,2}(undef,0,0), Array{Cdouble,2}(undef,0,0),Array{Cdouble,2}(undef,0,0), #
            Array{Cdouble,1}(undef,0), Array{Cdouble,2}(undef,0,0), Array{Cdouble,2}(undef,0,0), Array{Cdouble,3}(undef,0,0,0), Array{Cdouble,3}(undef,0,0,0), Array{Cdouble,3}(undef,0,0,0), Array{Cdouble,3}(undef,0,0,0),
            Array{Cdouble,3}(undef,0,0,0), Array{Cdouble,3}(undef,0,0,0), Array{Cdouble,3}(undef,0,0,0),
            false, false, false,
            0)
    end

    # ctor with known dimensions
    function wsKF(_mod::Int64, _meas::Int64,_time::Int64) #
        new(_mod,_meas,_time,
            Array{Cdouble,1}(undef,_mod), Array{Cdouble,1}(undef,_mod), Array{Cdouble,2}(undef,_mod,_mod), Array{Cdouble,2}(undef,_mod,_mod), Array{Cdouble,1}(undef,_meas),
            Array{Cdouble,1}(undef,_meas), Array{Cdouble,2}(undef,_mod,_mod), Array{Cdouble,2}(undef,_meas,_mod), Array{Cdouble,2}(undef,_mod,_meas),
            Array{Cdouble,2}(undef,_mod,_mod), Array{Cdouble,2}(undef,_mod,_mod),Array{Cdouble,2}(undef,_meas,_meas), #
            Array{Cdouble,1}(undef,_time), Array{Cdouble,2}(undef,_mod,_time), Array{Cdouble,2}(undef,_mod,_time), Array{Cdouble,3}(undef,_mod,_mod,_time), Array{Cdouble,3}(undef,_mod,_mod,_time), Array{Cdouble,3}(undef,_mod,_mod,_time), Array{Cdouble,3}(undef,_meas,_mod,_time),
            Array{Cdouble,3}(undef,_mod,_mod,_time), Array{Cdouble,3}(undef,_mod,_mod,_time), Array{Cdouble,3}(undef,_meas,_meas,_time),
            false, false, false,
            1)
    end

    # ctor: with known dimensions and time dependence of the models
    function wsKF(_mod::Int64, _meas::Int64,_time::Int64, _noise_model::Bool, _noise_data::Bool) #
        new(_mod,_meas,_time,
            Array{Cdouble,1}(undef,_mod), Array{Cdouble,1}(undef,_mod), Array{Cdouble,2}(undef,_mod,_mod), Array{Cdouble,2}(undef,_mod,_mod), Array{Cdouble,1}(undef,_meas),
            Array{Cdouble,1}(undef,_meas), Array{Cdouble,2}(undef,_mod,_mod), Array{Cdouble,2}(undef,_meas,_mod), Array{Cdouble,2}(undef,_mod,_meas),
            Array{Cdouble,2}(undef,_mod,_mod), Array{Cdouble,2}(undef,_mod,_mod),Array{Cdouble,2}(undef,_meas,_meas), #
            Array{Cdouble,1}(undef,_time), Array{Cdouble,2}(undef,_mod,_time), Array{Cdouble,2}(undef,_mod,_time), Array{Cdouble,3}(undef,_mod,_mod,_time), Array{Cdouble,3}(undef,_mod,_mod,_time), Array{Cdouble,3}(undef,_mod,_mod,_time), Array{Cdouble,3}(undef,_meas,_mod,_time),
            Array{Cdouble,3}(undef,_mod,_mod,_time), Array{Cdouble,3}(undef,_mod,_mod,_time), Array{Cdouble,3}(undef,_meas,_meas,_time),
            _noise_model, _noise_data, false,
            1)
    end

    # cptor
    function wsKF(ws::wsKF) #
        new(copy(ws.model_dim),copy(ws.meas_dim),copy(ws.n_samp),
            copy(ws.x_fil),copy(ws.x_pre),copy(ws.o_fil),copy(ws.o_pre),copy(ws.y_pre),
            copy(ws.y_t),copy(ws.F_ev),copy(ws.H_me),copy(ws.K_m),
            copy(ws.gam_w),copy(ws.sig_w),copy(ws.gam_v), #
            copy(ws.t_samp),copy(ws.x_fil_all),copy(ws.x_pre_all),copy(ws.o_fil_all),copy(ws.o_pre_all),copy(ws.F_ev_all),copy(ws.H_me_all),
            copy(ws.gam_w_all),copy(ws.sig_w_all),copy(ws.gam_v_all),
            copy(ws.time_dep_noise_model),copy(ws.time_dep_noise_data),copy(ws.time_dep_meas_op),
            copy(ws.ikf))
    end
end
