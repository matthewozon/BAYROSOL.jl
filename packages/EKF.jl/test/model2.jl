using SpecialMatrices
using Polynomials
using StatsBase

# let's consider of simple coupled model
# TODO: write the mathematical model here!

# constant of the model
  # model parameters
alpha = 0.98
A_ev  = alpha
  # state evolution (we strongly believe in the evolution model of the variables)
e_1 = 0.00000000001
Lx  = e_1
  # parameter evolution (we do not trust the evolution model of the parameters)
r_g   = 0.99#99
s_g   = 0.2
Lg    = sqrt(1.0-r_g^2)*s_g # s_g #

  # Cholesky factor of the model evolution covariance
L = zeros(model_dim,model_dim)
L[1,1] = Lx
L[2,2] = Lg
  # Cholesky factor of the noise covariance matrix
s_1 = 0.01 # 0.000002
Ly  = s_1


# evolution model: the state vector is organized as such [x_1 alpha]^T
# So, x[1] = x_1 and x[2] = alpha

# the known variables: x_fil_ (the current state), dt_ (only used for time continuous discretization but should be removed from the argument list), and L_ (the lower triangular Cholesky factor of the covariance matrix of the evolution model)
function EKF.apply_model!(x_fil_::Array{Cdouble,1},dt_::Cdouble,x_pre_::Array{Cdouble,1})
    # current state: x_fil
    # next state:    x_pre
    x_pre_[1] = x_fil_[2]*x_fil_[1]
    x_pre_[2] = r_g*x_fil_[2]

    # return the predicted state
    x_pre_
end


# Jacobian
function EKF.set_jacobian!(x_fil_::Array{Cdouble,1},dt_::Cdouble,F_ev_::Array{Cdouble,2})
    #
    F_ev_[1,1] = x_fil_[2]; F_ev_[1,2] = x_fil_[1]
    F_ev_[2,1] = 0.0;       F_ev_[2,2] = r_g
    F_ev_
end

function EKF.update_jacobian!(x_fil_::Array{Cdouble,1},dt_::Cdouble,F_ev_::Array{Cdouble,2})
    # for safety... but this is useless
    F_ev_ = set_jacobian!(x_fil_,dt_,F_ev_)
    F_ev_
end

function EKF.set_measurement_jacobian!(H_me_::Array{Cdouble,2})
    H_me_[1,1] = 1.0; H_me_[1,2] = 0.0;
    H_me_
end


function EKF.apply_perfect_measure!(x_pre_::Array{Cdouble,1},H_me_::Array{Cdouble,2})
    H_me_*x_pre_
end


  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # error models: only if time varying
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function EKF.update_model_covariance()
    # update the model covariance matrix
end

function EKF.update_data_covariance(diagPoisson::Array{Cdouble,1})
    # update the data covariance matrix
end
