using SpecialMatrices
using Polynomials
using StatsBase

# let's consider of simple coupled model
# TODO: write the mathematical model here!

# constant of the model
  # model parameters
alpha = 0.98
beta  = 0.97
gamma = 0.1
A_ev  = [alpha -gamma; gamma beta]
  # state evolution (we strongly believe in the evolution model of the variables)
e_1   = 0.00000000001
e_2   = 0.00000000001
Lx = [e_1 0.0; 0.0 e_2]
  # parameter evolution (we do not trust the evolution model of the parameters)
r_g   = 0.99
s_g   = 0.05
Lg    = sqrt(1.0-r_g^2)*s_g
r_p   = 0.99
s_a   = 10.1
s_b   = 10.1
tau   = 0.9
s_ab  = tau*s_a*s_b
Lp    = sqrt(1.0-r_p^2)*[s_a 0.0; tau*s_a s_b*sqrt(1.0-tau^2)]
  # Cholesky factor of the model evolution covariance
L = zeros(5,5)
L[1:2,1:2] = Lx
L[3,3]     = Lg
L[4:5,4:5] = Lp
  # Cholesky factor of the noise covariance matrix
s_1 = 0.000002
s_2 = 0.000002
Ly = [s_1 0.0; 0.0 s_2]


# evolution model: the state vector is organized as such [x_1 x_2 gamma alpha beta]^T
# So, x[1] = x_1, x[2] = x_2, x[3] = gamma, x[4] = alpha and x[5] = beta

# the known variables: x_fil_ (the current state), dt_ (only used for time continuous discretization but should be removed from the argument list), and L_ (the lower triangular Cholesky factor of the covariance matrix of the evolution model)
function EKF.apply_model!(x_fil_::Array{Cdouble,1},dt_::Cdouble,x_pre_::Array{Cdouble,1})
    # current state: x_fil
    # next state:    x_pre
    x_pre_[1] = x_fil_[4]*x_fil_[1] - x_fil_[3]*x_fil_[2]
    x_pre_[2] = x_fil_[3]*x_fil_[1] + x_fil_[5]*x_fil_[2]
    x_pre_[3] = r_g*x_fil_[3]
    x_pre_[4] = r_p*x_fil_[4]
    x_pre_[5] = r_p*x_fil_[5]


    # return the predicted state
    x_pre_
end


# Jacobian
function EKF.set_jacobian!(x_fil_::Array{Cdouble,1},dt_::Cdouble,F_ev_::Array{Cdouble,2})
    #
    F_ev_[1,1] = x_fil_[4]; F_ev_[1,2] = -x_fil_[3]; F_ev_[1,3] = -x_fil_[2]; F_ev_[1,4] = x_fil_[1]; F_ev_[1,5] = 0.0;

    F_ev_[2,1] = x_fil_[3]; F_ev_[2,2] = x_fil_[5];  F_ev_[2,3] = x_fil_[1];  F_ev_[2,4] = 0.0;       F_ev_[2,5] = x_fil_[2];

    F_ev_[3,1] = 0.0;       F_ev_[3,2] = 0.0;        F_ev_[3,3] = r_g;        F_ev_[3,4] = 0.0;       F_ev_[3,5] = 0.0;

    F_ev_[4,1] = 0.0;       F_ev_[4,2] = 0.0;        F_ev_[4,3] = 0.0;        F_ev_[4,4] = r_p;       F_ev_[4,5] = 0.0;

    F_ev_[5,1] = 0.0;       F_ev_[5,2] = 0.0;        F_ev_[5,3] = 0.0;        F_ev_[5,4] = 0.0;       F_ev_[5,5] = r_p;

    F_ev_
end

function EKF.update_jacobian!(x_fil_::Array{Cdouble,1},dt_::Cdouble,F_ev_::Array{Cdouble,2})
    # for safety... but this is useless
    F_ev_ = set_jacobian!(x_fil_,dt_,F_ev_)
    F_ev_
end

function EKF.set_measurement_jacobian!(H_me_::Array{Cdouble,2})
    H_me_[1,1] = 1.0; H_me_[1,2] = 0.0; H_me_[1,3] = 0.0; H_me_[1,4] = 0.0; H_me_[1,5] = 0.0;
    H_me_[2,1] = 0.0; H_me_[2,2] = 1.0; H_me_[2,3] = 0.0; H_me_[2,4] = 0.0; H_me_[2,5] = 0.0;
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
