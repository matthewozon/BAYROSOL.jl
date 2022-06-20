# This file is licensed under the MIT "Expat" License:

# Copyright (c) 2021: Matthew Ozon.

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


##################################################
### covariance and Cholesky factor: root model ###
##################################################

# compute the covaraince of correlated process (in diameter, not time) using a second order filter model
function covariance_psf_2nd(r_pol_::Cdouble,N_::Int64) #r_pol_::root of the filter, N_::length of the process
    # second order filter
    beta_ = poly([r_pol_, r_pol_])
    # impulse response of the filter
    impres_ = filt(1,reverse(beta_.a,dims=1),[1;zeros(2N_)]) # flit(b=[b0,b1,..., bN],a=[a0,a1,...aM],x) filter the sequence x using the filter of transfer function H(z^{-1}) = \frac{b0 + b1*z^{-1} + ... + bN*z^{-N}}{a0 + a1*z^{-1} + ... + aM*z^{-M}}
    # autocovariance of the impulse response
    xc_ = StatsBase.autocov(impres_,collect(0:2N_)) # size(impres_,1)-1
    sig_ = xc_[1]
    xcc_ = xc_[1:N_];
    # covariance
    # Gamma_atemp_ = Toeplitz([reverse(xcc_[2:end],dims=1); xcc_]);
    Gamma_atemp_ = Toeplitz(xcc_, xcc_);
    # return the covariance matrix and the autocovariance at index 0
    Gamma_atemp_, sig_
end

# compute the covariance of a correlated process with desired variances for each element
function covariance_process_2nd(r_pol_::Cdouble,D_tilde_::Array{Cdouble,1}) #r_pol_::root of the filter, d_::desired variance of the process
    # get the size of the
    N_ = length(D_tilde_)
    # compute a temporary covariance with the right shape
    Gamma_atemp_, sig_ = covariance_psf_2nd(r_pol_,N_)
    # normalize and scale the covariance with the right values
    D_bar_ =  sig_*ones(N_)
    D_bar_isqrt_ = diagm(sqrt.(D_bar_).^(-1));
    D_tilde_sqrt_ = diagm(sqrt.(D_tilde_));
    Gamma_atilde_ = D_tilde_sqrt_*D_bar_isqrt_*convert(Matrix{Cdouble},Gamma_atemp_)*D_bar_isqrt_*D_tilde_sqrt_;
    # return the covariance matrix
    Gamma_atilde_
end

# compute the asymptotic covriance of the source term of a vectorial first order stochastic process whose variances are defined by D_tilde_ and spacial coherence by the root r_pol_ of the second order filter
function covariance_2nd_order_space_1st_order_time(r_pol_::Cdouble,D_tilde_::Array{Cdouble,1},B_::Array{Cdouble,2})
    # get covariance of the process (asymptotic)
    Gamma_omega_ = covariance_process_2nd(r_pol_,D_tilde_)
    # compute the covariance of the source
    Gamma_omega_ - B_*Gamma_omega_*B_'
end

# draw time evolution sample with a first order time evolution model #TODO export
function space_covariance_chol(r_pol_::Cdouble,D_tilde_::Array{Cdouble,1},B_::Array{Cdouble,2})
    Gamma_noise_ = covariance_2nd_order_space_1st_order_time(r_pol_,D_tilde_,B_)
    # and its Cholesky decomposition
    L_noise_ = chol_PS_L(Gamma_noise_)
    # return
    L_noise_,Gamma_noise_
end

function space_covariance_chol(r_pol_::Cdouble,D_tilde_::Array{Cdouble,1})
    # get covariance of the process (asymptotic)
    Gamma_omega_ = covariance_process_2nd(r_pol_,D_tilde_)
    # and its Cholesky decomposition
    L_noise_ = chol_PS_L(Gamma_omega_)
    # return
    L_noise_,Gamma_omega_
end



###################################################
### covariance and Cholesky factor: array model ###
###################################################

# compute the covaraince of correlated process (in diameter, not time) using a user defined model
function covariance_psf_array(f_array::Array{Cdouble,1})
    # return the covariance matrix and the variance at index 0
    # Toeplitz([reverse(f_array[2:end],dims=1); f_array]), f_array[1]
    Toeplitz(f_array,f_array), f_array[1]
end

# compute the covariance of a correlated process with desired variances for each element
function covariance_process_array(f_array::Array{Cdouble,1},D_tilde_::Array{Cdouble,1})
    # get the size of the
    N_ = length(D_tilde_)
    # compute a temporary covariance with the right shape
    Gamma_atemp_, sig_ = covariance_psf_array(f_array)
    # normalize and scale the covariance with the right values
    D_bar_ =  sig_*ones(N_)
    D_bar_isqrt_ = diagm(sqrt.(D_bar_).^(-1));
    D_tilde_sqrt_ = diagm(sqrt.(D_tilde_));
    Gamma_atilde_ = D_tilde_sqrt_*D_bar_isqrt_*convert(Matrix{Cdouble},Gamma_atemp_)*D_bar_isqrt_*D_tilde_sqrt_;
    # return the covariance matrix
    Gamma_atilde_
end

# compute the asymptotic covriance of the source term of a vectorial first order stochastic process whose variances are defined by D_tilde_ and spacial coherence by the P.S.F. f_array
function covariance_psf_space_1st_order_time(f_array::Array{Cdouble,1},D_tilde_::Array{Cdouble,1},B_::Array{Cdouble,2})
    # get covariance of the process (asymptotic)
    Gamma_omega_ = covariance_process_array(f_array,D_tilde_)
    # compute the covariance of the source
    Gamma_omega_ - B_*Gamma_omega_*B_'
end

# draw time evolution sample with a first order time evolution model #TODO export
function space_covariance_array(f_array::Array{Cdouble,1},D_tilde_::Array{Cdouble,1},B_::Array{Cdouble,2})
    Gamma_noise_ = covariance_psf_space_1st_order_time(f_array,D_tilde_,B_)
    # and its Cholesky decomposition
    L_noise_ = chol_PS_L(Gamma_noise_)
    # return
    L_noise_,Gamma_noise_
end

function space_covariance_array(f_array::Array{Cdouble,1},D_tilde_::Array{Cdouble,1})
    # get covariance of the process (asymptotic)
    Gamma_omega_ = covariance_process_array(f_array,D_tilde_)
    # and its Cholesky decomposition
    L_noise_ = chol_PS_L(Gamma_omega_)
    # return
    L_noise_,Gamma_omega_
end


##############################################################
### covariance and Cholesky factor: range constraint model ###
##############################################################

# in this model, we assume that each variable x_i is probably in the range defined by the previous variable [(1-eps)*x_{i-1},(1+eps)*x_{i-1}]

function covariance_range_limit(mu_g::Cdouble,sig_g::Cdouble,eps::Cdouble,N_::Int64)
    # mu_g, sig_g: mean value and std of the first variable (0.0 and 1.0 seem to be a good start)
    # eps:         relative deviation between each consecutive variable
    # N_:          number of variables, size of the vector/matrix
    gamma_ = Array(Cdouble,N_,N_)
    for i in 1:N_
        gamma_[i,1:i] = (mu_g^2 + sig_g^2)*(1.0+eps^2).^(collect(1:i)) - mu_g^2
        gamma_[1:i-1,i] = gamma_[i,1:i-1]
    end
    gamma_,diag(gamma_)
end

# compute the covariance of a correlated process with desired variances for each element
function covariance_process_range_limit(mu_g::Cdouble,sig_g::Cdouble,eps::Cdouble,D_tilde_::Array{Cdouble,1})
    # get the size of the
    N_ = length(D_tilde_)
    # compute a temporary covariance with the right shape
    Gamma_atemp_, D_bar_ = covariance_range_limit(mu_g,sig_g,eps,N_)
    # normalize and scale the covariance with the right values
    D_bar_isqrt_ = diagm(sqrt.(D_bar_).^(-1));
    D_tilde_sqrt_ = diagm(sqrt.(D_tilde_));
    Gamma_atilde_ = D_tilde_sqrt_*D_bar_isqrt_*convert(Matrix{Cdouble},Gamma_atemp_)*D_bar_isqrt_*D_tilde_sqrt_;
    # return the covariance matrix
    Gamma_atilde_
end

# compute the asymptotic covriance of the source term of a vectorial first order stochastic process whose variances are defined by D_tilde_ and spacial coherence by the P.S.F. f_array
function covariance_range_limit_space_1st_order_time(mu_g::Cdouble,sig_g::Cdouble,eps::Cdouble,D_tilde_::Array{Cdouble,1},B_::Array{Cdouble,2})
    # get covariance of the process (asymptotic)
    Gamma_omega_ = covariance_process_range_limit(mu_g,sig_g,eps,D_tilde_)
    # compute the covariance of the source
    Gamma_omega_ - B_*Gamma_omega_*B_'
end

# draw time evolution sample with a first order time evolution model #TODO export
function space_covariance_range_limit(mu_g::Cdouble,sig_g::Cdouble,eps::Cdouble,D_tilde_::Array{Cdouble,1},B_::Array{Cdouble,2})
    Gamma_noise_ = covariance_range_limit_space_1st_order_time(mu_g,sig_g,eps,D_tilde_,B_)
    # and its Cholesky decomposition
    L_noise_ = chol_PS_L(Gamma_noise_)
    # return
    L_noise_,Gamma_noise_
end

function space_covariance_range_limit(mu_g::Cdouble,sig_g::Cdouble,eps::Cdouble,D_tilde_::Array{Cdouble,1})
    # get covariance of the process (asymptotic)
    Gamma_omega_ = covariance_process_range_limit(mu_g,sig_g,eps,D_tilde_)
    # and its Cholesky decomposition
    L_noise_ = chol_PS_L(Gamma_omega_)
    # return
    L_noise_,Gamma_omega_
end




function covariance_second_order_process(alpha::Cdouble,beta::Cdouble,sig1::Cdouble,sig2::Cdouble,sig12::Cdouble,var_vec::Array{Cdouble,1})
    # check if the process is convergent (but don't stop the creation of the matrix if it's not)
    if( !((alpha>=0.0) & (alpha^2>=(-4.0beta)) & (beta<=0.0) & ((alpha+beta)<1.0) & ((alpha^2+beta^2)<1.0)) )
        println("The stochastic process may not converge")
    end

    # create an empty matrix
    N = length(var_vec)
    gamma = Array(Cdouble,N,N)

    # initialize the main diagonal and the first upper and lower diagonal
    gamma[1,1] = sig1;  gamma[1,2] = sig12
    gamma[2,1] = sig12; gamma[2,2] = sig2;
    for i in 3:N
        gamma[i,i] = alpha^2*gamma[i-1,i-1] + beta^2*gamma[i-2,i-2] + (1.0-alpha^2+beta^2)*var_vec[i]
        # lower off diagonal
        gamma[i,i-1] = alpha*gamma[i-1,i-1] + beta*gamma[i-2,i-1]
        # upper off diagonal
        gamma[i-1,i] = gamma[i,i-1]
    end

    # fill the rest of the covariance matrix
    for i in 1:N
        for k in 2:(N-i)
            gamma[i+k,i] = alpha*gamma[i+k-1,i] + beta*gamma[i+k-2,i]
            gamma[i,i+k] = gamma[i+k,i]
        end
    end
    # gamma
    diagm(sqrt.(var_vec))*diagm(sqrt.(diag(gamma)).^(-1))*gamma*diagm(sqrt.(diag(gamma)).^(-1))*diagm(sqrt.(var_vec)) # with the desired variance for each element
end







#TODO: must add one from correlation matrix
