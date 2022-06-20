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


# estimate the variance of a transformed gaussian variable Y=f(X), X~N(mu,L*L')
function mean_and_varince_estimation(f::Function,mu::Union{Cdouble,Array{Cdouble,1}},L::Union{Cdouble,Array{Cdouble,2}},Ns::Int64)
    mu_y = 0.0
    gamma_y = 0.0
    if typeof(L)==Cdouble # scalar
        Y_samp = Array{Cdouble,1}(undef,Ns)
        for i in 1:Ns
            Y_samp[i] = f(L*randn()+mu)
        end
        mu_y = mean(Y_samp)
        gamma_y  = var(Y_samp)
    else # vector
        Nx = size(L,1)
        Y_samp = Array{Cdouble,2}(undef,Nx,Ns)
        for i in 1:Ns
            Y_samp[:,i] = f(L*randn(Nx)+mu)
        end
        mu_y = mean(Y_samp,2)
        gamma_y  = cov(Y_samp')
    end
    mu_y, gamma_y,Y_samp
end


# estimate the kth percentiles of a transformed gaussian variable Y=f(X), X~N(mu,sig^2)
function percentile_estimation(f::Function,mu::Cdouble,sig::Cdouble,Kth::Union{Int64,Array{Int64,1}}=[5;95],Ns::Int64=100000)
    Y_samp = Array{Cdouble,1}(undef,Ns)
    for i in 1:Ns
        Y_samp[i] = f(sig*randn()+mu)
    end
    sort!(Y_samp)
    idx = (Kth*Ns)/100.0
    if typeof(Kth)==Int # only one percentile
        idx_u = min(max(1,ceil.(Int64,idx)),Ns)
        idx_d = max(min(Ns,floor.(Int64,idx)),1)
        if idx_u==idx_d
            percentile = Y_samp[idx_u]
        else
            percentile = (idx_u-idx)*Y_samp[idx_d] +  (idx-idx_d)*Y_samp[idx_u]
        end
    else # vector
        idx_u = min(max(1*ones(Int64,length(Kth)),ceil.(Int64,idx)),Ns*ones(Int64,length(Kth)))
        idx_d = max(min(Ns*ones(Int64,length(Kth)),floor.(Int64,idx)),1*ones(Int64,length(Kth)))
        percentile = zeros(length(Kth))
        for i in 1:length(Kth)
            if idx_u[i]==idx_d[i]
                percentile[i] = Y_samp[idx_u[i]]
            else
                percentile[i] = (idx_u[i]-idx[i])*Y_samp[idx_d[i]] +  (idx[i]-idx_d[i])*Y_samp[idx_u[i]]
            end
        end
    end
    # mi = Y_samp[1]
    # ma = Y_samp[end]
    # histo(Y_samp,mi:(ma-mi)/20:ma,1)
    percentile
end

function percentile_estimation(f::Function,mu::Array{Cdouble,1},sig::Array{Cdouble,2},Kth::Union{Int64,Array{Int64,1}}=[5;95],Ns::Int64=100000) # this function is only suitable for independent variables
    # the function f must be a vectorial function (either scalar valued or vector valued)
    # mu is the vector of mean values
    # sig is the lower Cholesky factor (or square root) of a covariance matrix
    Nin = length(mu)
    N = length(f(mu))
    T_samp = Array{Cdouble,2}(undef,N,Ns)
    for i in 1:Ns
        T_samp[:,i] = f(sig*randn(Nin)+mu)
    end
    Y_samp = sort(T_samp,dims=2)
    idx = (Kth*Ns)/100.0
    if typeof(Kth)==Int # only one percentile
        idx_u = min(max(1,ceil.(Int64,idx)),Ns)
        idx_d = max(min(Ns,floor.(Int64,idx)),1)
        if idx_u==idx_d
            percentile = Y_samp[:,idx_u]
        else
            percentile = (idx_u-idx)*Y_samp[:,idx_d] +  (idx-idx_d)*Y_samp[:,idx_u]
        end
    else # vector
        idx_u = min(max(1*ones(Int64,length(Kth)),ceil.(Int64,idx)),Ns*ones(Int64,length(Kth)))
        idx_d = max(min(Ns*ones(Int64,length(Kth)),floor.(Int64,idx)),1*ones(Int64,length(Kth)))
        percentile = zeros(N,length(Kth))
        for i in 1:length(Kth)
            if idx_u[i]==idx_d[i]
                percentile[:,i] = Y_samp[:,idx_u[i]]
            else
                percentile[:,i] = (idx_u[i]-idx[i])*Y_samp[:,idx_d[i]] +  (idx[i]-idx_d[i])*Y_samp[:,idx_u[i]]
            end
        end
    end
    percentile
end


function histo(X::Array{Cdouble,1},R::StepRangeLen{Cdouble,Base.TwicePrecision{Cdouble},Base.TwicePrecision{Cdouble}},fig::Int64) # ,xerr::Union{Cdouble,Array{Cdouble,1}}=Union{},yerr::Union{Cdouble,Array{Cdouble,1}}=Union{})
    figure(fig)
    hist(X,R)
end




# compute the integral of sin^k over [0,pi]
function Ik(k::Int64)
    val = 0.0
    if k<0
        throw("Ik does not exist for negative integers")
    elseif k==0
        val = 1.0pi # just to make it real... pi is Irrational
    elseif k==1
        val = 2.0
    else
        val = ((k-1)/k)*Ik(k-2)
    end
    val
end

function Ik(k::Array{Int64,1})
    VAL = Array{Cdouble,1}(undef,length(k))
    for i in 1:length(k)
        try
            VAL[i] = Ik(k[i])
        catch(msgErr)
            println(msgErr)
            VAL[i] = NaN
        end
    end
    VAL
end



# compute an estimation of the integral of r^k exp(-0.5r^2) over [0,R]
function Hk(k::Int64,R::Cdouble)
    val = 0.0
    if k<0
        throw("Hk does not exist for negative integers")
    elseif k==0
        val = sqrt(pi/2.0)*erf(R/sqrt(2.0))
    elseif k==1
        val = 1.0-exp(-0.5*R^2)
    else
        if isinf(R)
            val = (k-1.0)*Hk(k-2,R)
        else
            val = (k-1.0)*Hk(k-2,R) - R^(k-1)*exp(-0.5*R^2)
        end
    end
    val
end

function Hk(k::Array{Int64,1},R::Cdouble)
    VAL = Array{Cdouble,1}(undef,length(k))
    for i in 1:length(k)
        try
            VAL[i] = Hk(k[i],R)
        catch(msgErr)
            println(msgErr)
            VAL[i] = NaN
        end
    end
    VAL
end

function Hk(k::Int64,R::Array{Cdouble,1})
    VAL = Array{Cdouble,1}(undef,length(R))
    for i in 1:length(R)
        try
            VAL[i] = Hk(k,R[i])
        catch(msgErr)
            println(msgErr)
            VAL[i] = NaN
        end
    end
    VAL
end


# compute the probability of the ball of radius R centered on 0_k for the multivariate gaussian distribution (0_k,I_k)
function multivariate_probability(k::Int64,R::Cdouble)
    p = -1.0
    if R<0.0
        throw("The radius of a ball must be a positive value")
    end
    if k<1
        throw("A probability density function must lie in a space with at least one dimension")
    elseif k==1
        p = erf(R/sqrt(2.0))
    elseif k==2
        p = 1.0-exp(-0.5*R^2)
    else # k>=3
        NormNorm = 2.0pi/sqrt((2.0pi)^k) # normalization factor and integration of the first angle over [0,2pi]
        IKK = Ik(collect(1:k-2))         # integration over the other angular variables
        p = NormNorm*Hk(k-1,R)*prod(IKK)
    end
    p
end




# conpute the variable change and its inverse for cartesian to spherical
function spherical_coordinate(N::Int64)
    g = Array{Function,1}(undef,N)
    g_inv = Array{Function,1}(undef,N)
    if N<1
        throw("A probability density function must lie in a space with at least one dimension")
    elseif N == 1 # nothing to do
        g[1] = (x::Cdouble->x); # radius
        g_inv[1] = (x::Cdouble->x); # radius_inv
    elseif N == 2 # polar coordinates (r,theta) \in [0,\infty)x[0,2pi) and (x,y) \in \mathbb{R}^2
        g[1] = (x::Array{Cdouble,1}->sqrt(x[1]^2+x[2]^2)); # radius
        g[2] = (x::Array{Cdouble,1}->atan2(x[2],x[1])); # phi
        g_inv[1] = (polar::Array{Cdouble,1}->polar[1]*cos(polar[2])); # coord_x
        g_inv[2] = (polar::Array{Cdouble,1}->polar[1]*sin(polar[2])); # coord_y
    elseif N == 3 # usual spherical transformation (r,theta,phi) \in [0,\infty)x[0,pi]x[0,2pi) and (x,y,z) \in \mathbb{R}^3
        g[1] = (x::Array{Cdouble,1}->sqrt(sum(x.^2))) # radius
        g[2] = (x::Array{Cdouble,1}->acos(x[3]/sqrt(sum(x.^2)))); # theta
        g[3] = (x::Array{Cdouble,1}->atan2(x[2],x[1])); # phi
        g_inv[1] = (spherical::Array{Cdouble,1}->spherical[1]*sin(spherical[2])*cos(spherical[3])); # coord_x
        g_inv[2] = (spherical::Array{Cdouble,1}->spherical[1]*sin(spherical[2])*sin(spherical[3])); # coord_y
        g_inv[3] = (spherical::Array{Cdouble,1}->spherical[1]*cos(spherical[2])); # coord_z
    else # hypergeometry (r,phi_1,ph_2,...,phi_{N-2},phi_{N-1}) \in [0,\infty)x[0,pi]^{N-2}x[0,2pi) and (x_1,x_2,...,x_N) \in \mathbb{R}^N
        # https://en.wikipedia.org/wiki/N-sphere#Spherical_coordinates
        # direct transformation
        # radius
        g[1] = x::Array{Cdouble,1}->sqrt(sum(x.^2))
        # angles
        for i in 1:N-2
            g[i+1] = x::Array{Cdouble,1}->acot(x[i]/sqrt(sum(x[i+1:end])))
        end
        g[N] = x::Array{Cdouble,1}->2.0*acot( (x[N-1]+sqrt(x[N]^2+x[N-1]^2))/x[N] )

        # inverse transformation
        g_inv[1] = spherical::Array{Cdouble,1}->spherical[1]*cos(spherical[2])
        for i in 1:N-2
            g_inv[i+1] = spherical::Array{Cdouble,1}->spherical[1]*prod(sin.(spherical[2:i+1]))*cos(spherical[i+2])
        end
        g_inv[N] = spherical::Array{Cdouble,1}->spherical[1]*prod(sin.(spherical[2:N]))
    end
    g,g_inv
end


# https://en.wikipedia.org/wiki/N-sphere#Generating_random_points
# uniform sampling on the (N-1)-sphere in R^N
function uniform_sample_unit_sphere(N::Int64)
    if N<2
        throw("Not sampling nothing")
    end
    X = randn(N)
    (1.0/sqrt(sum(X.^2)))*X
end

function uniform_sample_unit_sphere(N::Int64,Ns::Int64)
    S = Array{Cdouble,2}(undef,N,Ns)
    for i in 1:Ns
        S[:,i] = uniform_sample_unit_sphere(N)
    end
    S
end

# uniform sampling within the N-ball
function uniform_sample_unit_ball(N::Int64)
    if N<2
        throw("Not sampling nothing")
    end
    X = uniform_sample_unit_sphere(N)
    R = rand()
    R^(1.0/N)*X
end

function uniform_sample_unit_ball(N::Int64,Ns::Int64)
    B = Array{Cdouble,2}(undef,N,Ns)
    for i in 1:Ns
        B[:,i] = uniform_sample_unit_ball(N)
    end
    B
end


function region_high_probability(f::Array{Function,1},mu::Array{Cdouble,1},gamma::Symmetric,p::Cdouble=0.7)
    # get the number of dimensions
    N = length(f)
     # find the isocontour such that the probability of being in the enclosed volume is 0.7
    ra = 0.0
    while(multivariate_probability(N,ra)<p)
        ra = ra + 1.0
    end
    # refine the search by dicotomie
    minR = ra-1.0;
    maxR = ra;
    while((maxR-minR)>1.0e-10)
        if (multivariate_probability(N,0.5*(minR+maxR))>0.7)
            maxR = 0.5*(minR+maxR)
        else
            minR = 0.5*(minR+maxR)
        end
    end
    RA = 0.5*(minR+maxR)

    # new estimation
    S1 = RA*uniform_sample_unit_sphere(N,min(N*5000,500000)); # N=100 is the practical maximum for memory considerations

    # change of basis
    S = similar(S1)
    Fval,Fvect = eigen(gamma)
    P = Fvect*diagm(sqrt.(Fval))
    for i in 1:min(N*5000,500000)
        S[:,i] = P*S1[:,i]+mu
    end

    # change of variable
    X = similar(S)
    for i in 1:min(N*5000,500000)
        for n in 1:N
            X[n,i] = f[n](S[:,i])
        end
    end

    # return the levelset that define the border of the region
    X,S
end

# this function compute the corner' coordinates of the embedding N-box of the highest probability region
function percentile_estimation_ni(f::Array{Function,1},mu::Array{Cdouble,1},gamma::Symmetric,p::Cdouble=0.9) # better for non independent variables
    # get the number of dimensions
    N = length(f)

    # check a few dimensions
    if (length(mu)!=size(gamma,1))
        throw("In percentile_estimation_ni: the covariance and mean dimensions do not match!")
    end
    if (N!=length(mu))
        throw("In percentile_estimation_ni: the function only handles changes of variabele that are endomorphisms! Quick fix: increase the dimension of the smallest space and add trivial variable changes.")
    end

    # compute the contour of the region of probability p
    X,S = region_high_probability(f,mu,gamma,p)

    # new way to get the uncertainty... which does not account for the actual distribution, but I can't find anything better
    low_sig_new = dropdims(minimum(X,dims=2),dims=2)
    up_sig_new  = dropdims(maximum(X,dims=2),dims=2)
    [low_sig_new up_sig_new]
end
