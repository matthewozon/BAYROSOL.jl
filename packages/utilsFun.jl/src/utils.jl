#------------------------------------------------------------------------------
#
# This file is part of the utilsFun module which is licensed under the MIT "Expat" License:
#
# Copyright (C) 2019-2020,  Matthew Ozon.
#
#------------------------------------------------------------------------------



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# softmax (actually softPlus)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function softMax(x_::Cdouble,eps_::Cdouble=0.0,thX::Cdouble=33.275)
    if x_<-thX
        y_ = eps_*x_
    elseif x_>thX
        y_ = x_
    else
        y_ = (1.0-eps_)*log(1.0+exp(x_)) + eps_*x_
    end
    y_
end

function softMax(x_::Array{Cdouble,1},eps_::Cdouble=0.0,thX::Cdouble=33.275)
    idx1 = findall(q->q<-thX,x_)
    idx2 = findall(q->q>thX,x_)
    idx3 = findall(q->((q>=-thX) & (q<=thX)),x_)
    y_ = similar(x_)
    y_[idx1] = eps_*x_[idx1]
    y_[idx2] = x_[idx2]
    y_[idx3] = (1.0-eps_)*log.(1.0.+exp.(x_[idx3])) + eps_*x_[idx3]
    y_
end

function softMaxDeriv(x::Cdouble,eps_::Cdouble=0.0,thX::Cdouble=33.275)
    if x<-thX
        y = eps_
    elseif x>thX
        y = 1.0
    else
        y = (1.0-eps_)/(1.0+exp(-x)) + eps_
    end
    y
end

function softMaxDeriv(x_::Array{Cdouble,1},eps_::Cdouble=0.0,thX::Cdouble=33.275)
    idx1 = findall(q->q<-thX,x_)
    idx2 = findall(q->q>thX,x_)
    idx3 = findall(q->((q>=-thX) & (q<=thX)),x_)
    y_ = similar(x_)
    y_[idx1] .= eps_
    y_[idx2] .= 1.0
    y_[idx3] = (1.0-eps_)./(1.0.+exp.(-x_[idx3])) .+ eps_
    y_
end



function softMaxA(x_::Cdouble,alpha_::Cdouble=10.0,eps_::Cdouble=0.0,thX::Cdouble=33.275)
    (1.0/alpha_)*softMax(alpha_*x_,eps_,thX)
end

function softMaxA(x_::Array{Cdouble,1},alpha_::Cdouble=10.0,eps_::Cdouble=0.0,thX::Cdouble=33.275)
    (1.0/alpha_).*softMax(alpha_.*x_,eps_,thX)
end

function softMaxDerivA(x::Cdouble,alpha_::Cdouble=10.0,eps_::Cdouble=0.0,thX::Cdouble=33.275)
     softMaxDeriv(alpha_*x,eps_,thX)
end

function softMaxDerivA(x_::Array{Cdouble,1},alpha_::Cdouble=10.0,eps_::Cdouble=0.0,thX::Cdouble=33.275)
    softMaxDeriv(alpha_.*x_,eps_,thX)
end







# reciprocal function
function softMaxInv(x::Cdouble,eps_::Cdouble=1.0e-15)
    if x<0.0
        throw(DomainError())
    end
    if x<eps_
        y = log(x)
    elseif x>40.0
        y = x
    else
        y = log(exp(x)-1.0)
    end
    y
end


function softMaxInv(x::Array{Cdouble,1},eps_::Cdouble=1.0e-15)
    if any(x.<0.0)
        throw(DomainError())
    end
    idx1 = findall(q->q<eps_,x)
    idx2 = findall(q->q>40.0,x)
    idx3 = findall(q->((q>=eps_) & (q<=40.0)),x)
    y = similar(x)
    y[idx1] = log.(x[idx1])
    y[idx2] = x[idx2]
    y[idx3] = log.(exp.(x[idx3]).-1.0)
    y
end


function softMaxInvDeriv(x::Cdouble,eps_::Cdouble=1.0e-15)
    if x<0.0
        throw(DomainError())
    end
    if x<eps_
        y = 1.0/x
    elseif x>40.0
        y = 1.0
    else
        y = 1.0/(1.0-exp(x))
    end
    y
end


function softMaxInvDeriv(x::Array{Cdouble,1},eps_::Cdouble=1.0e-15)
    if any(x.<0.0)
        throw(DomainError())
    end
    idx1 = findall(q->q<eps_,x)
    idx2 = findall(q->q>40.0,x)
    idx3 = findall(q->((q>=eps_) & (q<=40.0)),x)
    y = similar(x)
    y[idx1] = 1.0./x[idx1]
    y[idx2] .= 1.0
    y[idx3] = 1.0./(1.0.-exp.(-x[idx3]))
    y
end




function softMaxInvA(x::Cdouble,alpha_::Cdouble=10.0,eps_::Cdouble=1.0e-15)
    (1.0/alpha_)*softMaxInv(alpha_*x,eps_)
end


function softMaxInvA(x::Array{Cdouble,1},alpha_::Cdouble=10.0,eps_::Cdouble=1.0e-15)
    (1.0/alpha_).*softMaxInv(alpha_.*x,eps_)
end


function softMaxInvDerivA(x::Cdouble,alpha_::Cdouble=10.0,eps_::Cdouble=1.0e-15)
    softMaxInvDeriv(alpha_*x,eps_)
end


function softMaxInvDerivA(x::Array{Cdouble,1},alpha_::Cdouble=10.0,eps_::Cdouble=1.0e-15)
    softMaxInvDeriv(alpha_.*x,eps_)
end



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# logistic function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function logistic(x::Cdouble,thX::Cdouble=33.275)
    if x<-thX
        y = 0.0
    elseif x>thX
        y = 1.0
    else
        y = 1.0/(1.0+exp(-x))
    end
    y
end

function logistic(x_::Array{Cdouble,1},thX::Cdouble=33.275)
    idx1 = findall(q->q<-thX,x_)
    idx2 = findall(q->q>thX,x_)
    idx3 = findall(q->((q>=-thX) & (q<=thX)),x_)
    y_ = similar(x_)
    y_[idx1] .= 0.0
    y_[idx2] .= 1.0
    y_[idx3] = 1.0./(1.0.+exp.(-x_[idx3]))
    y_
end

function logistic(x_::Cdouble,a_::Cdouble,b_::Cdouble,alpha_::Cdouble=1.0,thX::Cdouble=33.275)
    if (alpha_*x_)<-thX
        y = a_
    elseif (alpha_*x_)>thX
        y = b_
    else
        y = a_ + (b_-a_)/(1.0+(1.0/alpha_)*exp(-alpha_*x_))
    end
    y
end

function logistic(x_::Array{Cdouble,1},a_::Cdouble,b_::Cdouble,alpha_::Cdouble=1.0,thX::Cdouble=33.275)
    idx1 = findall(q->q<-thX,alpha_*x_)
    idx2 = findall(q->q>thX,alpha_*x_)
    idx3 = findall(q->((q>=-thX) & (q<=thX)),alpha_*x_)
    y_ = similar(x_)
    y_[idx1] .= a_
    y_[idx2] .= b_
    y_[idx3] = a_ .+ (b_-a_)./(1.0.+(1.0/alpha_)*exp.(-alpha_*x_[idx3]))
    y_
end






function logisticInv(x::Cdouble)
    if ((x<0.0) | (x>1.0))
        throw(DomainError())
    end
    log(x/(1.0-x))
end

function logisticInv(x_::Array{Cdouble,1})
    if ( any(x_.<0.0) | any(x_.>1.0))
        throw(DomainError())
    end
    log.(x_./(1.0.-x_))
end


function logisticInv(x::Cdouble,a_::Cdouble,b_::Cdouble,alpha_::Cdouble=1.0)
    if ((x<a_) | (x>b_))
        throw(DomainError())
    end
    (1.0/alpha_)*log((1.0/alpha_)*(x-a_)/(b_-x))
end

function logisticInv(x_::Array{Cdouble,1},a_::Cdouble,b_::Cdouble,alpha_::Cdouble=1.0)
    if ( any(x_.<a_) | any(x_.>b_))
        throw(DomainError())
    end
    (1.0/alpha_)*log.((1.0/alpha_)*(x_.-a_)./(b_.-x_))
end














function logisticDeriv(x::Cdouble,thX::Cdouble=33.275)
    if x<-thX
        y = 0.0
    elseif x>thX
        y = 0.0
    else
        y = exp(-x)/((1.0+exp(-x))^2)
    end
    y
end

function logisticDeriv(x_::Array{Cdouble,1},thX::Cdouble=33.275)
    idx1 = findall(q->q<-thX,x_)
    idx2 = findall(q->q>thX,x_)
    idx3 = findall(q->((q>=-thX) & (q<=thX)),x_)
    y_ = similar(x_)
    y_[idx1] .= 0.0
    y_[idx2] .= 0.0
    y_[idx3] = exp.(-x_[idx3])./((1.0.+exp.(-x_[idx3])).^2)
    y_
end


function logisticDeriv(x::Cdouble,a_::Cdouble,b_::Cdouble,alpha_::Cdouble=10.0,thX::Cdouble=33.275)
    if (alpha_*x)<-thX
        y = 0.0
    elseif (alpha_*x)>thX
        y = 0.0
    else
        y = (b_-a_)*exp(-alpha_*x)/((1.0+(1.0/alpha_)*exp(-alpha_*x))^2)
    end
    y
end

function logisticDeriv(x_::Array{Cdouble,1},a_::Cdouble,b_::Cdouble,alpha_::Cdouble=10.0,thX::Cdouble=33.275)
    idx1 = findall(q->q<-thX,alpha_*x_)
    idx2 = findall(q->q>thX,alpha_*x_)
    idx3 = findall(q->((q>=-thX) & (q<=thX)),x_)
    y_ = similar(x_)
    y_[idx1] .= 0.0
    y_[idx2] .= 0.0
    y_[idx3] = (b_-a_)*exp.(-alpha_*x_[idx3])./((1.0.+(1.0/alpha_)*exp.(-alpha_*x_[idx3])).^2)
    y_
end

# some function that can be used in the robust estimation framework: weighting and potental functions
function cauchy(x::Array{Cdouble,1},alpha_f::Cdouble;x0::Cdouble=0.0)
    ((x.-x0).^2)./(1.0 .+ ((x.-x0)/alpha_f).^2) # toto.^2 #
end
function cauchy_deriv(x::Array{Cdouble,1},alpha_f::Cdouble;x0::Cdouble=0.0)
    2.0*(x.-x0)./((1.0 .+ (x/alpha_f).^2).^2) # 2.0toto
end

function phi_hl(x::Array{Cdouble,1},alpha_f::Cdouble)
    log(1.0.+(x/alpha_f).^2)
end
function phi_hl_deriv(x::Array{Cdouble,1},alpha_f::Cdouble)
    2.0x./(alpha_f^2 .+ x.^2)
end 
