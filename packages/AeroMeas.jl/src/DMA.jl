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



# Some DMA model are implemented in this file:
#   - Gaussian model


################################################################
#                    Numerical integration                     #
################################################################

# a function to compute the numerical integration of the function f on the interval [a,b] with n evenly spaced samples
function riemann(f::Function, a::Real, b::Real, n::Int; method::String="right") # this should be in the utilsFun.jl package
  if method == "right"
      # meth(f,l,r) = f(r) * (r-l)
      xs = a .+ collect(0.:n)*(b-a)/n
      as = [f(r)*(r-l) for (l,r) in zip(xs[1:end-1], xs[2:end])]
    elseif method == "left"
        # meth(f,l,r) = f(l) * (r-l)
        xs = a .+ collect(0.:n)*(b-a)/n
        as = [f(l)*(r-l) for (l,r) in zip(xs[1:end-1], xs[2:end])]
    elseif method == "trapezoid"
        # meth(f,l,r) = (1/2) * (f(l) + f(r)) * (r-l)
        xs = a .+ collect(0.:n)*(b-a)/n
        as = [(1.0/2.0)*(f(l) + f(r))*(r-l) for (l,r) in zip(xs[1:end-1], xs[2:end])]
    elseif method == "simpsons"
        # meth(f,l,r) = (1.0/6.0) * (f(l) + 4.0*(f((l+r)/2.0)) + f(r))*(r-l)
        xs = a .+ collect(0.:n)*(b-a)/n
        as = [(1.0/6.0) * (f(l) + 4.0*(f((l+r)/2.0)) + f(r))*(r-l) for (l,r) in zip(xs[1:end-1], xs[2:end])]
    else
        throw(@sprintf "quadrature %s is not implemented" method)
  end
  sum(as)
end

################################################################
#                    Efficiency functions                      #
################################################################
# the triangular function should be the used one!!!!!! and the bandwidth of the each channel is given by dZ = (q_a/q_sh)*Z
function ch_eff(x0::Cdouble,xl::Cdouble,xh::Cdouble,x::Cdouble;model::String="gauss") # this should become the default function used later on
    # x0 is "the center of the measurement": it's the size at which we would measure is everything was perfect #REMARK: is the measurement device was "perfect", we would have eff_ch = 1 if x=x0, and 0 everywhere else.
    # xl and xh are the low and high cutoff sizes: x0\in[xl,xh]
    # x is the size for which we want the efficiency
    val = 0.0;
    if ((xl<=x0) & (x0<=xh) & (xl<xh))
        if model=="triangle" # this one should be used for the channel efficiency in terms of mobility
            if ((x>=xl) & (x<x0)) # the increasing segment of the triangle
                val = (x-xl)/(x0-xl);
            elseif ((x>=x0) & (x<xh))
                val = (xh-x)/(xh-x0);
            else
                val = 0.0 # there should be no need for this line, but it is more balanced that way
            end
        elseif model=="gauss"
            val = exp(-0.5*((x-x0)/(0.5*(xh-xl)))^2) # here I assume that x0 = 0.5(xl+xh) and that the std is 0.5(xh-xl)... it might not be exactly true though
        elseif model=="gate"
            if ((x<xh) & (xl<=x))
                val = 1.0
            end
        else # set the default model to gaussian
            val = exp(-0.5*((x-x0)/(0.5*(xh-xl)))^2)
        end
    else
        throw("Error in DMA.jl at eff_ch function: the limits are not set properly")
    end
    val
end

function ch_eff(x0::Cdouble,xl::Cdouble,xh::Cdouble,x::Array{Cdouble,1};model::String="gauss") # this should become the default function used later on
    # x0 is "the center of the measurement": it's the size at which we would measure is everything was perfect #REMARK: is the measurement device was "perfect", we would have eff_ch = 1 if x=x0, and 0 everywhere else.
    # xl and xh are the low and high cutoff sizes: x0\in[xl,xh]
    # x is the size for which we want the efficiency
    val = zeros(Cdouble,length(x));
    if ((xl<=x0) & (x0<=xh) & (xl<xh))
        if model=="triangle" # this one should be used for the channel efficiency in terms of mobility
            idx_i = (x.>=xl) .& (x.<x0); # the increasing segment of the triangle
            idx_d = (x.>=x0) .& (x.<xh);
            val[idx_i] = (x[idx_i].-xl)./(x0-xl);
            val[idx_d] = (xh.-x[idx_d])./(xh-x0);
        elseif model=="gauss"
            val = exp.(-0.5*((x.-x0)/(0.5*(xh-xl))).^2) # here I assume that x0 = 0.5(xl+xh) and that the std is 0.5(xh-xl)... it might not be exactly true though
        elseif model=="gate"
            idx = find(q->( ( q>=xl ) & ( q<xh ) ), x);
            val[idx] = 1.0;
        else # set the default model to gaussian
            val = exp.(-0.5*((x.-x0)/(0.5*(xh-xl))).^2)
        end
    else
        throw("Error in DMA.jl at eff_ch function: the limits are not set properly")
    end
    val
end


################################################################
#                             DMA                              #
################################################################

# separation efficiency of the DMA, considering the discretization bin is in a logarithmic scale
function chanel_efficiency(bin_center::Cdouble,cst_ratio::Cdouble,particle_size::Cdouble)
    # bin_center is the geometrical mean of the bin
    # cst_ratio is the constant ratio between two consecutive bin centers
    # particle_size is the size at which the detection efficiency is evaluated
    # assume a perfect efficiency over the chanel, and perfect cut offs (no overlap between measurements)
    # val = 0.0
    # if ( ( particle_size>=(bin_center*(cst_ratio^(-0.5))) ) & ( particle_size<(bin_center*(cst_ratio^(0.5))) ) )
    #     val = 1.0
    # end
    # val
    ch_eff(bin_center,bin_center/sqrt(cst_ratio),bin_center*sqrt(cst_ratio),particle_size;model="gate")
end

function chanel_efficiency(bin_center::Cdouble,cst_ratio::Cdouble,particle_size::Array{Cdouble,1})
    # bin_center is the geometrical mean of the bin
    # cst_ratio is the constant ratio between two consecutive bin centers
    # particle_size is the size at which the detection efficiency is evaluated
    # assume a perfect efficiency over the chanel, and perfect cut offs (no overlap between measurements)
    # val = zeros(length(particle_size))
    # idx = find(q->( ( q>=(bin_center*(cst_ratio^(-0.5))) ) & ( q<(bin_center*(cst_ratio^(0.5))) ) ), particle_size)
    # val[idx] = 1.0
    # val
    ch_eff(bin_center,bin_center/sqrt(cst_ratio),bin_center*sqrt(cst_ratio),particle_size;model="gate")
end


function chanel_efficiency_gaus(bin_center::Cdouble,cst_ratio::Cdouble,particle_size::Cdouble)
    # bin_center is the geometrical mean of the bin
    # cst_ratio is the constant ratio between two consecutive bin centers
    # particle_size is the size at which the detection efficiency is evaluated
    # assume a gaussian efficiency with overlap
    # sigg = 0.5*bin_center*(cst_ratio^(0.5)-cst_ratio^(-0.5))
    # exp(-0.5*((particle_size-bin_center)/sigg)^2)
    ch_eff(bin_center,bin_center/sqrt(cst_ratio),bin_center*sqrt(cst_ratio),particle_size;model="gauss")
end

function chanel_efficiency_gaus(bin_center::Cdouble,cst_ratio::Cdouble,particle_size::Array{Cdouble,1})
    # bin_center is the geometrical mean of the bin
    # cst_ratio is the constant ratio between two consecutive bin centers
    # particle_size is the size at which the detection efficiency is evaluated
    # assume a gaussian efficiency with overlap
    # sigg = 0.5*bin_center*(cst_ratio^(0.5)-cst_ratio^(-0.5))
    # exp(-0.5*((particle_size-bin_center)/sigg).^2)
    ch_eff(bin_center,bin_center/sqrt(cst_ratio),bin_center*sqrt(cst_ratio),particle_size;model="gauss")
end




# compute a measurement operator assuming that either the density of the efficiency kernel (or both) are piecewise constant and that the output number concentrations are collected on a logarithmic scale
# NOTE: it could be used in the Kalman Filter if we used the density evolution model
function DMA_gaus_OP(psd::Array{Cdouble,1},x_psd::Array{Cdouble,1},delta_psd::Array{Cdouble,1},bin_center::Array{Cdouble,1},cst_r::Cdouble)
    measOP = zeros(Cdouble,length(bin_center),length(psd));
    if ((length(psd)==length(x_psd)) & (length(x_psd)==length(delta_psd)))
        for i in 1:length(bin_center)
            measOP[i,:] = ch_eff(bin_center[i],bin_center[i]/sqrt(cst_r),bin_center[i]*sqrt(cst_r),x_psd;model="gauss").*delta_psd
        end
    end
    measOP
end






# the DMA at one instant in time
function DMA_gaus(psd::Array{Cdouble,1},particle_sizes::Array{Cdouble,1},delta_sizes::Array{Cdouble,1},bin_centers::Array{Cdouble,1},cst_ratio::Cdouble)
    # psd:            the particle size density
    # particle_sizes: the particle size at which the psd is sampled
    # delta_sizes:    the width of intervals for each sampling point
    # bin_centers:    the center of each chanels of the DMA
    # cst_ratio:      the constant ratio between two consecutive bin center

    # check the size of the 3 first argument: if they are not all of the same size, it is not possible to compute the number concentrations
    conc = similar(bin_centers)
    if ((length(psd)==length(particle_sizes)) & (length(particle_sizes)==length(delta_sizes)))
        for k in collect(1:length(bin_centers))
            conc[k] = sum(psd.*chanel_efficiency_gaus(bin_centers[k],cst_ratio,particle_sizes).*delta_sizes) # Reimann integration is good enough for this case, but it could be improved
        end
    else
        throw("The size of the arguments of the DMA function mismatches")
    end
    conc
end

function DMA_gate(psd::Array{Cdouble,1},particle_sizes::Array{Cdouble,1},delta_sizes::Array{Cdouble,1},bin_centers::Array{Cdouble,1},cst_ratio::Cdouble)
    # psd:            the particle size density
    # particle_sizes: the particle size at which the psd is sampled
    # delta_sizes:    the width of intervals for each sampling point
    # bin_centers:    the center of each chanels of the DMA
    # cst_ratio:      the constant ratio between two consecutive bin center

    # check the size of the 3 first argument: if they are not all of the same size, it is not possible to compute the number concentrations
    conc = similar(bin_centers)
    if ((length(psd)==length(particle_sizes)) & (length(particle_sizes)==length(delta_sizes)))
        for k in collect(1:length(bin_centers))
            function ff(d_::Cdouble)
                kk = findfirst(q->( (q<=(d_*sqrt(cst_ratio))) & (d_<(q*sqrt(cst_ratio)))  ),particle_sizes)
                psd[kk]*chanel_efficiency(bin_centers[k],cst_ratio,d_)
            end
            conc[k] = riemann(ff, bin_centers[k]/sqrt(cst_ratio), bin_centers[k]*sqrt(cst_ratio), 50, method="simpsons")#"trapezoid") #
            # conc[k] = sum(psd.*chanel_efficiency(bin_centers[k],cst_ratio,particle_sizes).*delta_sizes) # this is not precise enough: it induces numerical artifacts due to the integration method, not due to the model
        end
    else
        throw("The size of the arguments of the DMA function mismatches")
    end
    conc
end


function DMA_time_avg_gaus(PSD::Array{Cdouble,2},particle_sizes::Array{Cdouble,1},delta_sizes::Array{Cdouble,1},bin_centers::Array{Cdouble,1},cst_ratio::Cdouble)
    # psd:            the particle size density at different time
    # particle_sizes: the particle size at which the psd is sampled
    # delta_sizes:    the width of intervals for each sampling point
    # bin_centers:    the center of each chanels of the DMA
    # cst_ratio:      the constant ratio between two consecutive bin center
    nbin,Nt = size(PSD)
    DMA_gaus(dropdims(sum(PSD,dims=2),dims=2)/Nt,particle_sizes,delta_sizes,bin_centers,cst_ratio)
end


function DMA_time_avg_gate(PSD::Array{Cdouble,2},particle_sizes::Array{Cdouble,1},delta_sizes::Array{Cdouble,1},bin_centers::Array{Cdouble,1},cst_ratio::Cdouble)
    # psd:            the particle size density at different time
    # particle_sizes: the particle size at which the psd is sampled
    # delta_sizes:    the width of intervals for each sampling point
    # bin_centers:    the center of each chanels of the DMA
    # cst_ratio:      the constant ratio between two consecutive bin center
    nbin,Nt = size(PSD)
    DMA_gate(dropdims(sum(PSD,dims=2),dims=2)/Nt,particle_sizes,delta_sizes,bin_centers,cst_ratio)
end
