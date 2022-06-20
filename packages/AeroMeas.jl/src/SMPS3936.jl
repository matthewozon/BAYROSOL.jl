#------------------------------------------------------------------------------
#
# This file is part of the AeroMeas module which is licensed under the MIT "Expat" License:
#
# Copyright (C) 2021,  Matthew Ozon.
#
#------------------------------------------------------------------------------


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



######################################################################################
# this file (SMPS3936.jl) contains the implementation of the model of the SMPS-3936  #
# as described in the Instruction Manual provided by the manufacturer TSI            #
# Model 3936 Scanning Mobility Particle Sizer (SMPS) Spectrometer Instruction Manual #
######################################################################################

#####################################################################################
#        Let's try to model the lab's SMPS using the instruction manual             #
#####################################################################################

# cut off by inertial impactor
function impactor_eff(s::Union{Cdouble,Array{Cdouble,1}};s50::Cdouble=1.0e-6,delta50::Cdouble=0.1e-6)
    1.0./(1.0.+exp.((s.-s50)/delta50))
end

function impactor(u::Array{Cdouble,1},s::Array{Cdouble,1};s50::Cdouble=1.0e-6,delta50::Cdouble=0.1e-6)
    u.*impactor_eff(s;s50=s50,delta50=delta50)
end


# particle charging using a Kr85 neutralizer
function neutralizer_Kr_85(u_imp::Array{Cdouble,1},s::Array{Cdouble,1};T::Cdouble=293.0,Pr::Cdouble=Pr0,Nq::Int64=6)
    # create the set of mobilities coresponding to up to 6 elementary charges
    K = zeros(Cdouble,abs(Nq),length(s));
    if Nq>0
        R_charge = 1:Nq
    else
        R_charge = -1:-1:Nq
    end

    [K[abs(q),:]=mobility_from_size_and_charge(s,q;T=T,Pr=Pr) for q in R_charge]
    # probability for a particle of having q charges knowing its size
    Pqs = zeros(Cdouble,abs(Nq),length(s));
    [Pqs[abs(q),:]=P_q_charges_knowing_size(q,s,T=T) for q in R_charge]

    # create a set of particle mobility density
    u_q_k = zeros(Cdouble,abs(Nq),length(u_imp));
    [u_q_k[abs(q),:] = Pqs[abs(q),:].*u_imp for q in R_charge]
    # u_q_k[1,:] = u_imp; # get all particle to acquire one and only one elementary charge with 100% efficiency

    # return the estimated
    R_charge,K,u_q_k,Pqs
end

# size density through the analyser (at its outlet knowing the mobility densities at the inlet)
function DMA_size_density(u_q_k::Array{Cdouble,2},K::Array{Cdouble,2},s::Array{Cdouble,1},k_meas::Array{Cdouble,1},sig_k_meas::Array{Cdouble,1})
    u_dma_i_q_k = zeros(Cdouble,length(k_meas),size(K,1),size(K,2))
    Psi_dma = zeros(Cdouble,size(K,2));
    for i in 1:length(k_meas) # channel selection
        for q in 1:size(K,1)
            # Psi_dma = exp.(-400.0*((K[q,:].-k_meas[i])/k_meas[i]).^2) # set so that delt_k_i = (1/10)*k_i (flow ratio q_a/q_sh=1/10) q_a aerosol flow and q_sh sheath flow
            # Psi_dma = exp.(-0.5*((K[q,:].-k_meas[i])/sig_k_meas[i]).^2)
            idx_0 = (K[q,:].<(k_meas[i]-sig_k_meas[i])) .| (K[q,:].>(k_meas[i]+sig_k_meas[i]))
            idx_i = (K[q,:].>=(k_meas[i]-sig_k_meas[i])) .& (K[q,:].<k_meas[i]); # the increasing segment of the triangle
            idx_d = (K[q,:].>=k_meas[i]) .& (K[q,:].<(k_meas[i]+sig_k_meas[i]));
            Psi_dma[idx_0] .= 0.0;
            Psi_dma[idx_i] = (K[q,idx_i].-(k_meas[i]-sig_k_meas[i]))./sig_k_meas[i];
            Psi_dma[idx_d] = ((k_meas[i]+sig_k_meas[i]).-K[q,idx_d])./sig_k_meas[i];
            u_dma_i_q_k[i,q,:] = Psi_dma.*u_q_k[q,:] # this is for the set of mobility K[q,:]
        end
    end
    # now it's a bit more tricky because we have to make the ... wait a minute
    dropdims(sum(u_dma_i_q_k,dims=2),dims=2)
end

# complet model of the DMA from the inlet of the device to the outlet of the analyser and passing through the impactor and the neutralizer
function DMA_3010(u::Array{Cdouble,1},s::Array{Cdouble,1},k_meas::Array{Cdouble,1},sig_k_meas::Array{Cdouble,1};s50::Cdouble=1.0e-6,delta50::Cdouble=0.1e-6,T::Cdouble=293.0,Pr::Cdouble=Pr0,Nq::Int64=6)
    # cut off size by impactor
    u_imp = impactor(u,s;s50=s50,delta50=delta50)
    # charging the particle in the neutralizer
    R_charge,K,u_q_k = neutralizer_Kr_85(u_imp,s;T=T,Pr=Pr0,Nq=Nq)
    # mobility selection in the classifier: size density at the output of the mobility analyzer
    DMA_size_density(u_q_k,K,s,k_meas,sig_k_meas)
end



#####################################################################################
#                        Let's CPC ourselves out of this!                           #
#####################################################################################
function CPC_eff(s::Union{Cdouble,Array{Cdouble,1}}; s50::Cdouble=1.0e-8, delta50::Cdouble=1.0e-9)
    1.0./(1.0.+exp.(-(s.-s50)/delta50))
end


function CPC_density(u_i_s::Array{Cdouble,2},s::Array{Cdouble,1}; s50::Cdouble=1.0e-8, delta50::Cdouble=1.0e-9)
    u_cpc = similar(u_i_s);
    [u_cpc[i,:] = u_i_s[i,:].*CPC_eff(s;s50=s50,delta50=delta50) for i in 1:size(u_i_s,1)]
    u_cpc
end


# it is assumed, for now, that the density passed to this function is average over time or that it is time-constant over the integration period delta_t
function nb_count(u::Array{Cdouble,1},delta_s::Array{Cdouble,1},phi::Cdouble=1000.0,delta_t::Cdouble=30.0) # phi is the sample flow rate in [cm^{3} s^{-1}] assumed to be constant because it is in the instruction manual
    phi*delta_t*sum(u.*delta_s) # a simple centtral Reimann integration should be a good enough approximation in this case
end

#####################################################################################
#                           the complete SMPS 3936                                  #
#####################################################################################

function SMPS3936(u::Union{Array{Cdouble,1},Array{Cdouble,2}},s::Array{Cdouble,1},delta_s::Array{Cdouble,1},s_meas::Array{Cdouble,1},cst_r_meas::Cdouble,phi::Cdouble,dt::Cdouble;s50imp::Cdouble=1.0e-6,delta50imp::Cdouble=0.1e-6,s50cpc::Cdouble=1.0e-8, delta50cpc::Cdouble=1.0e-9,T::Cdouble=293.0,Pr::Cdouble=Pr0,Nq::Int64=-6,inv_data::Bool=true,Nq_inv::Int64=-3,q_a::Cdouble=0.3,q_sh::Cdouble=3.0)

    # create electric mobilities measured by the SMPS. Note that it is a good practice to have the range of sizes so that: s_meas \subset s
    if (Nq>0)
        k_meas = mobility_from_size_and_charge(s_meas,1;T=T,Pr=Pr)
    else
        k_meas = mobility_from_size_and_charge(s_meas,-1;T=T,Pr=Pr)
    end

    # mobility half bandwidth
    sig_k_meas = abs.(0.5*(q_a/q_sh)*k_meas);

    # create the time average of the density
    if (size(u,2)==1)
        u_time_avg = copy(u);
    else
        if ( (size(u,1)==length(s)) & (size(u,2)>1) )
            u_time_avg = dropdims(mean(u,dims=2),dims=2);
        else
            u_time_avg = dropdims(mean(u,dims=1),dims=1);
        end
    end

    # create the "monodisperse" size densities at the outlet of the DMA for each selected mobility
    u_dma = DMA_3010(u_time_avg,s,k_meas,sig_k_meas;s50=s50imp,delta50=delta50imp,T=T,Pr=Pr,Nq=Nq); # set the DMA so that it has a given size interval for each bin

    # for each channel, create the densities counted by the CPC
    u_cpc = CPC_density(u_dma,s; s50=s50cpc, delta50=delta50cpc)

    # count particles
    Ncount = similar(s_meas);
    Y = similar(s_meas);
    [Ncount[i] = nb_count(u_cpc[i,:],delta_s,phi,dt) for i in 1:size(u_cpc,1)]
    [Y[i] = (1.0/(phi*dt))*rand(Poisson(Ncount[i])) for i in 1:length(Ncount)]

    # if the user ask for data inversion
    if inv_data
        println("The data inversion is not implemented in this version: the simulated measurement are raw concentrations.")
        # # compute the mobility bin width
        # cst_r_meas = mean(s_meas[2:end]./s_meas[1:end-1]);
        # ds = s_meas*(sqrt(cst_r_meas)-1.0/sqrt(cst_r_meas)); # bin width in size
        # dk = abs.(dk_from_ds_and_charge(ds,s_meas,1;T=T,Pr=Pr)); # bin width in mobility
        # # dk = mobility_from_size_and_charge(s_meas/sqrt(cst_r_meas),1;T=T,Pr=Pr)-mobility_from_size_and_charge(s_meas*sqrt(cst_r_meas),1;T=T,Pr=Pr)
        # # compute the data inversion using Hoppel 1978 idea
        # Yinv = hoppel_1978_inversion(Y./dk,s_meas;T=T,Pr=Pr,Nq=Nq_inv,Niter=50)
        # Y = Yinv.*s_meas*(sqrt(cst_r_meas)-1.0/sqrt(cst_r_meas));
    end

    # return the corrupted number concentrations
    Y,u_dma,u_cpc
end


function SMPS3936_transfer_function(s::Array{Cdouble,1},s_meas::Array{Cdouble,1},cst_r_meas::Cdouble;s50imp::Cdouble=1.0e-6,delta50imp::Cdouble=0.1e-6,s50cpc::Cdouble=1.0e-8, delta50cpc::Cdouble=1.0e-9,T::Cdouble=293.0,Pr::Cdouble=Pr0,Nq::Int64=-6,q_a::Cdouble=0.3,q_sh::Cdouble=3.0)
    # create electric mobilities measured by the SMPS. Note that it is a good practice to have the range of sizes so that: s_meas \subset s
    if (Nq>0)
        k_meas = mobility_from_size_and_charge(s_meas,1;T=T,Pr=Pr)
    else
        k_meas = mobility_from_size_and_charge(s_meas,-1;T=T,Pr=Pr)
    end

    # mobility half bandwidth
    dk_meas = abs.(0.5*(q_a/q_sh)*k_meas);

    # create the "monodisperse" size densities at the outlet of the DMA for each selected mobility
    u_dma = DMA_3010(ones(Cdouble,length(s)),s,k_meas,dk_meas;s50=s50imp,delta50=delta50imp,T=T,Pr=Pr,Nq=Nq);

    # for each channel, create the densities counted by the CPC
    CPC_density(u_dma,s; s50=s50cpc, delta50=delta50cpc)
end
