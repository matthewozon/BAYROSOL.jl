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
# This file (charge_probability.jl) contains the implementation of the conditional   #
# probability of having a number of particles q knowing the size of the particle d,  #
# P(q|d) while being in an ion gaz mixture of with the same positive and negative    #
# ion concentration with mobilities Z_p (positive) and Z_m (negative) which can be   #
# different.                                                                         #
# It also contains some functions that can be used across different measurement      #
# devices, e.g. C_slip and mobility_from_size_and_charge.                            #
# The code is based on the paprs Wiedensohler 1988 (J. Aerosol Sci., Vol. 19, No. 3, #
# p. 387-389, and Hoppel 1978 (J. Aerosl Sci. Vol. 9, p. 41-54), and other papers.   #
######################################################################################


# cf Table 1 in Wiedensohler 1988 (An approxiamtion of the bipolar charge distribution for particles in the submicron size range)
AiN = [-26.3328 -2.3197 -0.0003 -2.3484 -44.4756;
       35.9044   0.6175 -0.1014  0.6044  79.3772;
       -21.4608  0.6201  0.3073  0.4800 -62.8900;
       7.0867   -0.1105 -0.3372  0.0013  26.4492;
       -1.3088  -0.1260  0.1023 -0.1544 -5.7480;
       0.1051    0.0297 -0.0105 0.0320 0.5049];

s_th_up = 1.0e-6   # [m], upper limit for the validity of AiN coefficient (approximation coef for the less than two charges on a particle)
s_th_down = 1.0e-9 # [m], lower limit for the validity of AiN coefficient

# a few physical quantities
e_charge = 1.60217733e-19;  # [coulomb],      elementary charge
e_cgs    = 4.8e-10          # [esu],          elementary charge in the cgs system
eps_0    = 8.854187817e-12; # [farad m^{-1}], dielectric constant of the air
kb       = 1.380658e-23;    # [J K^{-1}],     Boltzman constant
kb_cgs   = 1.380658e-16;    # [ergs K^{-1}],  Boltzman constant in the cgs system
Z_p0     = 1.35e-4;         # [m V s^{-1}],   positive ion electrical mobility
Z_m0     = 1.6e-4;          # [m V s^{-1}],   negative ion electrical mobility


# S_sutherland = 110.5 # [K], Sutherland temperature
Pr0 = 1.0e5         # [Pa], reference pressure
Pr0_cgs = 0.1*Pr0   # [Ba], reference pressure in cgs unit
Tr = 273.0          # [K], reference temperature
Mair = 28.96e-3     # [kg mol-1], mean molar mass of air
Mair_cgs = 28.96    # [g mol^{-1}], mean molar mass of air in cgs unit
Rg   = 8.314        # [J mol-1 K-1], universal gas constant
Rg_cgs = 8.314e7    # [erg mol^{-1} K^{-1}], universal gas constant in cgs unit
alpha_cun = 1.142   # Cunningham Slip coefficients
beta_cun  = 0.558   #
gamma_cun = 0.999   #
alpha_mil = 0.864
beta_mil  = 0.290
gamma_mil = 1.250



# compute the conditional probability P(q|s) (probability of q elementary charges knowing the size of the particle) # it is modified version in order to fit the physical reality
# this function is valid for particle in the size range [1,1000] nm # don't try to push your luck out of this range, it is a polynomial interpolation (it goes rogue out of it's validity range)
function P_q_charges_knowing_size(q::Int64,s::Union{Cdouble,Array{Cdouble,1}};T::Cdouble=293.0,Z_p::Cdouble=Z_p0,Z_m::Cdouble=Z_m0)
    if (typeof(s)==Array{Cdouble,1})
        val = zeros(Cdouble,length(s));
    else
        val = 0.0;
    end

    if (abs(q)<=2)
        if (abs(q)<=1)
            # formula valid for the range [1,1000] nm
            for i in 1:length(s)
                if ((s[i]>=1.0e-9) & (s[i]<=s_th_up))
                    tmp = 0.0;
                    for k in 1:6
                        tmp = tmp + AiN[k,q+3]*((log10(1.0e9s[i]))^(k-1))
                    end
                    val[i] = 10.0^tmp
                else
                    # cf Table 2 in Wiedensohler 1988 (An approxiamtion of the bipolar charge distribution for particles in the submicron size range)
                    if (q==0)
                        if (s[i]<=1.0e-9)
                            val[i] = 0.9909;
                        else
                            val[i] = 0.1236;
                        end
                    elseif (q==1)
                        if (s[i]<=1.0e-9)
                            val[i] = 0.0044;
                        else
                            val[i] = 0.1024;
                        end
                    else
                        if (s[i]<=1.0e-9)
                            val[i] = 0.0047;
                        else
                            val[i] = 0.1333;
                        end
                    end
                end
            end
        else
            # formula valid for the range [20,1000] nm
            for i in 1:length(s)
                if ((s[i]>=20.0e-9) & (s[i]<=s_th_up))
                    tmp = 0.0;
                    for k in 1:6
                        tmp = tmp + AiN[k,q+3]*((log10(1.0e9s[i]))^(k-1))
                    end
                    val[i] = 10.0^tmp
                else
                    if (q==2)
                        if (s[i]<=20.0e-9)
                            val[i] = 0.0001;
                        else
                            val[i] = 0.0759;
                        end
                    else
                        if (s[i]<=20.0e-9)
                            val[i] = 0.0001;
                        else
                            val[i] = 0.1286;
                        end
                    end
                end
            end
        end
    else
        # Gunn 1956: assuming equal concentration of positive and negative ions
        tmp1 = e_charge./sqrt.(4.0pi^2*eps_0*s*kb*T)
        tmp_exp_num = (q.-2.0pi*eps_0*s*kb*T*log(Z_p/Z_m)/(e_charge^2)).^2;
        tmp_exp_den = 4.0pi*eps_0*s*kb*T/(e_charge^2)
        val = tmp1.*exp.(-tmp_exp_num./tmp_exp_den);
    end
    # # Gunn 1956: assuming equal concentration of positive and negative ions
    # tmp1 = e_charge./sqrt.(4.0pi^2*eps_0*2.5*s*kb*T)
    # tmp_exp_num = (q-2.0pi*eps_0*2.5*s*kb*T*log(Z_p/Z_m)/(e_charge^2)).^2;
    # tmp_exp_den = 4.0pi*eps_0*2.5*s*kb*T/(e_charge^2)
    # val = tmp1.*exp.(-tmp_exp_num./tmp_exp_den);
    val
end


# Cunnigham slip coefficient
function C_slip(Kn::Union{Cdouble,Array{Cdouble,1}};alpha::Cdouble=alpha_cun,beta::Cdouble=beta_cun,gamma::Cdouble=gamma_cun)
    1.0 .+ Kn.*(alpha.+beta*exp.(-gamma./Kn))
end

function C_slip_deriv(Kn::Union{Cdouble,Array{Cdouble,1}};alpha::Cdouble=alpha_cun,beta::Cdouble=beta_cun,gamma::Cdouble=gamma_cun)
    alpha .+ beta*(1.0.+(gamma./Kn)).*exp.(-gamma./Kn)
end

# electric mobility knowing the size and charge of the particle
function mobility_from_size_and_charge(s::Union{Cdouble,Array{Cdouble,1}},q::Int64;T::Cdouble=293.0,Pr::Cdouble=Pr0)
    dyn_visc = 1.8e-5*(T/298.)^0.85                       # [kg.m^{-1}.s^{-1}] Dynamic viscosity of air
    l_gas=2.0*dyn_visc/(Pr*sqrt(8.0*Mair/(pi*Rg*T)))      # [m] Gas mean free path in air
    Kn = 2.0*l_gas./s
    q*e_charge*C_slip(Kn)./(3.0pi*dyn_visc*s)
end

function mobility_from_size_and_charge_deriv(s::Union{Cdouble,Array{Cdouble,1}},q::Int64;T::Cdouble=293.0,Pr::Cdouble=Pr0)
    dyn_visc = 1.8e-5*(T/298.)^0.85                       # [kg.m^{-1}.s^{-1}] Dynamic viscosity of air
    l_gas=2.0*dyn_visc/(Pr*sqrt(8.0*Mair/(pi*Rg*T)))      # [m] Gas mean free path in air
    Kn = 2.0*l_gas./s
    Kn_deriv = 2.0*l_gas*(-1.0./(s.^2));
    (q*e_charge/(3.0pi*dyn_visc)).*(1.0./(s.^2)).*(s.*Kn_deriv.*C_slip_deriv(Kn) - C_slip(Kn));
end

function dk_from_ds_and_charge(ds::Union{Cdouble,Array{Cdouble,1}},s::Union{Cdouble,Array{Cdouble,1}},q::Int64;T::Cdouble=293.0,Pr::Cdouble=Pr0)
    ds.*mobility_from_size_and_charge_deriv(s,q;T=T,Pr=Pr)
end


# size knowing the mobility and charge of the particle
function size_from_mobility_and_charge(k::Union{Cdouble,Array{Cdouble,1}},q::Int64;T::Cdouble=293.0,Pr::Cdouble=Pr0,alpha::Cdouble=alpha_cun,beta::Cdouble=beta_cun,gamma::Cdouble=gamma_cun)
    dyn_visc = 1.8e-5*(T/298.)^0.85                       # [kg.m^{-1}.s^{-1}] Dynamic viscosity of air
    l_gas=2.0*dyn_visc/(Pr*sqrt(8.0*Mair/(pi*Rg*T)))      # [m] Gas mean free path in air

    # use the root finding algorithm (Newtons' method, if necessary, change to a minimization problem) (it should be quite ok because it is a strictly increasing function for which we are looking for a root)
    y = 2.0*l_gas*3.0*pi*dyn_visc*k/(q*e_charge);
    function f(x::Union{Cdouble,Array{Cdouble,1}})
        x.*C_slip(x;alpha=alpha,beta=beta,gamma=gamma).-y
    end
    function fd(x::Union{Cdouble,Array{Cdouble,1}}) # the derivative of f
        1.0 .+ 2.0alpha*x .+ beta*(1.0.+2.0x).*exp.(-gamma./x)
    end

    # Newton's method seems to work!
    Kn_est = ones(Cdouble,length(k))
    for i in 1:20
        Kn_est[:] = Kn_est[:] - f(Kn_est)./fd(Kn_est)
    end
    2.0*l_gas./Kn_est
end
