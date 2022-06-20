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

# some constants for the calculation of the physical parameters
g    = 9.81                      # [m s-2], gravitational acceleration
Rg   = 8.3144598                 # [J mol-1 K-1], universal gas constant
Na   = 6.022140857e23            # [molec mol-1], Avogadro's number
Mair = 28.96e-3                  # [kg mol-1], mean molar mass of air
kb   = 1.38064852e-23            # [m2 kg s-2 K-1], Boltzmann constant
Cp   = 1012.0                    # [J kg-1 K-1], air specific heat at constant pressure

#WARNING: the coagulation coefficients are computed in SI units

# coagulation: OK
function coagulation!(dx_coag::Array{Cdouble,1},ws::AeroSys,x::Array{Cdouble})
    # loss
    # for i in 1:ws.nbin # loss due to collision with all particles
    #     dx_coag[i] = -x[i]*sum(ws.beta[i,1:end].*x)
    # end
    for i in 1:ws.nbin-1 # only loss from collision with bigger particles
        dx_coag[i] = -x[i]*sum(ws.beta[i,i+1:end].*x[i+1:end])
    end
    # gain
    for i in 1:ws.count_coag_index
        dx_coag[ws.Ic[1,i]] = dx_coag[ws.Ic[1,i]] + 0.5*x[ws.Ic[2,i]]*ws.beta[ws.Ic[2,i],ws.Ic[3,i]]*x[ws.Ic[3,i]]*ws.a_gain[i]
    end
    dx_coag[end] = 0.0
    dx_coag[:] = (ws.x0*ws.beta0*ws.t0)*dx_coag
end


function coagulation_loss!(dx_coag::Array{Cdouble,1},ws::AeroSys,x::Array{Cdouble})
    # loss
    for i in 1:ws.nbin-1 # only loss from collision with bigger particles
        dx_coag[i] = -x[i]*sum(ws.beta[i,i+1:end].*x[i+1:end])
    end
    dx_coag[end] = 0.0
    dx_coag[:] = (ws.x0*ws.beta0*ws.t0)*dx_coag
end



# jacobian of the coagulation term
function jacobian_coagulation!(F_ev::Array{Cdouble,2},ws::AeroSys,x::Array{Cdouble,1})
    # loss !nb op: 2*nbin^2
    # for s in 1:ws.nbin
    #     for k in 1:ws.nbin
    #         if s==k
    #             F_ev[s,k] = F_ev[s,k] - ws.beta[s,:]*x - ws.beta[s,k]*x[s]
    #         else
    #             F_ev[s,k] = F_ev[s,k] - ws.beta[s,k]*x[s]
    #         end
    #     end
    # end
    for s in 1:ws.nbin-1
        for k in s+1:ws.nbin
            F_ev[s,k] = F_ev[s,k] - (ws.beta[s,k]*x[s])[1]
        end
    end
    # gain
    for i in 1:ws.count_coag_indexj
       F_ev[ws.Icj[1,i],ws.Icj[2,i]] = F_ev[ws.Icj[1,i],ws.Icj[2,i]] + ws.a_gainj[i]*ws.beta[ws.Icj[2,i],ws.Icj[3,i]]*x[ws.Icj[3,i]]
    end
    F_ev[:,:] = (ws.x0*ws.beta0*ws.t0)*F_ev
end


# jacobian of the coagulation term
function jacobian_coagulation_loss!(F_ev::Array{Cdouble,2},ws::AeroSys,x::Array{Cdouble,1})
    # loss !nb op: 2*nbin^2
    for s in 1:ws.nbin-1
        for k in s+1:ws.nbin
            F_ev[s,k] = F_ev[s,k] - (ws.beta[s,k]*x[s])[1]
        end
    end
    F_ev[:,:] = (ws.x0*ws.beta0*ws.t0)*F_ev
end





# coagulation coefficient #WARNING: the coagulation coefficients are computed in SI units
function coagulation_coefficient!(ws::AeroSys,temperature::Cdouble=300.0, pressure::Cdouble=1.0e5,density::Cdouble=1400.0) # density [kg.m^-3]
    # Those parameters depend on the temperature, the pressure the molar mass of the air the diameter
    # and the particle mass which is assumed to be constant over time, at least for a few time steps.

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # The following parameters are well trusted, but could actually be part of the estimation
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # The Coagulation coefficient is calculated according to formula 13.56 in Seinfield and Pandis (2006), Page 603
    particle_mass = density*(pi/6.0)*(ws.d.^3)                                                   # [kg]
    dyn_visc = 1.8e-5*(temperature/298.)^0.85                                                    # [kg.m^{-1}.s^{-1}] Dynamic viscosity of air
    l_gas=2.0*dyn_visc/(pressure*sqrt(8.0*Mair/(pi*Rg*temperature)))                             # [m] Gas mean free path in air
    slip_correction = 1.0.+(2.0*l_gas./ws.d).*(1.257.+0.4*exp.(-1.1./(2.0*l_gas./ws.d)))            # [ ] Cunninghams slip correction factor (Seinfeld and Pandis eq 9.34)
    diffusivity = slip_correction*kb*temperature./(3.0*pi*dyn_visc*ws.d)                         # [m^2.s^{-1}] Diffusivity for the different particle sizes
    speed_p = sqrt.(8.0*kb*temperature./(pi*particle_mass))                                       # [m.s^{-1}] Speed of particles
    free_path_p = 8.0*diffusivity./(pi*speed_p)                                                  # [m] Particle mean free path (m)
    dist = (1.0./(3.0*ws.d.*free_path_p)).*((ws.d.+free_path_p).^3.0 .-(ws.d.^2.0.+free_path_p.^2.0).^(3.0/2.0)).-ws.d   # [m] mean distance from the center of a sphere reached by particles leaving the sphere's surface (m)

    # the coagulation coeficients
    beta_fs=Array{Cdouble,1}(undef,ws.nbin)                                                               # Fuchs-Sutugin correction coefficient
    for i in 1:ws.nbin
        beta_fs = 1.0./((ws.d.+ws.d[i])./(ws.d.+ws.d[i].+ 2.0*(dist.^2.0.+dist[i]^2.0).^0.5).+8.0*(diffusivity.+diffusivity[i])./(((speed_p.^2.0.+speed_p[i]^2.0).^0.5).*(ws.d.+ws.d[i]))) # Fuchs correction factor from Seinfeld and Pandis, 2006, p. 600
        ws.beta[i,:] = 2.0*pi*beta_fs.*(ws.d*diffusivity[i].+ws.d.*diffusivity.+ws.d[i]*diffusivity.+ws.d[i]*diffusivity[i]) # [m^3.s^{-1}] coagulation rates between two particles of all size combinations
    end
    ws.beta0 = mean(ws.beta) # [# m^3.s^{-1}] the mean colission frequency in [# m^3.s^{-1}] is given by 1.0e6*ws.beta0
    ws.beta[:,:] = (1.0/ws.beta0)*ws.beta # []
end





# shortcut for coagulation evolution computation
function init_coagulation_loop_indices(ws::AeroSys)
    volume_ = (pi/6.0)*(ws.d.^3)

    # indices pairs for coagulation
    count_coag_index_ = 0
    for s in 2:ws.nbin
        for i in 1:s
            for j in 1:s
                # volume of the new formed particle
                v_new = volume_[i]+volume_[j]

                # if the newly formed particle is in the good range, compute its contribution to the current bin (s)
                if ((v_new>=(volume_[s]/sqrt(ws.cst_v))) & (v_new<(volume_[s]*sqrt(ws.cst_v))))
                    count_coag_index_ = count_coag_index_ + 1
                end
                # if (ws.cst_v*v_new>volume_[s]) & (v_new<volume_[s]*ws.cst_v)
                #     count_coag_index_ = count_coag_index_ + 1
                # end
            end
        end
    end

    Ic_ = Array{Int64,2}(undef,3,count_coag_index_)
    a_gain_ = Array{Cdouble,1}(undef,count_coag_index_)
    k = 1
    for s in 2:ws.nbin
        for i in 1:s
            for j in 1:s
                # volume of the new formed particle
                v_new = volume_[i]+volume_[j]

                if ((v_new>=(volume_[s]/sqrt(ws.cst_v))) & (v_new<(volume_[s]*sqrt(ws.cst_v))))
                    Ic_[1,k] = s
                    Ic_[2,k] = i
                    Ic_[3,k] = j
                    if ((v_new>=(volume_[s]/sqrt(ws.cst_v))) & (v_new<=volume_[s]))
                        a_gain_[k] = 0.5*((volume_[s]/sqrt(ws.cst_v)) - v_new)/((volume_[s]/sqrt(ws.cst_v)) - volume_[s]) + 0.5
                    else
                        a_gain_[k] = 0.5 - 0.5*(v_new-(volume_[s]*sqrt(ws.cst_v)))/(volume_[s]-(volume_[s]*sqrt(ws.cst_v)))
                    end
                    a_gain_[k] = max(a_gain_[k],0.0)
                    k = k + 1
                end

                # if the newly formed particle is in the good range, compute its contribution to the current bin (s)
                # if (ws.cst_v*v_new>volume_[s]) & (v_new<volume_[s]*ws.cst_v)
                #     # index mapping
                #     Ic_[1,k] = s
                #     Ic_[2,k] = i
                #     Ic_[3,k] = j
                #     # relative contribution
                #     a_gain_[k] = max(0.0,(1.0/(ws.cst_v-1.0))*min((ws.cst_v*v_new/volume_[s])-1.0,ws.cst_v - (v_new/volume_[s])))
                #     k = k + 1
                # end
            end
        end
    end

    # idem for jacobian: coagulation gain
    count_coag_indexj_ = 0
    # count the number of indices
    for s in 1:ws.nbin
        for k in 1:ws.nbin
          # gain
            for i in 1:ws.nbin
                v_new = volume_[i]+volume_[k]
                if ((v_new>=(volume_[s]/sqrt(ws.cst_v))) & (v_new<(volume_[s]*sqrt(ws.cst_v))))
                    count_coag_indexj_ = count_coag_indexj_ + 1
                end
                # if  (v_new*ws.cst_v>volume_[s]) & (v_new<ws.cst_v*volume_[s])
                #     count_coag_indexj_ = count_coag_indexj_ + 1
                # end
            end
        end
    end
    # allocate
    Icj_ = Array{Int64,2}(undef,3,count_coag_indexj_)
    a_gainj_ = Array{Cdouble,1}(undef,count_coag_indexj_)
    # set the index list and the relative contribution
    j = 1
    for s in 1:ws.nbin
        for k in 1:ws.nbin
            # gain
            for i in 1:ws.nbin
                v_new = volume_[i]+volume_[k]
                if ((v_new>=(volume_[s]/sqrt(ws.cst_v))) & (v_new<(volume_[s]*sqrt(ws.cst_v))))
                    Icj_[1,j] = s
                    Icj_[2,j] = k # i
                    Icj_[3,j] = i # j
                    if ((v_new>=(volume_[s]/sqrt(ws.cst_v))) & (v_new<=volume_[s]))
                        a_gainj_[j] = 0.5*((volume_[s]/sqrt(ws.cst_v)) - v_new)/((volume_[s]/sqrt(ws.cst_v)) - volume_[s]) + 0.5
                    else
                        a_gainj_[j] = 0.5 - 0.5*(v_new-(volume_[s]*sqrt(ws.cst_v)))/(volume_[s]-(volume_[s]*sqrt(ws.cst_v)))
                    end
                    a_gainj_[j] = max(a_gainj_[j],0.0)
                    j = j + 1
                end
                # if (v_new*ws.cst_v>volume_[s]) & (v_new<ws.cst_v*volume_[s])
                #     # indices
                #     Icj_[1,j] = s
                #     Icj_[2,j] = k
                #     Icj_[3,j] = i
                #     # relative contribution
                #     a_gainj_[j] = max(0.0,(1.0/(ws.cst_v-1.0))*min((ws.cst_v*v_new/volume_[s])-1.0,ws.cst_v - (v_new/volume_[s])))
                #     # update iterator
                #     j = j + 1
                # end
            end
        end
    end
    count_coag_index_,Ic_,a_gain_,count_coag_indexj_,Icj_,a_gainj_
end


function init_coagulation_loop_indices!(ws::AeroSys)
    ws.count_coag_index,ws.Ic,ws.a_gain,ws.count_coag_indexj,ws.Icj,ws.a_gainj = init_coagulation_loop_indices(ws)
end
