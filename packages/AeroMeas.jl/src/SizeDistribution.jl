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



# some measurement devices are implemented in this file, such as:
#   - DMA
#   - CPC (the random number generator depends on the pakage Distributions.jl which can be installed via Pkg.add("Distributions"))
#   - DMPS

# load the DMA and CPC models
include("DMA.jl")
include("CPC.jl")

################################################################
#                            DMPS                              #
################################################################

# for a time invariant particle size density, it simulates a measurement made with a DMPS
function DMPS_gaus(psd::Array{Cdouble,1},particle_sizes::Array{Cdouble,1},delta_sizes::Array{Cdouble,1},bin_centers::Array{Cdouble,1},cst_ratio::Cdouble,volume::Cdouble)
    # psd:            the particle size density
    # particle_sizes: the size of particle corresponding to the PSD
    # delta_sizes:    the interval length over which the PSD is constant around each central size
    # bin_centers:    the centers of each bin of the DMPS
    # cst_ratio:      the constant ratio between two consecutive bins
    # volume:         the volume of the sample used for counting in the DMPS
    CPC(DMA_gaus(psd,particle_sizes,delta_sizes,bin_centers,cst_ratio),volume)
end

function DMPS_gate(psd::Array{Cdouble,1},particle_sizes::Array{Cdouble,1},delta_sizes::Array{Cdouble,1},bin_centers::Array{Cdouble,1},cst_ratio::Cdouble,volume::Cdouble)
    # psd:            the particle size density
    # particle_sizes: the size of particle corresponding to the PSD
    # delta_sizes:    the interval length over which the PSD is constant around each central size
    # bin_centers:    the centers of each bin of the DMPS
    # cst_ratio:      the constant ratio between two consecutive bins
    # volume:         the volume of the sample used for counting in the DMPS
    CPC(DMA_gate(psd,particle_sizes,delta_sizes,bin_centers,cst_ratio),volume)
end

# for a time varying particle size density, it simulates a measurement made with a DMPS
function DMPS_gaus(PSD::Array{Cdouble,2},particle_sizes::Array{Cdouble,1},delta_sizes::Array{Cdouble,1},bin_centers::Array{Cdouble,1},cst_ratio::Cdouble,volume::Cdouble)
    # psd:            the particle size density
    # particle_sizes: the size of particle corresponding to the PSD
    # delta_sizes:    the interval length over which the PSD is constant around each central size
    # bin_centers:    the centers of each bin of the DMPS
    # cst_ratio:      the constant ratio between two consecutive bins
    # volume:         the volume of the sample used for counting in the DMPS
    CPC(DMA_time_avg_gaus(PSD,particle_sizes,delta_sizes,bin_centers,cst_ratio),volume)
end

function DMPS_gate(PSD::Array{Cdouble,2},particle_sizes::Array{Cdouble,1},delta_sizes::Array{Cdouble,1},bin_centers::Array{Cdouble,1},cst_ratio::Cdouble,volume::Cdouble)
    # psd:            the particle size density
    # particle_sizes: the size of particle corresponding to the PSD
    # delta_sizes:    the interval length over which the PSD is constant around each central size
    # bin_centers:    the centers of each bin of the DMPS
    # cst_ratio:      the constant ratio between two consecutive bins
    # volume:         the volume of the sample used for counting in the DMPS
    CPC(DMA_time_avg_gate(PSD,particle_sizes,delta_sizes,bin_centers,cst_ratio),volume)
end


