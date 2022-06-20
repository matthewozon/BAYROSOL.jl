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



# Some CPC model are implemented in this file:
#   - Gaussian model
#   - Poisson model (the random number generator depends on the pakage Distributions.jl which can be installed via Pkg.add("Distributions"))
# in any case, the input of the devcie's model is the number concentration at the inlet and the volume that is used to do the counting.
# Beware! This volume may be quite different from the volume taken from the sample to measure. In my opinion, the volume maybe considered
# as an effective volume: the volume from which particles have actually been counted.


################################################################
#                             CPC                              #
################################################################

function CPC(conc::Cdouble,volume::Cdouble)
    # conc:   particle concentration at the inlet of the CPC
    # volume: volume of the sample used for counting in the CPC
    rand(Poisson(conc*volume))/volume
end

function CPC(conc::Array{Cdouble,1},volume::Cdouble)
    # conc:   particle concentration at the inlet of the CPC
    # volume: volume of the sample used for counting in the CPC
    conc_poi = similar(conc)
    for k in 1:length(conc)
        conc_poi[k] = rand(Poisson(conc[k]*volume))/volume
    end
    conc_poi
end

# some gaussian approximation of the counting process: valid if the number of particle conc*volume in the sample in the CPC is larger than 20
function CPC_gaus(conc::Cdouble,volume::Cdouble)
    # conc:   particle concentration at the inlet of the CPC
    # volume: volume of the sample used for counting in the CPC
    abs(conc + sqrt(conc/volume)*randn())
end

function CPC_gaus(conc::Array{Cdouble,1},volume::Cdouble)
    # conc:   particle concentration at the inlet of the CPC
    # volume: volume of the sample used for counting in the CPC
    conc_gaus = similar(conc)
    for k in 1:length(conc)
        conc_gaus[k] = abs(conc[k] + sqrt(conc[k]/volume)*randn())
    end
    conc_gaus
end
