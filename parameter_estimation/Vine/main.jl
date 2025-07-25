# This file is licensed under the MIT "Expat" License:

# Copyright (c) 2020: Matthew Ozon.

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



using PyPlot
# set matplotlib to allow fior the use of LaTeX formula in the graphs (side effect: the font is set to LaTeX's default)
rc("text", usetex=true)
using myPlot
using CSV
using DataFrames
using Printf
using LinearAlgebra
using Statistics
using XLSX

# NOTE: you need to run the data_simulation/nucleation_event/main.jl script to simulate the Ground Truth

##################################################
###    read data and init some variabeles      ###
##################################################
COAGULATION = true;
COAGULATION_GAIN = (true & COAGULATION)

# # loading the ground truth or not
# GT_loaded = false #WARNING need more memory for this!!!


# data with a different measurement model
input_folder = "../../data_vine/"
folder = "results"


# load data 
xf_exp_data = XLSX.readxlsx(joinpath(input_folder,"concentration_data_vine.xlsx")); 


df_data = DataFrame(XLSX.gettable(xf_exp_data["concentration"]))
raw_data = Matrix{Cdouble}(Matrix{Cdouble}(df_data)');

df_diameter = DataFrame(XLSX.gettable(xf_exp_data["centroid_diamter"]))
diameter_data = Vector{Cdouble}(df_diameter.values)
volume_data = (pi/6.0)*(diameter_data.^3)
cst_v_data = mean(volume_data[2:end]./volume_data[1:end-1])

df_time = DataFrame(XLSX.gettable(xf_exp_data["measurement_time"]))
t_samp = Vector{Cdouble}(df_time.values)

min_pc,max_pc = extrema(raw_data)
s = @sprintf "number concentration (%1.2e,%1.2e)" min_pc max_pc
displayLogData2D(1,t_samp/3600.0,diameter_data,raw_data,max(min_pc,0.001max_pc),max_pc,_title=s,_colorbar_label="concentration [cm\$^{-3}\$]")
tight_layout()

