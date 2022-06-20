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


############################################################
###                 display results                      ###
############################################################

x_fil_all = myWSKF.x_fil_all
x_pre_all = myWSKF.x_pre_all
o_fil_all = myWSKF.o_fil_all
o_pre_all = myWSKF.o_pre_all
gam_v_all = myWSKF.gam_v_all

idx_bin = [1; round(Int64,nbin/3); round(Int64,2nbin/3)]
# idx_time = [1; round(Int64,n_samp/2); n_samp]
idx_time = [1; round(Int64,n_samp/6); round(Int64,2n_samp/6); round(Int64,3n_samp/6); round(Int64,4n_samp/6); round(Int64,5n_samp/6); round(Int64,11n_samp/12); round(Int64,23n_samp/24); n_samp]



maxPSD = 1.0e-9maximum(y_all./delta)
minPSD = max(0.001maxPSD,1.0e-9minimum(y_all./delta))
maxPSDlog = maximum(y_all./log10(cst_r))
minPSDlog = max(0.001*maxPSDlog,minimum(y_all./log10(cst_r)))

titlePSD = @sprintf "particle size distribution (%1.2e,%1.2e)" 1.0e-9minimum(y_all./delta) 1.0e-9maximum(y_all./delta)
# cblabelPSD = "density [# cm\$^{-3}\$ nm\$^{-1}\$]"
cblabelPSD = "density [cm\$^{-3}\$ nm\$^{-1}\$]"
displayLogData2D(158,t_samp/3600.0,10e9diameter,1.0e-9y_all./delta,minPSD,maxPSD,_title=titlePSD,_colorbar_label=cblabelPSD)
savefig(string(folder,"data.png"))
savefig(string(folder,"data.pdf"))
if SAVE_PGF
    savefig(string(folder,"data.pgf"))
end
close("all")
titlePSDlog = @sprintf "log size distribution (%1.2e,%1.2e)" minimum(y_all./log10(cst_r)) maximum(y_all./log10(cst_r))
# cblabelPSDlog = "density [# cm\$^{-3}\$ log\$_{10}(\\frac{d_i}{d_{i-1}})^{-1}\$]"
cblabelPSDlog = "density [\$\\frac{\\mathrm{d} N}{\\mathrm{d} \\log_{10}(D_p)}\$]"
fig,ax,cbar = displayLogData2D(158,t_samp/3600.0,1.0e9diameter,y_all./log10(cst_r),minPSDlog,maxPSDlog,_title=titlePSDlog,_colorbar_label=cblabelPSDlog)
ylabel("diameter [nm]")
savefig(string(folder,"log_data.png"))
savefig(string(folder,"log_data.pdf"))
if SAVE_PGF
    savefig(string(folder,"log_data.pgf"))
end
close("all")


if GT_loaded
    maxPSD = maximum(PSD_simu)
    minPSD = max(0.001maxPSD,minimum(PSD_simu))
    maxPSDlog = maximum(1.0e9*PSD_simu.*delta_p/log10(cst_r_p))
    minPSDlog = max(0.001maxPSDlog,minimum(1.0e9PSD_simu.*delta_p/log10(cst_r_p)))
end

# displayPSD(159,t_samp/3600.0,diameter,x0*x_fil_all[R_psd,:])
titlePSD = @sprintf "particle size distribution (%1.2e,%1.2e)" 1.0e-9x0*minimum(x_fil_all[R_psd,:]./delta) 1.0e-9x0*maximum(x_fil_all[R_psd,:]./delta)
# cblabelPSD = "density [# cm\$^{-3}\$ nm\$^{-1}\$]"
cblabelPSD = "density [cm\$^{-3}\$ nm\$^{-1}\$]"
displayLogData2D(159,t_samp/3600.0,1.0e9diameter,1.0e-9x0*x_fil_all[R_psd,:]./delta,minPSD,maxPSD,_title=titlePSD,_colorbar_label=cblabelPSD)
ylabel("diameter [nm]")
savefig(string(folder,"reconstruction_fil.png"))
savefig(string(folder,"reconstruction_fil.pdf"))
if SAVE_PGF
    savefig(string(folder,"reconstruction_fil.pgf"))
end
close("all")

titlePSDlog = @sprintf "log size distribution (%1.2e,%1.2e)" x0*minimum(x_fil_all[R_psd,:]./log10(cst_r)) x0*maximum(x_fil_all[R_psd,:]./log10(cst_r))
# cblabelPSDlog = "density [# cm\$^{-3}\$ log\$_{10}(\\frac{d_i}{d_{i-1}})^{-1}\$]"
cblabelPSDlog = "density [\$\\frac{\\mathrm{d} N}{\\mathrm{d} \\log_{10}(D_p)}\$]"
displayLogData2D(159,t_samp/3600.0,1.0e9diameter,x0*x_fil_all[R_psd,:]./log10(cst_r),minPSDlog,maxPSDlog,_title=titlePSDlog,_colorbar_label=cblabelPSDlog)
ylabel("diameter [nm]")
savefig(string(folder,"log_reconstruction_fil.png"))
savefig(string(folder,"log_reconstruction_fil.pdf"))
if SAVE_PGF
    savefig(string(folder,"log_reconstruction_fil.pgf"))
end
close("all")

# displayPSD(753,t_samp/3600.0,diameter,x0*x_pre_all[R_psd,:])
titlePSD = @sprintf "particle size distribution (%1.2e,%1.2e)" 1.0e-9x0*minimum(x_pre_all[R_psd,:]./delta) 1.0e-9x0*maximum(x_pre_all[R_psd,:]./delta)
# cblabelPSD = "density [# cm\$^{-3}\$ nm\$^{-1}\$]"
cblabelPSD = "density [cm\$^{-3}\$ nm\$^{-1}\$]"
displayLogData2D(753,t_samp/3600.0,1.0e9diameter,1.0e-9x0*x_pre_all[R_psd,:]./delta,minPSD,maxPSD,_title=titlePSD,_colorbar_label=cblabelPSD)
ylabel("diameter [nm]")
savefig(string(folder,"reconstruction_pre.png"))
savefig(string(folder,"reconstruction_pre.pdf"))
if SAVE_PGF
    savefig(string(folder,"reconstruction_pre.pgf"))
end
close("all")

titlePSDlog = @sprintf "log size distribution (%1.2e,%1.2e)" x0*minimum(x_pre_all[R_psd,:]./log10(cst_r)) x0*maximum(x_pre_all[R_psd,:]./log10(cst_r))
# cblabelPSDlog = "density [# cm\$^{-3}\$ log\$_{10}(\\frac{d_i}{d_{i-1}})^{-1}\$]"
cblabelPSDlog = "density [\$\\frac{\\mathrm{d} N}{\\mathrm{d} \\log_{10}(D_p)}\$]"
displayLogData2D(753,t_samp/3600.0,1.0e9diameter,x0*x_pre_all[R_psd,:]./log10(cst_r),minPSDlog,maxPSDlog,_title=titlePSDlog,_colorbar_label=cblabelPSDlog)
ylabel("diameter [nm]")
savefig(string(folder,"log_reconstruction_pre.png"))
savefig(string(folder,"log_reconstruction_pre.pdf"))
if SAVE_PGF
    savefig(string(folder,"log_reconstruction_pre.pgf"))
end
close("all")

# displayPSD(123,t_samp/3600.0,diameter,x0*x_smo_all[R_psd,:])
titlePSD = @sprintf "particle size distribution (%1.2e,%1.2e)" 1.0e-9x0*minimum(x_smo_all[R_psd,:]./delta) 1.0e-9x0*maximum(x_smo_all[R_psd,:]./delta)
# cblabelPSD = "density [# cm\$^{-3}\$ nm\$^{-1}\$]"
cblabelPSD = "density [cm\$^{-3}\$ nm\$^{-1}\$]"
displayLogData2D(123,t_samp/3600.0,1.0e9diameter,1.0e-9x0*x_smo_all[R_psd,:]./delta,minPSD,maxPSD,_title=titlePSD,_colorbar_label=cblabelPSD)
ylabel("diameter [nm]")
savefig(string(folder,"reconstruction_smo.png"))
savefig(string(folder,"reconstruction_smo.pdf"))
if SAVE_PGF
    savefig(string(folder,"reconstruction_smo.pgf"))
end
close("all")


# displayPSD(123,t_samp/3600.0,diameter,x0*x_smo_all[R_psd,:])
titlePSDlog = @sprintf "log size distribution (%1.2e,%1.2e)" x0*minimum(x_smo_all[R_psd,:]./log10(cst_r)) x0*maximum(x_smo_all[R_psd,:]./log10(cst_r))
# cblabelPSDlog = "density [# cm\$^{-3}\$ log\$_{10}(\\frac{d_i}{d_{i-1}})^{-1}\$]"
cblabelPSDlog = "density [\$\\frac{\\mathrm{d} N}{\\mathrm{d} \\log_{10}(D_p)}\$]"
displayLogData2D(123,t_samp/3600.0,1.0e9diameter,x0*x_smo_all[R_psd,:]./log10(cst_r),minPSDlog,maxPSDlog,_title=titlePSDlog,_colorbar_label=cblabelPSDlog)
ylabel("diameter [nm]")
savefig(string(folder,"log_reconstruction_smo.png"))
savefig(string(folder,"log_reconstruction_smo.pdf"))
if SAVE_PGF
    savefig(string(folder,"log_reconstruction_smo.pgf"))
end
close("all")




minVarLoss = minimum(var_model[R_loss,R_loss])
maxVarLoss = maximum(var_model[R_loss,R_loss])
titleLoss = @sprintf "loss rate covariance model (%1.1e,%1.1e)" sqrt(minVarLoss) sqrt(maxVarLoss)
cblabelLoss = "loss rate [(s\$^{-1}\$)\$^{2}\$]"
fig,ax = imshowData(259,diameter,diameter,var_model[R_loss,R_loss],_norm=:Normalize,_vmin=0.0,_vmax=maxVarLoss,_edgecolors="face")
ylim(diameter[end],diameter[1])
yscale("log")
xscale("log")
title(titleLoss)
xlabel("diameter [nm]")
ylabel("diameter [nm]")
cbar=colorbar()
cbar.set_label(cblabelLoss)
cbar.formatter.set_powerlimits((0,0))
cbar.update_ticks()
savefig(string(folder,"covariance_model_loss.png"))
savefig(string(folder,"covariance_model_loss.pdf"))
if SAVE_PGF
    savefig(string(folder,"covariance_model_loss.pgf"))
end
close("all")








if GT_loaded
    # condensational growth rate
    figure(49)
    plot(tp,condensation_rate_all[1,:])
    s = @sprintf "growth rate"
    title(s)
    xlabel("time [h]")
    ylabel("condensation growth rate [m s\$^{-1}\$]")
    grid(true,which="major",ls="-")
    grid(true,which="minor",ls="-",alpha=0.5)
    xlim(t_samp[1]/3600.0,t_samp[end]/3600.0)
    savefig(string(folder,"cond_simu.png"))
    savefig(string(folder,"cond_simu.pdf"))
    if SAVE_PGF
        savefig(string(folder,"cond_simu.pgf"))
    end



    figure(43)
    plot(tp,nucleation_rate_all)
    s = @sprintf "nucleation rate"
    title(s)
    xlabel("time [h]")
    # ylabel("nucleation rate [# cm\$^{-3}\$ s\$^{-1}\$]")
    ylabel("nucleation rate [cm\$^{-3}\$ s\$^{-1}\$]")
    grid(true,which="major",ls="-")
    grid(true,which="minor",ls="-",alpha=0.5)
    xlim(t_samp[1]/3600.0,t_samp[end]/3600.0)
    savefig(string(folder,"nucleation_rate.png"))
    savefig(string(folder,"nucleation_rate.pdf"))
    if SAVE_PGF
        savefig(string(folder,"nucleation_rate.pgf"))
    end



    figure(44)
    loglog(dp*10.0^9,wall_rate_expected)
    s = @sprintf "wall losses"
    title(s)
    xlabel("diameter [nm]")
    ylabel("deposition rate [s\$^{-1}\$]")
    grid(true,which="major",ls="-")
    grid(true,which="minor",ls="-",alpha=0.5)
    xlim(diameter[1]*10.0^9,diameter[end]*10.0^9)
    ylim(8.0e-6,1.0e-3)
    savefig(string(folder,"wall_loss_rate.png"))
    savefig(string(folder,"wall_loss_rate.pdf"))
    if SAVE_PGF
        savefig(string(folder,"wall_loss_rate.pgf"))
    end
end

close("all")





# condensational growth rate
figure(50)
plot(t_samp/3600.0,GR0*CGR(x_fil_all[R_cond_init,:]))
fill_between(t_samp/3600.0,GR0*percentiles_cond_fil[:,1],GR0*percentiles_cond_fil[:,2],alpha=0.5,color=:blue)
plot(t_samp/3600.0,GR0*CGR(x_smo_all[R_cond_init,:]))
fill_between(t_samp/3600.0,GR0*percentiles_cond_smo[:,1],GR0*percentiles_cond_smo[:,2],alpha=0.5,color=:orange)
# plot(t_samp/3600.0,GR0*percentiles_cond_smo[:,1])
# plot(t_samp/3600.0,GR0*percentiles_cond_smo[:,2])
if GT_loaded
    plot(tp,condensation_rate_all[1,:],color=:green)
end
# s = @sprintf "evolution of the growth rate percentiles"
s = @sprintf "condensation growth rate"
title(s)
xlabel("time [h]")
ylabel("growth rate [m s\$^{-1}\$]")
grid(true,which="major",ls="-")
grid(true,which="minor",ls="-",alpha=0.5)
xlim(t_samp[1]/3600.0,t_samp[end]/3600.0)
ylim(-0.1e-12,40*2.78e-13)
legend(["filter","smoother","ground truth","EKF uncertainty","FIKS uncertainty"])
savefig(string(folder,"evolution_cond_percentiles_smo.png"))
savefig(string(folder,"evolution_cond_percentiles_smo.pdf"))
if SAVE_PGF
    savefig(string(folder,"evolution_cond_percentiles_smo.pgf"))
end





# wall loss rate
for k in 1:length(idx_time)
    figure(k)
    loglog(diameter*10.0^9,gamma0*wall_rate(x_fil_all[R_loss,idx_time[k]]))
    fill_between(diameter*10.0^9,gamma0*percentiles_wall_fil[:,idx_time[k],1],gamma0*percentiles_wall_fil[:,idx_time[k],2],alpha=0.5,color=:blue)
    loglog(diameter*10.0^9,gamma0*wall_rate(x_smo_all[R_loss,idx_time[k]]))
    fill_between(diameter*10.0^9,gamma0*percentiles_wall_smo[:,idx_time[k],1],gamma0*percentiles_wall_smo[:,idx_time[k],2],alpha=0.5,color=:orange)
    # loglog(diameter*10.0^9,gamma0*percentiles_wall_smo[:,idx_time[k],1])
    # loglog(diameter*10.0^9,gamma0*percentiles_wall_smo[:,idx_time[k],2])
    if GT_loaded
        loglog(dp*10.0^9,wall_rate_expected,color=:green)
    end
    # s = @sprintf "percentile of the wall loss at time %i" idx_time[k]
    local s = @sprintf "wall losses"
    title(s)
    xlabel("diameter [nm]")
    ylabel("deposition rate [s\$^{-1}\$]")
    grid(true,which="major",ls="-")
    grid(true,which="minor",ls="-",alpha=0.5)
    xlim(diameter[1]*10.0^9,diameter[end]*10.0^9)
    ylim(8.0e-6,1.0e-3)
    # legend(["smoother estimation","ground truth","uncertainties"])
    # legend(["smoother estimation","smoother: 15\$^{th}\$ percentile","smoother: 85\$^{th}\$ percentile", "expectation"])
    legend(["filter","smoother","ground truth","EKF uncertainty","FIKS uncertainty"])
    s = @sprintf "loss_percentiles_wrt_diameter_time_%i_smo.png" idx_time[k]
    savefig(string(folder,s))
    s = @sprintf "loss_percentiles_wrt_diameter_time_%i_smo.pdf" idx_time[k]
    savefig(string(folder,s))
    if SAVE_PGF
        s = @sprintf "loss_percentiles_wrt_diameter_time_%i_smo.pgf" idx_time[k]
        savefig(string(folder,s))
    end
end
close("all")



# nucleation rate
figure(54)
plot(t_samp/3600.0,J0*Nucleation_rate(x_fil_all[R_nuc_init,:]))
fill_between(t_samp/3600.0,J0*percentiles_nuc_fil[:,1],J0*percentiles_nuc_fil[:,2],alpha=0.5,color=:blue)
plot(t_samp/3600.0,J0*Nucleation_rate(x_smo_all[R_nuc_init,:]))
fill_between(t_samp/3600.0,J0*percentiles_nuc_smo[:,1],J0*percentiles_nuc_smo[:,2],alpha=0.5,color=:orange)
# plot(t_samp/3600.0,J0*percentiles_nuc_smo[:,1])
# plot(t_samp/3600.0,J0*percentiles_nuc_smo[:,2])
if GT_loaded
    plot(tp,nucleation_rate_all,color=:green)
end
s = @sprintf "nucleation rate"
title(s)
xlabel("time [h]")
# ylabel("nucleation rate [# cm\$^{-3}\$ s\$^{-1}\$]")
ylabel("nucleation rate [cm\$^{-3}\$ s\$^{-1}\$]")
grid(true,which="major",ls="-")
grid(true,which="minor",ls="-",alpha=0.5)
xlim(t_samp[1]/3600.0,t_samp[end]/3600.0)
ylim(-0.02,0.35)
# legend(["smoother estimation","ground truth","uncertainties"])
legend(["filter","smoother","ground truth","EKF uncertainty","FIKS uncertainty"])
savefig(string(folder,"evolution_nuc_percentiles_smo.png"))
savefig(string(folder,"evolution_nuc_percentiles_smo.pdf"))
if SAVE_PGF
    savefig(string(folder,"evolution_nuc_percentiles_smo.pgf"))
end



figure(56)
NN = 10
sample_nuc = Array{Cdouble,3}(undef,n_samp,2,NN)

sample_nuc[1,:,:] = 100.0*randn(2,NN) #  sig_nuc*randn(2,NN)

for i in 2:n_samp
    for k in 1:NN
        sample_nuc[i,:,k] = B_nuc_time*sample_nuc[i-1,:,k] + [100.0;0.0]*randn()
    end
end

plot(sample_nuc[:,1,:])




figure(57)
normFO = Array{Cdouble,1}(undef,n_samp)
normFN = Array{Cdouble,1}(undef,n_samp)
normFNO = Array{Cdouble,1}(undef,n_samp)
for i in 1:n_samp
    normFN[i] = norm(myWSKF.F_ev_all[R_psd,R_psd,i])
    normFO[i] = norm(myWSKF.F_ev_all[nbin+1:end,nbin+1:end,i])
    normFNO[i] = norm(myWSKF.F_ev_all[R_psd,nbin+1:end,i])
end
plot(normFN)
plot(normFO)
plot(normFNO)
legend(["\$\\|F_n\\|_2\$","\$\\|F_o\\|_2\$","\$\\|F_{no}\\|_2\$"])




idx_ttt = 301 # 120s
# idx_ttt = 61 # 600s
# idx_ttt_simu = 1201
idx_ttt_simu = 12001 # for FULL_SAVE simulated data
# for idx_ttt in 180:182
    figure()
    semilogx(1.0e9diameter,1.0e-9x0*x_fil_all[R_psd,idx_ttt]./delta)
    semilogx(1.0e9diameter,1.0e-9x0*x_smo_all[R_psd,idx_ttt]./delta)
    # semilogx(1.0e9diameter,1.0e-9y_all[:,idx_ttt]./delta)
    semilogx(1.0e9dp,PSD_simu[:,idx_ttt_simu])
    x_pre_up = x0*x_pre_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_pre_all[R_psd,R_psd,idx_ttt]))
    x_pre_do = x0*x_pre_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_pre_all[R_psd,R_psd,idx_ttt]))
    x_fil_up = x0*x_fil_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_fil_all[R_psd,R_psd,idx_ttt]))
    x_fil_do = x0*x_fil_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_fil_all[R_psd,R_psd,idx_ttt]))
    x_smo_up = x0*x_smo_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_smo_all[R_psd,R_psd,idx_ttt]))
    x_smo_do = x0*x_smo_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_smo_all[R_psd,R_psd,idx_ttt]))
    fill_between(1.0e9diameter,1.0e-9x_fil_do./delta,1.0e-9x_fil_up./delta,alpha=0.5,color=:blue)
    fill_between(1.0e9diameter,1.0e-9x_smo_do./delta,1.0e-9x_smo_up./delta,alpha=0.5,color=:orange)
    # fill_between(1.0e9diameter,1.0e-9x_pre_do./delta,1.0e-9x_pre_up./delta,alpha=0.5,color=:green)
    xlim(1.0e1,1.0e3)
    ylim(-2.0,60.0)
    xlabel("diameter [nm]")
    # ylabel("density [# cm\$^{-3}\$ nm\$^{-1}\$]")
    ylabel("density [cm\$^{-3}\$ nm\$^{-1}\$]")
    title("size distribution at time: 10h")
    legend(["filter", "smoother","ground truth","EKF uncertainty","FIKS uncertainty"])
    savefig(string(folder,"size_density_time_10h_smo.png"))
    savefig(string(folder,"size_density_time_10h_smo.pdf"))
    if SAVE_PGF
        savefig(string(folder,"size_density_time_10h_smo.pgf"))
    end
    # formatter.set_powerlimits((-1,2))
# end


figure()
semilogx(1.0e9diameter,x0*x_fil_all[R_psd,idx_ttt]./log10(cst_r))
semilogx(1.0e9diameter,x0*x_smo_all[R_psd,idx_ttt]./log10(cst_r))
# semilogx(1.0e9diameter,1.0e-9y_all[:,idx_ttt]./delta)
semilogx(1.0e9dp,1.0e9PSD_simu[:,idx_ttt_simu].*delta_p/log10(cst_r_p))
x_pre_up = x0*x_pre_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_pre_all[R_psd,R_psd,idx_ttt]))
x_pre_do = x0*x_pre_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_pre_all[R_psd,R_psd,idx_ttt]))
x_fil_up = x0*x_fil_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_fil_all[R_psd,R_psd,idx_ttt]))
x_fil_do = x0*x_fil_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_fil_all[R_psd,R_psd,idx_ttt]))
x_smo_up = x0*x_smo_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_smo_all[R_psd,R_psd,idx_ttt]))
x_smo_do = x0*x_smo_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_smo_all[R_psd,R_psd,idx_ttt]))
fill_between(1.0e9diameter,x_fil_do./log10(cst_r),x_fil_up./log10(cst_r),alpha=0.5,color=:blue)
fill_between(1.0e9diameter,x_smo_do./log10(cst_r),x_smo_up./log10(cst_r),alpha=0.5,color=:orange)
ax = gca();
# ax[:ticklabel_format](style="sci",axis="y",scilimits=(0,0))
ax.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
xlim(1.0e9diameter[1],1.0e9diameter[end])
ylim(-100.0)
xlabel("diameter [nm]")
# ylabel("density [# cm\$^{-3}\$ log\$_{10}(\\frac{d_i}{d_{i-1}})^{-1}\$]")
ylabel("density [\$\\frac{\\mathrm{d} N}{\\mathrm{d} \\log_{10}(D_p)}\$]")
title("size distribution at time: 10h")
legend(["filter", "smoother","ground truth","EKF uncertainty","FIKS uncertainty"])
savefig(string(folder,"log_size_density_time_10h_smo.png"))
savefig(string(folder,"log_size_density_time_10h_smo.pdf"))
if SAVE_PGF
    savefig(string(folder,"log_size_density_time_10h_smo.pgf"))
end










idx_ttt = 61 # 120s
# idx_ttt = 61 # 600s
# idx_ttt_simu = 1201
idx_ttt_simu = 2401 # for FULL_SAVE simulated data
# for idx_ttt in 180:182
    figure()
    semilogx(1.0e9diameter,1.0e-9x0*x_fil_all[R_psd,idx_ttt]./delta)
    semilogx(1.0e9diameter,1.0e-9x0*x_smo_all[R_psd,idx_ttt]./delta)
    # semilogx(1.0e9diameter,1.0e-9y_all[:,idx_ttt]./delta)
    semilogx(1.0e9dp,PSD_simu[:,idx_ttt_simu])
    x_pre_up = x0*x_pre_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_pre_all[R_psd,R_psd,idx_ttt]))
    x_pre_do = x0*x_pre_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_pre_all[R_psd,R_psd,idx_ttt]))
    x_fil_up = x0*x_fil_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_fil_all[R_psd,R_psd,idx_ttt]))
    x_fil_do = x0*x_fil_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_fil_all[R_psd,R_psd,idx_ttt]))
    x_smo_up = x0*x_smo_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_smo_all[R_psd,R_psd,idx_ttt]))
    x_smo_do = x0*x_smo_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_smo_all[R_psd,R_psd,idx_ttt]))
    fill_between(1.0e9diameter,1.0e-9x_fil_do./delta,1.0e-9x_fil_up./delta,alpha=0.5,color=:blue)
    fill_between(1.0e9diameter,1.0e-9x_smo_do./delta,1.0e-9x_smo_up./delta,alpha=0.5,color=:orange)
    # fill_between(1.0e9diameter,1.0e-9x_pre_do./delta,1.0e-9x_pre_up./delta,alpha=0.5,color=:green)
    xlim(1.0e1,1.0e3)
    ylim(-2.0,80.0)
    xlabel("diameter [nm]")
    # ylabel("density [# cm\$^{-3}\$ nm\$^{-1}\$]")
    ylabel("density [cm\$^{-3}\$ nm\$^{-1}\$]")
    title("size distribution at time: 2h")
    legend(["filter", "smoother","ground truth","EKF uncertainty","FIKS uncertainty"])
    savefig(string(folder,"size_density_time_2h_smo.png"))
    savefig(string(folder,"size_density_time_2h_smo.pdf"))
    if SAVE_PGF
        savefig(string(folder,"size_density_time_2h_smo.pgf"))
    end
    # formatter.set_powerlimits((-1,2))
# end





figure()
semilogx(1.0e9diameter,x0*x_fil_all[R_psd,idx_ttt]./log10(cst_r))
semilogx(1.0e9diameter,x0*x_smo_all[R_psd,idx_ttt]./log10(cst_r))
# semilogx(1.0e9diameter,1.0e-9y_all[:,idx_ttt]./delta)
semilogx(1.0e9dp,1.0e9PSD_simu[:,idx_ttt_simu].*delta_p/log10(cst_r_p))
x_pre_up = x0*x_pre_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_pre_all[R_psd,R_psd,idx_ttt]))
x_pre_do = x0*x_pre_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_pre_all[R_psd,R_psd,idx_ttt]))
x_fil_up = x0*x_fil_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_fil_all[R_psd,R_psd,idx_ttt]))
x_fil_do = x0*x_fil_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_fil_all[R_psd,R_psd,idx_ttt]))
x_smo_up = x0*x_smo_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_smo_all[R_psd,R_psd,idx_ttt]))
x_smo_do = x0*x_smo_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_smo_all[R_psd,R_psd,idx_ttt]))
fill_between(1.0e9diameter,x_fil_do./log10(cst_r),x_fil_up./log10(cst_r),alpha=0.5,color=:blue)
fill_between(1.0e9diameter,x_smo_do./log10(cst_r),x_smo_up./log10(cst_r),alpha=0.5,color=:orange)
ax = gca();
# ax[:ticklabel_format](style="sci",axis="y",scilimits=(0,0))
ax.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
xlim(1.0e9diameter[1],1.0e9diameter[end])
ylim(-100.0)
xlabel("diameter [nm]")
# ylabel("density [# cm\$^{-3}\$ log\$_{10}(\\frac{d_i}{d_{i-1}})^{-1}\$]")
ylabel("density [\$\\frac{\\mathrm{d} N}{\\mathrm{d} \\log_{10}(D_p)}\$]")
title("size distribution at time: 2h")
legend(["filter", "smoother","ground truth","EKF uncertainty","FIKS uncertainty"])
savefig(string(folder,"log_size_density_time_2h_smo.png"))
savefig(string(folder,"log_size_density_time_2h_smo.pdf"))
if SAVE_PGF
    savefig(string(folder,"log_size_density_time_2h_smo.pgf"))
end
