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

idx_bin = [1; round(Int64,nbin/3); round(Int64,2nbin/3)]
# idx_time = [1; round(Int64,n_samp/2); n_samp]
idx_time = [1; round(Int64,n_samp/6); round(Int64,2n_samp/6); round(Int64,3n_samp/6); round(Int64,4n_samp/6); round(Int64,5n_samp/6); round(Int64,11n_samp/12); round(Int64,23n_samp/24); n_samp]


# size distribution
maxPSD = maximum(1.0e-9y_all./delta) 
minPSD = max(0.001*maxPSD,1.0e-9minimum(1.0e-9y_all./delta))
maxPSDlog = maximum(y_all./log10(cst_r))
minPSDlog = max(0.001*maxPSDlog,minimum(y_all./log10(cst_r)))

titlePSD = @sprintf "particle size distribution (%1.2e,%1.2e)" 1.0e-9minimum(y_all./delta) 1.0e-9maximum(y_all./delta)
cblabelPSD = "density [cm\$^{-3}\$ nm\$^{-1}\$]"
fig,ax,cbar = displayLogData2D(158,t_samp/3600.0,1.0e9diameter,1.0e-9y_all./delta,minPSD,maxPSD,_title=titlePSD,_colorbar_label=cblabelPSD)
# yform  = ax[:yaxis][:set_major_formatter]
# stick = matplotlib[:ticker][:ScalarFormatter]
# yform(stick())
# ax[:ticklabel_format](style="sci",axis="y",scilimits=(0,0))
ylabel("diameter [nm]")
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
else
    maxPSD = maximum(1.0e-9x0*x_smo_all[R_psd,:]./delta) # maximum(1.0e-9y_all./delta) 
    minPSD = max(0.001*maxPSD,1.0e-9minimum(1.0e-9x0*x_smo_all[R_psd,:]./delta))
    maxPSDlog = x0*maximum(x_smo_all[R_psd,:]./log10(cst_r)) # maximum(y_all./log10(cst_r))
    minPSDlog = max(0.001*maxPSDlog,minimum(x0*x_smo_all[R_psd,:]./log10(cst_r)))
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



# covariance of the evolution model
minVarCond = minimum(var_model[R_cond,R_cond])
maxVarCond = maximum(var_model[R_cond,R_cond])
titleCond = @sprintf "condensation covariance model (%1.1e,%1.1e)" sqrt(minVarCond) sqrt(maxVarCond)
cblabelCond = "condensation rate [(m s\$^{-1}\$)\$^{2}\$]"
fig,ax = imshowData(258,1.0e9diameter,1.0e9diameter,var_model[R_cond,R_cond],_norm=:Normalize,_vmin=0.0minVarCond,_vmax=maxVarCond,_edgecolors="face")
ylim(1.0e9diameter[end],1.0e9diameter[1])
yscale("log")
xscale("log")
title(titleCond)
xlabel("diameter [nm]")
ylabel("diameter [nm]")
cbar=colorbar()
cbar.set_label(cblabelCond)
cbar.update_ticks()
savefig(string(folder,"covariance_model_condensation.png"))
savefig(string(folder,"covariance_model_condensation.pdf"))
if SAVE_PGF
    savefig(string(folder,"covariance_model_condensation.pgf"))
end



minVarLoss = minimum(var_model[R_loss,R_loss])
maxVarLoss = maximum(var_model[R_loss,R_loss])
titleLoss = @sprintf "loss rate covariance model (%1.1e,%1.1e)" sqrt(minVarLoss) sqrt(maxVarLoss)
cblabelLoss = "loss rate [(s\$^{-1}\$)\$^{2}\$]"
fig,ax = imshowData(259,1.0e9diameter,1.0e9diameter,var_model[R_loss,R_loss],_norm=:Normalize,_vmin=0.0minVarLoss,_vmax=maxVarLoss,_edgecolors="face")
ylim(1.0e9diameter[end],1.0e9diameter[1])
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


titleCond = "condensation rate"
cblabelCond = "growth rate [m s\$^{-1}\$]"
if GT_loaded
    minCond = max(5.0e-14,minimum(condensation_rate_all)) # max(1.0e-14,minimum(GR0*CGR(x_fil_all[R_cond,:]))) # max(1.0e-14,minimum(GR0*CGR(x_smo_all[R_cond,:])))
    maxCond = maximum(condensation_rate_all)              # maximum(GR0*CGR(x_fil_all[R_cond,:])) # maximum(GR0*CGR(x_smo_all[R_cond,:]))
else
    minCond,maxCond = extrema([GR0*CGR2D(x_fil_all[R_cond,:]); GR0*CGR2D(x_smo_all[R_cond,:])])
    minCond = max(5.0e-14,minCond) 
end
displayLogData2D(123,t_samp/3600.0,1.0e9diameter,GR0*CGR2D(x_fil_all[R_cond,:]),minCond,maxCond,_title=titleCond,_colorbar_label=cblabelCond)
ylabel("diameter [nm]")
savefig(string(folder,"cond_fil.png"))
savefig(string(folder,"cond_fil.pdf"))
if SAVE_PGF
    savefig(string(folder,"cond_fil.pgf"))
end
displayLogData2D(456,t_samp/3600.0,1.0e9diameter,GR0*CGR2D(x_smo_all[R_cond,:]),minCond,maxCond,_title=titleCond,_colorbar_label=cblabelCond)
ylabel("diameter [nm]")
savefig(string(folder,"cond_smo.png"))
savefig(string(folder,"cond_smo.pdf"))
if SAVE_PGF
    savefig(string(folder,"cond_smo.pgf"))
end
if GT_loaded
    displayLogData2D(789,tp,1.0e9dp,condensation_rate_all,minCond,maxCond,_title=titleCond,_colorbar_label=cblabelCond)
    ylabel("diameter [nm]")
    savefig(string(folder,"cond_gt.png"))
    savefig(string(folder,"cond_gt.pdf"))
    # if SAVE_PGF
    #     savefig(string(folder,"cond_gt.pgf"))
    # end

    minCond = 1.0e-1
    maxCond = 1.0
    displayLogData2D(753,tp,1.0e9dp,abs.(GR0*CGR2D(x_smo_all[R_cond,:])-condensation_rate_all)./condensation_rate_all,minCond,maxCond,_title=titleCond,_colorbar_label="relative error")
    ylabel("diameter [nm]")
    savefig(string(folder,"cond_relative_error_smo.png"))
    savefig(string(folder,"cond_relative_error_smo.pdf"))
    if SAVE_PGF
        savefig(string(folder,"cond_relative_error_smo.pgf"))
    end
    displayLogData2D(754,tp,1.0e9dp,abs.(GR0*CGR2D(x_fil_all[R_cond,:])-condensation_rate_all)./condensation_rate_all,minCond,maxCond,_title=titleCond,_colorbar_label="relative error")
    ylabel("diameter [nm]")
    savefig(string(folder,"cond_relative_error_fil.png"))
    savefig(string(folder,"cond_relative_error_fil.pdf"))
    if SAVE_PGF
        savefig(string(folder,"cond_relative_error_fil.pgf"))
    end
end
close("all")




if GT_loaded
    figure(1043)
    plot(tp,nucleation_rate_all)
    s = @sprintf "nucleation rate"
    title(s)
    xlabel("time [h]")
    ylabel("nucleation rate [cm\$^{-3}\$ s\$^{-1}\$]")
    grid(true,which="major",ls="-")
    grid(true,which="minor",ls="-",alpha=0.5)
    xlim(t_samp[1]/3600.0,t_samp[end]/3600.0)
    savefig(string(folder,"nucleation_rate.png"))
    savefig(string(folder,"nucleation_rate.pdf"))
    if SAVE_PGF
        savefig(string(folder,"nucleation_rate.pgf"))
    end
end




figure(1000)
# loglog(dp*10.0^9,wall_rate_expected)
loglog(1.0e9diameter,wall_rate_expected)
s = @sprintf "wall loss rate"
title(s)
xlabel("diameter [nm]")
ylabel("deposition rate [s\$^{-1}\$]")
grid(true,which="major",ls="-")
grid(true,which="minor",ls="-",alpha=0.5)
xlim(diameter[1]*10.0^9,diameter[end]*10.0^9)
ylim(5e-5,2e-3)
savefig(string(folder,"wall_loss_rate.png"))
savefig(string(folder,"wall_loss_rate.pdf"))
if SAVE_PGF
    savefig(string(folder,"wall_loss_rate.pgf"))
end

close("all")


o_fil_cond = Array{Cdouble}(undef,nbin,n_samp)
o_smo_cond = Array{Cdouble}(undef,nbin,n_samp)
o_fil_loss = Array{Cdouble}(undef,nbin,n_samp)
o_smo_loss = Array{Cdouble}(undef,nbin,n_samp)
for k in 1:n_samp
    o_fil_cond[:,k] = diag(o_fil_all[R_cond,R_cond,k])
    o_smo_cond[:,k] = diag(o_smo_all[R_cond,R_cond,k])
    o_fil_loss[:,k] = diag(o_fil_all[R_loss,R_loss,k])
    o_smo_loss[:,k] = diag(o_smo_all[R_loss,R_loss,k])
end
err_fil_cond = sqrt.(o_fil_cond) # 0.5GR0*(percentiles_cond_fil[:,:,2]-percentiles_cond_fil[:,:,1]);
displayLogData2D(123,t_samp/3600.0,1.0e9diameter,err_fil_cond,max(0.0,minimum(err_fil_cond)),maximum(err_fil_cond),_title=titleCond,_colorbar_label=cblabelCond)
ylabel("diameter [nm]")
savefig(string(folder,"cond_fil_std.png"))
savefig(string(folder,"cond_fil_std.pdf"))
if SAVE_PGF
    savefig(string(folder,"cond_fil_std.pgf"))
end

err_smo_cond = sqrt.(o_smo_cond) # 0.5GR0*(percentiles_cond_smo[:,:,2]-percentiles_cond_smo[:,:,1]);
displayLogData2D(456,t_samp/3600.0,1.0e9diameter,err_smo_cond,max(0.0,minimum(err_smo_cond)),maximum(err_smo_cond),_title=titleCond,_colorbar_label=cblabelCond)
ylabel("diameter [nm]")
savefig(string(folder,"cond_smo_std.png"))
savefig(string(folder,"cond_smo_std.pdf"))
if SAVE_PGF
    savefig(string(folder,"cond_smo_std.pgf"))
end
close("all")

for k in 1:length(idx_time)
    # condensational growth rate with percentiles
    figure(k)
    semilogx(1.0e9diameter,GR0*CGR(x_fil_all[R_cond,idx_time[k]])) # semilogx # loglog
    semilogx(1.0e9diameter,GR0*CGR(x_smo_all[R_cond,idx_time[k]])) # semilogx # loglog
    if GT_loaded
        semilogx(1.0e9dp,condensation_rate_all[:,idx_time[k]],color=:green)         # semilogx # loglog
    end
    fill_between(1.0e9diameter,GR0*percentiles_cond_fil[:,idx_time[k],1],GR0*percentiles_cond_fil[:,idx_time[k],2],alpha=0.5,color=:blue)
    fill_between(1.0e9diameter,GR0*percentiles_cond_smo[:,idx_time[k],1],GR0*percentiles_cond_smo[:,idx_time[k],2],alpha=0.5,color=:orange)
    local s = @sprintf "growth rate at %ih%imin" floor(Int64,idx_time[k]*dt/3600.0) floor(Int64,60.0*(idx_time[k]*dt/3600.0-floor(Int64,idx_time[k]*dt/3600.0)))
    title(s)
    xlabel("diameter [nm]")
    ylabel("growth rate [m s\$^{-1}\$]")
    grid(true,which="major",ls="-")
    grid(true,which="minor",ls="-",alpha=0.5)
    xlim(1.0e9diameter[1],1.0e9diameter[end])
    ylim(0.0,2.5e-12)
    if GT_loaded
        legend(["filter","smoother", "ground truth", "EKS uncertainty", "FIKS uncertainty"])
    else
        legend(["filter","smoother", "EKS uncertainty", "FIKS uncertainty"])
    end
    s = @sprintf "cond_wrt_diameter_time_%i.png" idx_time[k]
    savefig(string(folder,s))
    s = @sprintf "cond_wrt_diameter_time_%i.pdf" idx_time[k]
    savefig(string(folder,s))
    if SAVE_PGF
        s = @sprintf "cond_wrt_diameter_time_%i.pgf" idx_time[k]
        savefig(string(folder,s))
    end

    # wall loss rate with percentiles
    figure(1000+k)
    loglog(diameter*10.0^9,gamma0*wall_rate(x_fil_all[R_loss,idx_time[k]])) # wall_rate_fil_last)
    loglog(diameter*10.0^9,gamma0*wall_rate(x_smo_all[R_loss,idx_time[k]]))
    loglog(diameter*10.0^9,wall_rate_expected,color=:green)
    fill_between(diameter*10.0^9,gamma0*percentiles_wall_fil[:,idx_time[k],1],gamma0*percentiles_wall_fil[:,idx_time[k],2],alpha=0.5,color=:blue)
    fill_between(diameter*10.0^9,gamma0*percentiles_wall_smo[:,idx_time[k],1],gamma0*percentiles_wall_smo[:,idx_time[k],2],alpha=0.5,color=:orange)
    s = @sprintf "wall loss at %ih%imin" floor(Int64,idx_time[k]*dt/3600.0) floor(Int64,60.0*(idx_time[k]*dt/3600.0-floor(Int64,idx_time[k]*dt/3600.0)))
    title(s)
    xlabel("diameter [nm]")
    ylabel("deposition rate [s\$^{-1}\$]")
    grid(true,which="major",ls="-")
    grid(true,which="minor",ls="-",alpha=0.5)
    xlim(diameter[1]*10.0^9,diameter[end]*10.0^9)
    ylim(5e-5,2e-3)
    # legend(["filter","smoother", "ground truth"])
    legend(["filter","smoother", "ground truth", "EKF uncertainty","FIKS uncertainty"])
    s = @sprintf "loss_wrt_diameter_time_%i.png" idx_time[k]
    savefig(string(folder,s))
    s = @sprintf "loss_wrt_diameter_time_%i.pdf" idx_time[k]
    savefig(string(folder,s))
    if SAVE_PGF
        s = @sprintf "loss_wrt_diameter_time_%i.pgf" idx_time[k]
        savefig(string(folder,s))
    end
end
close("all")



estimation_cond_smo = GR0*CGR2D(x_smo_all[R_cond,:]);
estimation_cond_fil = GR0*CGR2D(x_fil_all[R_cond,:]);
# titleCond = "condensation rate"
# cblabelCond = "condensation rate [m s\$^{-1}\$]"
# minCond = max(5.0e-14,minimum(condensation_rate_all)) # max(1.0e-14,minimum(GR0*CGR(x_fil_all[R_cond,:]))) # max(1.0e-14,minimum(GR0*CGR(x_smo_all[R_cond,:])))
# maxCond = maximum(condensation_rate_all)              # maximum(GR0*CGR(x_fil_all[R_cond,:])) # maximum(GR0*CGR(x_smo_all[R_cond,:]))
# displayLogData2D(456,t_samp/3600.0,diameter,estimation_cond_smo,minCond,maxCond,_title=titleCond,_colorbar_label=cblabelCond)


close("all")



# nucleation with percentiles
figure(1043)
plot(t_samp/3600.0,J0*Nucleation_rate(x_fil_all[R_nuc_init,:]))
plot(t_samp/3600.0,J0*Nucleation_rate(x_smo_all[R_nuc_init,:]))
if GT_loaded
    plot(tp,nucleation_rate_all,color=:green)
    plot(tp,1.0e9condensation_rate_all[1,:].*PSD_simu[1,:],color=:darkgreen)
end
fill_between(t_samp/3600.0,J0*percentiles_nuc_fil[:,1],J0*percentiles_nuc_fil[:,2],alpha=0.5,color=:blue)
fill_between(t_samp/3600.0,J0*percentiles_nuc_smo[:,1],J0*percentiles_nuc_smo[:,2],alpha=0.5,color=:orange)
# legend(["filter","smoother", "ground truth"])
if GT_loaded
    legend(["filter", "smoother", "ground truth", "EKF uncertainty", "FIKS uncertainty"])
else
    legend(["filter", "smoother", "EKF uncertainty", "FIKS uncertainty"])
end
s = @sprintf "evolution of the nucleation rates"
title(s)
xlabel("time [h]")
ylabel("nucleation rate [cm\$^{-3}\$ s\$^{-1}\$]")
grid(true,which="major",ls="-")
grid(true,which="minor",ls="-",alpha=0.5)
xlim(t_samp[1]/3600.0,t_samp[end]/3600.0)
savefig(string(folder,"evolution_nuc.png"))
savefig(string(folder,"evolution_nuc.pdf"))
if SAVE_PGF
    savefig(string(folder,"evolution_nuc.pgf"))
end





idx_ttt = 31 # 120s
# idx_ttt = 61 # 600s
idx_ttt_simu = 31
# for idx_ttt in 180:182
    figure()
    semilogx(1.0e9diameter,1.0e-9x0*x_fil_all[R_psd,idx_ttt]./delta)
    semilogx(1.0e9diameter,1.0e-9x0*x_smo_all[R_psd,idx_ttt]./delta)
    # semilogx(1.0e9diameter,1.0e-9y_all[:,idx_ttt]./delta)
    if GT_loaded
        semilogx(1.0e9dp,PSD_simu[:,idx_ttt_simu])
    end
    x_pre_up = x0*x_pre_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_pre_all[R_psd,R_psd,idx_ttt]))
    x_pre_do = x0*x_pre_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_pre_all[R_psd,R_psd,idx_ttt]))
    x_fil_up = x0*x_fil_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_fil_all[R_psd,R_psd,idx_ttt]))
    x_fil_do = x0*x_fil_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_fil_all[R_psd,R_psd,idx_ttt]))
    x_smo_up = x0*x_smo_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_smo_all[R_psd,R_psd,idx_ttt]))
    x_smo_do = x0*x_smo_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_smo_all[R_psd,R_psd,idx_ttt]))
    fill_between(1.0e9diameter,1.0e-9x_fil_do./delta,1.0e-9x_fil_up./delta,alpha=0.5,color=:blue)
    fill_between(1.0e9diameter,1.0e-9x_smo_do./delta,1.0e-9x_smo_up./delta,alpha=0.5,color=:orange)
    ax = gca();
    ax.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
    xlim(1.0e9diameter[1],1.0e9diameter[end])
    ylim(-100.0)
    xlabel("diameter [nm]")
    ylabel("density [cm\$^{-3}\$ nm\$^{-1}\$]")
    title("size distribution at time: 1h")
    if GT_loaded
        legend(["filter", "smoother","ground truth","EKF uncertainty","FIKS uncertainty"])
    else
        legend(["filter", "smoother","EKF uncertainty","FIKS uncertainty"])
    end
    savefig(string(folder,"size_density_time_1h_smo.png"))
    savefig(string(folder,"size_density_time_1h_smo.pdf"))
    if SAVE_PGF
        savefig(string(folder,"size_density_time_1h_smo.pgf"))
    end
    # formatter.set_powerlimits((-1,2))
# end


figure()
semilogx(1.0e9diameter,x0*x_fil_all[R_psd,idx_ttt]./log10(cst_r))
semilogx(1.0e9diameter,x0*x_smo_all[R_psd,idx_ttt]./log10(cst_r))
# semilogx(1.0e9diameter,1.0e-9y_all[:,idx_ttt]./delta)
if GT_loaded
    semilogx(1.0e9dp,1.0e9PSD_simu[:,idx_ttt_simu].*delta_p/log10(cst_r_p))
end
x_pre_up = x0*x_pre_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_pre_all[R_psd,R_psd,idx_ttt]))
x_pre_do = x0*x_pre_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_pre_all[R_psd,R_psd,idx_ttt]))
x_fil_up = x0*x_fil_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_fil_all[R_psd,R_psd,idx_ttt]))
x_fil_do = x0*x_fil_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_fil_all[R_psd,R_psd,idx_ttt]))
x_smo_up = x0*x_smo_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_smo_all[R_psd,R_psd,idx_ttt]))
x_smo_do = x0*x_smo_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_smo_all[R_psd,R_psd,idx_ttt]))
fill_between(1.0e9diameter,x_fil_do./log10(cst_r),x_fil_up./log10(cst_r),alpha=0.5,color=:blue)
fill_between(1.0e9diameter,x_smo_do./log10(cst_r),x_smo_up./log10(cst_r),alpha=0.5,color=:orange)
ax = gca();
ax.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
xlim(1.0e9diameter[1],1.0e9diameter[end])
ylim(-200.0)
xlabel("diameter [nm]")
ylabel("density [\$\\frac{\\mathrm{d} N}{\\mathrm{d} \\log_{10}(D_p)}\$]")
title("size distribution at time: 1h")
if GT_loaded
    legend(["filter", "smoother","ground truth","EKF uncertainty","FIKS uncertainty"])
else
    legend(["filter", "smoother","EKF uncertainty","FIKS uncertainty"])
end
savefig(string(folder,"log_size_density_time_1h_smo.png"))
savefig(string(folder,"log_size_density_time_1h_smo.pdf"))
if SAVE_PGF
    savefig(string(folder,"log_size_density_time_1h_smo.pgf"))
end







idx_ttt = 46 # 120s
# idx_ttt = 61 # 600s
idx_ttt_simu = 46
# for idx_ttt in 180:182
    figure()
    semilogx(1.0e9diameter,1.0e-9x0*x_fil_all[R_psd,idx_ttt]./delta)
    semilogx(1.0e9diameter,1.0e-9x0*x_smo_all[R_psd,idx_ttt]./delta)
    # semilogx(1.0e9diameter,1.0e-9y_all[:,idx_ttt]./delta)
    if GT_loaded
        semilogx(1.0e9dp,PSD_simu[:,idx_ttt_simu])
    end
    x_pre_up = x0*x_pre_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_pre_all[R_psd,R_psd,idx_ttt]))
    x_pre_do = x0*x_pre_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_pre_all[R_psd,R_psd,idx_ttt]))
    x_fil_up = x0*x_fil_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_fil_all[R_psd,R_psd,idx_ttt]))
    x_fil_do = x0*x_fil_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_fil_all[R_psd,R_psd,idx_ttt]))
    x_smo_up = x0*x_smo_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_smo_all[R_psd,R_psd,idx_ttt]))
    x_smo_do = x0*x_smo_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_smo_all[R_psd,R_psd,idx_ttt]))
    fill_between(1.0e9diameter,1.0e-9x_fil_do./delta,1.0e-9x_fil_up./delta,alpha=0.5,color=:blue)
    fill_between(1.0e9diameter,1.0e-9x_smo_do./delta,1.0e-9x_smo_up./delta,alpha=0.5,color=:orange)
    ax = gca();
    ax.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
    xlim(1.0e9diameter[1],1.0e9diameter[end])
    ylim(-100.0)
    xlabel("diameter [nm]")
    ylabel("density [cm\$^{-3}\$ nm\$^{-1}\$]")
    title("size deistribution at time: 1h30min")
    if GT_loaded
        legend(["filter", "smoother","ground truth","EKF uncertainty","FIKS uncertainty"])
    else
        legend(["filter", "smoother","EKF uncertainty","FIKS uncertainty"])
    end
    savefig(string(folder,"size_density_time_1h30min_smo.png"))
    savefig(string(folder,"size_density_time_1h30min_smo.pdf"))
    if SAVE_PGF
        savefig(string(folder,"size_density_time_1h30min_smo.pgf"))
    end
    # formatter.set_powerlimits((-1,2))
# end


figure()
semilogx(1.0e9diameter,x0*x_fil_all[R_psd,idx_ttt]./log10(cst_r))
semilogx(1.0e9diameter,x0*x_smo_all[R_psd,idx_ttt]./log10(cst_r))
# semilogx(1.0e9diameter,1.0e-9y_all[:,idx_ttt]./delta)
if GT_loaded
    semilogx(1.0e9dp,1.0e9PSD_simu[:,idx_ttt_simu].*delta_p/log10(cst_r_p))
end
x_pre_up = x0*x_pre_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_pre_all[R_psd,R_psd,idx_ttt]))
x_pre_do = x0*x_pre_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_pre_all[R_psd,R_psd,idx_ttt]))
x_fil_up = x0*x_fil_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_fil_all[R_psd,R_psd,idx_ttt]))
x_fil_do = x0*x_fil_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_fil_all[R_psd,R_psd,idx_ttt]))
x_smo_up = x0*x_smo_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_smo_all[R_psd,R_psd,idx_ttt]))
x_smo_do = x0*x_smo_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_smo_all[R_psd,R_psd,idx_ttt]))
fill_between(1.0e9diameter,x_fil_do./log10(cst_r),x_fil_up./log10(cst_r),alpha=0.5,color=:blue)
fill_between(1.0e9diameter,x_smo_do./log10(cst_r),x_smo_up./log10(cst_r),alpha=0.5,color=:orange)
ax = gca();
ax.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
xlim(1.0e9diameter[1],1.0e9diameter[end])
ylim(-200.0)
xlabel("diameter [nm]")
ylabel("density [\$\\frac{\\mathrm{d} N}{\\mathrm{d} \\log_{10}(D_p)}\$]")
title("size distribution at time: 1h30min")
if GT_loaded
    legend(["filter", "smoother","ground truth","EKF uncertainty","FIKS uncertainty"])
else
    legend(["filter", "smoother","EKF uncertainty","FIKS uncertainty"])
end
savefig(string(folder,"log_size_density_time_1h30min_smo.png"))
savefig(string(folder,"log_size_density_time_1h30min_smo.pdf"))
if SAVE_PGF
    savefig(string(folder,"log_size_density_time_1h30min_smo.pgf"))
end









idx_ttt = 61 # 120s
# idx_ttt = 61 # 600s
idx_ttt_simu = 61
# for idx_ttt in 180:182
    figure()
    semilogx(1.0e9diameter,1.0e-9x0*x_fil_all[R_psd,idx_ttt]./delta)
    semilogx(1.0e9diameter,1.0e-9x0*x_smo_all[R_psd,idx_ttt]./delta)
    # semilogx(1.0e9diameter,1.0e-9y_all[:,idx_ttt]./delta)
    if GT_loaded
        semilogx(1.0e9dp,PSD_simu[:,idx_ttt_simu])
    end
    x_pre_up = x0*x_pre_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_pre_all[R_psd,R_psd,idx_ttt]))
    x_pre_do = x0*x_pre_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_pre_all[R_psd,R_psd,idx_ttt]))
    x_fil_up = x0*x_fil_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_fil_all[R_psd,R_psd,idx_ttt]))
    x_fil_do = x0*x_fil_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_fil_all[R_psd,R_psd,idx_ttt]))
    x_smo_up = x0*x_smo_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_smo_all[R_psd,R_psd,idx_ttt]))
    x_smo_do = x0*x_smo_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_smo_all[R_psd,R_psd,idx_ttt]))
    fill_between(1.0e9diameter,1.0e-9x_fil_do./delta,1.0e-9x_fil_up./delta,alpha=0.5,color=:blue)
    fill_between(1.0e9diameter,1.0e-9x_smo_do./delta,1.0e-9x_smo_up./delta,alpha=0.5,color=:orange)
    ax = gca();
    ax.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
    xlim(1.0e9diameter[1],1.0e9diameter[end])
    ylim(-100.0)
    xlabel("diameter [nm]")
    ylabel("density [cm\$^{-3}\$ nm\$^{-1}\$]")
    title("size distribution at time: 2h")
    if GT_loaded
        legend(["filter", "smoother","ground truth","EKF uncertainty","FIKS uncertainty"])
    else
        legend(["filter", "smoother","EKF uncertainty","FIKS uncertainty"])
    end
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
if GT_loaded
    semilogx(1.0e9dp,1.0e9PSD_simu[:,idx_ttt_simu].*delta_p/log10(cst_r_p))
end
x_pre_up = x0*x_pre_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_pre_all[R_psd,R_psd,idx_ttt]))
x_pre_do = x0*x_pre_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_pre_all[R_psd,R_psd,idx_ttt]))
x_fil_up = x0*x_fil_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_fil_all[R_psd,R_psd,idx_ttt]))
x_fil_do = x0*x_fil_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_fil_all[R_psd,R_psd,idx_ttt]))
x_smo_up = x0*x_smo_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_smo_all[R_psd,R_psd,idx_ttt]))
x_smo_do = x0*x_smo_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_smo_all[R_psd,R_psd,idx_ttt]))
fill_between(1.0e9diameter,x_fil_do./log10(cst_r),x_fil_up./log10(cst_r),alpha=0.5,color=:blue)
fill_between(1.0e9diameter,x_smo_do./log10(cst_r),x_smo_up./log10(cst_r),alpha=0.5,color=:orange)
ax = gca();
ax.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
xlim(1.0e9diameter[1],1.0e9diameter[end])
ylim(-200.0)
xlabel("diameter [nm]")
ylabel("density [\$\\frac{\\mathrm{d} N}{\\mathrm{d} \\log_{10}(D_p)}\$]")
title("size distribution at time: 2h")
if GT_loaded
    legend(["filter", "smoother","ground truth","EKF uncertainty","FIKS uncertainty"])
else
    legend(["filter", "smoother","EKF uncertainty","FIKS uncertainty"])
end
savefig(string(folder,"log_size_density_time_2h_smo.png"))
savefig(string(folder,"log_size_density_time_2h_smo.pdf"))
if SAVE_PGF
    savefig(string(folder,"log_size_density_time_2h_smo.pgf"))
end
