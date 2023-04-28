maxPSDlog = x0*maximum(x_smo_all[R_psd,:]./log10(cst_r))
minPSDlog = max(0.001maxPSDlog,x0*minimum(x_smo_all[R_psd,:]./log10(cst_r)))


VERTICAL = true
if VERTICAL
    _sub = 311
else
    _sub = 131
end
# density reconstruction
rc("ytick",color="black")
fig1, ax2, cb2, pcm2 =  displayLogData2D(1,t_samp[1:pad_factor+1:end]/3600.0,1.0e9diameter,x0*x_smo_all[R_psd,1:pad_factor+1:end]./log10(cst_r),minPSDlog,maxPSDlog,_title="",_colorbar_label="",_sub=_sub)
ax2.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
ax2.yaxis.set_minor_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
cb2.remove()
ylabel("diameter [nm]")
rc("ytick",color="white")
# cax2 = fig1.add_axes([0.58, .82, 0.03, 0.15])
cax2 = fig1.add_axes([0.135, .82, 0.03, 0.15])
cb22 = fig1.colorbar(pcm2, orientation="vertical", cax=cax2, shrink=0.6)
cb22.set_label("distribution [cm\$^{-3}\$]", color="white")
cb22.ax.yaxis.set_tick_params(color="white")
cb22.outline.set_edgecolor("white")




if VERTICAL
    _sub = 312
else
    _sub = 132
end
# nucleation rate
rc("ytick",color="black")
ax5 = subplot(_sub)
plot(t_samp/3600.0,J0*Nucleation_rate(x_smo_all[R_nuc_init,:]),color=:darkorange)
fill_between(t_samp/3600.0,J0*percentiles_nuc_smo[:,1],J0*percentiles_nuc_smo[:,2],alpha=0.5,color=:darkorange)
plot(J_time/3600.0,J_est)
nuc_th = 5.0
nuc_rel = 0.3
J_est_low = (1.0-nuc_rel)*J_est.-nuc_th
J_est_low = J_est_low.*(J_est_low.>0.0)
J_est_up  = (1.0+nuc_rel)*J_est.+nuc_th
fill_between(J_time/3600.0,J_est_low,J_est_up,alpha=0.5)
xlabel("time [h]")
ylabel("nucleation rate [cm\$^{-3}\$ s\$^{-1}\$]")
grid(true,which="major",ls="-")
grid(true,which="minor",ls="-",alpha=0.5)
xlim(t_samp[1]/3600.0,t_samp[end]/3600.0)
ylim(-5.0,170.0)
legend(["FIKS estimate","PSM estimate","FIKS uncertainty","PSM uncertainty"])
# ax5.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.1f"))
# ax5.yaxis.set_ticklabels([""; "0"; ""; "1"; "";  "2"; "";  "3"],color="black") #WARNING: check that it's actually the right values




if VERTICAL
    _sub = 313
else
    _sub = 133
end
# growth rate
ax3 = subplot(_sub)
rc("ytick",color="black")

idx_k = 7
idx_time_gr = findfirst(t_samp.>=(3600.0*0+16*60.0))
semilogx(1.0e9diameter,GR0*CGR(x_smo_all[R_cond,idx_time_gr])/2.77778e-13,color=:orange)
idx_gr_max = findfirst(GR_dp.>=3.2e-9)-1
semilogx(1.0e9GR_dp[1:idx_gr_max],GR_est[1:idx_gr_max],color=:green)
fill_between(1.0e9diameter,GR0*percentiles_cond_smo[:,idx_time_gr,1]/2.77778e-13,GR0*percentiles_cond_smo[:,idx_time_gr,2]/2.77778e-13,alpha=0.5,color=:orange)
xlabel("diameter [nm]")
ylabel("growth rate [nm h\$^{-1}\$]")
grid(true,which="major",ls="-")
grid(true,which="minor",ls="-",alpha=0.5)
xlim(1.0e9diameter[1],1.0e9diameter[end])
legend(["FIKS estimate","INSIDE estimate","FIKS uncertainty"])
#
# if true
#     idx_k = 7
#     idx_time_gr = findfirst(t_samp.>=(3600.0*2+30*60.0))
#     semilogx(1.0e9diameter,GR0*CGR(x_smo_all[R_cond,idx_time_gr])/2.77778e-13,color=:orange) # semilogx # loglog
#     idx_gr_max = findfirst(GR_dp.>=4.3e-9)-1
#     semilogx(1.0e9GR_dp[1:idx_gr_max],GR_est[1:idx_gr_max],color=:green)
#     fill_between(1.0e9diameter,GR0*percentiles_cond_smo[:,idx_time_gr,1]/2.77778e-13,GR0*percentiles_cond_smo[:,idx_time_gr,2]/2.77778e-13,alpha=0.5,color=:orange)
# else
#     idx_k = 8
#     idx_time_gr = findfirst(t_samp.>=(3600.0*3+30*60.0))
#     semilogx(1.0e9diameter,GR0*CGR(x_smo_all[R_cond,idx_time_gr])/2.77778e-13,color=:orange)
#     idx_gr_max = findfirst(GR_dp_2.>=6.2e-9)-1
#     semilogx(1.0e9GR_dp_2[1:idx_gr_max],GR_est_2[1:idx_gr_max],color=:green)
#     fill_between(1.0e9diameter,GR0*percentiles_cond_smo[:,idx_time_gr,1]/2.77778e-13,GR0*percentiles_cond_smo[:,idx_time_gr,2]/2.77778e-13,alpha=0.5,color=:orange)
# end
# xlabel("diameter [nm]")
# ylabel("growth rate [nm h\$^{-1}\$]")
# grid(true,which="major",ls="-")
# grid(true,which="minor",ls="-",alpha=0.5)
# xlim(1.0e9diameter[1],1.0e9diameter[end])
# legend(["FIKS estimate","INSIDE estimate","FIKS uncertainty"])
#











fig1.set_figwidth(5.3333) # (5.6)
fig1.set_figheight(9.0) # (10.0)

tight_layout(pad=1.0, w_pad=0.5, h_pad=0.5)
tight_layout(pad=1.0, w_pad=0.5, h_pad=0.5)


ax3.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(-0.085, 2.19), textcoords="axes fraction", color="black")
ax3.annotate("c)", xy=(3, 1),  xycoords="data", xytext=(-0.085, 0.99), textcoords="axes fraction", color="black")
ax3.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-0.085, 3.39), textcoords="axes fraction", color="black")


# savefig(string(folder,"1906_03_estimation_weak_short.pdf"))
# savefig(string(folder,"1906_03_estimation_weak_short.png"))

# savefig(string(folder,"1906_03_estimation_weak_long.pdf"))
# savefig(string(folder,"1906_03_estimation_weak_long.png"))

# savefig(string(folder,"1906_03_estimation_strong_very_short.pdf"))
# savefig(string(folder,"1906_03_estimation_strong_very_short.png"))

# savefig(string(folder,"1906_03_estimation_strong_short.pdf"))
# savefig(string(folder,"1906_03_estimation_strong_short.png"))

# savefig(string(folder,"1906_03_estimation_strong_long.pdf"))
# savefig(string(folder,"1906_03_estimation_strong_long.png"))

# savefig(string(folder,"1906_03_estimation_strong_shorter_long.pdf"))
# savefig(string(folder,"1906_03_estimation_strong_shorter_long.png"))






# titleCond = "condensation rate"
# cblabelCond = "growth rate [nm h\$^{-1}\$]"
# minCond = max(5.0e-14,minimum(GR0*CGR2D(x_smo_all[R_cond,:])))/2.77778e-13 # max(1.0e-14,minimum(GR0*CGR(x_fil_all[R_cond,:]))) # max(1.0e-14,minimum(GR0*CGR(x_smo_all[R_cond,:])))
# maxCond = maximum(GR0*CGR2D(x_smo_all[R_cond,:]))/2.77778e-13
# displayLogData2D(456,t_samp[1:pad_factor+1:end]/3600.0,1.0e9diameter,GR0*CGR2D(x_smo_all[R_cond,1:pad_factor+1:end])/2.77778e-13,minCond,maxCond,_title=titleCond,_colorbar_label=cblabelCond)
# ylabel("diameter [nm]")
# savefig(string(folder,"1906_03_estimation_weak_short_cond_smo.png"))
# savefig(string(folder,"1906_03_estimation_weak_short_cond_smo.pdf"))
