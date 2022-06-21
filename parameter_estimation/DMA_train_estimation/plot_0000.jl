VERTICAL = false
if VERTICAL
    _sub = 323
else
    _sub = 233
end

# growth rate contour plot
ax3 = subplot(_sub)
rc("ytick",color="black")
minCond = max(5.0e-14,minimum(GR0*CGR2D(x_smo_all[R_cond,:])))/2.77778e-13
maxCond = maximum(GR0*CGR2D(x_smo_all[R_cond,:]))/2.77778e-13
fig1, ax3, pcm3 = imshowData(1,t_samp[1:pad_factor+1:end]/3600.0,1.0e9diameter,GR0*CGR2D(x_smo_all[R_cond,1:pad_factor+1:end])/2.77778e-13,_norm=:Normalize,_vmin=0.0,_vmax=maxCond,_edgecolors="face",_sub=_sub)
yscale("log")
xlabel("time [h]",fontsize=12)
ylabel("diameter [nm]",fontsize=12)
ax3.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
ax3.yaxis.set_minor_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
ax3.yaxis.set_tick_params(color="black")
rc("ytick",color="white")
# cax3 = fig1.add_axes([0.11, .48, 0.03, 0.15])
# cb33 = fig1.colorbar(pcm3, orientation="vertical", cax=cax3, shrink=0.6)
if VERTICAL
    cax3 = fig1.add_axes([0.11, .48, 0.03, 0.15])
    cb33 = fig1.colorbar(pcm3, orientation="vertical", cax=cax3, shrink=0.6)
else
    cax3 = fig1.add_axes([0.71, .75, 0.01, 0.2])
    cb33 = fig1.colorbar(pcm3, orientation="vertical", cax=cax3, shrink=0.6)
end
cb33.set_label("growth rate [nm h\$^{-1}\$]", color="white")
cb33.ax.yaxis.set_tick_params(color="white")
cb33.outline.set_edgecolor("white")
# ax3.yaxis.set_ticklabels(["2"; "3"; "4"; "5"; "6";  "7"; "8";  "9"; "10"],color="black")



if VERTICAL
    _sub = 321
else
    _sub = 231
end
# density simulation
rc("ytick",color="black")
cst_r_simu = mean(dp[2:end]./dp[1:end-1]);
delta_simu = dp*(sqrt(cst_r_simu)-1.0/sqrt(cst_r_simu));
psd_log_simu = 1.0e9psd_simu.*delta_simu/log10(cst_r_simu);
maxPSDlog = maximum(psd_log_simu)
minPSDlog = max(0.001maxPSDlog,minimum(psd_log_simu))
fig1, ax1, cb1, pcm1 =  displayLogData2D(1,tp[1:50:end],1.0e9dp[idx_p_simu],psd_log_simu[idx_p_simu,1:50:end],minPSDlog,maxPSDlog,_title="",_colorbar_label="",_sub=_sub)
ax3.get_shared_y_axes().join(ax3,ax1)
ax1.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-0.08, 0.98), textcoords="axes fraction")
ax1.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
ax1.yaxis.set_minor_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
# ax1.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(0.5, 3.5), textcoords="axes fraction")
# A = matplotlib.ticker.ScalarFormatter()
# A.set_powerlimits((0,0))
# ax1.ticklabel_format(self, *, axis='both', style='', scilimits=None, useOffset=None, useLocale=None, useMathText=None)
cb1.remove()
# xscale("log")
# xlabel("diameter [nm]")
# yscale("linear")
ylabel("diameter [nm]",fontsize=12)
xlabel("time [h]",fontsize=12)
# title("")

# fig1.set_figwidth(5.6) # (5.625)
# fig1.set_figheight(9.0) # (10.0)

rc("ytick",color="white")
# cax1 = fig1.add_axes([0.11, .82, 0.03, 0.15])
# cb11 = fig1.colorbar(pcm1, orientation="vertical", cax=cax1, shrink=0.6)
if VERTICAL
    cax1 = fig1.add_axes([0.11, .82, 0.03, 0.15])
    cb11 = fig1.colorbar(pcm1, orientation="vertical", cax=cax1, shrink=0.6)
else
    cax1 = fig1.add_axes([0.045, .62, 0.01, 0.33])
    cb11 = fig1.colorbar(pcm1, orientation="vertical", cax=cax1, shrink=0.6)
end
cb11.set_label("distribution [cm\$^{-3}\$]", color="white")
cb11.ax.yaxis.set_tick_params(color="white")
cb11.outline.set_edgecolor("white")
# cb1.update_ticks()
# tight_layout()



if VERTICAL
    _sub = 322
else
    _sub = 232
end
# density reconstruction
rc("ytick",color="black")
fig1, ax2, cb2, pcm2 =  displayLogData2D(1,t_samp[1:pad_factor+1:end]/3600.0,1.0e9diameter,x0*x_smo_all[R_psd,1:pad_factor+1:end]./log10(cst_r),minPSDlog,maxPSDlog,_title="",_colorbar_label="",_sub=_sub)
ax2.get_shared_y_axes().join(ax3,ax2)
ax2.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
ax2.yaxis.set_minor_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
cb2.remove()
ylabel("diameter [nm]",fontsize=12)
xlabel("time [h]",fontsize=12)
rc("ytick",color="white")
# cax2 = fig1.add_axes([0.58, .82, 0.03, 0.15])
# cax2 = fig1.add_axes([0.6, .82, 0.03, 0.15])
# cb22 = fig1.colorbar(pcm2, orientation="vertical", cax=cax2, shrink=0.6)
if VERTICAL
    cax2 = fig1.add_axes([0.6, .82, 0.03, 0.15])
    cb22 = fig1.colorbar(pcm2, orientation="vertical", cax=cax2, shrink=0.6)
else
    cax2 = fig1.add_axes([0.375, .62, 0.01, 0.33]) # [0.67, .40, 0.01, 0.5])
    cb22 = fig1.colorbar(pcm2, orientation="vertical", cax=cax2, shrink=0.6)
end
cb22.set_label("distribution [cm\$^{-3}\$]", color="white")
cb22.ax.yaxis.set_tick_params(color="white")
cb22.outline.set_edgecolor("white")





if VERTICAL
    _sub = 324
else
    _sub = 234
end
# growth rate: slice
idx_time = [1; round(Int64,n_samp/6); round(Int64,3n_samp/12); round(Int64,2n_samp/6); round(Int64,5n_samp/12); round(Int64,3n_samp/6); round(Int64,4n_samp/6); round(Int64,5n_samp/6); round(Int64,11n_samp/12); round(Int64,23n_samp/24); n_samp]
rc("ytick",color="black")
ax4 = subplot(_sub)
k = 9 # 8
semilogx(1.0e9diameter,GR0*CGR(x_smo_all[R_cond,idx_time[k]])/2.77778e-13,color=:orange) # semilogx # loglog
if GT_loaded
    idx_simu = findfirst(3600.0tp.>=t_samp[idx_time[k]])
    semilogx(1.0e9dp,condensation_rate_all[:,idx_simu]/2.77778e-13,color=:green)         # semilogx # loglog
end
fill_between(1.0e9diameter,GR0*percentiles_cond_smo[:,idx_time[k],1]/2.77778e-13,GR0*percentiles_cond_smo[:,idx_time[k],2]/2.77778e-13,alpha=0.5,color=:orange)
xlabel("diameter [nm]",fontsize=12)
ylabel("growth rate [nm h\$^{-1}\$]",fontsize=12)
grid(true,which="major",ls="-")
grid(true,which="minor",ls="-",alpha=0.5)
xlim(1.0e9diameter[1],1.0e9diameter[end])
ylim(-0.05,5.05)
if GT_loaded
    legend(["FIKS estimate", "ground truth", "FIKS uncertainty"],fontsize=12)
else
    legend(["FIKS estimate", "FIKS uncertainty"],fontsize=12)
end
ax4.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
ax4.yaxis.set_minor_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
ax4.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
ax4.xaxis.set_minor_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))



if VERTICAL
    _sub = 325
else
    _sub = 235
end
# nucleation rate
rc("ytick",color="black")
ax5 = subplot(_sub)
plot(t_samp/3600.0,J0*Nucleation_rate(x_smo_all[R_nuc_init,:]),color=:darkorange)
fill_between(t_samp/3600.0,J0*percentiles_nuc_smo[:,1],J0*percentiles_nuc_smo[:,2],alpha=0.5,color=:darkorange)
if GT_loaded
    plot(tp,1.0e9condensation_rate_all[idx_s_0,:].*psd_simu[idx_s_0,:],color=:green)
    # plot(tp,nucleation_rate_all,color=:green)
else
    plot(J_time/3600.0,J_est)
end
xlabel("time [h]",fontsize=12)
ylabel("nucleation rate [cm\$^{-3}\$ s\$^{-1}\$]",fontsize=12)
grid(true,which="major",ls="-")
grid(true,which="minor",ls="-",alpha=0.5)
xlim(t_samp[1]/3600.0,t_samp[end]/3600.0)
# ylim(-0.2,55.0)
if GT_loaded
    legend(["smoother","ground truth","FIKS uncertainty"],fontsize=12)
else
    # legend(["filter","smoother","PSM: 1.7 nm","EKF uncertainty","FIKS uncertainty"]) # ,"particle flux"
    legend(["FIKS estimate","PSM estimate","FIKS uncertainty"],fontsize=12)
end
ax5.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.1f"))
ax5.yaxis.set_ticklabels([""; "0"; ""; "1"; "";  "2"; "";  "3"],color="black") #WARNING: check that it's actually the right values


if VERTICAL
    _sub = 326
else
    _sub = 236
end
# size distribution slice
rc("ytick",color="black")
ax6 = subplot(_sub)
idx_ttt = idx_time[k] #  241
semilogx(1.0e9diameter,x0*x_smo_all[R_psd,idx_ttt]./log10(cst_r),color=:orange)
if GT_loaded
    idx_t_simu = findfirst(3600.0tp.>t_samp[idx_ttt])
    semilogx(1.0e9dp[1:10:end],psd_log_simu[1:10:end,idx_t_simu],color=:green)
    # scatter(dp[idx_p_simu],psd_log_simu[idx_p_simu,idx_t_simu],color=:green)
    # psd_log_mean_plot = zeros(Cdouble,nbin)
    # for i in 1:nbin
    #     psd_log_mean_plot[i] = mean(psd_log_simu[idx_p_simu[i]-24:idx_p_simu[i]+23,idx_t_simu])
    # end
    # scatter(1.0e9dp[idx_p_simu],psd_log_mean_plot,color=:green)
end
x_smo_up = (x0*x_smo_all[R_psd,idx_ttt] + x0*sqrt.(diag(o_smo_all[R_psd,R_psd,idx_ttt])))./log10(cst_r)
x_smo_do = (x0*x_smo_all[R_psd,idx_ttt] - x0*sqrt.(diag(o_smo_all[R_psd,R_psd,idx_ttt])))./log10(cst_r)
fill_between(1.0e9diameter,x_smo_do,x_smo_up,alpha=0.5,color=:orange)
# ax = gca();
xlim(1.0e9diameter[1],1.0e9diameter[end])
ylim(-1.0,2.0e4)
xlabel("diameter [nm]",fontsize=12)
# ylabel("density \$\\frac{dN}{d\\log D_p}\$ [cm\$^{-3}\$]")
ylabel("distribution [cm\$^{-3}\$]",fontsize=12)
grid(true,which="major",ls="-")
grid(true,which="minor",ls="-",alpha=0.5)
if GT_loaded
    # legend(["FIKS estimate","simulation","simulation mean","FIKS uncertainty"])
    legend(["FIKS estimate","simulation","FIKS uncertainty"],fontsize=12)
else
    legend(["FIKS estimate","FIKS uncertainty"],fontsize=12)
end
ax6.ticklabel_format(style="scientific",axis="y",scilimits=(0,0))
# ax6.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
# ax6.yaxis.set_minor_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
ax6.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
ax6.xaxis.set_minor_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))





# I don't know why!
# cb22.outline.set_edgecolor("white")

# fig1.set_figwidth(5.6) # (5.625)
# fig1.set_figheight(9.0) # (10.0)
if VERTICAL
    fig1.set_figwidth(5.6) # (5.6)
    fig1.set_figheight(9.0) # (10.0)
else
    fig1.set_figwidth(15.0) # (5.6)
    # fig1.set_figheight(6.0) # (10.0)
    fig1.set_figheight(7.0) # (10.0)
end

# tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

tight_layout(pad=1.0, w_pad=0.5, h_pad=0.5)
tight_layout(pad=1.0, w_pad=0.5, h_pad=0.5)

if VERTICAL
    ax5.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-0.18, 3.5), textcoords="axes fraction", color="black")
    ax5.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(1.07, 3.5), textcoords="axes fraction", color="black")
    ax5.annotate("c)", xy=(3, 1),  xycoords="data", xytext=(-0.18, 2.24), textcoords="axes fraction", color="black")
    ax5.annotate("d)", xy=(3, 1),  xycoords="data", xytext=(1.07, 2.24), textcoords="axes fraction", color="black")
    ax5.annotate("e)", xy=(3, 1),  xycoords="data", xytext=(-0.18, 0.98), textcoords="axes fraction", color="black")
    ax5.annotate("f)", xy=(3, 1),  xycoords="data", xytext=(1.07, 0.98), textcoords="axes fraction", color="black")
else
    ax5.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-1.25, 2.22), textcoords="axes fraction", color="black",fontsize=15)
    ax5.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(-0.10, 2.22), textcoords="axes fraction", color="black",fontsize=15)
    ax5.annotate("c)", xy=(3, 1),  xycoords="data", xytext=(1.03, 2.22), textcoords="axes fraction", color="black",fontsize=15)
    ax5.annotate("d)", xy=(3, 1),  xycoords="data", xytext=(-1.25, 1.0), textcoords="axes fraction", color="black",fontsize=15)
    ax5.annotate("e)", xy=(3, 1),  xycoords="data", xytext=(-0.10, 1.0), textcoords="axes fraction", color="black",fontsize=15)
    ax5.annotate("f)", xy=(3, 1),  xycoords="data", xytext=(1.03, 1.0), textcoords="axes fraction", color="black",fontsize=15)

    # time stamps
    s = @sprintf "time %ih%imin" floor(t_samp[idx_time[k]]/3600) 60*(t_samp[idx_time[k]]/3600.0-floor(t_samp[idx_time[k]]/3600))
    ax5.annotate(s, xy=(3, 1),  xycoords="data", xytext=(-1.1, 0.875), textcoords="axes fraction", color="black",fontsize=12)
    s = @sprintf "time %ih%imin" floor(t_samp[idx_ttt]/3600) 60*(t_samp[idx_ttt]/3600.0-floor(t_samp[idx_ttt]/3600))
    ax6.annotate(s, xy=(3, 1),  xycoords="data", xytext=(0.1, 0.875), textcoords="axes fraction", color="black",fontsize=12)
end

# savefig(string(folder,"0000_00_simulation_estimation.pdf"))
# savefig(string(folder,"0000_00_simulation_estimation.png"))

# savefig(string(folder,"0000_04_simulation_estimation.pdf"))
# savefig(string(folder,"0000_04_simulation_estimation.png"))






if FLAG_0000_04
    rc("ytick",color="black")
    min_H_full,max_H_full = extrema(H_DMATRAIN_FULL)
    s = @sprintf "kernels (%1.2e,%1.2e)" min_H_full max_H_full
    _sub=211
    fig2, ax1, cb1, pcm1 =  displayLogData2D(2,1.0e9diameter_model_full[idx_s_0:idx_s_z]/sqrt(cst_r_raw),collect(0.5:6.5),H_DMATRAIN_FULL[:,idx_s_0:idx_s_z],max(min_H_full,0.001max_H_full),max_H_full,_title=s,_colorbar_label="transfer function [cm\$^{3}\$]",_sub=_sub)
    cb1.remove()
    xscale("log")
    yscale("linear")
    ylabel("channel")
    title("")
    # tight_layout()
    ax1.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-0.08, 0.98), textcoords="axes fraction")
    fig2.set_figwidth(5.6) # (5.625)
    fig2.set_figheight(6.0) # (10.0)
    xlabel("diameter [nm]")
    # xlabel("")
    rc("ytick",color="white")
    cax1 = fig2.add_axes([0.12, .711, 0.03, 0.25])
    cb11 = fig2.colorbar(pcm1, orientation="vertical", cax=cax1, shrink=0.6)
    cb11.set_label("transfer function [cm\$^{3}\$]", color="white")
    cb11.ax.yaxis.set_tick_params(color="white")
    cb11.outline.set_edgecolor("white")
    tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
    # cb1.update_ticks()

    # density reconstruction
    rc("ytick",color="black")
    _sub = 212
    fig2, ax2, cb2, pcm2 =  displayLogData2D(2,t_samp[1:pad_factor+1:end]/3600.0,1.0e9diameter,x0*x_smo_all[R_psd,1:pad_factor+1:end]./log10(cst_r),minPSDlog,maxPSDlog,_title="",_colorbar_label="",_sub=_sub)
    ax2.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
    ax2.yaxis.set_minor_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
    cb2.remove()
    ylabel("diameter [nm]")
    rc("ytick",color="white")
    # cax2 = fig1.add_axes([0.58, .82, 0.03, 0.15])
    # cax2 = fig2.add_axes([0.6, .82, 0.03, 0.15])
    cax2 = fig2.add_axes([0.12, .21+0.01, 0.03, 0.25])
    cb22 = fig2.colorbar(pcm2, orientation="vertical", cax=cax2, shrink=0.6)
    cb22.set_label("distribution [cm\$^{-3}\$]", color="white")
    cb22.ax.yaxis.set_tick_params(color="white")
    cb22.outline.set_edgecolor("white")
    tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
    rc("ytick",color="black")
    ax1.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(-0.08, -0.2), textcoords="axes fraction")

    # savefig(string(folder,"0000_04_design_transfer_function.pdf"))
    # savefig(string(folder,"0000_04_design_transfer_function.png"))
end
