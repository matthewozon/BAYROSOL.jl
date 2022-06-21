# plotting transfer function in one figure #WARNING: annotation good for 1952 model

rc("ytick",color="black")
min_H_full,max_H_full = extrema(H_DMATRAIN_FULL)
s = @sprintf "kernels (%1.2e,%1.2e)" min_H_full max_H_full
# if FLAG_0000_04
#     _sub=211
# else
_sub=311
# end
fig1, ax1, cb1, pcm1 =  displayLogData2D(20,1.0e9diameter_model_full[idx_s_0:idx_s_z]/sqrt(cst_r_raw),collect(0.5:6.5),H_DMATRAIN_FULL[:,idx_s_0:idx_s_z],max(min_H_full,0.001max_H_full),max_H_full,_title=s,_colorbar_label="transfer function [cm\$^{3}\$]",_sub=_sub)
cb1.remove()
xscale("log")
yscale("linear")
ylabel("channel")
title("")
# tight_layout()
ax1.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-0.08, 0.98), textcoords="axes fraction")
# if FLAG_0000_04
#     fig1.set_figwidth(5.6) # (5.625)
#     fig1.set_figheight(6.0) # (10.0)
# else
fig1.set_figwidth(5.6) # (5.625)
fig1.set_figheight(9.0) # (10.0)
# end
xlabel("diameter [nm]")
# xlabel("")
rc("ytick",color="white")
# if FLAG_0000_04
#     cax1 = fig1.add_axes([0.12, .711, 0.03, 0.25])
# else
cax1 = fig1.add_axes([0.12, .815, 0.03, 0.15])
# end
cb11 = fig1.colorbar(pcm1, orientation="vertical", cax=cax1, shrink=0.6)
cb11.set_label("transfer function [cm\$^{3}\$]", color="white")
cb11.ax.yaxis.set_tick_params(color="white")
cb11.outline.set_edgecolor("white")
# cb1.update_ticks()

rc("ytick",color="black")
min_H,max_H = extrema(H_DMATRAIN)
s = @sprintf "kernels (%1.2e,%1.2e)" min_H max_H
# if FLAG_0000_04
#     _sub=212
# else
_sub=312
# end
fig2, ax2, cb2, pcm2 =  displayLogData2D(20,1.0e9diameter_model/sqrt(cst_r),collect(0.5:6.5),H_DMATRAIN,max(min_H_full,0.001max_H_full),max_H_full,_title=s,_colorbar_label="volume [cm\$^{3}\$]",_sub=_sub)
cb2.remove()
xscale("log")
yscale("linear")
xlabel("diameter [nm]")
# xlabel("")
ylabel("channel")
title("")
# tight_layout()
tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
ax2.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(-0.08, 0.98), textcoords="axes fraction")
rc("ytick",color="white")
# if FLAG_0000_04
#     cax2 = fig2.add_axes([0.12, .21, 0.03, 0.25])
# else
cax2 = fig2.add_axes([0.12, .475, 0.03, 0.15])
# end
cb22 = fig2.colorbar(pcm2, orientation="vertical", cax=cax2, shrink=0.6)
cb22.set_label("transfer function [cm\$^{3}\$]", color="white")
cb22.ax.yaxis.set_tick_params(color="white")
cb22.outline.set_edgecolor("white")

# if !FLAG_0000_04
rc("ytick",color="black")
# ax3 = subplot(313)
ax3 = subplot(3,1,3)
# plot(1.0e9diameter_model_full,H_DMATRAIN_FULL[2,:],color=:green)
semilogx(1.0e9diameter_model_full,H_DMATRAIN_FULL[2,:],color=:green)
scatter(1.0e9diameter_model,H_DMATRAIN[2,:],color=:darkorange)
fill_between(1.0e9diameter_model_full,0.8*H_DMATRAIN_FULL[2,:]-dH_FULL[2,:],1.2*H_DMATRAIN_FULL[2,:]+dH_FULL[2,:],alpha=0.5,color=:darkorange)
xlim(1.7,3.5)
# xlabel("diameter [nm]")
# ylabel("kernel gain [cm\$^{3}\$]")
# title("kernel (channel 2)")
legend(["fine model";"coarse model";"model uncertainty"])
ax2.annotate("c)", xy=(3, 1),  xycoords="data", xytext=(-0.08, -0.175), textcoords="axes fraction")
ax3.text(1.72, 0.27, "kernel function", rotation=0)
# ax3.set_aspect(2.75)
ax3.text(2.2, -0.04, "diameter [nm]", rotation=0)
# end
# tight_layout()
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
# ax3.text(2.5, 0.1, "text 0", rotation=0)



savefig(string(folder,"transfer_function.png"))
savefig(string(folder,"transfer_function.pdf"))
