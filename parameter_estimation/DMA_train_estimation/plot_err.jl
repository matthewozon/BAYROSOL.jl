# figure(); imshow(1.0e-9ERR_DISC'./delta); colorbar()
# figure(); imshow(1.0e-9ERR_LOSS'./delta); colorbar()
# figure(); imshow(1.0e-9ERR_GROW'./delta); colorbar()

ERR_DISC = abs.(ERR_DISC)
ERR_LOSS = abs.(ERR_LOSS)
ERR_GROW = abs.(ERR_GROW)
maxDISCerr = x0*maximum(ERR_DISC'./log10(cst_r))
minDISCerr = max(0.001maxDISCerr,minimum(ERR_DISC'./log10(cst_r)))

maxLOSSerr = x0*maximum(ERR_LOSS'./log10(cst_r))
minLOSSerr = max(0.001maxLOSSerr,minimum(ERR_LOSS'./log10(cst_r)))

maxGROWerr = x0*maximum(ERR_GROW'./log10(cst_r))
minGROWerr = max(0.001maxGROWerr,minimum(ERR_GROW'./log10(cst_r)))

maxERR = max(maxDISCerr,maxLOSSerr,maxGROWerr)
minERR = min(minDISCerr,minLOSSerr,minGROWerr)

VERTICAL = true
if VERTICAL
    _sub = 311
else
    _sub = 131
end
# density reconstruction
rc("ytick",color="black")
fig1, ax2, cb2, pcm2 =  displayLogData2D(2,t_samp[1:pad_factor+1:end]/3600.0,1.0e9diameter,ERR_DISC'./log10(cst_r),minERR,maxERR,_title="",_colorbar_label="",_sub=_sub)
ax2.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
ax2.yaxis.set_minor_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
cb2.remove()
ylabel("diameter [nm]")
rc("ytick",color="white")
# cax2 = fig1.add_axes([0.58, .82, 0.03, 0.15])
# cb22 = fig1.colorbar(pcm2, orientation="vertical", cax=cax2, shrink=0.6)
# cax2 = fig1.add_axes([0.12, .82, 0.03, 0.15])
if VERTICAL
    cax2 = fig1.add_axes([0.12, .82, 0.03, 0.15])
    cb22 = fig1.colorbar(pcm2, orientation="vertical", cax=cax2, shrink=0.6)
else
    cax2 = fig1.add_axes([0.05, .40, 0.01, 0.5])
    cb22 = fig1.colorbar(pcm2, orientation="vertical", cax=cax2, shrink=0.6)
end
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
fig1, ax5, cb5, pcm5 =  displayLogData2D(2,t_samp[1:pad_factor+1:end]/3600.0,1.0e9diameter,ERR_LOSS'./log10(cst_r),minERR,maxERR,_title="",_colorbar_label="",_sub=_sub)
ax5.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
ax5.yaxis.set_minor_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
cb5.remove()
ylabel("diameter [nm]")
rc("ytick",color="white")
# cax2 = fig1.add_axes([0.58, .82, 0.03, 0.15])
# cb22 = fig1.colorbar(pcm2, orientation="vertical", cax=cax2, shrink=0.6)
# cax2 = fig1.add_axes([0.12, .82, 0.03, 0.15])
if VERTICAL
    cax5 = fig1.add_axes([0.12, .495, 0.03, 0.15])
    cb55 = fig1.colorbar(pcm5, orientation="vertical", cax=cax5, shrink=0.6)
else
    cax5 = fig1.add_axes([0.05, .40, 0.01, 0.5])
    cb55 = fig1.colorbar(pcm5, orientation="vertical", cax=cax5, shrink=0.6)
end
cb55.set_label("distribution [cm\$^{-3}\$]", color="white")
cb55.ax.yaxis.set_tick_params(color="white")
cb55.outline.set_edgecolor("white")

if VERTICAL
    _sub = 313
else
    _sub = 133
end
# nucleation rate
rc("ytick",color="black")
fig1, ax6, cb6, pcm6 =  displayLogData2D(2,t_samp[1:pad_factor+1:end]/3600.0,1.0e9diameter,ERR_GROW'./log10(cst_r),minERR,maxERR,_title="",_colorbar_label="",_sub=_sub)
ax6.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
ax6.yaxis.set_minor_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
cb6.remove()
ylabel("diameter [nm]")
rc("ytick",color="white")
# cax2 = fig1.add_axes([0.58, .82, 0.03, 0.15])
# cb22 = fig1.colorbar(pcm2, orientation="vertical", cax=cax2, shrink=0.6)
# cax2 = fig1.add_axes([0.12, .82, 0.03, 0.15])
if VERTICAL
    cax6 = fig1.add_axes([0.12, .17, 0.03, 0.15])
    cb66 = fig1.colorbar(pcm6, orientation="vertical", cax=cax6, shrink=0.6)
else
    cax6 = fig1.add_axes([0.05, .40, 0.01, 0.5])
    cb66 = fig1.colorbar(pcm6, orientation="vertical", cax=cax6, shrink=0.6)
end
cb66.set_label("distribution [cm\$^{-3}\$]", color="white")
cb66.ax.yaxis.set_tick_params(color="white")
cb66.outline.set_edgecolor("white")






# fig1.set_figwidth(5.3333) # (5.6)
# fig1.set_figheight(9.0) # (10.0)
if VERTICAL
    fig1.set_figwidth(5.3333) # (5.6)
    fig1.set_figheight(9.0) # (10.0)
else
    fig1.set_figwidth(15.0) # (5.6)
    fig1.set_figheight(3.3333) # (10.0)
end

# tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

tight_layout(pad=1.0, w_pad=0.5, h_pad=0.5)
tight_layout(pad=1.0, w_pad=0.5, h_pad=0.5)


# ax5.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-0.085, 2.19), textcoords="axes fraction", color="black")
# ax5.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(-0.085, 0.99), textcoords="axes fraction", color="black")
# ax5.annotate("c)", xy=(3, 1),  xycoords="data", xytext=(-0.085, -0.22), textcoords="axes fraction", color="black")
# rc("ytick",color="black")
# if VERTICAL
#     ax5.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-0.085, 2.19), textcoords="axes fraction", color="black")
#     ax5.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(-0.085, 0.99), textcoords="axes fraction", color="black")
#     ax5.annotate("c)", xy=(3, 1),  xycoords="data", xytext=(-0.085, -0.22), textcoords="axes fraction", color="black")
# else
#     ax5.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-1.21, 0.99), textcoords="axes fraction", color="black")
#     ax5.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(-0.085, 0.99), textcoords="axes fraction", color="black")
#     ax5.annotate("c)", xy=(3, 1),  xycoords="data", xytext=(1.06, 0.99), textcoords="axes fraction", color="black")
# end


# savefig(string(folder,"1952_02_gde_err.pdf"))
# savefig(string(folder,"1952_02_gde_err.png"))

# savefig(string(folder,"1802_01_gde_err.pdf"))
# savefig(string(folder,"1802_01_gde_err.png"))

# savefig(string(folder,"1906_03_gde_err_strong_long.pdf"))
# savefig(string(folder,"1906_03_gde_err_strong_long.png"))

# savefig(string(folder,"1906_03_gde_err_strong_shorter_long.pdf"))
# savefig(string(folder,"1906_03_gde_err_strong_shorter_long.png"))

# savefig(string(folder,"0000_00_gde_err.pdf"))
# savefig(string(folder,"0000_00_gde_err.png"))

# savefig(string(folder,"0000_04_gde_err.pdf"))
# savefig(string(folder,"0000_04_gde_err.png"))
