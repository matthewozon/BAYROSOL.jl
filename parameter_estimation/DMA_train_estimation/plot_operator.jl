
min_H_full,max_H_full = extrema(H_DMATRAIN_FULL)
s = @sprintf "kernels (%1.2e,%1.2e)" min_H_full max_H_full
displayLogData2D(20,1.0e9diameter_model_full,1.0e9diameter_data,H_DMATRAIN_FULL,max(min_H_full,0.001max_H_full),max_H_full,_title=s,_colorbar_label="volume [cm\$^{3}\$]")
xscale("log")
xlabel("diameter [nm]")
ylabel("channel center [nm]")
# savefig(string(folder,"H_DMATRAIN_FULL.png"))
# savefig(string(folder,"H_DMATRAIN_FULL.pdf"))


min_H,max_H = extrema(H_DMATRAIN)
s = @sprintf "kernels (%1.2e,%1.2e)" min_H max_H
displayLogData2D(2,1.0e9diameter_model,1.0e9diameter_data,H_DMATRAIN,max(min_H,0.001max_H),max_H_full,_title=s,_colorbar_label="volume [cm\$^{3}\$]")
xscale("log")
xlabel("diameter [nm]")
ylabel("channel center [nm]")
# savefig(string(folder,"H_DMATRAIN.png"))
# savefig(string(folder,"H_DMATRAIN.pdf"))

min_dH,max_dH = extrema(dH)
s = @sprintf "kernels (%1.2e,%1.2e)" min_dH max_dH
displayLogData2D(23,1.0e9diameter_model_full,1.0e9diameter_data,dH_FULL,max(min_dH,0.001max_dH),max_H_full,_title=s,_colorbar_label="volume [cm\$^{3}\$]")
xscale("log")
xlabel("diameter [nm]")
ylabel("channel center [nm]")
# savefig(string(folder,"dH_DMATRAIN_FULL.png"))
# savefig(string(folder,"dH_DMATRAIN_FULL.pdf"))


min_dH,max_dH = extrema(dH)
s = @sprintf "kernels (%1.2e,%1.2e)" min_dH max_dH
displayLogData2D(21,1.0e9diameter_model,1.0e9diameter_data,dH_SMALL,max(min_dH,0.001max_dH),max_H_full,_title=s,_colorbar_label="volume [cm\$^{3}\$]")
xscale("log")
xlabel("diameter [nm]")
ylabel("channel center [nm]")
# savefig(string(folder,"dH_DMATRAIN.png"))
# savefig(string(folder,"dH_DMATRAIN.pdf"))

figure();
semilogx(1.0e9diameter_model_full,H_DMATRAIN_FULL[1,:],color=:green)
scatter(1.0e9diameter_model,H_DMATRAIN[1,:],color=:darkorange)
fill_between(1.0e9diameter_model_full,0.8*H_DMATRAIN_FULL[1,:]-dH_FULL[1,:],1.2*H_DMATRAIN_FULL[1,:]+dH_FULL[1,:],alpha=0.5,color=:darkorange)
xlim(1.65,2.25)
xlabel("diameter [nm]")
ylabel("gain [cm\$^{3}\$]")
title("kernel (channel 1)")
legend(["fine model";"coarse model";"model uncertainty"])
# savefig(string(folder,"channel_1_uncertainty.png"))
# savefig(string(folder,"channel_1_uncertainty.pdf"))

figure();
semilogx(1.0e9diameter_model_full,H_DMATRAIN_FULL[2,:],color=:green)
# scatter(1.0e9diameter_model*sqrt(cst_r),H_DMATRAIN[2,:],color=:darkorange)
scatter(1.0e9diameter_model,H_DMATRAIN[2,:],color=:darkorange)
fill_between(1.0e9diameter_model_full,0.8*H_DMATRAIN_FULL[2,:]-dH_FULL[2,:],1.2*H_DMATRAIN_FULL[2,:]+dH_FULL[2,:],alpha=0.5,color=:darkorange)
xlim(1.7,3.5)
xlabel("diameter [nm]")
ylabel("gain [cm\$^{3}\$]")
title("kernel (channel 2)")
legend(["fine model";"coarse model";"model uncertainty"])
# savefig(string(folder,"channel_2_uncertainty.png"))
# savefig(string(folder,"channel_2_uncertainty.pdf"))


figure();
semilogx(1.0e9diameter_model_full,H_DMATRAIN_FULL[4,:],color=:green)
# scatter(1.0e9diameter_model*sqrt(cst_r),H_DMATRAIN[4,:],color=:darkorange)
scatter(1.0e9diameter_model,H_DMATRAIN[4,:],color=:darkorange)
fill_between(1.0e9diameter_model_full,0.8*H_DMATRAIN_FULL[4,:]-dH_FULL[4,:],1.2*H_DMATRAIN_FULL[4,:]+dH_FULL[4,:],alpha=0.5,color=:darkorange)
xlim(2.7,4.0)
xlabel("diameter [nm]")
ylabel("gain [cm\$^{3}\$]")
title("kernel (channel 4)")
legend(["fine model";"coarse model";"model uncertainty"])
# savefig(string(folder,"channel_4_uncertainty.png"))
# savefig(string(folder,"channel_4_uncertainty.pdf"))



figure();
semilogx(1.0e9diameter_model_full,H_DMATRAIN_FULL[6,:],color=:green)
# scatter(1.0e9diameter_model*sqrt(cst_r),H_DMATRAIN[6,:],color=:darkorange)
scatter(1.0e9diameter_model,H_DMATRAIN[6,:],color=:darkorange)
fill_between(1.0e9diameter_model_full,0.8*H_DMATRAIN_FULL[6,:]-dH_FULL[6,:],1.2*H_DMATRAIN_FULL[6,:]+dH_FULL[6,:],alpha=0.5,color=:darkorange)
xlim(5.0,7.5)
xlabel("diameter [nm]")
ylabel("gain [cm\$^{3}\$]")
title("kernel (channel 6)")
legend(["fine model";"coarse model";"model uncertainty"])
# savefig(string(folder,"channel_6_uncertainty.png"))
# savefig(string(folder,"channel_6_uncertainty.pdf"))
