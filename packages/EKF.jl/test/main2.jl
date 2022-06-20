using PyPlot



##################################################
####   load the Kalman Filter implementation   ###
##################################################

using EKF
# create a workspace for the Kalman filter
model_dim = 2
meas_dim  = 1
n_samp    = 50 #0
myWSKF = wsKF(model_dim,meas_dim,n_samp,false,false) # neither the model covariance nor the noise covariance matrices are time varying
myWSKF.t_samp = collect(1.0:n_samp)
t = myWSKF.t_samp




##################################################
###              load the model                ###
##################################################

include("model2.jl")


##################################################
###           simulate acquisition             ###
##################################################

X_samp = Array{Cdouble,2}(undef,meas_dim,n_samp)
Y_samp = Array{Cdouble,2}(undef,meas_dim,n_samp)

X_samp[:,1] .= 1.0
Y_samp[:,1] .= X_samp[1] + Ly*randn()
alpha_GT = Array{Cdouble,1}(undef,n_samp)
alpha_GT[1] = 0.99
dA = -0.004
if false
    for i in 2:n_samp
        if i<20 #0
            # X_samp[:,i] = 0.999*X_samp[:,i-1]
            X_samp[:,i] = 0.99*X_samp[:,i-1]
            alpha_GT[i] = 0.99
        else
            # X_samp[:,i] = A_ev*X_samp[:,i-1]
            X_samp[:,i] = 0.9*X_samp[:,i-1]
            alpha_GT[i] = 0.9
        end
        Y_samp[:,i] = X_samp[:,i] + Ly*randn(1)
    end
else
    for i in 2:n_samp
        alpha_GT[i] = alpha_GT[i-1] + dA
        X_samp[:,i] = alpha_GT[i-1]*X_samp[:,i-1]
        Y_samp[:,i] = X_samp[:,i] + Ly*randn(1)
    end
end


##################################################
###           parameter estimation             ###
##################################################

# initial state
x00       = Array{Cdouble,1}(undef,model_dim)
gam00_vec = Array{Cdouble,1}(undef,model_dim)
var_model = Array{Cdouble,2}(undef,model_dim,model_dim)
L_model   = Array{Cdouble,2}(undef,model_dim,model_dim)
var_noise = Array{Cdouble,2}(undef,meas_dim,meas_dim)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# initialization of the Kalman Filter
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# how much do we trust the model
var_model = L*L'
L_model   = L

#how much do we trust the data
var_noise[:,:] .= Ly*Ly'

# it is important to initialize the measurement model before starting so that we can compute the initial guess
myWSKF.H_me = set_measurement_jacobian!(myWSKF.H_me) # time independent, must be initialized # TODO: move it to the initializing fuction

# draw the initial state
x00[1:meas_dim] = Y_samp[:,1] # or H_me'*y
# x00[2]          = 0.1+Lg*randn()
x00[2]          = 0.8+Lg*randn()

# I don't know how to explain those values
# gam00_vec = [s_1^2, 0.2^2]
gam00_vec = [2.0Ly^2, 0.05^2]

############################################################
###            fixed interval smoother                  ####
############################################################
t0 = 1.0 # time normalization
myWSKF,x_smo_all,o_smo_all = KF_FIS(x00,gam00_vec,var_model,var_noise,t0,myWSKF,Y_samp,n_samp)

# create aliases
x_pre_all = myWSKF.x_pre_all
o_pre_all = myWSKF.o_pre_all
x_fil_all = myWSKF.x_fil_all
o_fil_all = myWSKF.o_fil_all

initColor = "darkcyan"
predColor = "blue" # "#0080d0"
filtColor = "orange" #"#ff7060"
smooColor = "darkgreen" # "#60d060"
dataColor = "darkmagenta" # "cyan" #"darkcyan" # "#ff00ff"
simuColor = "darkred" #"#ff0000"





# explanation plot
#######################################
###       simulation and data       ###
#######################################
figure(1,figsize=(7,4.4))
subplot(121)
plot(t,X_samp[1,:],color=simuColor)
plot(t,Y_samp[1,:],color=dataColor)
title("Variable \$x\$")
xlim(-1.0,t[end]+1.0)
ylim(-0.05,1.05)
legend(["simulation","data"])

subplot(122)
plot(t,alpha_GT,color=simuColor)
title("Parameter \$\\alpha\$")
xlim(-1.0,t[end]+1.0)
ylim(0.69,1.02)
savefig("simulation_and_data.pdf")
savefig("simulation_and_data.png")
close("all")








#######################################
###          initial guess          ###
#######################################
figure(1,figsize=(7,4.4))
subplot(121)
scatter(0.0,x00[1],color=initColor)
plot(t,X_samp[1,:],color=simuColor)
plot(t,Y_samp[1,:],color=dataColor)
title("Variable \$x\$")
xlim(-1.0,t[end]+1.0)
ylim(-0.05,1.05)
legend(["simulation","data","initial guess"])

subplot(122)
scatter(0.0,x00[2],color=initColor)
plot(t,alpha_GT,color=simuColor)
title("Parameter \$\\alpha\$")
xlim(-1.0,t[end]+1.0)
ylim(0.69,1.02)
savefig("KF_state_initial_guess.pdf")
savefig("KF_state_initial_guess.png")


# std of the estimates
figure(2,figsize=(7,4.4))
subplot(121)
scatter(0.0,sqrt(gam00_vec[1]),color=initColor)
title("Variable std \$\\sigma_{x}\$")
xlim(-1.0,t[end]+1.0)
ylim(-0.001,0.055)
legend(["inital guess"])

subplot(122)
scatter(0.0,sqrt(gam00_vec[2]),color=initColor)
title("Parameter std \$\\sigma_{\\alpha}\$")
xlim(-1.0,t[end]+1.0)
ylim(0.005,0.091)
legend(["inital guess"])
savefig("KF_state_std_initial_guess.pdf")
savefig("KF_state_std_initial_guess.png")
close("all")









#######################################
###   first iteration: prediction   ###
#######################################
figure(1,figsize=(7,4.4))
subplot(121)
scatter(0.0,x00[1],color=initColor)
scatter(t[1],x_pre_all[1,1],color=predColor)
plot(t,X_samp[1,:],color=simuColor)
plot(t,Y_samp[1,:],color=dataColor)
title("Variable \$x\$")
xlim(-1.0,t[end]+1.0)
ylim(-0.05,1.05)
legend(["simulation","data","initial guess","prediction"])

subplot(122)
scatter(0.0,x00[2],color=initColor)
scatter(t[1],x_pre_all[2,1],color=predColor)
plot(t,alpha_GT,color=simuColor)
title("Parameter \$\\alpha\$")
xlim(-1.0,t[end]+1.0)
ylim(0.69,1.02)
savefig("KF_state_first_prediction.pdf")
savefig("KF_state_first_prediction.png")


# std of the estimates
figure(2,figsize=(7,4.4))
subplot(121)
scatter(0.0,sqrt(gam00_vec[1]),color=initColor)
scatter(t[1],sqrt(o_pre_all[1,1,1]),color=predColor)
title("Variable std \$\\sigma_{x}\$")
xlim(-1.0,t[end]+1.0)
ylim(-0.001,0.055)
legend(["inital guess","prediction"])

subplot(122)
scatter(0.0,sqrt(gam00_vec[2]),color=initColor)
scatter(t[1],sqrt(o_pre_all[2,2,1]),color=predColor)
title("Parameter std \$\\sigma_{\\alpha}\$")
xlim(-1.0,t[end]+1.0)
ylim(0.005,0.091)
savefig("KF_state_std_first_prediction.pdf")
savefig("KF_state_std_first_prediction.png")
close("all")








########################################
###first iteration: prediction+filter###
########################################
figure(1,figsize=(7,4.4))
subplot(121)
scatter(0.0,x00[1],color=initColor)
scatter(t[1],x_pre_all[1,1],color=predColor)
scatter(t[1],x_fil_all[1,1],color=filtColor)
plot(t,X_samp[1,:],color=simuColor)
plot(t,Y_samp[1,:],color=dataColor)
title("Variable \$x\$")
xlim(-1.0,t[end]+1.0)
ylim(-0.05,1.05)
legend(["simulation","data","initial guess","prediction","filter"])

subplot(122)
scatter(0.0,x00[2],color=initColor)
scatter(t[1],x_pre_all[2,1],color=predColor)
scatter(t[1],x_fil_all[2,1],color=filtColor)
plot(t,alpha_GT,color=simuColor)
title("Parameter \$\\alpha\$")
xlim(-1.0,t[end]+1.0)
ylim(0.69,1.02)
savefig("KF_state_first_estimation.pdf")
savefig("KF_state_first_estimation.png")


# std of the estimates
figure(2,figsize=(7,4.4))
subplot(121)
scatter(0.0,sqrt(gam00_vec[1]),color=initColor)
scatter(t[1],sqrt(o_pre_all[1,1,1]),color=predColor)
scatter(t[1],sqrt(o_fil_all[1,1,1]),color=filtColor)
title("Variable std \$\\sigma_{x}\$")
xlim(-1.0,t[end]+1.0)
ylim(-0.001,0.055)
legend(["inital guess","prediction","filter"])

subplot(122)
scatter(0.0,sqrt(gam00_vec[2]),color=initColor)
scatter(t[1],sqrt(o_pre_all[2,2,1]),color=predColor)
scatter(t[1],sqrt(o_fil_all[2,2,1]),color=filtColor)
title("Parameter std \$\\sigma_{\\alpha}\$")
xlim(-1.0,t[end]+1.0)
ylim(0.005,0.091)
savefig("KF_state_std_first_estimation.pdf")
savefig("KF_state_std_first_estimation.png")
close("all")













#######################################
###          the full KF            ###
#######################################
figure(1,figsize=(7,4.4))
subplot(121)
scatter(0.0,x00[1],color=initColor)
plot(t,X_samp[1,:],color=simuColor)
plot(t,Y_samp[1,:],color=dataColor)
plot(t,x_pre_all[1,:],color=predColor)
plot(t,x_fil_all[1,:],color=filtColor)
title("Variable \$x\$")
xlim(-1.0,t[end]+1.0)
ylim(-0.05,1.05)
legend(["simulation","data","prediction","filter","initial guess"])

subplot(122)
scatter(0.0,x00[2],color=initColor)
plot(t,x_pre_all[2,:],color=predColor)
plot(t,x_fil_all[2,:],color=filtColor)
plot(t,alpha_GT,color=simuColor)
title("Parameter \$\\alpha\$")
xlim(-1.0,t[end]+1.0)
ylim(0.69,1.02)
savefig("KF_state_estimation.pdf")
savefig("KF_state_estimation.png")


# std of the estimates
figure(2,figsize=(7,4.4))
subplot(121)
scatter(0.0,sqrt(gam00_vec[1]),color=initColor)
plot(t,sqrt.(o_pre_all[1,1,:]),color=predColor)
plot(t,sqrt.(o_fil_all[1,1,:]),color=filtColor)
title("Variable std \$\\sigma_{x}\$")
xlim(-1.0,t[end]+1.0)
ylim(-0.001,0.055)
legend(["prediction","filter","inital guess"])

subplot(122)
scatter(0.0,sqrt(gam00_vec[2]),color=initColor)
plot(t,sqrt.(o_pre_all[2,2,:]),color=predColor)
plot(t,sqrt.(o_fil_all[2,2,:]),color=filtColor)
title("Parameter std \$\\sigma_{\\alpha}\$")
xlim(-1.0,t[end]+1.0)
ylim(0.005,0.091)
savefig("KF_state_std_estimation.pdf")
savefig("KF_state_std_estimation.png")
close("all")






if false
    # the full FIKS
    figure(1,figsize=(7,4.4))
    subplot(121)
    scatter(0.0,x00[1],color=initColor)
    plot(t,X_samp[1,:],color=simuColor)
    plot(t,Y_samp[1,:],color=dataColor)
    plot(t,x_pre_all[1,:],color=predColor)
    plot(t,x_fil_all[1,:],color=filtColor)
    plot(t,x_smo_all[1,:],color=smooColor)
    title("Variable \$x\$")
    xlim(-1.0,t[end]+1.0)
    ylim(-0.05,1.05)
    legend(["simulation","data","prediction","filter","smoother","initial guess"])

    subplot(122)
    scatter(0.0,x00[2],color=initColor)
    plot(t,x_pre_all[2,:],color=predColor)
    plot(t,x_fil_all[2,:],color=filtColor)
    plot(t,x_smo_all[2,:],color=smooColor)
    plot(t,alpha_GT,color=simuColor)
    title("Parameter \$\\alpha\$")
    xlim(-1.0,t[end]+1.0)
    ylim(0.69,1.02)
    savefig("FIKS_state_estimation.pdf")
    savefig("FIKS_state_estimation.png")


    # std of the estimates
    figure(2,figsize=(7,4.4))
    subplot(121)
    scatter(0.0,sqrt(gam00_vec[1]),color=initColor)
    plot(t,sqrt.(o_pre_all[1,1,:]),color=predColor)
    plot(t,sqrt.(o_fil_all[1,1,:]),color=filtColor)
    plot(t,sqrt.(o_smo_all[1,1,:]),color=smooColor)
    title("Variable std \$\\sigma_{x}\$")
    xlim(-1.0,t[end]+1.0)
    ylim(-0.001,0.055)
    legend(["prediction","filter","inital guess"])

    subplot(122)
    scatter(0.0,sqrt(gam00_vec[2]),color=initColor)
    plot(t,sqrt.(o_pre_all[2,2,:]),color=predColor)
    plot(t,sqrt.(o_fil_all[2,2,:]),color=filtColor)
    plot(t,sqrt.(o_smo_all[2,2,:]),color=smooColor)
    title("Parameter std \$\\sigma_{\\alpha}\$")
    xlim(-1.0,t[end]+1.0)
    ylim(0.005,0.091)
    savefig("FIKS_state_std_estimation.pdf")
    savefig("FIKS_state_std_estimation.png")
    close("all")
end
