using PyPlot



##################################################
####   load the Kalman Filter implementation   ###
##################################################

using EKF
# create a workspace for the Kalman filter
model_dim = 5
meas_dim  = 2
n_samp    = 500
myWSKF = wsKF(model_dim,meas_dim,n_samp,false,false) # neither the model covariance nor the noise covariance matrices are time varying
myWSKF.t_samp = collect(1.0:n_samp)




##################################################
###              load the model                ###
##################################################

include("model.jl")


##################################################
###           simulate acquisition             ###
##################################################

X_samp = Array{Cdouble,2}(undef,meas_dim,n_samp)
Y_samp = Array{Cdouble,2}(undef,meas_dim,n_samp)

X_samp[:,1] = [1.0 2.0]
Y_samp[:,1] = X_samp[:,1] + Ly*randn(2)
for i in 2:n_samp
    X_samp[:,i] = A_ev*X_samp[:,i-1]
    Y_samp[:,i] = X_samp[:,i] + Ly*randn(2)
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
var_noise = Ly*Ly'

# it is important to initialize the measurement model before starting so that we can compute the initial guess
myWSKF.H_me = set_measurement_jacobian!(myWSKF.H_me) # time independent, must be initialized # TODO: move it to the initializing fuction

# draw the initial state
x00[1:meas_dim] = Y_samp[:,1] # or H_me'*y
x00[3]          = 0.1+Lg*randn()
x00[4:5]        = 0.5 .+ Lp*randn(2)

# I don't know how to explain those values
gam00_vec = [s_1^2, s_2^2, 1.1^2, 2.1^2, 1.1^2]

############################################################
###             fix interval smoother                   ####
############################################################
t0 = 1.0 # time normalization
# myWSKF,x_smo_all,o_smo_all = KF_FIS(x00,gam00_vec,var_model,L_model,var_noise,t0,myWSKF, Y_samp,n_samp)
myWSKF,x_smo_all,o_smo_all = KF_FIS(x00,gam00_vec,var_model,var_noise,t0,myWSKF,Y_samp,n_samp)

# create aliases
x_pre_all = myWSKF.x_pre_all
o_pre_all = myWSKF.o_pre_all
x_fil_all = myWSKF.x_fil_all
o_fil_all = myWSKF.o_fil_all

# plot the results
figure(1)
plot(x_pre_all[1,:])
plot(x_fil_all[1,:])
plot(x_smo_all[1,:])
plot(X_samp[1,:])
plot(Y_samp[1,:])
title("First variable")
legend(["prediction","filter","smoother","simulation","data"])


figure(2)
plot(x_pre_all[2,:])
plot(x_fil_all[2,:])
plot(x_smo_all[2,:])
plot(X_samp[2,:])
plot(Y_samp[2,:])
title("Second variable")
legend(["prediction","filter","smoother","simulation","data"])


figure(3)
plot(x_pre_all[4,:])
plot(x_fil_all[4,:])
plot(x_smo_all[4,:])
plot(alpha*ones(n_samp))
title("parameter \$\\alpha\$")
legend(["prediction","filter","smoother","simulation"])


figure(4)
plot(x_pre_all[5,:])
plot(x_fil_all[5,:])
plot(x_smo_all[5,:])
plot(beta*ones(n_samp))
title("parameter \$\\beta\$")
legend(["prediction","filter","smoother","simulation"])


figure(5)
plot(x_pre_all[3,:])
plot(x_fil_all[3,:])
plot(x_smo_all[3,:])
plot(gamma*ones(n_samp))
title("parameter \$\\gamma\$")
legend(["prediction","filter","smoother","simulation"])




# std of the estimates
figure(11)
plot(sqrt.(o_pre_all[1,1,:]))
plot(sqrt.(o_fil_all[1,1,:]))
plot(sqrt.(o_smo_all[1,1,:]))
title("std of the first variable")
legend(["prediction","filter","smoother"])

figure(12)
plot(sqrt.(o_pre_all[2,2,:]))
plot(sqrt.(o_fil_all[2,2,:]))
plot(sqrt.(o_smo_all[2,2,:]))
title("std of the second variable")
legend(["prediction","filter","smoother"])

figure(13)
plot(sqrt.(o_pre_all[4,4,:]))
plot(sqrt.(o_fil_all[4,4,:]))
plot(sqrt.(o_smo_all[4,4,:]))
title("std of the parameter \$\\alpha\$")
legend(["prediction","filter","smoother"])

figure(14)
plot(sqrt.(o_pre_all[5,5,:]))
plot(sqrt.(o_fil_all[5,5,:]))
plot(sqrt.(o_smo_all[5,5,:]))
title("std of the parameter \$\\beta\$")
legend(["prediction","filter","smoother"])

figure(15)
plot(sqrt.(o_pre_all[3,3,:]))
plot(sqrt.(o_fil_all[3,3,:]))
plot(sqrt.(o_smo_all[3,3,:]))
title("std of the parameter \$\\gamma\$")
legend(["prediction","filter","smoother"])
