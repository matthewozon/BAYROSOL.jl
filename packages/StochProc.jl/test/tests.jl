using PyPlot     # for plotting
# using StochProc  # stochastic process functions
using Printf
using Statistics # for std
using LinearAlgebra # for diagm

# define some test functions

# test first order process
function test_1st_order(r_time::Cdouble,s_proc::Cdouble,Nt::Int64)
    Ns = 1000
    X_sample           = Array{Cdouble,2}(undef,Nt,Ns)
    X_sample_corr      = Array{Cdouble,2}(undef,Nt,Ns)
    X_sample_corr_time = Array{Cdouble,2}(undef,Nt,Ns)
    for n in 1:Ns
        x0 = s_proc*sqrt(1.0-r_time^2)randn()
        # the basic 1st order
        X_sample[:,n]           = StochProc.SP_1T(r_time,s_proc,Nt,x0)
        # the asymptotic correction
        X_sample_corr[:,n]      = StochProc.SP_1TC(r_time,s_proc,Nt,x0)
        # the time varying correction
        X_sample_corr_time[:,n] = StochProc.SP_1TCV(r_time,s_proc,Nt,x0)
    end

    # plot the samples
    fig, (axs1,axs2,axs3) = subplots(1,3)
    if Ns<30
        axs1.plot(X_sample)
        axs1.set_title("process")
        axs2.plot(X_sample_corr)
        axs2.set_title("process asymptotical corrected")
        axs3.plot(X_sample_corr_time)
        axs3.set_title("process corrected")
    else
        axs1.plot(X_sample[:,1:30])
        axs1.set_title("process")
        axs2.plot(X_sample_corr[:,1:30])
        axs2.set_title("process asymptotical corrected")
        axs3.plot(X_sample_corr_time[:,1:30])
        axs3.set_title("process corrected")
    end
    # compute the variance
    s_sample           = std(X_sample,dims=2)
    s_sample_corr      = std(X_sample_corr,dims=2)
    s_sample_corr_time = std(X_sample_corr_time,dims=2)
    # plot the standard deviations with the expected ones
    figure(2)
    plot(s_proc*ones(Nt))
    plot(s_sample)
    plot(s_sample_corr)
    plot(s_sample_corr_time)
    legend(["expected","process","asymptotic correction","time varying correction"])
    title("Standard deviation")
end

# run test
# test_1st_order(0.995,1.0,200)





# test second order process

function test_2nd_order(r_1::Cdouble,r_2::Cdouble,s_proc::Cdouble,Nt::Int64)
    Ns = 1000
    X_sample               = Array{Cdouble,2}(undef,Nt,Ns)
    X_sample_corr          = Array{Cdouble,2}(undef,Nt,Ns)
    X_sample_corr_time     = Array{Cdouble,2}(undef,Nt,Ns)
    X_sample_corr_time_mat = Array{Cdouble,2}(undef,Nt,Ns)
    for n in 1:Ns
        x0 = s_proc*sqrt(1.0 -(r_1+r_2)^2 -(r_1*r_2)^2 +2.0*r_1*r_2*((r_1+r_2)^2)/(1.0+r_1*r_2) )randn()
        x1 = s_proc*sqrt(1.0 -(r_1+r_2)^2 -(r_1*r_2)^2 +2.0*r_1*r_2*((r_1+r_2)^2)/(1.0+r_1*r_2) )randn()
        # the basic 1st order
        X_sample[:,n]               = StochProc.SP_2T(r_1,r_2,s_proc,Nt,x0,x1)
        # the asymptotic correction
        X_sample_corr[:,n]          = StochProc.SP_2TC(r_1,r_2,s_proc,Nt,x0,x1)
        # the time varying correction
        X_sample_corr_time[:,n]     = StochProc.SP_2TCV(r_1,r_2,s_proc,Nt,x0,x1)
        X_sample_corr_time_mat[:,n] = StochProc.SP_2TCV_mat(r_1,r_2,s_proc,Nt,x0,x1)
    end

    # plot the samples
    fig, (axs1,axs2,axs3,axs4) = subplots(1,4)
    if Ns<20
        axs1.plot(X_sample)
        axs1.set_title("process")
        axs2.plot(X_sample_corr)
        axs2.set_title("process asymptotical corrected")
        axs3.plot(X_sample_corr_time)
        axs3.set_title("process corrected")
        axs4.plot(X_sample_corr_time_mat)
        axs4.set_title("process corrected mat")
    else
        axs1.plot(X_sample[:,1:20])
        axs1.set_title("process")
        axs2.plot(X_sample_corr[:,1:20])
        axs2.set_title("process asymptotical corrected")
        axs3.plot(X_sample_corr_time[:,1:20])
        axs3.set_title("process corrected")
        axs4.plot(X_sample_corr_time_mat[:,1:20])
        axs4.set_title("process corrected mat")
    end
    # compute the variance
    s_sample               = std(X_sample,dims=2)
    s_sample_corr          = std(X_sample_corr,dims=2)
    s_sample_corr_time     = std(X_sample_corr_time,dims=2)
    s_sample_corr_time_mat = std(X_sample_corr_time_mat,dims=2)
    # plot the standard deviations with the expected ones
    figure(2)
    plot(s_proc*ones(Nt))
    plot(s_sample)
    plot(s_sample_corr)
    plot(s_sample_corr_time)
    plot(s_sample_corr_time_mat)
    legend(["expected std","process","asymptotic correction","time varying correction","time varying correction mat"])
    title("Standard deviation")
end


# run test
# test_2nd_order(0.9,0.9,1.0,500)



# test first order process for vector series
function test_1st_vector(r_time::Cdouble, Nt::Int64=3000)
    # define the expected variance of each element in the vector
    D_tilde = exp.(-collect(range(0.0,5.0,length=19)))

    # R_time = 0.9collect(0.99:0.0005:0.999)
    R_time = 0.99collect(0.99:0.0005:0.999)
    B = diagm(r_time*ones(19))

    L_source, Gamma_source = StochProc.space_covariance_chol(0.92,D_tilde,B);
    X = StochProc.SP_1T_V(R_time,L_source,Nt,1.0randn(19));
    figure(3); plot(X');
    figure(1); imshow(Gamma_source); colorbar();
    figure(2); imshow(cov(X[:,1000:end],dims=2)); colorbar()


    #
    B = diagm(r_time*ones(19))
    L_source, Gamma_source = StochProc.space_covariance_chol(0.92,D_tilde,B);
    X = StochProc.SP_1TC_V(r_time,L_source,Nt,1.0randn(19));
    figure(6); plot(X');
    figure(4); imshow(Gamma_source); colorbar();
    figure(5); imshow(cov(X[:,1000:end],dims=2)); colorbar()



    #
    X = StochProc.SP_1TC_V(R_time,L_source,Nt,1.0randn(19));
    figure(9); plot(X');
    figure(7); imshow(Gamma_source); colorbar();
    figure(8); imshow(cov(X[:,1000:end],dims=2)); colorbar()


    #
    B = zeros(19,19)
    b = [2.0r_time  -r_time^2; 1.0 0.0]
    B[1,1] = r_time
    B[2:3,2:3] = b
    B[4:5,4:5] = b
    B[6:7,6:7] = b
    B[8:9,8:9] = b
    B[10:11,10:11] = b
    B[12:13,12:13] = b
    B[14:15,14:15] = b
    B[16:17,16:17] = b
    B[18:19,18:19] = b
    # B =  diagm(0.99*ones(19)) - diagm(-1 => 0.001ones(18)) #
    B = diagm(R_time) #
    L_source, Gamma_source = StochProc.space_covariance_chol(0.92,D_tilde,B);
    eigvals(Gamma_source - B*Gamma_source*B')
    # X = SP_1TC_V(B,L_source,3000,1.0randn(19)); figure(12); plot(X'); figure(10); imshow(Gamma_source); colorbar(); figure(11); imshow(cov(X[:,1000:end],2)); colorbar()
    X = StochProc.SP_1TC_V(B,L_source,Nt,sqrt.(D_tilde).*randn(19));
    figure(12); plot(X');
    figure(10); imshow(Gamma_source); colorbar();
    figure(11); imshow(cov(X[:,1000:end],dims=2)); colorbar()
    X
end


# run test
# test_1st_vector(0.99, 3000)



function test_covariance(r_pol_::Cdouble,Nr::Int64=50,Ns::Int64=20)
    # define the variance of each element of the vector
    # D_tilde = 1.0+0.95cos(collect(linspace(0.0,4pi,Nr)))
    # D_tilde = 1.0+1.0cos(collect(linspace(0.0,4pi,Nr)))
    D_tilde = 1.0.+0.5cos.(collect(range(0.0,4pi,length=Nr)))
    # D_tilde = D_tilde.-minimum(D_tilde)
    D_tilde = exp.(-collect(range(0.0,5.0,length=Nr)))
    # for a few time process, compute the covariance matrices
    fig = 1
    for r_time in [0.0; 0.5; 0.9; 0.99]
        # time process
        B = r_time*Matrix{Cdouble}(I,Nr,Nr)
        # compute the covariance and its cholesky factor
        println(r_time)
        L_source, Gamma_source = StochProc.space_covariance_chol(r_pol_,D_tilde,B)
        # plot the result
        figure(fig)
        imshow(Gamma_source)
        ti = @sprintf "covaraince, root of the time process %f" r_time
        title(ti)
        colorbar()
        # plot the diagonal of the covariance matrix
        figure(fig+1)
        plot(sqrt.(diag(Gamma_source)))
        # draw a sample from the distribution
        X = Array{Cdouble,2}(undef,Nr,Ns)
        Y = Array{Cdouble,2}(undef,Nr,Ns)
        for n in 1:Ns
            X[:,n] = randn(Nr)
            Y[:,n] = L_source*X[:,n]
        end
        if Ns<21
            figure(fig+2)
            plot(X)
            title("samples")
            figure(fig+3)
            plot(Y)
            title("correlated samples")
            figure(fig+4)
            plot(std(Y,dims=2))
            title("estimation of the std")
            # figure(fig+5)
            # plot(sum(full(L_source),2))
        else
            figure(fig+2)
            plot(X[:,1:20])
            title("samples")
            figure(fig+3)
            plot(Y[:,1:20])
            title("correlated samples")
            figure(fig+4)
            plot(std(Y,dims=2))
            title("estimation of the std")
        end
        fig += 10
    end

end

# run test
# test_covariance(0.9,50,200)



function test_iter_scalar(Nt::Int64=200)
    # scalar first order
    ord = 1
    dim = 1
    r = 0.9
    B = r
    V = 1.0
    ws1stS = StochProc.linSP(ord,dim,r,V)
    X = Array{Cdouble,1}(undef,Nt)
    X[1] = ws1stS.L_source*randn()
    for t in 2:Nt
        X[t] = StochProc.iterator(X[t-1],ws1stS)
    end
    figure(1)
    plot(X)
    title(@sprintf "1\$^{st}\$ order process, r = %f" r )

    # scalar second order
    ord = 2
    dim = 1
    r = 0.9
    B = [2.0r -r^2; 1.0 0.0]
    V = 1.0
    ws2ndS = StochProc.linSP(ord,dim,B,V)
    X = Array{Cdouble,2}(undef,ord,Nt)
    X[:,1] = ws2ndS.L_source*randn(ord)
    for t in 2:Nt
        X[:,t] = StochProc.iterator(X[:,t-1],ws2ndS)
    end
    figure(2)
    plot(X[1,:])
    title(@sprintf "2\$^{nd}\$ order process, r = %f" r )

    # scalar third order
    ord = 3
    dim = 1
    r = 0.9
    B = [3.0r -3.0r^2 r^3; 1.0 0.0 0.0; 0.0 1.0 0.0]
    V = 1.0
    ws3rdS = StochProc.linSP(ord,dim,B,V)
    X = Array{Cdouble,2}(undef,ord,Nt)
    X[:,1] = ws3rdS.L_source*randn(ord)
    for t in 2:Nt
        X[:,t] = StochProc.iterator(X[:,t-1],ws3rdS)
    end
    figure(3)
    plot(X[1,:])
    title(@sprintf "3\$^{rd}\$ order process, r = %f" r )
end

# run test_iter_scalar()


function test_iter_vector(Nt::Int64=2000)
    # vector first order
    ord = 1
    dim = 10
    r = 0.9
    B = r*Matrix{Cdouble}(I,dim,dim)
    r_pol = 0.95
    D_tilde = ones(dim)
    G_proc = StochProc.covariance_process_2nd(r_pol,D_tilde)
    V = StochProc.covariance_2nd_order_space_1st_order_time(r_pol,D_tilde,B)
    ws1stV = StochProc.linSP(ord,dim,r,V)
    X = Array{Cdouble,2}(undef,dim,Nt)
    X[:,1] = ws1stV.L_source*randn(dim)
    for t in 2:Nt
        X[:,t] = StochProc.iterator(X[:,t-1],ws1stV)
    end
    figure(1)
    plot(X[:,Nt-500:end]')
    title(@sprintf "1\$^{st}\$ order process, r = %f" r )

    figure(2)
    imshow(G_proc)
    title("Expected covariance of the process")
    colorbar()

    figure(3)
    imshow(cov(X[:,500:end],dims=2))
    title("Estimated covariance of the process")
    colorbar()



    # vector second order
    ord = 2
    dim = 10
    # time evolution model
    r = 0.9
    B_sub = [2.0r -r^2; 1.0 0.0]
    B = StochProc.timeBlockMatrix(B_sub,dim)
    # covariance model
    r_pol = 0.95
    D_tilde = ones(dim)

    G_source = StochProc.covariance_process_2nd(r_pol,D_tilde)
    G_proc_expend = zeros(ord*dim,ord*dim)
    G_proc_expend[1:ord:end,1:ord:end] = G_source
    for t in 2:Nt
        G_proc_expend = B*G_proc_expend*B'
        G_proc_expend[1:ord:end,1:ord:end] = G_proc_expend[1:ord:end,1:ord:end] + G_source
    end
    G_proc = G_proc_expend[1:ord:end,1:ord:end]


    # create a linear stochastic process
    ws2ndV = StochProc.linSP(ord,dim,B,G_source)
    X = Array{Cdouble,2}(undef,dim*ord,Nt)
    X[:,1] = StochProc.initGaussRep(0.0,ws2ndV.L_source,ws2ndV.procOrd)
    for t in 2:Nt
        X[:,t] = StochProc.iterator(X[:,t-1],ws2ndV)
    end

    figure(11)
    plot(X[1:ws2ndV.procOrd:end,Nt-500:end]')
    title(@sprintf "2\$^{nd}\$ order process, r = %f" r )

    figure(12)
    imshow(G_proc_expend)
    title("Expected covariance of the process")
    colorbar()

    figure(13)
    imshow(cov(X[1:ws2ndV.procOrd:end,500:end],dims=2))
    title("Estimated covariance of the process")
    colorbar()

    figure(14)
    imshow(G_source)
    title("covariance of the source")
    colorbar()
end




# another test
function toto(N::Int64,p::Cdouble=0.7) # at most 25, more dimensions will make too unstable (it's much better to keep it smaller than 5 dimensions)
    # current estimation of the uncertainty ranges
    low_sig = log(1.0+exp(-1.0))*collect(1.0:N)
    up_sig = log(1.0+exp(1.0))*collect(1.0:N)
    estimate_ = log(1.0+exp(0.0))*collect(1.0:N)

    # find the isocontour such that the probability of being in the enclosed volume is 0.7
    ra = 0.0
    while(StochProc.multivariate_probability(N,ra)<p)
        ra = ra + 1.0
    end
    # refine the search by dicotomie
    minR = ra-1.0;
    maxR = ra;
    while((maxR-minR)>1.0e-10)
        if (StochProc.multivariate_probability(N,0.5*(minR+maxR))>0.7)
            maxR = 0.5*(minR+maxR)
        else
            minR = 0.5*(minR+maxR)
        end
    end
    RA = 0.5*(minR+maxR)
    # new estimation
    S = RA*StochProc.uniform_sample_unit_sphere(N,min(N*5000,500000)); # N=100 is the practical maximum for memory considerations
    X = similar(S)
    X[1,:] = log.(1.0.+exp.(S[1,:]))
    for i in 2:N
        X[i,:] = log.(1.0.+exp.(S[i,:])) + X[i-1,:]
    end
    low_sig_new = dropdims(minimum(X,dims=2),dims=2)
    up_sig_new  = dropdims(maximum(X,dims=2),dims=2)

    # plot
    figure(456); plot(collect(1.0:N),estimate_)
    fill_between(collect(1.0:N),low_sig,up_sig,alpha=0.5)
    fill_between(collect(1.0:N),low_sig_new,up_sig_new,alpha=0.5)
    RA
end


function tata(p::Cdouble=0.7)
    # get the number of dimensions
    N = 2

    # create a transformation
    f = Array{Function,1}(undef,N)
    for n in 1:N
        f[n] = x::Array{Cdouble,1}->sum(log.(1.0.+exp.(x[1:n])))
    end

    # statistics of the original gaussian distribution
    mu = -ones(2)
    gamma = [4.0 -1.0; -1.0 1.0] # [1.0 -1.99; -1.99 4.0]

    # current estimation of the uncertainty ranges
    low_sig  = Array{Cdouble,1}(undef,N)
    up_sig   = Array{Cdouble,1}(undef,N)
    mean_est = Array{Cdouble,1}(undef,N)
    for i in 1:N
        low_sig[i]  = f[i](mu.-sqrt.(diag(gamma)))
        up_sig[i]   = f[i](mu.+sqrt.(diag(gamma)))
        mean_est[i] = f[i](mu)
    end

    # compute the contour of the region of probability p
    X,S = StochProc.region_high_probability(f,mu,Symmetric(gamma),p)

    # new way to get the uncertainty... which does not account for the actual distribution, but I can't find anything better
    low_sig_new = dropdims(minimum(X,dims=2),dims=2)
    up_sig_new  = dropdims(maximum(X,dims=2),dims=2)

    # plot
    figure(456); plot(collect(1.0:N),mean_est)
    fill_between(collect(1.0:N),low_sig,up_sig,alpha=0.5)
    fill_between(collect(1.0:N),low_sig_new,up_sig_new,alpha=0.5)
    figure(10); scatter(S[1,:],S[2,:]); scatter(X[1,:],X[2,:])
    minX,maxX = extrema(S[1,:])
    minY,maxY = extrema(S[2,:])
    plot([minX; maxX; maxX; minX; minX],[minY; minY; maxY; maxY; minY]) # ,color=:darkorange)
    minX,maxX = extrema(X[1,:])
    minY,maxY = extrema(X[2,:])
    plot([minX; maxX; maxX; minX; minX],[minY; minY; maxY; maxY; minY]) # ,color=:darkorange)
    title("highest probability regions")
    legend(["overestimation X,Y", "overestimation U,V","X,Y trust region","U,V trust region"])
    S,X
end

function tutu(p::Cdouble=0.7)
    # get the number of dimensions
    N = 2

    # create a transformation
    f = Array{Function,1}(undef,N)
    for n in 1:N
        f[n] = x::Array{Cdouble,1}->sum(log.(1.0.+exp.(x[1:n])))
    end

    # statistics of the original gaussian distribution
    mu = -ones(2)
    gamma = [1.0 -1.99; -1.99 4.0]

    # current estimation of the uncertainty ranges
    low_sig  = Array{Cdouble,1}(undef,N)
    up_sig   = Array{Cdouble,1}(undef,N)
    mean_est = Array{Cdouble,1}(undef,N)
    for i in 1:N
        low_sig[i]  = f[i](mu.-sqrt.(diag(gamma)))
        up_sig[i]   = f[i](mu.+sqrt.(diag(gamma)))
        mean_est[i] = f[i](mu)
    end

    # compute the contour of the region of probability p
    sig_new = StochProc.percentile_estimation_ni(f,mu,Symmetric(gamma),p)
    low_sig_new = sig_new[:,1] # squeeze(sig_new[1,:],1)
    up_sig_new  = sig_new[:,2] # squeeze(sig_new[2,:],1)


    # plot
    figure(456); plot(collect(1.0:N),mean_est)
    fill_between(collect(1.0:N),low_sig,up_sig,alpha=0.5)
    fill_between(collect(1.0:N),low_sig_new,up_sig_new,alpha=0.5)
    # figure(10); scatter(S[1,:],S[2,:]); scatter(X[1,:],X[2,:])
    # S,X
    sig_new
end


