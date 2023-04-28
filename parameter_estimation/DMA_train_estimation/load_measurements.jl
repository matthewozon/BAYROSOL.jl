##
## load data
##

if FLAG_1952_02
    df_y_all = CSV.File(string(input_folder,"cts_1952_02_dmatrain.csv"); header=true)|> DataFrame
    y_all = Matrix{Cdouble}(Matrix{Cdouble}(df_y_all)');
end


if FLAG_0000_00
    df_y_all = CSV.File(string(input_folder,"cts_0000_00_dmatrain.csv"); header=false)|> DataFrame
    y_all = Matrix{Cdouble}(Matrix{Cdouble}(df_y_all)');
end
if FLAG_0000_01
    df_y_all = CSV.File(string(input_folder,"cts_0000_01_dmatrain.csv"); header=false)|> DataFrame
    y_all = Matrix{Cdouble}(Matrix{Cdouble}(df_y_all)');
end
if FLAG_0000_02
    df_y_all = CSV.File(string(input_folder,"cts_0000_02_dmatrain.csv"); header=false)|> DataFrame
    y_all = Matrix{Cdouble}(Matrix{Cdouble}(df_y_all)');
end
if FLAG_0000_03
    df_y_all = CSV.File(string(input_folder,"cts_0000_03_dmatrain.csv"); header=false)|> DataFrame
    y_all = Matrix{Cdouble}(Matrix{Cdouble}(df_y_all)');
end
if FLAG_0000_04
    df_y_all = CSV.File(string(input_folder,"cts_0000_04_dmatrain.csv"); header=false)|> DataFrame
    y_all = Matrix{Cdouble}(Matrix{Cdouble}(df_y_all)');
end

if FLAG_1802_01
    df_y_all = CSV.File(string(input_folder,"cts_1802_01_dmatrain.csv"); header=true)|> DataFrame
    y_all = Matrix{Cdouble}(Matrix{Cdouble}(df_y_all)');
    # get rid of the outlier. You may comment this out to see the difference. Their is on pixel that was saturated for one reason or another during the acquisition, and it is not well handled by the Gaussian aproximation underlying the KF filter, so it must beignored. The following two lines do the trick, but it could be done automatically using a median filter along the time direction.
    idx_pair = findfirst(y_all.>1000)
    y_all[idx_pair] = 0.5*(y_all[idx_pair.I[1],idx_pair.I[2]-1]+y_all[idx_pair.I[1],idx_pair.I[2]+1])
end
if FLAG_1906_03
    df_y_all = CSV.File(string(input_folder,"cts_1906_03_dmatrain.csv"); header=true)|> DataFrame
    y_all = Matrix{Cdouble}(Matrix{Cdouble}(df_y_all)');
end

if (!GT_loaded)
    if FLAG_1802_01
        J_psm = CSV.File(string(input_folder,"J_1802.csv"); header=true)|> DataFrame; # smooth PSM estimation
        J_psm = CSV.File(string(input_folder,"J_1802_new.csv"); header=true)|> DataFrame; # much more bumpy
    end
    if FLAG_1952_02
        J_psm = CSV.File(string(input_folder,"J_1952.csv"); header=true)|> DataFrame;
    end
    if FLAG_1906_03
        J_psm = CSV.File(string(input_folder,"J_1906_03.csv"); header=true)|> DataFrame;
    end

    # J_psm = DataFrame(CSV.File(string(input_folder,"J_1916.csv"); header=true));
    J_TIME = J_psm[:,1];
    J_time = zeros(Cdouble,length(J_TIME));
    for i in 1:length(J_TIME)
        # J_time[i] = DateTime(J_TIME[i][1:end-7],"yyyy-mm-dd HH:MM:SS").instant.periods.value/1000.0
        J_time[i] = DateTime(J_TIME[i][1:19],"yyyy-mm-dd HH:MM:SS").instant.periods.value/1000.0
    end
    # J_time = J_time .- J_time[1];
    J_est = J_psm[:,2];
end

if GR_EST_1952_02
    GR_inside = CSV.File(string(input_folder,"GR_INSIDE_1952_02_2h20min.csv"); header=true)|> DataFrame;
    GR_dp = 1.0e-9GR_inside[:,1]
    GR_est = GR_inside[:,2]
end

if GR_EST_1802_01
    GR_inside = CSV.File(string(input_folder,"GR_INSIDE_1802_01_2h30min.csv"); header=true)|> DataFrame;
    GR_dp = 1.0e-9GR_inside[:,1]
    GR_est = GR_inside[:,2]

    GR_inside_2 = CSV.File(string(input_folder,"GR_INSIDE_1802_01_3h30min.csv"); header=true)|> DataFrame;
    GR_dp_2 = 1.0e-9GR_inside[:,1]
    GR_est_2 = GR_inside[:,2]
end

if GR_EST_1906_03
    GR_inside = CSV.File(string(input_folder,"GR_INSIDE_1906_03_0h16min.csv"); header=true)|> DataFrame;
    GR_dp = 1.0e-9GR_inside[:,1]
    GR_est = GR_inside[:,2]
end

if FLAG_0000_00
    df_diam_data = CSV.File(string(input_folder,"diam_0000_00_dmatrain.csv"); header=false)|> DataFrame;
    diameter_data = Matrix{Cdouble}(df_diam_data)[1,:]
end
if FLAG_0000_01
    df_diam_data = CSV.File(string(input_folder,"diam_0000_01_dmatrain.csv"); header=false)|> DataFrame;
    diameter_data = Matrix{Cdouble}(df_diam_data)[1,:]
end
if FLAG_0000_02
    df_diam_data = CSV.File(string(input_folder,"diam_0000_02_dmatrain.csv"); header=false)|> DataFrame;
    diameter_data = Matrix{Cdouble}(df_diam_data)[1,:]
end
if FLAG_0000_03
    df_diam_data = CSV.File(string(input_folder,"diam_0000_03_dmatrain.csv"); header=false)|> DataFrame;
    diameter_data = Matrix{Cdouble}(df_diam_data)[1,:]
end
if FLAG_0000_04
    df_diam_data = CSV.File(string(input_folder,"diam_0000_04_dmatrain.csv"); header=false)|> DataFrame;
    diameter_data = Matrix{Cdouble}(df_diam_data)[1,:]
end
if FLAG_1952_02
    df_diam_data = CSV.File(string(input_folder,"diam_1952_02_dmatrain.csv"); header=true)|> DataFrame;
    diameter_data = 1.0e-9*Matrix{Cdouble}(df_diam_data)[1,:]
end
if FLAG_1802_01
    df_diam_data = CSV.File(string(input_folder,"diam_1802_01_dmatrain.csv"); header=true)|> DataFrame;
    diameter_data = 1.0e-9*Matrix{Cdouble}(df_diam_data)[1,:]
end

if FLAG_1906_03
    df_diam_data = CSV.File(string(input_folder,"diam_1906_03_dmatrain.csv"); header=true)|> DataFrame;
    diameter_data = 1.0e-9*Matrix{Cdouble}(df_diam_data)[1,:]
end

volume_data = (pi/6.0)*(diameter_data.^3)
N = length(diameter_data)


if FLAG_0000_00
    df_T_samp = CSV.File(string(input_folder,"time_0000_00_dmatrain.csv"); header=false)|> DataFrame;
    T_samp = Matrix{Cdouble}(df_T_samp)[1,:];
end
if FLAG_0000_01
    df_T_samp = CSV.File(string(input_folder,"time_0000_01_dmatrain.csv"); header=false)|> DataFrame;
    T_samp = Matrix{Cdouble}(df_T_samp)[1,:];
end
if FLAG_0000_02
    df_T_samp = CSV.File(string(input_folder,"time_0000_02_dmatrain.csv"); header=false)|> DataFrame;
    T_samp = Matrix{Cdouble}(df_T_samp)[1,:];
end
if FLAG_0000_03
    df_T_samp = CSV.File(string(input_folder,"time_0000_03_dmatrain.csv"); header=false)|> DataFrame;
    T_samp = Matrix{Cdouble}(df_T_samp)[1,:];
end
if FLAG_0000_04
    df_T_samp = CSV.File(string(input_folder,"time_0000_04_dmatrain.csv"); header=true)|> DataFrame
    T_samp = dropdims(Matrix(df_T_samp),dims=2)
end
if FLAG_1952_02
    df_T_samp = CSV.File(string(input_folder,"time_1952_02_dmatrain.csv"); header=true)|> DataFrame
    T_samp = dropdims(Matrix(df_T_samp),dims=2)
end
if FLAG_1802_01
    df_T_samp = CSV.File(string(input_folder,"time_1802_01_dmatrain.csv"); header=true)|> DataFrame
    T_samp = dropdims(Matrix(df_T_samp),dims=2)
end
if FLAG_1906_03
    df_T_samp = CSV.File(string(input_folder,"time_1906_03_dmatrain.csv"); header=true)|> DataFrame
    T_samp = dropdims(Matrix(df_T_samp),dims=2)
end

n_samp = length(T_samp)
t_samp = zeros(Cdouble,n_samp);
if (FLAG_1952_02 | FLAG_1802_01 | FLAG_1906_03)
    for i in 1:n_samp
        t_samp[i] = DateTime(T_samp[i][1:end-4],"yyyy-mm-dd HH:MM:SS").instant.periods.value/1000.0
    end
    J_time = J_time .- t_samp[1];
    t_samp = t_samp .- t_samp[1];
else
    # for 0000_00
    if (FLAG_0000_00 | FLAG_0000_01 | FLAG_0000_02 | FLAG_0000_03)
        t_samp = 3600.0*T_samp
    end
    # for 0000 04
    if FLAG_0000_04
        for i in 1:n_samp
            t_samp[i] = T_samp[i].instant.periods.value/1000.0
        end
    end
    t_samp = t_samp .- t_samp[1];
end



if PLOT_MORE
    min_count,max_count = extrema(y_all);
    s = @sprintf "raw counts (%1.2f,%1.2f)" min_count max_count
    fig1, ax1 =  displayLogData2D(1,t_samp/3600.0,1.0e9diameter_data,y_all,max(0.01max_count,min_count),max_count,_title=s,_colorbar_label="numbers []")

    ylabel("diameter [nm]")
    savefig(string(folder,"raw_counts.png"))
    savefig(string(folder,"raw_counts.pdf"))
end

##
## pad the data to emulate a btter time resolution of the acquisition
##
dt = median(t_samp[2:end]-t_samp[1:end-1])
if need_padding
    # padding of the data so that the algorithm's timestep length is in the range [30,60]s
    dt_data = dt
    n_samp_data = n_samp
    pad_factor = floor(Int64,0.5dt_data/30.0)
    if FLAG_1906_03
        pad_factor = floor(Int64,0.5dt_data/1.0)
    end
    # allocate a new time array
    n_samp = (n_samp_data-1)*(pad_factor+1)+1 # (n_samp_data)*(pad_factor+1)+1 # 
    t_samp_start = t_samp[1];
    t_samp = t_samp_start .+ (dt_data/(pad_factor+1))*collect(0:1:((n_samp_data-1)*(pad_factor+1)));
    # t_samp = t_samp_start .+ (dt_data/(pad_factor+1))*collect(0:1:((n_samp_data)*(pad_factor+1)));
    dt = dt_data/(pad_factor+1)
    # allocate a new data array
    y_all_data = copy(y_all)
    y_all = zeros(Cdouble,length(diameter_data),n_samp)
    y_all[:,1:pad_factor+1:end] = y_all_data
end

##
## load the ground truth if it is available
##
if GT_loaded

    if FLAG_0000_04
        df_tp                 = CSV.File(string(input_folder,"time_hour_parameter.csv"); delim=",", header=false)|> DataFrame;
        tp                    = dropdims(Matrix{Cdouble}(df_tp),dims=2); # for 0000 04
    else
        df_tp                 = CSV.File(string(input_folder,"time_hour_parameter.csv"); delim=",", header=false)|> DataFrame;
        tp                    = Matrix{Cdouble}(df_tp)[1,:]
    end
    if FLAG_0000_00
        df_dp                 = CSV.File(string(input_folder,"measurement_operator_discretization_0000_00.csv"); delim=",", header=false)|> DataFrame;
        dp                    = Matrix{Cdouble}(df_dp)[1,:]
    end
    if FLAG_0000_01
        df_dp                 = CSV.File(string(input_folder,"measurement_operator_discretization_0000_01.csv"); delim=",", header=false)|> DataFrame;
        dp                    = Matrix{Cdouble}(df_dp)[1,:]
    end
    if FLAG_0000_02
        df_dp                 = CSV.File(string(input_folder,"measurement_operator_discretization_0000_02.csv"); delim=",", header=false)|> DataFrame;
        dp                    = Matrix{Cdouble}(df_dp)[1,:]
    end
    if FLAG_0000_03
        df_dp                 = CSV.File(string(input_folder,"measurement_operator_discretization_0000_03.csv"); delim=",", header=false)|> DataFrame;
        dp                    = Matrix{Cdouble}(df_dp)[1,:]
    end
    if FLAG_0000_04
        df_dp                 = CSV.File(string(input_folder,"measurement_operator_discretization_0000_04.csv"); delim=",", header=false)|> DataFrame;
        dp                    = dropdims(Matrix{Cdouble}(df_dp),dims=2);
    end
    df_GR_GT                  = CSV.File(string(input_folder,"condensation_rate.csv"); delim=",", header=false)|> DataFrame;
    condensation_rate_all     = Matrix{Cdouble}(Matrix{Cdouble}(df_GR_GT)');
    if FLAG_0000_04
        df_J_GT               = CSV.File(string(input_folder,"J_0000.csv"); delim=",", header=true) |> DataFrame
        nucleation_rate_all   = df_J_GT[:,2]; # for 0000 04
    else
        df_J_GT               = CSV.File(string(input_folder,"J_0000.csv"); delim=",", header=false) |> DataFrame
        nucleation_rate_all   = Matrix{Cdouble}(df_J_GT)[1,:];
    end
    if FLAG_0000_04
        df_wr                 = CSV.File(string(input_folder,"wall_rate_constant.csv"); delim=",", header=false)|> DataFrame;
        wall_rate_expected    = df_wr[:,1];  # for 0000 04
    else
        df_wr                 = CSV.File(string(input_folder,"wall_rate_constant.csv"); delim=",", header=false)|> DataFrame;
        wall_rate_expected    = Matrix{Cdouble}(df_wr)[1,:];
    end
    df_psd_simu           = CSV.File(string(input_folder,"simulated_psd.csv"); delim=",", header=false) |> DataFrame
    psd_simu              = Matrix{Cdouble}(Matrix{Cdouble}(df_psd_simu)');
    df_psd_log_simu       = CSV.File(string(input_folder,"simulated_log_psd.csv"); delim=",", header=false) |> DataFrame
    psd_log_simu          = Matrix{Cdouble}(Matrix{Cdouble}(df_psd_log_simu)');
    maxPSD = maximum(psd_simu)
    minPSD = max(0.001maxPSD,minimum(psd_simu))
    if PLOT_MORE
        titlePSD = @sprintf "particle size distribution (%1.2e,%1.2e)" minimum(psd_simu) maximum(psd_simu)
        cblabelPSD = "density [cm\$^{-3}\$ nm\$^{-1}\$]"
        displayLogData2D(159,tp[1:50:end],1.0e9dp[1:10:end],psd_simu[1:10:end,1:50:end],minPSD,maxPSD,_title=titlePSD,_colorbar_label=cblabelPSD)
        ylabel("diameter [nm]")
        savefig(string(folder,"psd_gt.png"))
        savefig(string(folder,"psd_gt.pdf"))
        close("all")
    end
end

# if some data have been save in matlab format, these functions  might be handy
# ##################################################
# #        transform matlab timestamps             #
# ##################################################
# const MATLAB_EPOCH_MS   = Dates.DateTime(-0001,12,31)
# const MATLAB_EPOCH_DATE = Dates.Date(-0001,12,31)
# function dateFromMatlabInteger(intday::Int64)
#     Date(Dates.UTD(intday+Dates.value(MATLAB_EPOCH_DATE)))
# end
#
# function datetimeFromMatlabInteger(figday::Cdouble)
#     intday = floor(Int64,figday)
#     the_date = dateFromMatlabInteger(intday)
#     hh_tot = 24*(figday-intday)
#     mm_tot = 24*60*(figday-intday)
#     ss_tot = 24*60*60*(figday-intday)
#     hh = floor(Int64,hh_tot)
#     mm = floor(Int64,mm_tot-60*hh)
#     ss = floor(Int64,ss_tot-60*60*hh-60*mm)
#     ms = floor(Int64,1000.0*((ss_tot-60*60*hh-60*mm)-ss))
#     the_time = Time(hh,mm,ss,ms)
#     DateTime(the_date,the_time)
# end
