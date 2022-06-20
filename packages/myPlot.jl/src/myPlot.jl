module myPlot

using PyPlot
using Printf

export imshowData, displayPSD, displayGR, displayLogData2D, displayLogData2DDownSample
export imshowDataPolar, displayCov, setVerticalColorbar


# display in polar coordinates
function imshowDataPolar(figNum::Int64,r::Array{Cdouble,1},θ::Array{Cdouble,1},H::Array{Cdouble,2};subFig::Int64=111,cb_label::String="gain [a.u.]",cb_ax_loc::Tuple{Cdouble,Cdouble}=(0.0,0.0),cb_ax_size::Tuple{Cdouble,Cdouble}=(0.02,0.25),rlabel_position::Cdouble=-22.5)
    fig = figure(figNum) # ,figsize=[9,5])
    ax = subplot(subFig,polar=true)
    ax.set_rlabel_position(rlabel_position)
    # ax.set_rticks([0.97, 0.98, 0.99, 1.0])
    pcm = ax.pcolormesh(θ,r,H,edgecolors="face")
    ylim(r[1],r[end])
    cax = fig.add_axes([cb_ax_loc[1], cb_ax_loc[2], cb_ax_size[1], cb_ax_size[2]])
    cb = fig.colorbar(pcm, orientation="vertical", cax=cax, shrink=0.6)
    cb.set_label(cb_label, fontsize=10)

    fig,ax,pcm,cax,cb
end



# display data on the actual grid (taking into account the pixel coordinates)
function imshowData(figNum::Int64,_t::Array{Cdouble,1},_h::Array{Cdouble,1},Z::Array{Cdouble,2};_norm=:NoNorm,_vmin=0.0,_vmax=1.0,_edgecolors="face",_shading="None", _sub=111)
    # get the focus on the figure or create a figure
    fig=figure(figNum)

    # get the current axis
    ax=subplot(_sub)


    # display image if possible
    if (length(_h),length(_t))==size(Z)
        Z_display = [Z zeros(Cdouble,length(_h)); zeros(Cdouble,1,length(_t)+1)];
        t_display = [_t; _t[end]+(_t[end]-_t[end-1])];
        h_display = [_h; _h[end]+(_h[end]-_h[end-1])];
        if (_norm==:NoNorm)
            pcm = pcolormesh(repeat(t_display,1,length(h_display))',repeat(h_display,1,length(t_display)),Z_display,edgecolors=_edgecolors,shading=_shading)
        elseif (_norm==:LogNorm)
            pcm = pcolormesh(repeat(t_display,1,length(h_display))',repeat(h_display,1,length(t_display)),Z_display,edgecolors=_edgecolors,shading=_shading,norm=matplotlib.colors.LogNorm(vmin=_vmin,vmax=_vmax))
        else
            pcm = pcolormesh(repeat(t_display,1,length(h_display))',repeat(h_display,1,length(t_display)),Z_display,edgecolors=_edgecolors,shading=_shading,norm=matplotlib.colors.Normalize(vmin=_vmin,vmax=_vmax)) #shading="gouraud"
        end
    elseif (length(_t),length(_h))==size(Z)
        Z_display = [Z zeros(Cdouble,length(_t)); zeros(Cdouble,1,length(_h)+1)]
        t_display = [_t; _t[end]+(_t[end]-_t[end-1])]
        h_display = [_h; _h[end]+(_h[end]-_h[end-1])]
        if (_norm==:NoNorm)
            pcm = pcolormesh(repeat(t_display,1,length(h_display))',repeat(h_display,1,length(t_display)),Z_display',edgecolors=_edgecolors,shading=_shading)
        elseif (_norm==:LogNorm)
            pcm = pcolormesh(repeat(t_display,1,length(h_display))',repeat(h_display,1,length(t_display)),Z_display',edgecolors=_edgecolors,shading=_shading,norm=matplotlib.colors.LogNorm(vmin=_vmin,vmax=_vmax))
        else
            pcm = pcolormesh(repeat(t_display,1,length(h_display))',repeat(h_display,1,length(t_display)),Z_display',edgecolors=_edgecolors,shading=_shading,norm=matplotlib.colors.Normalize(vmin=_vmin,vmax=_vmax)) #shading="gouraud"
        end
    else
        throw(@sprintf "Cannot display the image because of a size mismatch: (%i,%i)!=%i,%i" size(Z,1) size(Z,2) length(_h) length(_t))
    end

    #return axis handler
    fig,ax,pcm
end


function displayPSD(figNum::Int64,t::Array{Cdouble,1},d::Array{Cdouble,1},PSD::Array{Cdouble,2}, _sub=111)
    PSD_r=PSD.*(PSD.>0.0).+1.0e-16
    min_PSD,max_PSD = extrema(PSD)
    fig,ax,pcm = imshowData(figNum,t,d,log.(PSD_r),_norm=:Normalize,_vmin=0.0,_vmax=log(max_PSD),_edgecolors="face",_shading="None", _sub=_sub)
    yscale("log")
    # unfortunately the following commented code does not work properly in julia 0.4.5, but it is fixed later, so it will be available at some point
    # yform  = ax[:yaxis][:set_major_formatter]
    # stick = matplotlib[:ticker][:ScalarFormatter]
    # yform(stick())
    # ax[:ticklabel_format](style="sci",axis="y",scilimits=(0,0))
    s = @sprintf "size distribution (%1.2e,%1.2e)" min_PSD max_PSD
    title(s)
    xlabel("time [h]")
    ylabel("diameter [m]")
    cbar=colorbar()
    cbar.set_label("concentration [log(# cm\$^{-3}\$)]")
    cbar.formatter.set_powerlimits((-1,2))
    cbar.update_ticks()
    fig,ax,cbar,pcm
end

function displayLogData2D(figNum::Int64,t::Array{Cdouble,1},d::Array{Cdouble,1},PSD::Array{Cdouble,2},min_PSD::Cdouble,max_PSD::Cdouble;_title::AbstractString="size distribution",_xlabel::AbstractString="time [h]",_ylabel::AbstractString="diameter [m]",_colorbar_label::AbstractString="concentration [log(# cm\$^{-3}\$)]", _sub=111)
    PSD_r=PSD.*(PSD.>0.0).+1.0e-16
    fig,ax,pcm = imshowData(figNum,t,d,PSD_r,_norm=:LogNorm,_vmin=min_PSD.+1.0e-16,_vmax=max_PSD,_edgecolors="face",_shading="None", _sub=_sub)
    yscale("log")
    title(_title)
    xlabel(_xlabel)
    ylabel(_ylabel)
    cbar=colorbar()
    cbar.set_label(_colorbar_label)
    cbar.update_ticks()
    fig,ax,cbar,pcm
end

function displayLogData2DDownSample(figNum::Int64,t::Array{Cdouble,1},d::Array{Cdouble,1},PSD::Array{Cdouble,2},min_PSD::Cdouble,max_PSD::Cdouble;_title::AbstractString="size distribution",_xlabel::AbstractString="time [h]",_ylabel::AbstractString="diameter [m]",_colorbar_label::AbstractString="concentration [log(# cm\$^{-3}\$)]",lmax::Int64=400,cmax::Int64=400, _sub=111)
    Z_display = copy(PSD)
    t_display = copy(t)
    d_display = copy(d)
    # down sample
    if (size(Z_display,1)>lmax) & (size(Z_display,2)>cmax)
        # down sample until it reaches a displayable size
        ker = (1.0/16.0)*[1.0 2.0 1.0; 2.0 4.0 2.0; 1.0 2.0 1.0]
        while (size(Z_display,1)>lmax) & (size(Z_display,2)>cmax)
            Z_display = conv2(ker,Z_display)[2:end-1,2:end-1]
            Z_display = copy(Z_display[1:2:end,1:2:end])
            t_display = copy(t_display[1:2:end])
            d_display = copy(d_display[1:2:end])
        end
    end
    if (size(Z_display,1)>lmax) & (size(Z_display,2)<=cmax)
        # down sample in the first dimension
        ker = (1.0/4.0)*[1.0 2.0 1.0]
        while (size(Z_display,1)>lmax)
            Z_display = conv2(ker',Z_display)[2:end-1,:]
            Z_display = copy(Z_display[1:2:end,:])
            if (length(t)==size(PSD,1))
                t_display = copy(t_display[1:2:end])
            else
                d_display = copy(d_display[1:2:end])
            end
        end
    end
    if (size(Z_display,1)<=lmax) & (size(Z_display,2)>cmax)
        # down sample in the second dimension
        ker = (1.0/4.0)*[1.0 2.0 1.0]
        while (size(Z_display,2)>cmax)
            Z_display = conv2(ker,Z_display)[:,2:end-1]
            Z_display = copy(Z_display[:,1:2:end])
            if (length(t)==size(PSD,2))
                t_display = copy(t_display[1:2:end])
            else
                d_display = copy(d_display[1:2:end])
            end
        end
    end

    Z_display=Z_display.*(Z_display.>0.0).+1.0e-16
    fig,ax,pcm = imshowData(figNum,t_display,d_display,Z_display,_norm=:LogNorm,_vmin=min_PSD.+1.0e-16,_vmax=max_PSD,_edgecolors="face",_shading="None", _sub=_sub)
    yscale("log")
    title(_title)
    xlabel(_xlabel)
    ylabel(_ylabel)
    cbar=colorbar()
    cbar.set_label(_colorbar_label)
    cbar.update_ticks()
    fig,ax,cbar,pcm
end

function displayGR(figNum::Int64,t::Array{Cdouble,1},d::Array{Cdouble,1},GR_::Array{Cdouble,2}, _sub=111)
    GR_r=GR_.*(GR_.>0.0)
    min_GR,max_GR = extrema(GR_)
    imshowData(figNum,t,d,GR_r,_norm=:Normalize,_vmin=0.0,_vmax=max_GR,_edgecolors="face",_shading="None", _sub=_sub)
    yscale("log")
    s = @sprintf "growth rate (%1.2e,%1.2e)" min_GR max_GR
    title(s)
    xlabel("time [h]")
    ylabel("diameter [m]")
    cbar=colorbar()
    cbar.set_label("growth rate [m.s\$^{-1}\$]")
    cbar.formatter.set_powerlimits((-1,2))
    cbar.update_ticks()
end

function displayCov(figNum::Int64,r::Array{Cdouble,1},Γ::Array{Cdouble,2};_sub::Int64=111,fontsize_ticks::Int64=14)
    min_Γ,max_Γ = extrema(Γ);
    fig,ax,pcm = imshowData(figNum,r,r,Γ;_norm=:NoNorm,_vmin=min_Γ,_vmax=max_Γ,_edgecolors="face",_shading="None", _sub=_sub);
    xticks(fontsize=fontsize_ticks)
    yticks(fontsize=fontsize_ticks)
    fig,ax,pcm
end

function setVerticalColorbar(fig::Figure,pcm::PyPlot.PyObject,x::Cdouble,y::Cdouble,dx::Cdouble,dy::Cdouble,slabel::String;fontsize_label::Int64=10,fontsize_ticks::Int64=10,color::String="white",_power_lim::Bool=true)
    rc("ytick",color=color)
    cax = fig.add_axes([x, y, dx, dy])
    cb = fig.colorbar(pcm, orientation="vertical", cax=cax)
    cb.set_label(slabel, fontsize=fontsize_label, color=color) # 
    cb.ax.yaxis.set_tick_params(color=color)
    cb.ax.tick_params(labelsize=fontsize_ticks)
    cb.outline.set_edgecolor(color)
    if _power_lim
        cb.formatter.set_powerlimits((-1,2))
        cb.ax.yaxis.offsetText.set_size(fontsize_ticks)
    end
    cb.update_ticks()
    rc("ytick",color="black") # back to black
    cb
end

end
