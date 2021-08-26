#------------------------------------------------------------------------------
# script to plot MT 1D anisotropic forward modeling results
#------------------------------------------------------------------------------

using Plots
using Printf
#------------------------------------------------------------------------------

include("readMT3DResp.jl")

# First, read the forward response file
respfile  = "commemi3d2_RhoPhs.resp"
dataInfo, resp = readMT3DResp(respfile)

# Second, select the forward response
freqNo = 2  # The second frequency.
subIdx = findall(dataInfo.freqID .== freqNo)


# Third, prepare for plotting
perValueString = @sprintf("%g", 1 / dataInfo.freqArray[freqNo])
titleString = "Period = " * perValueString * " sec"

# choose a plotting backend: PlotlyJS, GR, PyPlot, etc
# plotlyjs()
# gr()
# pyplot()

# pre-set some plotting attributes
labelx = "Position [km]"
labelfontsize  = 16
legendfontsize = 15
tickfontsize   = 13
linewidth = 2

pos = dataInfo.rxLoc[:, 2]/1000


# rhoxy
plt1 = plot(
            pos,
            resp[subIdx, 3],
            legend = false,
            seriestype = :scatter,
            color = :blue,
            yscale=:log10,
            xlabel = labelx,
            ylabel = "RhoXY [Ohm]",
            title = titleString,
            minorticks=true,
            gridlinewidth=0.3,
            foreground_color_grid=:gray,
            gridstyle=:dash,
            )


# phsxy
plt2 = plot(
            pos,
            resp[subIdx, 4],
            legend = false,
            seriestype = :scatter,
            color = :blue,
            xlabel = labelx,
            ylabel = "Phase XY [degree]",
            title = titleString,
            minorticks=true,
            gridlinewidth=0.3,
            foreground_color_grid=:gray,
            gridstyle=:dash,
            )



# rhoyx
plt3 = plot(
            pos,
            resp[subIdx, 5],
            legend = false,
            seriestype = :scatter,
            color = :blue,
            yscale=:log10,
            xlabel = labelx,
            ylabel = "RhoYX [Ohm]",
            title = titleString,
            minorticks=true,
            gridlinewidth=0.3,
            foreground_color_grid=:gray,
            gridstyle=:dash,
            )


# phsyx
plt4 = plot(
            pos,
            resp[subIdx, 6],
            legend = false,
            seriestype = :scatter,
            color = :blue,
            xlabel = labelx,
            ylabel = "Phase YX [degree]",
            title = titleString,
            minorticks=true,
            gridlinewidth=0.5,
            foreground_color_grid=:gray,
            gridstyle=:dash,
            )


#plot(plt1, plt2, plt3, plt4, layout = (2, 2))
#
savefig(plt1, "rhoxy.pdf")
savefig(plt2, "phsxy.pdf")
savefig(plt3, "rhoyx.pdf")
savefig(plt4, "phsyx.pdf")
#
