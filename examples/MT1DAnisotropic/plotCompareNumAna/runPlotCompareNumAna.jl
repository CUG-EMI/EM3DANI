#------------------------------------------------------------------------------
# script to plot MT 1D anisotropic forward modeling results
#------------------------------------------------------------------------------

using Plots
using Printf
#------------------------------------------------------------------------------

include("readMT1DAniResp.jl")
include("readMT3DResp.jl")

# First, read the forward response file
anaRespfile       = "Pek1DModel_31Freqs.resp"
numRespfile_total = "31Freqs_RhoPhs_Total_MUMPS.resp"
numRespfile_bg1   = "31Freqs_RhoPhs_BGM1_MUMPS.resp"
numRespfile_bg2   = "31Freqs_RhoPhs_BGM2_MUMPS.resp"

anaResp = readMT1DAniResp(anaRespfile)
dataInfo, numResp_t = readMT3DResp(numRespfile_total)
(), numResp_bg1   = readMT3DResp(numRespfile_bg1)
(), numResp_bg2   = readMT3DResp(numRespfile_bg2)


# Second, select the forward response
periods = 1 ./ anaResp.freqs

rhoxy0 =  anaResp.appRho[:, 3:3]
phsxy0 =  anaResp.appRho[:, 4:4]
rhoyx0 =  anaResp.appRho[:, 5:5]
phsyx0 =  anaResp.appRho[:, 6:6]

siteNo = Int( round( size(dataInfo.rxLoc, 1)/2 ) )
subRxID = findall(dataInfo.rxID .== siteNo)

rhoxy_t =  numResp_t[subRxID, 3:3]
phsxy_t =  numResp_t[subRxID, 4:4]
rhoyx_t =  numResp_t[subRxID, 5:5]
phsyx_t =  numResp_t[subRxID, 6:6]

rhoxy_bg1 =  numResp_bg1[subRxID, 3:3]
phsxy_bg1 =  numResp_bg1[subRxID, 4:4]
rhoyx_bg1 =  numResp_bg1[subRxID, 5:5]
phsyx_bg1 =  numResp_bg1[subRxID, 6:6]

rhoxy_bg2 =  numResp_bg2[subRxID, 3:3]
phsxy_bg2 =  numResp_bg2[subRxID, 4:4]
rhoyx_bg2 =  numResp_bg2[subRxID, 5:5]
phsyx_bg2 =  numResp_bg2[subRxID, 6:6]


# Third, prepare for plotting

# choose a plotting backend: PlotlyJS, GR, PyPlot, etc
# plotlyjs()
# gr()
pyplot()

# pre-set some plotting attributes
labelfontsize  = 16
legendfontsize = 15
tickfontsize   = 13
linewidth = 2

plt1 = plot(
            periods,
            [rhoxy0 rhoxy_t rhoyx0 rhoyx_t],
            label = ["Analytic, XY"  "FV, XY"  "Analytic, YX"  "FV, YX"],
            seriestype = [:line :scatter :line :scatter],
            color = [:blue :blue :red :red],
            markersize = 6,
            scale=:log10,
            xlabel = "Period (s)",
            ylabel = "App.Res (Ohm)",
            labelfontsize  = labelfontsize,
            legendfontsize = legendfontsize,
            tickfontsize   = tickfontsize,
            minorticks=true,
            grid = :off
            )


plt2 = plot(
            periods,
            [phsxy0 phsxy_t phsyx0 phsyx_t],
            legend = false,
            seriestype = [:line :scatter :line :scatter],
            color = [:blue :blue :red :red],
            markersize = 6,
            xscale=:log10,
            xlabel = "Period (s)",
            ylabel = "Phase (degree)",
            labelfontsize  = labelfontsize,
            legendfontsize = legendfontsize,
            tickfontsize   = tickfontsize,
            minorticks=true,
            grid = :off
            )


rhoxy_t_e = abs.(rhoxy_t-rhoxy0) ./ abs.(rhoxy0) * 100
rhoyx_t_e = abs.(rhoyx_t-rhoyx0) ./ abs.(rhoyx0) * 100
rhoxy_bg1_e = abs.(rhoxy_bg1-rhoxy0) ./ abs.(rhoxy0) * 100
rhoyx_bg1_e = abs.(rhoyx_bg1-rhoyx0) ./ abs.(rhoyx0) * 100
rhoxy_bg2_e = abs.(rhoxy_bg2-rhoxy0) ./ abs.(rhoxy0) * 100
rhoyx_bg2_e = abs.(rhoyx_bg2-rhoyx0) ./ abs.(rhoyx0) * 100

plt3_rhoxye = plot(
            periods,
            [rhoxy_t_e  rhoxy_bg1_e  rhoxy_bg2_e],
            label = ["Total" "BGM1"  "BGM2"],
            linestyle = [:solid :dash :dot],
            linewidth = linewidth,
            color = [:blue :blue :blue],
            xscale=:log10,
            xlabel = "Period (s)",
            ylabel = "App.Res relative error (%)",
            labelfontsize  = labelfontsize,
            legendfontsize = legendfontsize,
            tickfontsize   = tickfontsize,
            minorticks=true,
            gridlinewidth=0.3,
            foreground_color_grid=RGB(0.8, 0.8, 0.8),
            gridstyle=:dash,
            )


plt4_rhoyxe = plot(
            periods,
            [rhoyx_t_e  rhoyx_bg1_e  rhoyx_bg2_e],
            label = ["Total" "BGM1"  "BGM2"],
            linestyle = [:solid :dash :dot],
            linewidth = linewidth,
            color = [:red :red :red],
            xscale=:log10,
            xlabel = "Period (s)",
            ylabel = "App.Res relative error (%)",
            labelfontsize  = labelfontsize,
            legendfontsize = legendfontsize,
            tickfontsize   = tickfontsize,
            minorticks=true,
            gridlinewidth=0.3,
            foreground_color_grid=RGB(0.8, 0.8, 0.8),
            gridstyle=:dash,
            )


phsxy_t_e = phsxy_t-phsxy0
phsyx_t_e = phsyx_t-phsyx0
phsxy_bg1_e = phsxy_bg1-phsxy0
phsyx_bg1_e = phsyx_bg1-phsyx0
phsxy_bg2_e = phsxy_bg2-phsxy0
phsyx_bg2_e = phsyx_bg2-phsyx0

plt5_phsxye = plot(
            periods,
            [phsxy_t_e  phsxy_bg1_e  phsxy_bg2_e],
            legend = false,
            linestyle = [:solid :dash :dot],
            linewidth = linewidth,
            color = [:blue :blue :blue],
            xscale=:log10,
            xlabel = "Period (s)",
            ylabel = "Phase difference (degree)",
            labelfontsize  = labelfontsize,
            legendfontsize = legendfontsize,
            tickfontsize   = tickfontsize,
            minorticks=true,
            gridlinewidth=0.3,
            foreground_color_grid=RGB(0.8, 0.8, 0.8),
            gridstyle=:dash
            )



plt6_phsyxe = plot(
            periods,
            [phsyx_t_e  phsyx_bg1_e  phsyx_bg2_e],
            legend = false,
            linestyle = [:solid :dash :dot],
            linewidth = linewidth,
            color = [:red :red :red],
            xscale=:log10,
            xlabel = "Period (s)",
            ylabel = "Phase difference (degree)",
            labelfontsize  = labelfontsize,
            legendfontsize = legendfontsize,
            tickfontsize   = tickfontsize,
            minorticks=true,
            gridlinewidth=0.3,
            foreground_color_grid=RGB(0.8, 0.8, 0.8),
            gridstyle=:dash
            )
# plot(plt1, plt2, plt3, plt4, layout = (2, 2))
#
savefig(plt1, "rho.pdf")
savefig(plt2, "phs.pdf")
savefig(plt3_rhoxye, "rhoxy_e.pdf")
savefig(plt4_rhoyxe, "rhoyx_e.pdf")
savefig(plt5_phsxye, "phsxy_e.pdf")
savefig(plt6_phsyxe, "phsyx_e.pdf")
#


#=
savefig(plt1, "rho.eps")
savefig(plt2, "phs.eps")
savefig(plt3_rhoxye, "rhoxy_e.eps")
savefig(plt4_rhoyxe, "rhoyx_e.eps")
savefig(plt5_phsxye, "phsxy_e.eps")
savefig(plt6_phsyxe, "phsyx_e.eps")
=#
