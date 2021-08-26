#------------------------------------------------------------------------------
# script to plot CSEM forward modeling results
#------------------------------------------------------------------------------

using Plots
using Printf
using DSP   # providing unwrap function
#------------------------------------------------------------------------------

include("readCSEMResp.jl")
include("getAmpPhs.jl")


# First, read the forward response file
anaRespfile = "canonical1D_emmod.resp"
numRespfile = "0.25Hz_inline_broadside_MUMPS.resp"

dataInfo, anaResp = readCSEMResp(anaRespfile)
(), numResp = readCSEMResp(numRespfile)

# Second, select the forward response
# ------------------------ inline ------------------------- #
freqID = 1
txID   = 1
offset = dataInfo.rxLoc[:, 2]
txRange = [1000, 8000]
rxID = findall(x -> x>=txRange[1] && x<=txRange[2], offset)
gap=2;  rxID=rxID[1:gap:end]


# Ey
Ey0_amp, Ey0_phs = getAmpPhs(anaResp.Ey[rxID, txID, freqID])
Ey1_amp, Ey1_phs = getAmpPhs(numResp.Ey[rxID, txID, freqID])

# Bx
Bx0_amp, Bx0_phs = getAmpPhs(anaResp.Bx[rxID, txID, freqID])
Bx1_amp, Bx1_phs = getAmpPhs(numResp.Bx[rxID, txID, freqID])

# Ez
Ez0_amp, Ez0_phs = getAmpPhs(anaResp.Ez[rxID, txID, freqID])
Ez1_amp, Ez1_phs = getAmpPhs(numResp.Ez[rxID, txID, freqID])
# --------------------------------------------------------- #

# ----------------------- broadside ----------------------- #
freqID = 1
txID   = 2
offset = dataInfo.rxLoc[:, 2]
txRange = [1000, 8000]
rxID = findall(x -> x>=txRange[1] && x<=txRange[2], offset)
gap=2;  rxID=rxID[1:gap:end]


# Ex
Ex0_amp, Ex0_phs = getAmpPhs(anaResp.Ex[rxID, txID, freqID])
Ex1_amp, Ex1_phs = getAmpPhs(numResp.Ex[rxID, txID, freqID])

# By
By0_amp, By0_phs = getAmpPhs(anaResp.By[rxID, txID, freqID])
By1_amp, By1_phs = getAmpPhs(numResp.By[rxID, txID, freqID])

# Bz
Bz0_amp, Bz0_phs = getAmpPhs(anaResp.Bz[rxID, txID, freqID])
Bz1_amp, Bz1_phs = getAmpPhs(numResp.Bz[rxID, txID, freqID])
# --------------------------------------------------------- #


# Third, prepare for plotting

# choose a plotting backend: PlotlyJS, GR, PyPlot, etc
# plotlyjs()
# gr()
pyplot()

# pre-set some plotting attributes
labelfontsize  = 16
legendfontsize = 10
tickfontsize   = 13
linewidth = 1.5

label = ["Ana. Ex",  "MFV Ex",  "Ana. Ey",  "MFV Ey", "Ana. Ez",  "MFV Ez",
         "Ana. Bx",  "MFV Bx",  "Ana. By",  "MFV By", "Ana. Bz",  "MFV Bz"]

c1 = RGB(0.5, 0, 0)
c2 = RGB(1.0, 0.3125, 0)
c3 = RGB(0.5625, 1.0, 0.4375)
c4 = RGB(0.125, 1.0, 0.875)
c5 = RGB(0, 0.375, 1.0)
c6 = RGB(0, 0, 0.5625)

x = offset[rxID]/1000


#------------ 1st, Amplitude ------------#
plt1 = plot(
            x,
            log10.([Ex0_amp Ex1_amp Ey0_amp Ey1_amp Ez0_amp Ez1_amp Bx0_amp Bx1_amp By0_amp By1_amp Bz0_amp Bz1_amp]),
            label = reshape(label, 1, 12),
            #legend = false,
            seriestype = [:line :scatter],
            # palette = cgrad(:jet),
            color = [c1 c1 c2 c2 c3 c3 c4 c4 c5 c5 c6 c6],
            markerstrokecolor = [c1 c1 c2 c2 c3 c3 c4 c4 c5 c5 c6 c6],
            linewidth = linewidth,
            markersize = 5,
            xlabel = "Offset (km)",
            ylabel = "log10 (Amplitude (V/m or T))",
            labelfontsize  = labelfontsize,
            legendfontsize = legendfontsize,
            tickfontsize   = tickfontsize,
            minorticks=true,
            grid = :off,
            ylims = (-20, -10),
            xticks = collect(1:8),
            yticks = collect(-20:2:-10)
            )



#-------------- 2nd, Phase --------------#
plt2 = plot(
            x,
            [Ex0_phs Ex1_phs Ey0_phs Ey1_phs Ez0_phs Ez1_phs Bx0_phs Bx1_phs By0_phs By1_phs Bz0_phs Bz1_phs],
            legend = false,
            seriestype = [:line :scatter],
            color = [c1 c1 c2 c2 c3 c3 c4 c4 c5 c5 c6 c6],
            markerstrokecolor = [c1 c1 c2 c2 c3 c3 c4 c4 c5 c5 c6 c6],
            linewidth = linewidth,
            markersize = 5,
            xlabel = "Offset (km)",
            ylabel = "Phase (degree)",
            labelfontsize  = labelfontsize,
            legendfontsize = legendfontsize,
            tickfontsize   = tickfontsize,
            minorticks=true,
            grid = :off,
            xticks = collect(1:8),
            )



#------------ 3rd, Amplitude error ------------#
Ex_amp_e = abs.(Ex1_amp - Ex0_amp) ./ abs.(Ex0_amp) * 100
Ey_amp_e = abs.(Ey1_amp - Ey0_amp) ./ abs.(Ey0_amp) * 100
Ez_amp_e = abs.(Ez1_amp - Ez0_amp) ./ abs.(Ez0_amp) * 100
Bx_amp_e = abs.(Bx1_amp - Bx0_amp) ./ abs.(Bx0_amp) * 100
By_amp_e = abs.(By1_amp - By0_amp) ./ abs.(By0_amp) * 100
Bz_amp_e = abs.(Bz1_amp - Bz0_amp) ./ abs.(Bz0_amp) * 100

plt3 = plot(
            x,
            [Ex_amp_e Ey_amp_e Ez_amp_e Bx_amp_e By_amp_e Bz_amp_e],
            legend = false,
            color = [c1 c2 c3 c4 c5 c6],
            linewidth = linewidth,
            xlabel = "Offset (km)",
            ylabel = "Amplitude relative error (%)",
            labelfontsize  = labelfontsize,
            legendfontsize = legendfontsize,
            tickfontsize   = tickfontsize,
            minorticks=true,
            gridlinewidth=0.2,
            #foreground_color_grid=:gray,
            foreground_color_grid=RGB(0.8, 0.8, 0.8),
            gridstyle=:dash,
            xticks = collect(1:8),
            )



#------------ 4th, Phase error ------------#
Ex_phs_e = Ex1_phs - Ex0_phs
Ey_phs_e = Ey1_phs - Ey0_phs
Ez_phs_e = Ez1_phs - Ez0_phs
Bx_phs_e = Bx1_phs - Bx0_phs
By_phs_e = By1_phs - By0_phs
Bz_phs_e = Bz1_phs - Bz0_phs

plt4 = plot(
            x,
            [Ex_phs_e Ey_phs_e Ez_phs_e Bx_phs_e By_phs_e Bz_phs_e],
            legend = false,
            color = [c1 c2 c3 c4 c5 c6],
            linewidth = linewidth,
            xlabel = "Offset (km)",
            ylabel = "Phase difference (degree)",
            labelfontsize  = labelfontsize,
            legendfontsize = legendfontsize,
            tickfontsize   = tickfontsize,
            minorticks=true,
            gridlinewidth=0.2,
            #foreground_color_grid=:gray,
            foreground_color_grid=RGB(0.8, 0.8, 0.8),
            gridstyle=:dash,
            xticks = collect(1:8),
            )



savefig(plt1, "amp.pdf")
savefig(plt2, "phs.pdf")
savefig(plt3, "ampe.pdf")
savefig(plt4, "phse.pdf")

#
savefig(plt1, "amp.eps")
savefig(plt2, "phs.eps")
savefig(plt3, "ampe.eps")
savefig(plt4, "phse.eps")
#
