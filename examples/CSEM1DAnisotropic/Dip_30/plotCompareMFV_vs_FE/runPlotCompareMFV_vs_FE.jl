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
anaRespfile = "canonical1D_FE.resp"
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

label = ["FE Ex",  "MFV Ex",  "FE Ey",  "MFV Ey", "FE Ez",  "MFV Ez",
         "FE Bx",  "MFV Bx",  "FE By",  "MFV By", "FE Bz",  "MFV Bz"]

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


savefig(plt1, "amp.pdf")
savefig(plt2, "phs.pdf")

#
savefig(plt1, "amp.eps")
savefig(plt2, "phs.eps")
#
