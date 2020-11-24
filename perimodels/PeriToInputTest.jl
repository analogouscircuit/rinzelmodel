module PeriToInputTest

using Plots; pyplot()
include("/home/dahlbom/research/rinzelmodel/perimodels/PeriToInput.jl")
using .PeriToInput: ZilanyParams, sigtospikedata, spikedatatochannels, channelstoinputfuncs
include("/home/dahlbom/research/AdaptiveSieves/util/DahlPlot.jl")
using .DahlPlot: stagger_plot

f0 = 155.0
fs = 100e3
db = 90.0
params = ZilanyParams(fs, db, (10, 7, 3), (125.0, 6000.0, 29), "human", 1)
cfs_sd = exp.(range(log(125.0), stop=log(3000.0), length=20))
t = collect(range(0.0, stop=0.1, step=1/params.fs))
sig = max.( 0.0, sin.(2π*f0 .* t) )

(cfs_an, spikes) =  sigtospikedata(sig, params) |> spikedatatochannels
g_syns = channelstoinputfuncs(cfs_sd, cfs_an, spikes; g_E=4.5)
inputs = [f.(t) for f ∈ g_syns]

p_stagger = stagger_plot(t, inputs; offset=2.0, xlims=(0.0, t[end]))
p_sum = plot(t, sum(inputs))

l = @layout [a{0.8h}; b]
p = plot(p_stagger, p_sum, layout=l)

end
