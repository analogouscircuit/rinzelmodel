module VThetaFilterPeri


using AuditoryFilters
using Plots

include("/home/dahlbom/research/AdaptiveSieves/util/StimGen.jl")
include("/home/dahlbom/research/AdaptiveSieves/util/DahlPlot.jl")
include("/home/dahlbom/research/AdaptiveSieves/util/EdgePitchUtils.jl")
include("/home/dahlbom/research/AdaptiveSieves/core/PitchModelPeaks.jl")
include("/home/dahlbom/research/AdaptiveSieves/scripts/ChannelShiftUtils.jl")
include("/home/dahlbom/research/rinzelmodel/vthetamodel/VTheta.jl")
include("/home/dahlbom/research/rinzelmodel/util/InputFunc.jl")
include("/home/dahlbom/research/rinzelmodel/util/SpikeUtils.jl")

using .ChannelShiftUtils

# -------------------- Parameters --------------------
# Periphery parameters
fs = 50e3
num_chan = 50
f_lo = 100.0
fb = make_erb_filterbank(fs, num_chan, f_lo)
cfs = fb.center_frequencies
cfs_out = exp.( range( log(50.0), log(5000), length=25 ) )
println(cfs_out)

# Stimulus parameters
f0 = 200.0
dur = 0.15
dur_used = 0.05
mth = 4
mtp = 1.05
mistunings = [1.0 for _ ∈ 1:12]
mistunings[mth] = mtp
stim = StimGen.hc(dur, f0, fs; mistunings=mistunings)

# Model Parameters
skew = 0.0
bump = :gauss
σ0 = 5.0e-4
vtparams = VTheta.VThetaParams()


# -------------------- Main Script --------------------
# Filter and process channels
chans = filt(fb, stim)
idx_max = Int(round(dur_used*fs))
println("max idx: $idx_max")
chans = chans[1:idx_max,:]
ts = collect(0:idx_max-1) ./ fs

# Half-wave rectify
chans = map( x -> x > 0.0 ? x : 0.0, chans)
chans_hwr = copy(chans)
chans_hwr ./= maximum(chans_hwr)

p0 = DahlPlot.stagger_plot_mat(ts, chans; offset=0.5, rev=true)

# Autocorrelate
# chans = ac_channels(chans)
# chans_ac = copy(chans)
# peaks_cf_ac = Float64[]
# for n ∈ 1:size(chans_ac)[2]
# 	peaks = PitchModelPeaks.findpeaks(ts, chans_ac[:,n], 1/f0, fs, 1)
# 	push!(peaks_cf_ac, peaks[1])
# end

# Run each channel through V-Θ model
isis = []
v0 = 0.1
θ0 = 0.2
factor = 100.0
ts .*= 1000.0
which_chan = 38
# v = zeros(size(chans))
# th = zeros(size(chans))
for n ∈ which_chan:which_chan#size(chans)[2]
	global v
	global th
	ifunc(x) = InputFunc.interpds(x, ts, factor .* chans[:,n])
	spikes, v, th = VTheta.vthetarun(v0, θ0, dur_used*1000, 1000.0/fs, ifunc, vtparams)
	push!(isis, SpikeUtils.isi(spikes))
end

p1 = plot(v)
plot!(p1, th)
p2 = plot(ts, factor .* chans[:,which_chan])
p = plot(p1, p2, layout=(2,1))

end
