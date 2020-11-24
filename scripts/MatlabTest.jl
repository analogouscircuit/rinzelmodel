module MatlabTest

# Packages
using Plots
using PyCall
using DataFrames

# PyCall Imports
co = pyimport("cochlea"); setdb = co.set_dbspl

# Local Modules
include("/home/dahlbom/research/AdaptiveSieves/util/PeakPick.jl")
using .PeakPick: findpeaknear, peakinterp
include("/home/dahlbom/research/rinzelmodel/scripts/MatlabModel.jl")
using .MatlabModel: traces 
include("/home/dahlbom/research/rinzelmodel/util/StimsSamp.jl")
using .StimsSamp: nepstim, sinstim
include("/home/dahlbom/research/rinzelmodel/util/SpikeUtils.jl")

# -------------------- Generate Stimulus --------------------
dur_on = 10.0 
dur_off = 0.05
dur_tot = dur_on + 2dur_off
fs = 50e3
f0 = 800.0
db = 70.0
stim = nepstim(f0, fs, dur_on, dur_off; kind=:hp) 
# stim = sinstim(f0, fs, dur_on, dur_off)
# for k ∈ 2:8
# 	k == 4 ? p = 1.08 : p = 1.0
# 	stim .+= sinstim(k*f0*p, fs, dur_on, dur_off)
# end
stim = setdb(stim, db)

# -------------------- Run through Matlab Model --------------------
@time sd_data = traces(stim, fs)
cfs = unique(sd_data[!, :cfs])

# sd_fs = sd_data[sd_data.cfs .== cf, :trace]
t = collect(0.0:1/fs:dur_tot);
p = plot(t, sd_data[sd_data.cfs .== 100.0, :trace][1])

p_chan = []
spikes_chan = []
isis_chan = []
isis_ao_chan = []
for cf ∈ cfs
	ps =[]
	spikes = []
	for (k, trace) ∈ enumerate(sd_data[sd_data.cfs .== cf, :trace])
		push!(spikes, SpikeUtils.findspikes(trace, 0.0, dur_tot; scale=:s))
		push!(ps, plot(t, trace, legend=false))
	end
	push!(p_chan, ps)
	push!(spikes_chan, spikes)
	push!(isis_chan, vcat([ SpikeUtils.isi(spks) for spks ∈ spikes]...))
	push!(isis_ao_chan, vcat([ SpikeUtils.isi_ao(spks) for spks ∈ spikes]...))
end
isis_summary = vcat( isis_chan... )
isis_ao_summary = vcat( isis_ao_chan... )

@time begin
	τ_min = 0.002
	τ_max = 0.012
	hist = histogram(isis_summary, bins=range(τ_min, τ_max, step=0.0001))
	hist_ao = histogram(isis_ao_summary, bins=range(τ_min, τ_max, step=0.0001))
	τs = collect(τ_min:0.00001:τ_max);
	smoothed = SpikeUtils.smoothhistfunc(isis_summary).(τs)
	smoothed_ao = SpikeUtils.smoothhistfunc(isis_ao_summary).(τs)
	maxval = maximum(smoothed)
	# i = findpeaknear(τs, smoothed, 1/f0; tol=0.2)
	# peak = peakinterp(τs[i-1:i+1], smoothed[i-1:i+1])
	maxval_ao = maximum(smoothed_ao)
	# i = findpeaknear(τs, smoothed_ao, 1/f0)
	# peak_ao = peakinterp(τs[i-1:i+1], smoothed_ao[i-1:i+1])
	hist_smooth = plot(τs, smoothed)#, title="$(1/peak)")
	hist_ao_smooth = plot(τs, smoothed_ao)#, title="$(1/peak_ao)")
	for k ∈ 1:8
		plot!(hist_smooth, [k/f0, k/f0], [0.0, 1.05*maxval], ls=:dash, color=:black)
		# plot!(hist_smooth, [peak, peak], [0.0, 1.05*maxval], color=:red)
		plot!(hist_ao_smooth, [k/f0, k/f0], [0.0, 1.05*maxval_ao], ls=:dash, color=:black)
		# plot!(hist_ao_smooth, [peak_ao, peak_ao], [0.0, 1.05*maxval_ao], color=:red)
	end
end

l = @layout [a b ; c d]
p = plot(hist, hist_smooth, hist_ao, hist_ao_smooth; layout=l, size=(1200,900))
# println("Estimated Pitch (first order): $(1/peak)")
# println("Estimated Pitch (all order):   $(1/peak_ao)")

# spikes_all = vcat(vcat(spikes_chan...)...) |> sort
# f0s = 190.0:0.05:210.0
# strengths = []
# for f0 ∈ f0s
# 	push!(strengths, SpikeUtils.vectorstrength(spikes_all, f0)[1])
# end
# p2 = plot(f0s, strengths)
# i = findpeaknear(f0s, strengths, 200.0)
# peak_vs = peakinterp(f0s[i-1:i+1], strengths[i-1:i+1])
# println("Estimate Pitch (vs): $peak_vs")
# plot!(p2, [peak_vs, peak_vs], [0, 1.0])

end
