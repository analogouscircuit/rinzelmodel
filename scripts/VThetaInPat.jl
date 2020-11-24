module VThetaInPat

include("/home/dahlbom/research/rinzelmodel/vthetamodel/VTheta.jl")
include("/home/dahlbom/research/rinzelmodel/util/SpikeUtils.jl")
include("/home/dahlbom/research/rinzelmodel/util/InputFunc.jl")
using Plots
using Plots.Measures
using Statistics: mean
using .VTheta: vthetarun
using .SpikeUtils: isi, isi_ao, vectorstrength
using .InputFunc: sinstim

params = VTheta.params_default

# -------------------- functions --------------------
w_gauss_full(cf_sd, cf_an; σ=2) = exp( - abs( log(cf_an) - log(cf_sd) )^2 / σ^2 )

function w_gauss_half(cf_sd, cf_an; σ=2)
	if cf_an >= cf_sd
		return w_gauss_full(cf_sd, cf_an; σ=σ) 
	else
		return 0.0
	end
end

# -------------------- Parameters --------------------
f0 = 200.0
num_h = 3
mistuned_h = 2
mt = 1.00
amp = 1.0 
σ=2.0
dur = 250.0
fs = 100e3
dt = 1/fs * 1000.0
bins = range(0.0, 50.0, step=0.025)
freqs_in = [ f0*k for k ∈ 1:num_h ]
freqs_in[mistuned_h] *= mt
cfs_sd = exp.( collect(range(log(100), log(maximum(freqs_in)), length=6)) )


# -------------------- Generate Inputs --------------------
channels = []
for freq ∈ freqs_in
	push!(channels, sinstim(freq/1000.0; ϕ=0.0, hwr=true))
end

inputs = []
plts_in = []
t = collect(range(0, dur-dt, step=dt))
w(cf_sd, cf_an) = w_gauss_full(cf_sd, cf_an; σ=σ)
for cf_sd ∈ cfs_sd
	input(t) = sum([amp * w(cf_sd, cf_an) * channels[k](t) for (k, cf_an) ∈ enumerate(freqs_in)])
	push!(inputs, input)	
	push!(plts_in, plot(t, input.(t), legend=:none))
end
p1 = plot(plts_in..., layout=(length(plts_in), 1), legend=:none)

## Test plot
# t = collect(range(0.0, 100.0, step=0.01))
# p = plot(t, inputs[1])
# for (k, input) ∈ enumerate(inputs[2:end])
# 	plot!(p, t, input.(t) .+ 5*k)
# end

# -------------------- Run through VΘ Model --------------------
isis_chan = []
plts_chan = []
hists_chan = [] 
for infunc ∈ inputs
	v0 = 0.1
	θ0 = 0.2
	(spikes, v, θ) = vthetarun(v0, θ0, dur, dt, infunc, params)
	println(length(spikes))
	isis = isi(spikes)
	push!(isis_chan, isis)
	push!(hists_chan, histogram(isis, bins=bins))
	push!(plts_chan, plot(t, v, legend=:none))
end
p2 = plot(plts_chan..., layout=(length(plts_chan), 1), legend=:none)
p3 = plot(hists_chan..., layout=(length(hists_chan), 1), legend=:none)

p = plot(p1, p2, p3,  layout=(1,3), legend=:none, size=(1800, 800))

isis_all = vcat(isis_chan...)
p_hist_all = histogram(isis_all, bins=bins)



end
