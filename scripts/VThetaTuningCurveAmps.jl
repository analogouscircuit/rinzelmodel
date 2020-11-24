module VThetaTuningCurve

include("/home/dahlbom/research/rinzelmodel/vthetamodel/VTheta.jl")
include("/home/dahlbom/research/rinzelmodel/util/SpikeUtils.jl")
using Plots
using Statistics: mean
using .VTheta: vthetarun
using .SpikeUtils: isi

params = VTheta.params_default


sinusoids(t, f0, p) = max(0.0, sin(2π*f0*t) + sin(2π*f0*p*t))

mistunings = collect(range(0.88, stop=1.12, step=0.005))
amps = collect(range(0.2, 0.5, step=0.01))
# amps=amps[25]
dt = 1/50000*1000
dur = 2.0*1000
f0 = 200.0/1000

spikes_all = []
isis_all = []
v_all = []
θ_all = []
plots = []
for (k, amp) ∈ enumerate(amps)
	println("$k of $(length(amps))...")
	ests = []
	for mt ∈ mistunings
		v0 = 0.1
		θ0 = 0.3
		infunc(t) = amp*sinusoids(t, f0, mt)
		(spikes, v, θ) = vthetarun(v0, θ0, dur, dt, infunc, params)
		isis = isi(spikes)
		isis1 = filter(x -> 0.9*5.0 < x < (1/.9)*5.0, isis)
		est = 1000.0/mean(isis1)
		push!(ests, est)
		# push!(spikes_all, spikes)
		# push!(θ_all, θ)
		# push!(v_all, v)
		# push!(isis_all, isis)
	end
	push!(plots, plot(mistunings, ests, legend=:none, title="$amp"))
end

p = plot(plots...)

end
