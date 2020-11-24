module HRStochTuningCurve

using Plots; pyplot()
using DifferentialEquations
using Statistics

include("/home/dahlbom/research/rinzelmodel/util/SpikeUtils.jl")
include("/home/dahlbom/research/rinzelmodel/util/Inputs.jl")
include("/home/dahlbom/research/rinzelmodel/sdmodels/HRDirectInput.jl")
using .SpikeUtils: findspikes, isi, phasefromspikes
using .Inputs: ramp, bump, bumps, hc
using .HRDirectInput: hr, σ_hr

################################################################################
# Solve problem
################################################################################

# input = t -> ramp(t, 10.0, 10.2, 0.0, 1200.0) 
# input = t -> bumps(t, [10.0, 10.15], 0.05, 1200.0)
num_h = 6
mistuned_h = 4
amps = [1.0 for k ∈ 1:num_h] .* 1350.0/num_h
ϕ = [0.0 for k ∈ 1:num_h]
mts = [1.0 for k ∈ 1:num_h]
mts[mistuned_h] = 1.035
mistunings = [0.88:0.01:1.12...]
f0 = 200.0
τ0 = 1000.0/f0 	# fundamental period in ms
σ = 5.0

u0 = [-63.641, 0.431]
tspan = (0.0, 1000.0)

ests = []
sols_mt = []
spikes_mt = []
isis_mt = []
for (k, mt) ∈ enumerate(mistunings)
	println("_______________ $k of $(length(mistunings)) _______________")
	mts = [1.0 for k ∈ 1:num_h]
	mts[mistuned_h] = mt
	input = t -> hc(t, f0, amps, ϕ, mts)
	p = [input, σ] 
	prob = SDEProblem(hr, σ_hr, u0, tspan, p)
	sol = solve(prob)
	spikes = findspikes(sol)
	isis = isi(spikes)
	hbw = 0.5 	# half bin width in ms
	isis_1 = filter( x -> τ0-hbw < x < τ0+hbw, isis)
	isis_2 = filter( x -> 2*τ0-hbw < x < 2*τ0+hbw, isis)
	println(length(isis_1) > 0 ? mean(isis_1) : "cycles skipped")
	est = 1000.0/(mean(vcat(isis_1, isis_2 .- 5.0)))
	push!(ests, est)
	push!(sols_mt, sol)
	push!(spikes_mt, spikes)
	push!(isis_mt, isis)
end
plt = plot(mistunings, ests, legend=:none)

end
