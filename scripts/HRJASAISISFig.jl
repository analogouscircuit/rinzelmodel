module HRJASAISISFig

using Plots
using Plots.Measures
using DifferentialEquations
using Statistics
using KernelDensity

include("/home/dahlbom/research/rinzelmodel/util/SpikeUtils.jl")
include("/home/dahlbom/research/rinzelmodel/util/InputFunc.jl")
include("/home/dahlbom/research/rinzelmodel/sdmodels/HRDirectInput.jl")
using .SpikeUtils: findspikes, isi, isi_ao, phasefromspikes
using .InputFunc: ramp, bump, bumps, hc
using .HRDirectInput: hr

################################################################################
# Functions
################################################################################
gaussbump(t, t0; σ=0.5) = exp( - (t - t0) ^ 2 / ( 2 * σ^2 ) )

function histconvolve(t, points;
					  kernel= (t, t0) -> gaussbump(t, t0; σ=0.5))
	sum( [kernel(t, t0) for t0 ∈ points] )
end

function kdemean(density)
	expectation = 0.0
	for k ∈ 1:length(density.x)-1
		expectation += density.x[k] *
			   (density.x[k+1] - density.x[k]) *
			   (0.5*(density.density[k+1]+density.density[k])) 
	end
	return expectation
end

################################################################################
# Solve problem
################################################################################

# input = t -> ramp(t, 10.0, 10.2, 0.0, 1200.0) 
# input = t -> bumps(t, [10.0, 10.15], 0.05, 1200.0)
num_h = 3
mistuned_h = 1
amps = [1.0 for k ∈ 1:num_h] .* 1400.0/num_h
ϕ = [0.0 for k ∈ 1:num_h]
mts = [1.0 for k ∈ 1:num_h]
mts[mistuned_h] = 0.9 
#mistunings = collect(1.0:0.02:1.075)
mistunings = collect(1.0:0.02:1.12)
f0 = 200.0/1000.0
df = (f0 - mts[mistuned_h]*f0)/2
τ0 = 1000.0/f0 	# fundamental period in ms
bw = 0.07
bw_kde = 0.01

u0 = [-63.641, 0.431]
tspan = (0.0, 200.0)

isis_all = []

for mistuning ∈ mistunings
	mts = [1.0 for k ∈ 1:num_h]
	mts[mistuned_h] = mistuning
	global input = t -> hc(t, f0, amps, ϕ, mts)
	p = [input] 
	prob = ODEProblem(hr, u0, tspan, p)
	global sol = solve(prob, AutoTsit5(Rosenbrock23()), reltol=1e-5, abstol=1e-5)
	spikes = findspikes(sol)
	# ϕs = phasefromspikes(spikes, f0)
	isis = isi(spikes)
	push!(isis_all, isis)
	hbw = 0.5 	# half bin width in ms
	isis1 = filter( x -> τ0-hbw < x < τ0+hbw, isis)
	println(length(isis1) > 1 ? 1000.0/mean(isis1[2:end]) : "no spikes or skipped cycles")
end

################################################################################
# Plot results
################################################################################

plts = []
#bins = range(4.0, 6.0, step=0.07)
bins = range(4.0, 6.0, step=bw)
titles = ["200 Hz + 400 Hz",
		  "200 Hz + 410 Hz",
		  "200 Hz + 420 Hz",
		  "200 Hz + 430 Hz"]
ts = collect(range(4.0, 6.0, length=1000))
for (k, isis) ∈ enumerate(isis_all)
	mts = [1.0 for k ∈ 1:num_h]
	mts[mistuned_h] = mistunings[k]
	# println(mts)
	input = t -> hc(t, f0, [1.0 for k ∈ 1:num_h] ./ num_h, [π/2 for k ∈ 1:num_h], mts)
	vals = input.(ts)
	p = plot(ts, vals, xlim=(4.0, 6.0))
	#p = plot([1/f0, 1/f0], [0, 1.1], lc=:red, ls=:dash)
	plot!(p, [1/f0, 1/f0], [0, 1.1], lc=:red, ls=:dash)
	plot!(p, [1/f0, 1/f0] .* (1.0/mistunings[k]), [0, 1.1], lc=:red, ls=:dash)
	# plot!(p, ([1/f0, 1/f0] .+ ([1/f0, 1/f0] .* (1.0/mistunings[k]))) ./ 2,
	#  	  [0, 1.1], lc=:red, ls=:dashdot)
	histogram!(p, isis,
				  bins=bins,
				  fillcolor=:gray,
				  legend=:none,
				  normalize=:probability,
				  xtickfontsize=14,
				  xrotation=45,
				  ytickfontsize=14,
				  title="$(mistunings[k])",
				  # title=titles[k],
				  ylims=(0.0, 1.1),
				  xlabel="time (ms)",
				  ylabel="probability",
				  guidefontsize=16,
				  framestyle=:box)
	isis_filt = filter(x -> 4.25 <= x < 5.75, isis)
	density = kde(isis_filt, bandwidth=0.04)
	plot!(p, density.x, density.density ./ maximum(density.density), lc=:red)
	println("--------------------$(mistunings[k])--------------------")
	max_lin = ts[argmax(vals)]
	max_den = density.x[argmax(density.density)]
	mean_den = kdemean(density)
	println("Linear Superposition:\t$max_lin")
	println("Maximum of density:\t$max_den")
	println("Mean of density:\t$mean_den")
	println()
	push!(plts, p)
end
p = plot(plts..., size=(700, 700), margin=5mm)


end
