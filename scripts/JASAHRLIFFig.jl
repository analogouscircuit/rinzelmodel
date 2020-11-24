module JASAHRLIFFig

using Plots
using Plots.Measures
using DifferentialEquations
using Statistics

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


################################################################################
# Solve problem
################################################################################

# input = t -> ramp(t, 10.0, 10.2, 0.0, 1200.0) 
# input = t -> bumps(t, [10.0, 10.15], 0.05, 1200.0)
num_h = 2
mistuned_h = 2
amps = [1.0 for k ∈ 1:num_h] .* 1200.0/num_h
ϕ = [0.0 for k ∈ 1:num_h]
mts = [1.0 for k ∈ 1:num_h]
mts[mistuned_h] = 0.9 
#mistunings = collect(1.0:0.02:1.075)
mistunings = collect(1.0:0.025:1.08)
f0 = 200.0/1000.0
df = (f0 - mts[mistuned_h]*f0)/2
τ0 = 1000.0/f0 	# fundamental period in ms
bw = 0.07

## Meng-Rinzel 
u0 = [-63.641, 0.431]
tspan = (0.0, 200.0)
isis_mr = []
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
	push!(isis_mr, isis)
end

## LIF
fs = 100e3
dt = 1.0/fs
ts = collect(range(0.0, 0.400, step=1.0/fs))
amps = [1.0 for k ∈ 1:num_h] .* (15.0/num_h) * 1.0e-9
τ_m = 0.006
R_m = 10.0e6
V_thresh = -50.0e-3
V_reset  = -65.0e-3
E_leak   = -65.0e-3
isis_lif = []
for mistuning ∈ mistunings
	global us = zeros(size(ts))
	us[1]    = -65.0e-3
	mts = [1.0 for k ∈ 1:num_h]
	mts[mistuned_h] = mistuning
	global input = t -> hc(t, f0*1000.0, amps, ϕ, mts; scale=:s)
	spikes = []
	for k ∈ 2:length(ts)
		us[k] = us[k-1] + (dt/τ_m) * (E_leak - us[k-1] + R_m * input(ts[k]))
		if us[k] > V_thresh
			us[k] = -65.0e-3
			push!(spikes, ts[k] )
		end
	end
	println("Number of spikes: $(length(spikes))")
	isis = isi(spikes)
	push!(isis_lif, isis .* 1000.0 )
end

################################################################################
# Plot results
################################################################################

bins = range(4.0, 6.0, step=0.07)
#bins = range(0.0, 20.0, step=bw)
titles = ["200 Hz + 400 Hz",
		  "200 Hz + 410 Hz",
		  "200 Hz + 420 Hz",
		  "200 Hz + 430 Hz"]
plts_lif = []
for (k, isis) ∈ enumerate(isis_lif)
	global plts_lif
	# mts = [1.0 for k ∈ 1:num_h]
	# mts[mistuned_h] = mistunings[k]
	# println(mts)
	# input = t -> hc(t, f0, [0.5 for k ∈ 1:num_h] ./ num_h, [π/2 for k ∈ 1:num_h], mts)
	# vals = input.(ts)
	#p = plot(ts, vals, xlim=(4.0, 6.0))
	p = plot([1/f0, 1/f0], [0, 1.1], lc=:red, ls=:dash)
	plot!(p, [1/f0, 1/f0] .* (1.0/mistunings[k]), [0, 1.1], lc=:red, ls=:dash)
	# plot!(p, ([1/f0, 1/f0] .+ ([1/f0, 1/f0] .* (1.0/mistunings[k]))) ./ 2,
	# 	  [0, 1.1], lc=:red, ls=:dashdot)
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
	push!(plts_lif, p)
end
p_lif = plot(plts_lif..., size=(700, 700), margin=5mm)

plts_mr = []
for (k, isis) ∈ enumerate(isis_mr)
	global plts_mr
	# mts = [1.0 for k ∈ 1:num_h]
	# mts[mistuned_h] = mistunings[k]
	# println(mts)
	# input = t -> hc(t, f0, [0.5 for k ∈ 1:num_h] ./ num_h, [π/2 for k ∈ 1:num_h], mts)
	# vals = input.(ts)
	#p = plot(ts, vals, xlim=(4.0, 6.0))
	p = plot([1/f0, 1/f0], [0, 1.1], lc=:red, ls=:dash)
	plot!(p, [1/f0, 1/f0] .* (1.0/mistunings[k]), [0, 1.1], lc=:red, ls=:dash)
	# plot!(p, ([1/f0, 1/f0] .+ ([1/f0, 1/f0] .* (1.0/mistunings[k]))) ./ 2,
	# 	  [0, 1.1], lc=:red, ls=:dashdot)
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
	push!(plts_mr, p)
end
p_mr = plot(plts_mr..., size=(700, 700), margin=5mm)

p = plot(p_lif, p_mr)

end
