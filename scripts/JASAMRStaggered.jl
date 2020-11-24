module JASAMRStaggered

using Plots
using Plots.Measures
using DifferentialEquations
using Statistics
using KernelDensity
using Formatting

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

num_h = 6
mistuned_h = 3
num_grouped = 3
amps = [1.0 for k ∈ 1:num_h] .* 1400.0/num_grouped
ϕ = [0.0 for k ∈ 1:num_h]
mts = [1.0 for k ∈ 1:num_h]
mts[mistuned_h] = 0.9 
#mistunings = collect(1.0:0.02:1.075)
mistunings = collect(1.0:0.025:1.11)
f0 = 200.0/1000.0
df = (f0 - mts[mistuned_h]*f0)/2
τ0 = 1000.0/f0 	# fundamental period in ms
bw = 0.07



## Meng-Rinzel 
u0 = [-63.641, 0.431]
tspan = (0.0, 300.0)
isis_mr = []
for (i, mistuning) ∈ enumerate(mistunings)
	mts = [1.0 for k ∈ 1:num_h]
	mts[mistuned_h] = mistuning
	harmonics = []
	for k ∈ 1:num_h
		h = t -> hc(t, f0*k*mts[k], amps[k], ϕ[k], [1.0])
		push!(harmonics, h)
	end
	isis_all = []
	for k ∈ 1:num_h-2
		input = t -> harmonics[k](t) + harmonics[k+1](t) + harmonics[k+2](t)
		p = [input] 
		prob = ODEProblem(hr, u0, tspan, p)
		global sol = solve(prob, AutoTsit5(Rosenbrock23()), reltol=1e-5, abstol=1e-5)
		spikes = findspikes(sol)
		# ϕs = phasefromspikes(spikes, f0)
		isis = isi(spikes)
		push!(isis_all, isis)
	end
	isis_collected = vcat(isis_all...)
	push!(isis_mr, isis_collected)
end


################################################################################
# Plot results
################################################################################

# bins = range(4.0, 6.0, step=0.07)
bins = range(0.0, 20.0, step=bw)
titles = ["200 Hz + 400 Hz",
		  "200 Hz + 410 Hz",
		  "200 Hz + 420 Hz",
		  "200 Hz + 430 Hz"]

plts_mr = []
ts = collect(range(4.0, 6.0, length=1000))
for (k, isis) ∈ enumerate(isis_mr)
	global plts_mr
	mts = [1.0 for k ∈ 1:num_h]
	mts[mistuned_h] = mistunings[k]
	input = t -> hc(t, f0, [1.0 for k ∈ 1:num_h] ./ num_h, [π/2 for k ∈ 1:num_h], mts)
	vals = input.(ts) 	# summed inputs for plotting
	isis_filt = filter(x -> 4.25 <= x < 5.75, isis)
	max_lin = ts[argmax(vals)]
	mean_spikes = mean(isis_filt)
	p = plot(ts, vals, xlim=(4.0, 6.0), lw=2.0, la=0.9)
	plot!(p, [mean_spikes, mean_spikes], [0, 1.1], lc=:red, ls=:dash, lw=1.75)
	plot!(p, [max_lin, max_lin], [0, 1.1], lc=:black, ls=:dot, lw=1.75)
	#legend = length(isis_mr) == k ? :outertopright : :none
	legend = :none 
	histogram!(p, isis,
				  bins=bins,
				  fillcolor=:gray,
				  fillalpha=1.00,
				  legend=legend,
				  normalize=:probability,
				  xtickfontsize=14,
				  xrotation=45,
				  ytickfontsize=14,
				  title=format("Mistuning: {:.1f}%",(mistunings[k]-1.0)*100.0),
				  # title=titles[k],
				  ylims=(0.0, 1.1),
				  xlabel="time (ms)",
				  ylabel="probability",
				  guidefontsize=16,
				  framestyle=:box)
	println("--------------------$(mistunings[k])--------------------")
	println("Linear Superposition:\t$max_lin")
	println("Mean of spikes:\t\t$mean_spikes")
	println()
	push!(plts_mr, p)
end

# hack for legend
for (k, isis) ∈ enumerate([isis_mr[end]])
	global plts_mr
	mts = [1.0 for k ∈ 1:num_h]
	mts[mistuned_h] = mistunings[k]
	input = t -> hc(t, f0, [1.0 for k ∈ 1:num_h] ./ num_h, [π/2 for k ∈ 1:num_h], mts)
	vals = input.(ts) 	# summed inputs for plotting
	isis_filt = filter(x -> 4.25 <= x < 5.75, isis)
	max_lin = ts[argmax(vals)]
	mean_spikes = mean(isis_filt)
	p = plot(ts, vals, xlim=(4.0, 6.0), lw=2.5, label="sum of harmonics")
	plot!(p, [max_lin, max_lin], [0, 1.1], lc=:black, ls=:dot, lw=1.75, label="max of sum")
	plot!(p, [mean_spikes, mean_spikes], [0, 1.1], lc=:red, ls=:dash, lw=1.75, label="mean of spikes")
	legend=:topleft
	histogram!(p, isis,
			      label="ISIs",
				  bins=bins,
				  fillcolor=:gray,
				  fillalpha=0.75,
				  legend=legend,
				  legendfontsize=12,
				  normalize=:probability,
				  # title=titles[k],
				  ylims=(2.0, 3.0),
				  guidefontsize=16,
				  grid=false,
				  framestyle=:none)
	println("--------------------$(mistunings[k])--------------------")
	println("Linear Superposition:\t$max_lin")
	println("Mean of spikes:\t\t$mean_spikes")
	println()
	push!(plts_mr, p)
end

p_mr = plot(plts_mr..., size=(900, 650), margin=5mm)

p = plot(p_mr)

end
