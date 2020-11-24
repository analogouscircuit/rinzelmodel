module HRDataAnalysis

using JLD
using Plots
using Statistics
using KernelDensity
using Interpolations
using Distributions
using Measures


include("/home/dahlbom/research/rinzelmodel/util/SpikeUtils.jl")

# -------------------- functions --------------------

function mymean(x, density)
	@assert length(x) == length(density)
	expectation = 0.0
	for k ∈ 1:length(x)-1
		expectation += x[k] *
			   (x[k+1] - x[k]) *
			   (0.5*(density[k+1]+density[k])) 
	end
	expectation
end

function gkern(xs, x0, σ)
	exp.( - ((xs .- x0) .^ 2) ./ (2 * σ ^ 2) )
end

function fastconv(τs::Array{Float64, 1}, points::Array{Float64, 1}, σ::Float64)
	dt = τs[2] - τs[1] 	# assumes equal spacing
	n = length(τs)
	offset = Int(round(3σ/dt))
	out = zeros(length(τs))
	for p ∈ points
		idx = searchsortedfirst(τs, p)
		i_lo = max(1, idx-offset)
		i_hi = min(n, idx+offset)
		out[i_lo:i_hi] += gkern(τs[i_lo:i_hi], p, σ)
	end
	out
end



# -------------------- script --------------------
file = "/home/dahlbom/research/rinzelmodel/data/"*"f0-200.0 dB-80.0 mth-4 bi-true width-3.0 numtrials-5.jld"

data = load(file)
spikes_all = data["spikes_all"]
bw = 0.10

plots = []
ests_f0 = []
ests_all = []
ests_peak_all = []

isis_all = []
τs = collect(0.0:(1/1000e3):0.01) .* 1000.0

for (k,mt) ∈ enumerate(keys(spikes_all))
	global isis_all
	global ests_all
	ests_local = []
	ests_peak_local = []
	isis = Float64[]
	for spikes ∈ spikes_all[mt]
		#@time isis_new = SpikeUtils.isi_ao(1000.0 .* spikes)
		@time isis_new = SpikeUtils.isi(1000.0 .* spikes)
		isis_new = filter( x -> 4.0 < x < 6.0, isis_new)
		density = kde(isis_new, bandwidth=bw)
		expectation = 0.0
		for k ∈ 1:length(density.x)-1
			expectation += density.x[k] *
				   (density.x[k+1] - density.x[k]) *
				   (0.5*(density.density[k+1]+density.density[k])) 
		end
		push!(ests_local, 1000.0/expectation)
		vec = fastconv(τs, isis_new, bw)
		idx = argmax(vec)
		push!(ests_peak_local, 1000.0/τs[idx])
		isis = [isis; isis_new]
	end
	@time density = kde(isis, bandwidth=bw)
	expectation = 0.0
	for k ∈ 1:length(density.x)-1
		expectation += density.x[k] *
			   (density.x[k+1] - density.x[k]) *
			   (0.5*(density.density[k+1]+density.density[k])) 
	end
	push!(ests_f0, (mt, 1000.0/expectation))
	p = plot(density.x, density.density, xlims=(3.5, 6.5), title="$mt")
	push!(plots, p)
	push!(isis_all, (mt, isis))
	push!(ests_all, (mt, ests_local))
	push!(ests_peak_all, (mt, ests_peak_local))
end

isis_all = sort(isis_all, by=first)
isis_all = [x[2] for x ∈ isis_all]
ests_all = sort(ests_all, by=first)
ests_all = [x[2] for x ∈ ests_all]
ests_peak_all = sort(ests_peak_all, by=first)
ests_peak_all = [x[2] for x ∈ ests_peak_all]


ests_μ = [mean(est) for est ∈ ests_all]
ests_std = [std(est) for est ∈ ests_all]
ests_peak_μ = [mean(est) for est ∈ ests_peak_all]
ests_peak_std = [std(est) for est ∈ ests_peak_all]

ests_f0 = sort(ests_f0, by=first)

@time p1 = plot(plots...)
p2 = plot([x[1] for x ∈ ests_f0], [x[2] for x ∈ ests_f0], grid=true)#ylims=(198.0, 202.0))

idx=3
p3 = plot([x[1] for x ∈ ests_f0[3:end-2]], ests_μ[3:end-2],
		  xlabel="Mistuning Factor",
		  ylabel="Pitch (Hz)",
		  xtickfontsize=14,
		  ytickfontsize=14,
		  yerr=ests_std,
		  legend=:none,
		  size=(700,500))

# plts = []
# σ = 0.01
# ests_peak = []
# for isis ∈ isis_all
# 	global plts
# 	vec = fastconv(τs, isis, σ)
# 	p = plot(τs, vec, xlims=(3.5, 6.5))
# 	push!(plts, p)
# 	idx = argmax(vec)
# 	push!(ests_peak, 1000.0/τs[idx])
# end
# p4 = plot(plts[3:end-3]...)
mts = [x[1] for x ∈ ests_f0]
p5 = plot(mts[2:end-0], ests_peak_μ[2:end-0],
		  yerr=ests_peak_std,
		  xlabel="Mistuning Factor",
		  ylabel="Pitch (Hz)",
		  guidefontsize=16,
		  xtickfontsize=14,
		  ytickfontsize=14,
		  legend=:none,
		  framestyle=:box,
		  size=(700,550),
		  margin=10mm)


end
