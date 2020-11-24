module SpikeUtils

using Roots
using DiffEqBase

export findspikes, isi, phasefromspikes

"""
	findspikes(sol; dt=0.01, mint = 0.75)

Find spikes via zero-crossings from below.  dt is the step size for the search.
mint is the minimum time between spikes (important for stochastic scenarios).

NOTE: Time-scale is in ms, not s.
"""
function findspikes(sol::DiffEqBase.AbstractODESolution; dt=0.01, mint=0.75, thresh=0.0)
	t = sol.t[1]
	t0s = [] 
	while t < sol.t[end] - dt
		if sol(t)[1] < thresh && sol(t+dt)[1] >= thresh	
			push!(t0s, (t, t+dt) )
			t += mint
		else	
			t += dt
		end
	end
	# D(f) = t -> ForwardDiff.derivative(f, float(t))
	# [find_zero((t -> sol(t)[1], D(t -> sol(t)[1])), t0, Roots.Newton()) for t0 ∈ t0s]
	[find_zero(t -> sol(t)[1] - thresh, t0) for t0 ∈ t0s]
end

function findspikes(sol::Function, t0, t1; scale=:s, dt=nothing, mint=nothing, thresh=0.0)
	if scale == :s
		dt   == nothing && (dt = 0.01/1000.0)
		mint == nothing && (mint = 0.75/1000.0)
	elseif scale == :ms
		dt   == nothing && (dt = 0.01)
		mint == nothing && (mint = 0.75)
	else
		error("Time scale unrecognized: choose :s or :ms")
	end
	t = t0 
	t0s = [] 
	while t < t1 - dt
		if sol(t) < thresh && sol(t+dt) >= thresh	
			push!(t0s, (t, t+dt) )
			t += mint
		else	
			t += dt
		end
	end
	[find_zero(t -> sol(t) - thresh, t0) for t0 ∈ t0s]
end

"""
	isi(spikes)

Get the interspike times from an array of spike times.
"""
function isi(spikes)
	spikes[2:end] .- spikes[1:end-1]
end


"""
	isi_ao(spikes)

Get the interspike times from an array of spike times (all order)
"""
function isi_ao(spikes)
	isis = Float64[]
	for k ∈ 1:length(spikes)
		for j ∈ k+1:length(spikes)
			push!(isis, spikes[j] - spikes[k])
		end
	end
	isis
end


"""
	phasefromspikes(spikes, f_ref)

Get the firing phase for each spike given a reference frequency.
"""
function phasefromspikes(spikes, f_ref; scale=:s)
	τ0 = 1 / f_ref
	(scale == :ms) && (τ0 *= 1000.0) 
	[(2π/τ0) * (spike % τ0) for spike ∈ spikes]
end

"""
"""
function vectorstrength(spikes, f_ref; scale=:s)	
	ϕs = phasefromspikes(spikes, f_ref; scale=scale)
	vec = exp.( im * ϕs ) |> sum
	return abs(vec)/length(spikes), angle(vec)
end

"""
"""
function vsweights(spikes, edges)
	binvals = [ (edges[i] + edges[i+1])/2.0 for i ∈ 1:length(edges)-1]
	weights = zeros(length(binvals))
	for (k, τ) ∈ enumerate(binvals)
		weights[k] = vectorstrength(spikes, 1/τ)[1]
	end
	return weights
end

function histcounts(vals, bins)
	# This is ugly -- make nice later
	num_edges = length(bins) + 1
	num_bins = num_edges-1
	counts = zeros(num_bins)
	vals = sort(vals)
	i = 1
	for val ∈ vals
		while (i < num_edges) && !( edges[i] <= val < edges[i+1] )
			i += 1
		end
		if i == num_edges
			println("aborting after $i tries!")
			return counts
		else
			counts[i] += 1.0
		end
	end
	return counts
end

# function histcounts(vals, edgess)
# 	# This is ugly -- make nice later
# 	num_edges = length(edges)
# 	num_bins = num_edges-1
# 	counts = zeros(num_bins)
# 	vals = sort(vals)
# 	i = 1
# 	for val ∈ vals
# 		while (i < num_edges) && !( edges[i] <= val < edges[i+1] )
# 			i += 1
# 		end
# 		if i == num_edges
# 			println("aborting after $i tries!")
# 			return counts
# 		else
# 			counts[i] += 1.0
# 		end
# 	end
# 	return counts
# end


"""
	smoothhistfunc(points; kernel = (x, x0) -> exp(-(x-x0)^2/(2 * (0.1)^2)) )
"""
function smoothhistfunc(points; kernel = (x, x0) -> exp(-(x-x0)^2/(2 * (0.1)^2)) )
	x -> sum( kernel.(x, points) )		
end


"""
	smoothhistfunc(points; σ = 0.1)
"""
function smoothhistfunc(points; scale=:s, σ = nothing)
	if scale == :s
		(σ == nothing) && (σ = 0.1/1000.0)
	elseif scale == :ms
		(σ == nothing) && (σ = 0.1)
	else
		error("Scale unrecognized: choose :s or :ms")
	end
	kernel(x, x0) = exp(-(x-x0)^2/(2*σ^2)) 
	x -> sum( kernel.(x, points) )		
end

"""
	kdsvector(points; scale=:s, σ=nothing)
"""
function kdsvector(ts, points; scale=:s, σ=nothing)
	if scale == :s
		(σ == nothing) && (σ = 0.1/1000.0)
	elseif scale == :ms
		(σ == nothing) && (σ = 0.1)
	else
		error("Scale unrecognized: choose :s or :ms")
	end
	kernel(x, x0) = exp.(- (x .- x0).^2 ./ (2*σ^2)) 
	out = zeros(length(ts))
	for t0 ∈ points
		out .+= kernel(ts, t0)
	end
	out
end


end
