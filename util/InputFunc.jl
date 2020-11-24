module InputFunc


"""
	interpds(t, ds)	

Take sampled data and interpolate at any point. (Just linear interpolation here.)
"""
function interpds(t, ts, ds)	
	if t < ts[1] || t > ts[end-1]
		return 0.0
	end
	i1 = searchsortedfirst(ts, t)
	dt = ts[i1+1] - ts[i1]
	c2 = (t - ts[i1]) / dt 
	(1.0 - c2) * ds[i1] + c2 * ds[i1+1]
end


"""
	ramp(t, t0=0.0, t1=10.0, i0=0.0, i1=1200.0)

Current ramp.  Specify start of ramp, end, and values of current before and after ramp.
"""
function ramp(t, t0=0.0, t1=10.0, i0=0.0, i1=1200.0)
	if t <= t0
		return i0
	elseif t > t1
		return i1
	else
		return i0 + ( (t - t0) / (t1 - t0) ) * (i1 - i0)
	end
end


"""
	bump(t, t0, σ, mag)

Single bump at specified time with specified width and magnitude.
"""
bump(t, t0, σ, mag) = mag * exp( - ((t - t0) ^ 2) / (2 * σ ^2 ) )


"""
	bumps(t, t0s, σ=0.25, mag=1200.0)

Gaussian bumps at specified times, all with specified width and magnitude.
"""
bumps(t, t0s, σ=0.25, mag=1200.0) = sum( [bump(t, t0, σ, mag) for t0 ∈ t0s] )


"""
"""
function sinstim(f0; ϕ=0.0, hwr=true)
	if hwr
		return t -> max(0.0, sin(2π*f0*t + ϕ))
	else
		return t -> sin(2π*f0*t + ϕ)
	end
end

"""
	hc(t, f0, a, ϕ, mts)

Half-wave rectified harmonic complex with specified partial amplitudes, phases and mistunings.

... seems like this should go in InputSamp ...
"""
function hc(t, f0, a, ϕ, mts; scale=:s)
	if scale == :ms
		f0 /= 1000.0
	end
	val = sum( [a[k] * sin(2π*f0*mt*k*t + ϕ[k]) for (k, mt) ∈ enumerate(mts)] )
	max(val, 0.0)
end


end
