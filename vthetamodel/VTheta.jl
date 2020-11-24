module VTheta


# Stucts and Types
struct VThetaParams
	v_r::Real
	v_reset::Real
	Δθ::Real
	a::Real
	b::Real
	c::Real
	τθ::Real
end

# Default parameters
VThetaParams() = VThetaParams(0.1, 0.0, 0.3, 0.08, 4.9, 0.53, 2.0)
params_default = VThetaParams()

# ODE functions
f(v, p) = p.a + exp( p.b * (v - p.c) ) 
dV(v, i, p) = -v + p.v_r + i	
dθ(θ, v, p) = -(θ - f(v, p))/p.τθ

# Simulation Function
function vthetarun(v0, θ0, dur, dt, ifunc::Function, p::VThetaParams)
	ts = collect(0:(Int(round(dur/dt))-1)) .* dt
	n = length(ts)
	vs = Array{Float64, 1}(undef, n); vs[1] = v0
	θs = Array{Float64, 1}(undef, n); θs[1] = θ0
	spikes = Float64[]
	for (k, t) ∈ enumerate(ts[2:end])
		vs[k+1] = vs[k]	+ dt * dV(vs[k], ifunc(t), p)
		θs[k+1] = θs[k] + dt * dθ(θs[k], vs[k], p)
		if vs[k+1] >= θs[k+1]
			vs[k+1] = p.v_reset
			θs[k+1] += p.Δθ
			push!(spikes, t)
		end
	end
	return (spikes, vs, θs)
end

end
