module HRIPD

using Plots
using Plots.PlotMeasures
using DifferentialEquations
using Statistics

include("/home/dahlbom/research/rinzelmodel/util/SpikeUtils.jl")
include("/home/dahlbom/research/rinzelmodel/util/InputFunc.jl")
#include("/home/dahlbom/research/rinzelmodel/sdmodels/HRDirectInput.jl")
using .SpikeUtils: findspikes, isi, isi_ao, phasefromspikes
using .InputFunc: ramp, bump, bumps, hc
#using .HRDirectInput: hr

################################################################################
# Experiment Parameters
################################################################################
f0 = 200.0
factors = 1.0:9.0 
num_cycles = 15 
dur = 1000.0 * (num_cycles/f0)
amp = 600.0
ϕs = collect(range(-π, π, length=50)) 
ϕs_normed = ϕs ./ 2π

################################################################################
# Differential equation (copied here for tweaking without messing master file)
################################################################################

function hr(du, u, p, t)
	# Parameters
	## Constants
	### Maximal Conductances, Reversal Potentials and Membrane Capacitance
	g_Na  = 1000.0;		e_Na = 55.0
	g_h   = 20.0; 		e_h = -43.0
	g_L   = 2.0; 		e_L = -65.0
	g_KLT = 200.0; 		e_K = -70.0
	g_KHT = 150.0;

	c_m = 12.0

	### Frozen gating parameters
	w_0 = 0.511
	h_0 = 0.445
	z_0 = 0.662
	n_0 = 0.0077
	p_0 = 0.0011
	r_0 = 0.147

	### Other gating parameters
	a = 0.9
	b = ( a - w_0 ) / h_0

	## Functions
	### Time Constants
	τ_w(v) = 100.0 / ( 6.0 * exp( (v + 60.0) / 6.0 ) + 16.0 * exp( - (v + 60.0) / 45.0 ) ) + 1.5
	τ_h(v) = 100.0 / ( 7.0 * exp( (v + 60.0) / 11.0 ) + 10.0 * exp( -(v + 60.0) / 25.0 ) ) + 0.6
	τ_U(v) = min( τ_w(v), τ_h(v) )

	### Steady-state gating functions
	m_inf(v) = 1.0 / ( 1.0 + exp( -(v + 38.0) / 7.0 ) )
	h_inf(v) = 1.0 / ( 1.0 + exp( (v + 65.0) / 6.0 ) )
	w_inf(v) = ( 1.0 + exp( -(v + 48.0) / 6.0 ) ) ^ (-0.25)
	u_inf(v) = b * ( h_inf(v) + b * ( a - w_inf(v) ) ) / ( a * ( 1.0 + b ^ 2 ) )

	### Input function is first given parameter
	i = p[1]
	factor = p[2]

	du[1] = - (2.0 / c_m) * ( g_Na * (a/b) * (m_inf(u[1]) ^ 3.0) * u[2] * (u[1] - e_Na) +
							   g_KLT * ((a - a*u[2]) ^ 4) * z_0 * (u[1] - e_K) +
							   g_KHT * ( 0.85 * n_0^2 + 0.15 * p_0 ) * (u[1] - e_K) +
							   g_h * r_0 * (u[1] - e_h) +
							   g_L * (u[1] - e_L)
							   ) + i(t)/c_m
	du[2] = factor * ( u_inf(u[1]) - u[2] ) / τ_U(u[1])
	# du[1] = - ((2.0 / c_m) * ( g_Na * (a/b) * (m_inf(u[1]) ^ 3.0) * u[2] * (u[1] - e_Na) +
	# 						   g_KLT * ((a - a*u[2]) ^ 4) * z_0 * (u[1] - e_K) +
	# 						   g_KHT * ( 0.85 * n_0^2 + 0.15 * p_0 ) * (u[1] - e_K) +
	# 						   g_h * r_0 * (u[1] - e_h) +
	# 						   g_L * (u[1] - e_L)
	# 						   ) + i(t)/c_m) * 1000.0
	# du[2] = 1000.0 * ( u_inf(u[1]) - u[2] ) / τ_U(u[1])
end

################################################################################
# Solve problem and generate plots
################################################################################
inputs = []
for ϕ ∈ ϕs
	input = t -> hc(t, f0, [amp], [0.0], [1.0]) + hc(t, f0, [amp], [ϕ], [1.0])
	push!(inputs, input)
end

# u0 = [-63.641, 0.431]
u0 = [-60.0, 0.380]
tspan = (0.0, dur)
# sols = []
plots = []
for factor ∈ factors
	spikes_per_cycle = []
	for (k, input) ∈ enumerate(inputs)
		println("Starting $k of $(length(inputs))...")
		p = [] 
		push!(p, input)
		push!(p, factor)
		prob = ODEProblem(hr, u0, tspan, p)
		sol = solve(prob, AutoTsit5(Rosenbrock23()), reltol=1e-5, abstol=1e-5)
		spikes = findspikes(sol)
		push!(spikes_per_cycle, length(spikes)/num_cycles)
		# push!(sols, sol)
	end
	plt = plot(ϕs_normed, spikes_per_cycle; 
			 xticks=-0.5:0.1:0.5,
			 lw=1.5,
			 legend=:none, 
			 title="$factor ($(f0)Hz, $(amp)pA)",
			 xlabel="Normalized Phase Difference",
			 ylabel="spikes per cycle", 
			 titlefontsize=8)
	push!(plots, plt)
end

################################################################################
# Create final plot 
################################################################################
p = plot(plots...;
		 margin=10mm,
		 titlefontsize=14,
		 size=(1200,900),
		 layout=(3,3))

# v_trace = [v[1] for v ∈ sol.u]
# t = collect(range(sol.t[input_num], sol.t[end], step=0.1))
# p1 = plot(sol.t, v_trace)
# for spike ∈ spikes
# 	plot!(p1, [spike, spike], [-60.0, 30.0], color=:black, ls=:dash)
# end
# p2 = plot(t, inputs[input_num])
# p = plot(p1, p2, layout=(2,1), size=(600,900))

end
