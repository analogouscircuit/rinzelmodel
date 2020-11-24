module HRTuningCurve

using Plots
using DifferentialEquations
using Roots
using ForwardDiff
using Statistics
using Formatting: format

include("/home/dahlbom/research/rinzelmodel/util/SpikeUtils.jl")
include("/home/dahlbom/research/rinzelmodel/util/InputFunc.jl")
# include("/home/dahlbom/research/rinzelmodel/sdmodels/HRDirectInput.jl")
using .SpikeUtils: findspikes, isi, phasefromspikes
using .InputFunc: ramp, bump, bumps, hc
# using .HRDirectInput: hr

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
# Solve problem
################################################################################
#maxamps = 900.0:50.0:1600.0
maxamps = [1350.0]
#maxamps = 1345.0:1.0:1355.0
num_h = 2
amps = [1.0 for k ∈ 1:num_h] .* 1200.0/num_h
ϕ = [0.0 for k ∈ 1:num_h]
mtp = 2
mistunings = [0.88:0.005:1.12...]
f0 = 200.0
factor = 3.0
τ0 = 1000.0/f0 	# fundamental period in ms
# σ = 0.1 	# noise for input -- DON'T DO THIS -- STEP SIZE ISSUES -- USE PROPER SDE

# u0 = [-63.641, 0.431]
u0 = [-60.0, 0.380]
tspan = (0.0, 500.0)

plots = []
hists = []
sols_a = []
spikes_a = []
isis_a = []

ests = []
strengths = []
for (k,mt) ∈ enumerate(mistunings)
	println(" $k of $(length(mistunings))")
	mts = [1.0 for k ∈ 1:num_h]
	mts[mtp] = mt
	input = t -> hc(t, f0, amps, ϕ, mts)
	p = [input, factor]
	prob = ODEProblem(hr, u0, tspan, p)
	sol = solve(prob, AutoTsit5(Rosenbrock23()), reltol=1e-4, abstol=1e-4)
	spikes = findspikes(sol)
	isis = isi(spikes)[2:end]
	hbw = 0.5 	# half bin width in ms
	isis_1 = filter( x -> τ0-hbw < x < τ0+hbw, isis)
	est = length(isis_1) > 1 ? est = 1000.0/(mean(isis_1)) : 200.0
	# strength = SpikeUtils.vectorstrength(spikes, est; scale=:ms)[1]
	bin_vals = SpikeUtils.histcounts(isis_1, range(0.0, stop=10.0, step=0.1))
	bin_vals ./= length(isis_1)
	strength = maximum(bin_vals)
	hist = histogram(isis_1; 
					 bins=range(0.0, 10.0, step=0.1),
					 normalize=:probability,
					 legend=:none,
					 ylims=(0.0, 1.0))
	annotate!(hist, 6.5, 0.8, text("n: $(length(isis_1))", 8, :left))
	annotate!(hist, 6.5, 0.65, text("mt: $mt", 8, :left))
	annotate!(hist, 6.5, 0.5, text(format("f0: {:.1f}", est), 8, :left))
	annotate!(hist, 6.5, 0.35, text(format("s: {:.1f}", strength), 8, :left))
	push!(hists, hist)
	push!(ests, est)
	push!(strengths, strength)
	push!(isis_a, isis_1)
	push!(sols_a, sol)
	push!(spikes_a, spikes)
end

p1 = plot(mistunings, ests, legend=:none, title="$factor", xlabel="mistuning", ylabel="f0 estimate")
p2 = plot(mistunings, strengths, xlabel="mistuning", ylabel="strength")
p1 = plot(p1, p2, layout=(2,1))
p3 = plot(hists...)
l = @layout [ a{0.3w} b ]
p = plot(p1, p3, layout=l, size=(1800, 800))

end
