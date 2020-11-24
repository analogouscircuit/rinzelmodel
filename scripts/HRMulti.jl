module HRMulti
# Huang-Rinzel-Tuning-Curve-Multiple-Slope-Detector-Units

using Plots
using DifferentialEquations
using Statistics: mean
using Formatting: format
using DataFrames

include("/home/dahlbom/research/rinzelmodel/util/SpikeUtils.jl")
include("/home/dahlbom/research/rinzelmodel/util/InputFunc.jl")
include("/home/dahlbom/research/AdaptiveSieves/util/PeakPick.jl")
# include("/home/dahlbom/research/rinzelmodel/sdmodels/HRDirectInput.jl")
# using .HRDirectInput: hr
using .SpikeUtils: findspikes, isi, isi_ao, phasefromspikes, smoothhistfunc
using .InputFunc: ramp, bump, bumps, hc
using .PeakPick: findpeaknear, peakinterp

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
end


################################################################################
# Solve problem
################################################################################
num_h = 8
num_sd_units = num_h - 1 
amp = 600.0 
mtp = 4 # mistuned partial
mt = 1.1
# mistunings = [0.88:0.005:1.12...]
f0 = 100.0
τ0 = 1000.0/f0 # in ms
factor = 3.0
bins = range(4.0, 16.0, step=0.1)
bins_all = range(0.0, 50.0, step=0.1)
tol = 0.5

# u0 = [-63.641, 0.431]
u0 = [-60.0, 0.380]
tspan = (0.0, 500.0)


isis_sd = []
spikes_sd = []
τs = collect(range(0.0, 50.0, step=0.01))

for sdi ∈ 1:num_sd_units 
	mt1 = sdi == mtp ? mt : 1.0
	mt2 = (sdi+1) == mtp ? mt : 1.0
	input1 = t -> hc(t, f0*sdi*mt1, [600.0], [0.0], [1.0]; scale=:ms)
	input2 = t -> hc(t, f0*(sdi+1)*mt2, [600.0], [0.0], [1.0]; scale=:ms)
	input(t) = input1(t) + input2(t) 
	p = [input, factor]
	prob = ODEProblem(hr, u0, tspan, p)
	sol = solve(prob, AutoTsit5(Rosenbrock23()), reltol=1e-4, abstol=1e-4)
	spikes = findspikes(sol)
	#isis = isi(spikes)
	isis = isi_ao(spikes)
	push!(spikes_sd, spikes)
	push!(isis_sd, isis)
end


hists = []
for k ∈ 1:num_sd_units
	hist = histogram(isis_sd[k], bins=bins_all, legend=:none)
	f = smoothhistfunc(isis_sd[k]; σ=0.05)
	plot!(hist, τs, f)
	for q ∈ 1:4
		isis_filt = filter( x -> q*10.0-tol < x < q*10.0+tol, isis_sd[k])
		est_avg = length(isis_filt) > 0 ? mean(isis_filt) : 0.0
		f_vals = f.(τs)
		i = findpeaknear(τs, f_vals, q*τ0; tol=0.2)
		est_pp = peakinterp(τs[i-1:i+1], f_vals[i-1:i+1])
		annotate!(hist, 10.0*q+0.5, 10.0, text(format("Avg: {:.2f}", q*1000.0/est_avg), 8, :left))
		annotate!(hist, 10.0*q+0.5, 12.0, text(format("PP:  {:.2f}", q*1000.0/est_pp), 8, :left))
	end
	push!(hists, hist)
end


isis_all = vcat(isis_sd...)
isis_all_filt = filter( x -> 10.0-tol < x < 10.0+tol, isis_all)
est = mean(isis_all_filt)
est = 1000.0/est
p_each = plot(hists..., size=(1500, 900))
p_hist_chan_sum = histogram(isis_all, bins=bins_all, legend=:none, size=(1700, 900))
annotate!(p_hist_chan_sum, 13.0, 20.0, text(format("{:.2f}", est), 8, :left))
isis_spikes_all = isi_ao(vcat(spikes_sd...))
isis_all_filt = filter( x -> 10.0-tol < x < 10.0+tol, isis_spikes_all)
est = mean(isis_all_filt)
est = 1000.0/est
p_hist_all = histogram(isis_spikes_all, bins=bins_all, legend=:none, size=(1700, 900))
annotate!(p_hist_all, 30.0, 150.0, text(format("{:.2f}", est), 14, :left))
p_all = plot(p_hist_chan_sum, p_hist_all, layout=(1,2), size=(1800, 900), legend=:none)


end
