module HRBiHistTC
# Huang-Rinzel-Tuning-Curve-Multiple-Slope-Detector-Units

using Plots
using Plots.Measures
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
mistunings = [0.89:0.02:1.11...]
f0 = 200.0
factor = 3.0
bins = range(4.0, 16.0, step=0.1)
bins_all = range(0.0, 50.0, step=0.1)
τ0 = 1000.0/f0 # in ms
tol = 0.5

u0 = [-60.0, 0.380]
tspan = (0.0, 500.0)
τs = collect(range(0.0, 50.0, step=0.01))

# Solve the ODE and get spike times and ISIs
isis_lo = []
isis_hi = []
spikes_lo = []
spikes_hi = []
for (k, mt) ∈ enumerate(mistunings)
	println("$k of $(length(mistunings))")
	# Lower harmonic mistuned
	in_lo = t -> hc(t, f0*mt, [amp], [0.0], [1.0]; scale=:ms) +
				 hc(t, 2*f0, [amp], [0.0], [1.0]; scale=:ms)
	p_lo = [in_lo, factor]
	prob_lo = ODEProblem(hr, u0, tspan, p_lo)
	sol_lo = solve(prob_lo, AutoTsit5(Rosenbrock23()), reltol=1e-4, abstol=1e-4)
	spks_lo = findspikes(sol_lo)
	isi_lo = isi_ao(spks_lo)
	push!(isis_lo, isi_lo)
	push!(spikes_lo, spks_lo)
	# Upper harmonic mistuned
	in_hi = t -> hc(t, f0, [amp], [0.0], [1.0]; scale=:ms) + 
				 hc(t, 2*f0*mt, [amp], [0.0], [1.0]; scale=:ms)
	p_hi = [in_hi, factor]
	prob_hi = ODEProblem(hr, u0, tspan, p_hi)
	sol_hi = solve(prob_hi, AutoTsit5(Rosenbrock23()), reltol=1e-4, abstol=1e-4)
	spks_hi = findspikes(sol_hi)
	isi_hi = isi_ao(spks_hi)
	push!(isis_hi, isi_hi)
	push!(spikes_hi, spks_hi)

end


# Generate histograms and peak estimates
hists = []
y1 = 0.5
y2 = 0.9
ests_lo_pp = []
ests_lo_avg = []
ests_hi_pp = []
ests_hi_avg = []
for k ∈ 1:length(mistunings)
	# lo mistuned
	hist_lo	= histogram(isis_lo[k], bins=bins_all, legend=:none, normalize=:pdf) 
	f_lo = smoothhistfunc(isis_lo[k]; σ=0.05)
	f_lo_vals = f_lo.(τs)
	f_lo_vals ./= maximum(f_lo_vals)
	plot!(hist_lo, τs, f_lo_vals)
	# hi mistuned
	hist_hi = histogram(isis_hi[k], bins=bins_all, legend=:none, normalize=:pdf) 
	f_hi = smoothhistfunc(isis_hi[k]; σ=0.05)
	f_hi_vals = f_hi.(τs)
	f_hi_vals ./= maximum(f_hi_vals)
	plot!(hist_hi, τs, f_hi_vals)
	for q ∈ 1:4
		# lo mistuned
		isis_lo_filt = filter( x -> q*τ0-tol < x < q*τ0+tol, isis_lo[k])
		est_lo_avg = length(isis_lo_filt) > 0 ? mean(isis_lo_filt) : 0.0
		i = findpeaknear(τs, f_lo_vals, q*τ0; tol=0.2)
		if i == 0
			est_lo_pp = 0.0
		else
			est_lo_pp = peakinterp(τs[i-1:i+1], f_lo_vals[i-1:i+1])
		end
		annotate!(hist_lo, τ0*q+0.5, y1,
				  text(format("Avg: {:.2f}", q*1000.0/est_lo_avg), 8, :left))
		annotate!(hist_lo, τ0*q+0.5, y2,
				  text(format("PP:  {:.2f}", q*1000.0/est_lo_pp), 8, :left))
		# hi mistuned
		isis_hi_filt = filter( x -> q*τ0-tol < x < q*τ0+tol, isis_hi[k])
		est_hi_avg = length(isis_hi_filt) > 0 ? mean(isis_hi_filt) : 0.0
		i = findpeaknear(τs, f_hi_vals, q*τ0; tol=0.2)
		if i == 0
			est_hi_pp = 0.0
		else
			est_hi_pp = peakinterp(τs[i-1:i+1], f_hi_vals[i-1:i+1])
		end
		annotate!(hist_hi, τ0*q+0.5, y1,
				  text(format("Avg: {:.2f}", q*1000.0/est_hi_avg), 8, :left))
		annotate!(hist_hi, τ0*q+0.5, y2,
				  text(format("PP:  {:.2f}", q*1000.0/est_hi_pp), 8, :left))
		if q == 1
			push!(ests_lo_pp, est_lo_pp)
			push!(ests_lo_avg, est_lo_avg)
			push!(ests_hi_pp, est_hi_pp)
			push!(ests_hi_avg, est_hi_avg)
		end
	end
	hist = plot(hist_lo, hist_hi, layout=(2,1), legend=:none, margin=5mm, title="$(mistunings[k])")
	push!(hists, hist)
end

p_hist = plot(hists..., size=(1800, 900))

p_lo = plot(mistunings, 1000.0 ./ ests_lo_pp, label="peak pick", title="Lower mistuned")
plot!(p_lo, mistunings, 1000.0 ./ ests_lo_avg, ls=:dash, lc=:black, label="average")
p_hi = plot(mistunings, 1000.0 ./ ests_hi_pp, label="peak pick", title="Upper mistuned")
plot!(p_hi, mistunings, 1000.0 ./ ests_hi_avg, ls=:dash, lc=:black, label="average")
p_curves = plot(p_lo, p_hi, layout=(1,2), size=(1200,600))


end
