module HRInPat

using Plots
using Plots.Measures
using Formatting
using DifferentialEquations
using Statistics: mean

include("/home/dahlbom/research/rinzelmodel/util/SpikeUtils.jl")
include("/home/dahlbom/research/rinzelmodel/util/InputFunc.jl")
include("/home/dahlbom/research/rinzelmodel/sdmodels/HRDirectInput.jl")
include("/home/dahlbom/research/AdaptiveSieves/util/PeakPick.jl")
using .SpikeUtils: isi, isi_ao, vectorstrength, findspikes, smoothhistfunc
using .InputFunc: sinstim
using .HRDirectInput: hr
using .PeakPick: findpeaknear, peakinterp


# -------------------- functions --------------------
w_gauss_full(cf_sd, cf_an; σ=2) = exp( - abs( log(cf_an) - log(cf_sd) )^2 / σ^2 )

function w_gauss_half(cf_sd, cf_an; σ=2)
	if cf_an >= cf_sd
		return w_gauss_full(cf_sd, cf_an; σ=σ) 
	else
		return 0.0
	end
end


# -------------------- Parameters --------------------
f0 = 200.0
num_h = 7
mistuned_h = 4
mt = 1.04
amp = 1200.0 
σ=1.0
dur = 250.0
fs = 100e3
dt = 1/fs * 1000.0
bin_step = 0.025
bin_lo = 4.0
bin_hi = 6.0
smooth_step = 0.005
bins_hist = range(bin_lo, bin_hi, step=bin_step)
bins_sm = range(bin_lo, bin_hi, step=smooth_step)
freqs_in = [ f0*k for k ∈ 1:num_h ]
freqs_in[mistuned_h] *= mt
cfs_sd = exp.( collect(range(log(100), log(maximum(freqs_in)), length=6)) )


# -------------------- Generate Inputs --------------------
channels = []
for freq ∈ freqs_in
	push!(channels, sinstim(freq/1000.0; ϕ=0.0, hwr=true))
end

inputs = []
plts_in = []
t = collect(range(0, dur-dt, step=dt))
w(cf_sd, cf_an) = w_gauss_full(cf_sd, cf_an; σ=σ)
for cf_sd ∈ cfs_sd
	input(t) = sum([amp * w(cf_sd, cf_an) * channels[k](t) for (k, cf_an) ∈ enumerate(freqs_in)])
	push!(inputs, input)	
	push!(plts_in, plot(t, input.(t), legend=:none, ylabel=format("{:.2f}", cf_sd)))
end
p_inputs = plot(reverse(plts_in)..., layout=(length(plts_in), 1), legend=:none)


# -------------------- SD Unit Solve --------------------
isis_chan = []
sols_chan = []
plts_sols = []
plts_hists = []
u0 = [-63.641, 0.431]
tspan = (0.0, dur)
for input ∈ inputs
	p = [input] 
	prob = ODEProblem(hr, u0, tspan, p)
	sol = solve(prob, AutoTsit5(Rosenbrock23()), reltol=1e-5, abstol=1e-5)
	spikes = findspikes(sol)
	isis = isi(spikes)
	smoothed = smoothhistfunc(isis; scale=:ms).(bins_sm)
	i = findpeaknear(bins_sm, smoothed, bins_sm[argmax(smoothed)])
	τ_peak = peakinterp(bins_sm[i-1:i+1], smoothed[i-1:i+1])
	push!(isis_chan, isis)
	push!(plts_sols, plot(sol, legend=:none, xlabel=""))
	h_temp = histogram(isis, bins=bins_hist, legend=:none)
	annotate!(h_temp, 5.5, 5.0, text("$(1000.0/τ_peak)", 8))
	push!(plts_hists, h_temp)
end
p_hists = plot(reverse(plts_hists)..., legend=:none, layout=(length(plts_hists), 1))
p_sols = plot(reverse(plts_sols)..., layout=(length(plts_sols), 1), legend=:none)
p_stacked = plot(p_inputs, p_sols, p_hists, layout=(1,3))
isis_all = vcat(isis_chan...)
p_hist_all = histogram(isis_all, bins=bins_hist, legend=:none, 
					   title="width: $σ, amp: $amp, f0: $f0, num: $mistuned_h, %: $mt",
					   titlefontsize=9)
smoothed = smoothhistfunc(isis_all; scale=:ms).(bins_sm)
i = findpeaknear(bins_sm, smoothed, bins_sm[argmax(smoothed)])
τ_peak = peakinterp(bins_sm[i-1:i+1], smoothed[i-1:i+1])
annotate!(p_hist_all, 5.5, 20.0, text("$(1000.0/τ_peak)", 8))
l = @layout [a{0.8w} b]
p = plot(p_stacked, p_hist_all,
		 layout=l, 
		 legend=:none,
		 size=(1800, 600))



end
