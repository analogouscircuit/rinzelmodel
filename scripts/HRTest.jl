module HRTest

using Plots
using Plots.Measures
using DifferentialEquations
using Statistics

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


################################################################################
# Solve problem
################################################################################

# input = t -> ramp(t, 10.0, 10.2, 0.0, 1200.0) 
# input = t -> bumps(t, [10.0, 10.15], 0.05, 1200.0)
num_h = 2
mistuned_h = 2
amps = [1.0 for k ∈ 1:num_h] .* 1200.0/num_h
ϕ = [0.0 for k ∈ 1:num_h]
mts = [1.0 for k ∈ 1:num_h]
mts[mistuned_h] = 0.9 
f0 = 200.0/1000.0
df = (f0 - mts[mistuned_h]*f0)/2
τ0 = 1000.0/f0 	# fundamental period in ms
input = t -> hc(t, f0, amps, ϕ, mts)

u0 = [-63.641, 0.431]
tspan = (0.0, 200.0)
p = [input] 
prob = ODEProblem(hr, u0, tspan, p)
sol = solve(prob, AutoTsit5(Rosenbrock23()), reltol=1e-5, abstol=1e-5)

spikes = findspikes(sol)
ϕs = phasefromspikes(spikes, f0)
isis = isi(spikes)
hbw = 0.5 	# half bin width in ms
isis1 = filter( x -> τ0-hbw < x < τ0+hbw, isis)
println(length(isis1) > 1 ? 1000.0/mean(isis1[2:end]) : "no spikes or skipped cycles")

################################################################################
# Plot results
################################################################################

p_in = plot(sol.t, input, xlims=(sol.t[1], sol.t[end]), color=:red, legend=:none, la=0.4, title="Input")
p_out = plot(sol, vars=(1), ylims=(-75.0, 50.0), legend=:none, title="SD Output")
l1 = @layout[ a{0.7h}; b ]
p1 = plot(p_out, p_in, layout=l1, margin=5mm)

p2 = histogram(isis, bins=(range(2.5, 7.5, step=0.05)),
	  		   legend=:none, title="ISI Histogram")

l2 = @layout[ a; b]
p = plot(p1, p2, layout=l2, size=(800,800), margin=5mm)


# plt1 = plot(sol, vars=(1), ylims=(-75.0, 50.0), legend=:none, title="SD Output")
# for spike ∈ spikes
# 	plot!(plt1, [spike, spike], [-75.0, 50.0], ls=:dash, color=:black, label=:none)
# end
# plt2 = plot(sol.t, input, xlims=(sol.t[1], sol.t[end]), color=:red, legend=:none, la=0.4)
# plot!(plt2, sol.t, t -> abs(1000.0 * sin(2π*df*4/1000.0 * t)), la=0.75)
# for (k, spike) ∈ enumerate(spikes)
# 	plot!(plt2, [spike, spike], 1000.0 .* [0.0, ϕs[k]], color=:blue)
# end
# l = @layout [a{0.9h}; b]
# p1 = plot(plt1, plt2, layout=l)
# 
# hist_rad = histogram(ϕs, proj=:polar, nbins=360, legend=:none, title="Phase Histogram")
# #hist_isi = histogram(isis, bins=(range(minimum(isis), maximum(isis) + 1.0, step=0.01)),
# hist_isi = histogram(isis, bins=(range(0.0, 10.0, step=0.01)),
# 					 legend=:none, title="ISI Histogram")
# l = @layout [ a{0.8w} b ]
# p2 = plot(hist_isi, hist_rad, layout=l)
# 
# l = @layout [ a{0.65h};  b ]
# p = plot(p1, p2, layout=l)
# 
# p_isi = plot(isis, line=:stem, marker=:circ)
# plot!(p_isi, 1:length(isis), [5.0 for k ∈ 1:length(isis)], ls=:dash, color=:black)

end
