module HRArrayTest

using DifferentialEquations
using Plots 
include("/home/dahlbom/research/rinzelmodel/perimodels/PeriToInput.jl")
using .PeriToInput: MeddisParams, ZilanyParams, sigtospikedata, spikedatatochannels, channelstosdinputfuncs, channelstofuncs
include("/home/dahlbom/research/rinzelmodel/sdmodels/HR.jl")
using .HR: hr
include("/home/dahlbom/research/rinzelmodel/util/SpikeUtils.jl")
using .SpikeUtils: findspikes, isi, phasefromspikes
include("/home/dahlbom/research/AdaptiveSieves/util/StimGen.jl")
using .StimGen: edgestim


println("Setting up test signal and parameters...")
## Set up signal and run through periphery, get SD inputs
f0 = 200.0
τ0 = 1000.0 / f0
fs = 50e3 	# must be greater that 100k for Zilany
db = 80.0
dur_sig = 0.10
dur_pad = 0.050
dur = 2*dur_pad + dur_sig
#params = ZilanyParams(fs, db, (10, 7, 3), (125.0, 6000.0, 29), "human", 1)
cfs = exp.( range(log(125), stop=log(6000.0), length=29) )
cfs_sd = exp.(range(log(125.0), stop=log(3000.0), length=20))
@time params = MeddisParams(fs, db, (200, 0, 0), cfs, 1)
### signal
#### out-of-phase sinusoids
# sig = max.( 0.0, sin.(2π*f0 .* t) )
# sig .+= max.( 0.0, sin.(2π*2*f0 .* t) )
# sig = max.( 0.0, sin.(2π*3*f0 .* t) )
# sig .+= max.( 0.0, sin.(2π*4*f0 .* t) )
# sig .+= max.( 0.0, sin.(2π*5*f0 .* t) )
# pad = zeros(Int(round(dur_pad*fs)))
# sig = [pad; sig; pad]
sig = PeriToInput.teststim(f0, fs, dur_sig, dur_pad)
for k ∈ 2:8
	k == 4 ? p = 1.2 : p = 1.0
	sig .+= PeriToInput.teststim(f0*k*p, fs, dur_sig, dur_pad)
end
t = collect(0:length(sig)-1) ./ fs
#### noise edge pitch
# fe = 800.0
# sig = edgestim(fe, fs, dur; kind=:hp) 

println("Running periphery...")
@time begin
	(cfs_an, spikes) =  sigtospikedata(sig, params) |> spikedatatochannels
	g_syns = channelstosdinputfuncs(cfs_sd, cfs_an, spikes; g_E=1.5)
	an_funcs = channelstofuncs(spikes)
	inputs = [f.(t) for f ∈ g_syns]
	an_outputs = [f.(t) for f ∈ an_funcs]
end


println("Simulating SD units...")
@time begin
	# Solve for each channel
	#u0 = [-63.641, 0.431]
	u0 = [-60, 0.38]
	tspan = (0.0, t[end]*1000.0)
	sols_sd = []
	spikes_sd = []
	isis_sd = []
	for input ∈ g_syns
		g(t) = input(t/1000.0)
		p = [g] 
		prob = ODEProblem(hr, u0, tspan, p)
		sol = solve(prob, AutoTsit5(Rosenbrock23()), reltol=1e-5, abstol=1e-5)
		#sol = solve(prob)
		spikes = findspikes(sol; thresh=-20.0)
		isis = isi(spikes)
		push!(sols_sd, sol)
		push!(spikes_sd, spikes)
		push!(isis_sd, isis)
	end
end

println("Setting up plots...")
using Plots.Measures
ps = []
@time for sol ∈ sols_sd
	push!(ps, plot(sol, legend=:none))
end
p_sols = plot(ps...)

# Auditory Nerve output
ps = []
for an ∈ an_outputs[end:-2:1]
	push!(ps, plot(t, an, legend=false, ylims=(0,30), margin=0mm))
end
p_an = plot(ps..., layout=(length(ps), 1))

# Inputs to SD Units
ps = []
for sdin ∈ inputs[end:-2:1]
	push!(ps, plot(t, sdin, legend=false, ylims=(0,120), margin=0mm))
end
p_in = plot(ps..., layout=(length(ps), 1))

# Output of SD Units 
ps = []
for sdout ∈ [sol(t .* 1000.0)[1,:] for sol ∈ sols_sd[end:-2:1]]
	push!(ps, plot(t, sdout, legend=false, ylims=(-60.0, 25.0), margin=0mm))
end
p_out = plot(ps..., layout=(length(ps), 1))

# Put it all together
p = plot(p_an, p_in, p_out, layout = @layout [a b c])

isis_ao = vcat(isis_sd...)
hist = histogram(isis_ao, bins=(range(0.0, stop=min(5*τ0, maximum(isis_ao)), step=0.1)))


end
