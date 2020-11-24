module VThetaTuningCurve

include("/home/dahlbom/research/rinzelmodel/vthetamodel/VTheta.jl")
include("/home/dahlbom/research/rinzelmodel/util/SpikeUtils.jl")
using Plots
using Plots.Measures
using Statistics: mean
using .VTheta: vthetarun
using .SpikeUtils: isi, isi_ao, vectorstrength

params = VTheta.params_default


sinusoids(t, f0, p) = max(0.0, sin(2π*f0*t) + sin(2π*f0*p*t))

mistunings = collect(range(0.88, stop=1.12, step=0.01))
amps = collect(range(0.2, 0.5, step=0.01))
amp = 0.4
dt = 1/50000*1000
dur = 2.0*1000
f0 = 200.0/1000

spikes_all = []
isis_all = []
v_all = []
θ_all = []
plots = []

ests = []
hists = []
hists_ao = []
bins = collect(4.25:0.025:5.75)
bins_ao = collect(0.0:0.1:50.0)

for mt ∈ mistunings
	v0 = 0.1
	θ0 = 0.3
	infunc(t) = amp*sinusoids(t, f0, mt)
	(spikes, v, θ) = vthetarun(v0, θ0, dur, dt, infunc, params)
	#vs = vectorstrength(spikes)[1]
	isis = isi(spikes)
	isis_ao = isi_ao(spikes)
	isis1 = filter(x -> 0.9*5.0 < x < (1/.9)*5.0, isis)
	est = 1000.0/mean(isis1)
	push!(ests, est)
	p_vs = plot([0,0], [0,1.0], 
				line=:stem, 
				lc=:red, 
				legend=:none, 
				xlims=(-0.5, 0.5),
				ylims=(0.0, 1.0), 
				xticks=[],
				yticks=[],
				yaxis=false,
				lw=5.0)
	l = @layout [ a{0.9w} b ]
	hist_ao = histogram(isis_ao, 
						bins=bins_ao, 
						legend=:none,
						ylims=(0.0, 0.25),
						title="$mt",
						titlefontsize=18,
						normalize=:probability)
	push!(hists_ao, plot(hist_ao, p_vs, layout=l))
	l = @layout [ a{0.9w} b ]
	hist = histogram(isis,
					 bins=bins,
					 legend=:none,
					 ylims=(0.0, 405.0),
					 title="$mt",
					 titlefontsize=8)
	push!(hists, hist)
	#push!(hists, plot(hist, p_vs, layout=l))
	# push!(spikes_all, spikes)
	# push!(θ_all, θ)
	# push!(v_all, v)
	# push!(isis_all, isis)
end

p = plot(mistunings, ests, legend=:none, title="$amp", size=(1800,900), margin=5mm)
p_hist = plot(hists..., size=(1800,900))
p_hist_ao = plot(hists_ao..., size=(1800,900))


end
