module MatlabTuningCurve

# Packages
using Plots
using JLD
using PyCall
using DataFrames
using Statistics: mean, std

# PyCall Imports
co = pyimport("cochlea"); setdb = co.set_dbspl

# Local Modules
include("/home/dahlbom/research/AdaptiveSieves/util/PeakPick.jl")
using .PeakPick: findpeaknear, peakinterp
include("/home/dahlbom/research/rinzelmodel/scripts/MatlabModel.jl")
using .MatlabModel: traces 
include("/home/dahlbom/research/rinzelmodel/util/StimsSamp.jl")
using .StimsSamp: nepstim, sinstim
include("/home/dahlbom/research/rinzelmodel/util/SpikeUtils.jl")


# -------------------- Stimulus and Model Parameters --------------------
println("Generating stimuli...")
fs = 50e3
dur_on = 0.5
dur_off = 0.05
dur_tot = dur_on + 2dur_off
num_h = 12 
mistunings = collect(0.88:0.01:1.12)
bin_size = 2.5e-5

f0s = [200.0]
dbs = [80.0]
bis = [true, false]
widths = [0.5, 1.0, 2.0, 4.0]
mistuned_hs = [3]
trialsperstims = [5]

# Generate data for different parameter settings
for (f0, db, mistuned_h, bi, width, trialsperstim) ∈ Iterators.product(f0s,
																	   dbs,
																	   mistuned_hs,
																	   bis,
																	   widths,
																	   trialsperstims)

	title = "f0-$f0 dB-$db mth-$mistuned_h bi-$bi width-$width numtrials-$trialsperstim"

	# -------------------- Generate stimuli --------------------
	stims = []
	for mt ∈ mistunings
		# stim = sinstim(f0, fs, dur_on, dur_off)
		for k ∈ 1:num_h
			k == mistuned_h ? p = mt : p = 1.0
			if k == 1
				stim = sinstim(k*f0*p, fs, dur_on, dur_off) 
			else
				stim .+= sinstim(k*f0*p, fs, dur_on, dur_off)
			end
		end
		stim = setdb(stim, db)
		push!(stims, stim)
	end

	# -------------------- Run through Matlab Model --------------------
	println(title)
	println("Running through model and calculating pitches...")
	ests = []
	ests_μ = []
	# isis_mts = []
	spikes_all = Dict()
	# isis_all = Dict()
	for (k, stim) ∈ enumerate(stims)
		spikes_all[mistunings[k]] = []
		# isis_all[mistunings[k]] = []
		println("Starting mistuning $(mistunings[k]), $k of $(length(mistunings))...")
		ests_local = []
		print("Of $trialsperstim: ")
		for q ∈ 1:trialsperstim
			print(" $q ")
			sd_data = traces(stim, fs; bi=bi, width=width); # run Matlab model and get SD voltage traces
			cfs = unique(sd_data[!, :cfs])
			# isis_chan = []
			spikes_chan = []
			for cf ∈ cfs
				ps =[]
				spikes = []
				for (k, trace) ∈ enumerate(sd_data[sd_data.cfs .== cf, :trace])
					push!(spikes, SpikeUtils.findspikes(trace, 0.0, dur_tot; scale=:s))
				end
				# push!(isis_chan, vcat([ SpikeUtils.isi(spks) for spks ∈ spikes]...))
				push!(spikes_chan, spikes)
			end
			spikes_summary = vcat( vcat(spikes_chan...)... )
			# isis_summary = vcat( isis_chan... )
			push!(spikes_all[mistunings[k]], spikes_summary)
			# push!(isis_all[mistunings[k]], isis_summary)
			# τs = collect(0.0:0.00001:0.01);
			# bins = (τs[1:end-1] .+ τs[2:end]) ./ 2
			# smoothed = SpikeUtils.smoothhistfunc(isis_summary; σ=bin_size).(bins)
			# maxval = maximum(smoothed)
			# i = findpeaknear(bins, smoothed, 1.0/f0; tol=0.02)
			# if i != 0
			# 	peak = peakinterp(bins[i-1:i+1], smoothed[i-1:i+1])
			# 	push!(ests_local, 1.0/peak)
			# else
			# 	push!(ests_local, f0) 
			# end
		end
		# push!(ests, ests_local)
		# push!(ests_μ, mean(ests_local))
		println()
	end

	# errs = [std(e) for e ∈ ests]
	# p = plot(mistunings, ests_μ, yerror=errs,
	# 		 title=title, legend=:none, size=(1800,900))

	save("/home/dahlbom/research/rinzelmodel/data/"*title*".jld",
		 # "mistunings", mistunings,
		 "spikes_all", spikes_all,
		 # "isis_all", isis_all,
		 # "ests_μ", ests_μ,
		 # "errs", errs,
		 "title", title,
		 "dur", dur_on)
end


end
