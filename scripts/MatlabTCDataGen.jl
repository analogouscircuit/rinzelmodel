module MatlabTuningCurve

# Packages
using Plots
using PyCall
using DataFrames
using JLD
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


# -------------------- Stimulus Parameters --------------------
println("Generating stimuli...")
dur_on = 0.25
dur_off = 0.05
dur_tot = dur_on + 2dur_off
fs = 50e3
f0 = 200.0
# stim = edgestim(f0, fs, dur_on; kind=:hp)
trialsperstim = 5
mistunings = collect(0.88:0.01:1.12)
mistuned_h = 4
num_h = 8

mistuned_hs = 3:4
widths = 0.25:0.25:2
bis = [true false]
dbs = [60.0, 80.0] 
combos = Iterators.product(mistuned_hs, widths, bis, dbs)
num_combos = length(combos)
count = 0
for (mistuned_h, width, bi, db) ∈ combos
	global count += 1
	println("-------------------- $count of $num_combos --------------------")
	# Set up data frame
	spikedata = DataFrame(f0 = Float64[], 
						  fs = Float64[], 
						  dur_on = Float64[], 
						  dur_tot = Float64[],
						  db = Float64[], 
						  numh = Int64[], 
						  mth = Int64[], 
						  mistuning = Float64[], 
						  cf = Float64[], 
						  spikes = Array{Array{Float64,1},1}(undef,0))

	# -------------------- Run through Matlab Model --------------------
	println("Running through model and calculating pitches...")
	# ests = []
	# ests_μ = []
	#ests_w = []
	#ests_μ_w = []
	for (k, mt) ∈ enumerate(mistunings)
		# Generate Stimulus
		for q ∈ 1:num_h
			q == mistuned_h ? p = mt : p = 1.0
			if q == 1
				stim = sinstim(q*f0*p, fs, dur_on, dur_off) 
			else
				stim .+= sinstim(q*f0*p, fs, dur_on, dur_off)
			end
		end
		stim = setdb(stim, db)
		println("Starting mistuning $(mistunings[k]), $k of $(length(mistunings))...")
		ests_local = []
		#ests_local_w = []
		print("Of $trialsperstim: ")
		for q ∈ 1:trialsperstim
			print(" $q ")
			sd_data = traces(stim, fs; bi=bi, width=width); 	# run Matlab model and get SD voltage traces
			cfs = unique(sd_data[!, :cfs])
			isis_chan = []
			spikes_chan = []
			for cf ∈ cfs
				ps =[]
				spikes = []
				for trace ∈ sd_data[sd_data.cfs .== cf, :trace]
					newspikes = SpikeUtils.findspikes(trace, 0.0, dur_tot; scale=:s)
					push!(spikes, newspikes)
					push!(spikedata, (f0, 
									  fs, 
									  dur_on,
									  dur_off,
									  db,
									  num_h,
									  mistuned_h,
									  mistunings[k],
									  cf,
									  newspikes))
				end
				push!(isis_chan, vcat([ SpikeUtils.isi(spks) for spks ∈ spikes]...))
				push!(spikes_chan, spikes)
			end
			# spikes_summary = vcat( vcat(spikes_chan...)... )
			# isis_summary = vcat( isis_chan... )
			# τs = collect(0.0:0.00001:0.01);
			# bins = (τs[1:end-1] .+ τs[2:end]) ./ 2
			# smoothed = SpikeUtils.smoothhistfunc(isis_summary).(bins)
			# weights = SpikeUtils.vsweights(spikes_summary, τs)
			# #smoothed_w = smoothed .* weights
			# maxval = maximum(smoothed)
			# #maxval_w = maximum(smoothed_w)
			# i = findpeaknear(bins, smoothed, 0.005)
			# #i_w = findpeaknear(bins, smoothed_w, 0.005)
			# peak = peakinterp(bins[i-1:i+1], smoothed[i-1:i+1])
			# #peak_w = peakinterp(bins[i_w-1:i_w+1], smoothed_w[i_w-1:i_w+1])
			# push!(ests_local, 1.0/peak)
			# #push!(ests_local_w, 1.0/peak_w)
		end
		# push!(ests, ests_local)
		# push!(ests_μ, mean(ests_local))
		# # push!(ests_w, ests_local_w)
		# # push!(ests_μ_w, mean(ests_local_w))
		println()
	end

	name = "mtpartial_$mistuned_h-width_$width-bi_$bi-db_$db"
	println("Saving data...")
	save("/home/dahlbom/research/rinzelmodel/data/$name.jld",
		 "spikedata", spikedata,
		 "mistunings", mistunings)
	println("Saved.")
	#errs = [std(e) for e ∈ ests]
	# errs_w = [std(e) for e ∈ ests_w]
	#p = plot(mistunings, ests_μ, yerror=errs)
	# p_w = plot(mistunings, ests_μ_w, yerror=errs_w)
	end

end
