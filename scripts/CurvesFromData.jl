module CurvesFromData

# -------------------- Load Packages and Modules --------------------
using Plots
using DataFrames
using DataFramesMeta
using Statistics: mean, std
using JLD
include("/home/dahlbom/research/rinzelmodel/util/SpikeUtils.jl")
using .SpikeUtils: isi, isi_ao

# -------------------- Load Trial Data --------------------
data_home = "/home/dahlbom/research/rinzelmodel/data/"
mistuned_hs = 3:4
widths = 0.25:0.25:2
bis = [true false]
dbs = [60.0, 80.0] 
combos = Iterators.product(mistuned_hs, widths, bis, dbs)
names = []
for (mistuned_h, width, bi, db) ∈ combos
	push!(names, data_home*"mtpartial_$mistuned_h-width_$width-bi_$bi-db_$db.jld")
end

# -------------------- Parameters for Data Processing --------------------
f0 = 200.0 	# get this from data later
τ0 = 1/f0
tol = 0.0005
bins = collect(range(0.0, 30.0, step=0.1))

# -------------------- Process Data --------------------
curves = []
for filename ∈ names
	global curves
	spikedata = load(filename, "spikedata") 
	mistunings = @with(spikedata, :mistuning) |> unique
	cfs = @with(spikedata, :cf) 			  |> unique
	ests = []
	stds = []
	#hists = []
	for mt ∈ mistunings
		rows_mt = @where(spikedata, :mistuning .== mt)
		isis = Float64[]
		for cf ∈ cfs
			rows_cf = @where(rows_mt, :cf .== cf)
			spike_set = rows_cf[!, :spikes]
			isis_cf = Float64[]
			for spikes ∈ spike_set
				isis = isi(spikes)
				vcat(isis_cf, isis)
			end
			vcat(isis, isis_cf)
		end
		isis_filt = filter( x -> τ0-tol < x < τ0+tol, isis)
		est = mean(isis_filt)
		err = std(1.0 ./ isis_filt)
		push!(ests, 1/est)
		push!(stds, err)
		#push!(hists, histogram(isis, bins=bins))
	end
	title = split(filename, "/")[end]
	push!(curves, plot(mistunings, ests,
					   yerror=stds,
					   title = title,
					   titlefontsize = 7,
					   legend=:none))
end




end
