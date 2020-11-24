module PeriToInput

using Plots
using DataFrames
using DataFramesMeta
import Pandas
using PyCall

co = pyimport("cochlea")
coex = pyimport("cochlea.external")
th = pyimport("thorns")
matwrap = pyimport("matlab_wrapper")

################################################################################
# Types
################################################################################
struct ZilanyParams
	fs::Real
	db::Real
	anf_num::Tuple{Int, Int, Int}
	cf::Tuple{Real, Real, Int}
	species::String
	seed::Int
end

struct MeddisParams
	fs::Real
	db::Real
	anf_num::Tuple{Int, Int, Int}
	cf::Union{AbstractFloat, Tuple{Real, Real, Int}, Array{Float64, 1}}
	seed::Int
end

################################################################################
# Functions
################################################################################
function sigtospikedata(sig, p::MeddisParams)
 	if length(p.cf) == 3
 		cf = collect( exp.( range(log(cf[1]), log(cf[2]), length=cf[3]) ) )
 	else
 		cf = p.cf
 	end
	sig_db = co.set_dbspl(sig, p.db)
	sd = coex.run_matlab_auditory_periphery(sig_db,
											 p.fs,
											 anf_num=p.anf_num,
											 cf=cf,
											 seed=p.seed)
	sd_df = sd |> Pandas.DataFrame |> Pandas.reset_index |> DataFrames.DataFrame
	select!(sd_df, Not(:index))
end


function sigtospikedata(sig, p::ZilanyParams)
	sig_db = co.set_dbspl(sig, p.db)
	sd  = co.run_zilany2014(sig_db,
							p.fs,
							anf_num=p.anf_num,
							cf=p.cf,
							species=p.species,
							seed=p.seed)
	sd |> Pandas.DataFrame |> DataFrames.DataFrame
end

function spikedatatochannels(sd)
	cfs = unique( @with(sd, :cf) )
	spikes = [Float64[] for cf ∈ cfs]
	for (k,cf) ∈ enumerate(cfs)
		spikes[k] = vcat( @with(sd, sd[ :cf .== cf, ^(:spikes)])... ) |> sort
	end
	return(cfs, spikes)
end

function channelstofuncs(spikes)
	[ t -> g_AN(t, spks) for spks ∈ spikes ]
end

function channelstosdinputfuncs(cfs_sd, cfs_an, spikes; g_E=1.5)
	[ t -> g_syn(t, cf_sd, cfs_an, spikes; g_E=g_E) for cf_sd ∈ cfs_sd ]
end

function ω(f_an, f_sd)
	σ = 2.0 	#octaves -- hardcode this in equation after verification
	if f_an >= f_sd
		return exp( - abs( log(f_an) - log(f_sd) ) ^ 2 / σ ^ 2 )
	end
	return 0.0
end

function η(t; τ_E = 0.07e-3)
	if t < 0.0
		return 0.0
	end
	return (t / τ_E) * exp( 1.0 - t / τ_E)
end

function g_AN(t, spikes; dt=1.05e-3, g_E=1.5)
	# NOTE: paper gives g_E=1.5, but the code uses 34 -- big difference!
	# NOTE: matlab code uses 15 * τ_E, where τ_E=0.07e-3 -- .75e-3 seems
	# 		sufficient
	i_lo = searchsortedfirst(spikes, t-dt)
	i_hi = searchsortedlast(spikes, t)
	g_E * sum(η.( t .- spikes[i_lo:i_hi] ))
end

function g_syn(t, cf_sd, cfs_an, spikes_an; g_E=34)
	out = 0.0
	for (k,cf_an) ∈ enumerate(cfs_an)
		if cf_an < cf_sd
			continue
		end
		out += ω(cf_an, cf_sd) * g_AN(t, spikes_an[k]; g_E=g_E)
	end
	out
end


end
