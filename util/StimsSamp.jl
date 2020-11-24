module StimsSamp

using Formatting

include("/home/dahlbom/research/AdaptiveSieves/util/StimGen.jl")
using .StimGen: edgestim 		# Edge pitch


# -------------------- data types --------------------
struct MTHParams
	fs::Float64
	f0::Float64
	dur_on::Float64
	dur_off::Float64
	numharmonics::Int64
	mistunedharmonic::Int64
	mistuning::Float64
end

struct NEPParams
	fs::Float64
	fe::Float64
	dur_on::Float64
	dur_off::Float64
	kind::Symbol 	# :hp or :lp
end

# -------------------- test wrangling --------------------
function params_to_string(params::Union{MTHParams, NEPParams})
	if typeof(params) <: MTHParams
		name = "mth"
	elseif typeof(params) <: NEPParams
		name = "nep"
	else
		error("Unrecognized stimulus type")
	end
	for f ∈ fieldnames(typeof(params))
		name *= format("_{}", getfield(params, f))
	end
	return replace(name, "." => "-")
end

# -------------------- stimulus generation --------------------

"""
	nepstim(fe, fs, dur_on, dur_off; kind=:hp, rampdur=0.005)

Edge pitch stimulus
"""
function nepstim(fe, fs, dur_on, dur_off; kind=:hp, rampdur=0.005)
	nep = edgestim(fe, fs, dur_on; kind=kind)	
	silence = zeros( Int(round(fs*dur_off)) )
	ramp = sinramp(rampdur, fs)
	nep[1:length(ramp)] .*= ramp
	nep[end-length(ramp)+1:end] .*= ramp[end:-1:1]
	[ silence; nep; silence ]
end


"""
	sinstim(f0, fs, dur_on, dur_off; rampdur=0.005)

Generate a single sinuoid stimulus in the manner of Huang and Rinzel's (2016) paper. Unit height. Appropriate scaling must be performed separating (e.g. into dB SPL).   
"""
function sinstim(f0, fs, dur_on, dur_off; rampdur=0.005)
	silence = zeros( Int(round(fs*dur_off)) )
	t_on = collect( range( 0, dur_on, step=1/fs) )
	sinusoid = sin.( 2π*f0 .* t_on)
	ramp = sinramp(rampdur, fs)
	sinusoid[1:length(ramp)] .*= ramp
	sinusoid[end-length(ramp)+1:end] .*= ramp[end:-1:1]
	[ silence; sinusoid; silence ]
end


"""
	sinramp(dur, fs; kind=:on)

Sine ramp.
"""
function sinramp(dur, fs; kind=:on)
	t = collect( range(0, dur, step=1/fs) )
	f0 = 1.0 / (dur * 4.0)
	ramp = sin.( 2π*f0 .* t )
	if kind != :on
		ramp = ramp[end:-1:1]
	end
	return ramp
end


end
