module StimForPeriNEP

using Plots
using WAV
using PyCall

include("../../AdaptiveSieves/util/StimGen.jl")

co = pyimport("cochlea")

fs = 50e3
fes = [150.0 * 2 ^ k for k ∈ 0:4]
dur = 0.5
silence1 = 0.05
silence2 = 0.05
kinds = [:hp, :lp]
db = 60.0
stimdir = "/home/dahlbom/research/rinzelmodel/stimuli/"

for kind ∈ kinds
	for fe ∈ fes
		silence1_n = zeros( Int(round(silence1*fs)) )
		silence2_n = zeros( Int(round(silence2*fs)) ) 
		stim = StimGen.edgestim(fe, fs, dur; kind=kind)
		stim_db = co.set_dbspl(stim, db)
		stim_db = [silence1_n; stim_db; silence2_n]
		name = "nep_$(kind)_$(Int(round(fe)))hz_$(Int(round(db)))db.wav"
		wavwrite(stim_db, stimdir * name, Fs=fs)
	end
end


end
