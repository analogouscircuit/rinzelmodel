module StimForPeriMTH

using Plots
using WAV
using PyCall

include("../../AdaptiveSieves/util/StimGen.jl")

co = pyimport("cochlea")

fs = 50e3
f0 = 155.0
dur = 0.1
silence1 = 0.05
silence2 = 0.05
mt_range = range(0.88, stop=1.12, length=25) 
num_h = 12
mt_num = 4
mts   = [1.0 for k ∈ 1:num_h]
amps  = [1.0 for k ∈ 1:num_h]
phase = [0.0 for k ∈ 1:num_h]
db = 80.0
stimdir = "/home/dahlbom/research/rinzelmodel/stimuli/"


for mt ∈ mt_range
	silence1_n = zeros( Int(round(silence1*fs)) )
	silence2_n = zeros( Int(round(silence2*fs)) ) 
	mts = [1.0 for k ∈ 1:num_h]
	mts[mt_num] = mt
	stim = StimGen.hc(dur, f0, fs; amps=amps,
					  			   mistunings=mts,
								   phase=phase)
	stim_db = co.set_dbspl(stim, db)
	stim_db = [silence1_n; stim_db; silence2_n]
	name = "mth_$(f0)_$(Int(round(mt_num)))_$(mt).wav"
	wavwrite(stim_db, stimdir * name, Fs=fs)
end


end
