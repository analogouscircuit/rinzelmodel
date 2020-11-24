module MatlabModel

using MATLAB
using DataFrames

"""
	traces(sig, fs)

Get the cfs and voltage traces of the SD units, using the Huang and Rinzel model.
"""
function traces(sig, fs;
				bi=true,
				width=2.0,
				bflist= 50.0 .* 2.0 .^ collect(0.0:0.25:7.0),
				scriptpath = "/home/dahlbom/research/pitchmodelhrmat")
	mxcall(:addpath, 0, scriptpath)
	cfs, traces = mxcall(:AN_SD_from_sig, 2, mxarray(sig), fs, width, bi)
	tracedata = DataFrame(cfs = Float64[], trace = Any[])
	for set ∈ traces
		for (k, cf) ∈ enumerate(cfs)
			push!(tracedata, (cf, interpf(set[k][:,1], set[k][:,2])))
		end
	end
	return tracedata 
end


function lininterp(x0, x, y)
	(x0 <= x[1])   && return x[1]
	(x0 >= x[end]) && return x[end]
	i1 = searchsortedfirst(x, x0)
	i1 == length(x) && return x[end]
	k = (x0 - x[i1]) / (x[i1+1] - x[i1])	
	y[i1] + k * (y[i1+1] - y[i1])
end


function interpf(x, y)
	x0 -> lininterp(x0, x, y)
end


end
