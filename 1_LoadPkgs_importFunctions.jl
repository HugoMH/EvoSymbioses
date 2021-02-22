using BenchmarkTools # @btime

using MultivariateStats
using Distributions
# using Gadfly
# using Combinatorics # combinations
using Random # shuffle!, rand!
using Plots
using StatsPlots
using Colors
using DataFrames
using Dates
using StatsBase # Weights
using Distributed
# using HypothesisTests # OneSampleTTest
using JLD
# using GLM # ccdf(FDist(1, df_residual(mm)), abs2(tt)))

# import Base.prod ; function prod(X::Array{Any,1}) 1.0 end
import Base.keys, Base.getindex, Base.length, Base.setindex!, Base.println
× = * # helps to copy-past formulas from word
Float = Float64
Value = Union{String,Symbol,Number}

function flat(arr)
    rst = Any[]
    T = Vector{Type}()
    grep(v) = for x in v
        if isa(x, Array) | isa(x, Tuple) grep(x) else push!(rst, x)  ; push!(T,typeof(x)) end
    end
    grep(arr)
    if length(unique(T))==1
        rst = convert(Vector{unique(T)[1]},rst)
    end
    rst
end
