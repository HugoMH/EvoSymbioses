if FASTMATH== "@fastmath"
    FASTMATHstart = "@fastmath(begin"
    FASTMATHend = "end)"
elseif FASTMATH== ""
    FASTMATHstart = ""
    FASTMATHend = ""
else
    error("FASTMATH must be either \"@fastmath\" or \"\" but it was \"$FASTMATH\"")
end

if INBOUND== "@inbounds"
    INBOUNDstart = "@inbounds(begin"
    INBOUNDend = "end)"
elseif INBOUND== ""
    INBOUNDstart = ""
    INBOUNDend = ""
else
    error("INBOUND must be either \"@inbounds\" or \"\" but it was \"$INBOUND\"")
end

if NO_INLINE == "@inline / @noinline"
    NOINLINE = "@noinline"
    INLINE = "@inline"
elseif NO_INLINE == ""
    NOINLINE = ""
    INLINE = ""
else
    error("NO_INLINE must be either \"@inline / @noinline\" or \"\" but it was \"$NO_INLINE\"")
end

#############################
##### General functions #####
#############################

function FastLinearModel(y::Vector{Float})
    n = length(y)
    x = collect(1:n)
    x_Mx = x .- (n/2+0.5)
    My = mean(y)
    slop = sum((x_Mx).*(y.-My)) / sum(x_Mx.^2)
    # y_hat = My .+ slop .* (x_Mx)
    # VarSlop = sum((y.-y_hat).^2) / ((n-2)*sum(x_Mx.^2))
    # tobs = (slop^2) / VarSlop
    # return([slop, ccdf(FDist(1, n-2), tobs)])
    return slop
end

function insertsort!(x::AbstractArray)
  # from https://github.com/Dawny33/SearchSortAlgos.jl/blob/master/src/Sorting.jl
  for index = 2:length(x)
    current  = x[index]
    position = index
    # For sublist sorting
    while position > 1 && x[position - 1] > current
      x[position] = x[position - 1]
      position -= 1
    end
    x[position] = current
  end
  return x
end

function getindex(X::Tuple,ii::Union{UnitRange,Vector{Int},NTuple{N,Int} where N})
    out=Vector{Union{typeof.(X)...}}(undef,length(ii))
    for (iout,i) in enumerate(ii)
        out[iout] = X[i]
    end
    return(out)
end
# function aA_(x::Tuple) x[[1,3,5]] end
# function aT_(x::Tuple) x[(1,3,5)] end
#
# @btime for _ in 1:1e3 aA_((1,3,4,5,6)) end
# # 267.024 μs (4000 allocations: 312.50 KiB)
# @btime for _ in 1:1e3 aT_((1,3,4,5,6)) end
# # 21.459 μs (1000 allocations: 109.38 KiB)

function getindex(X::Array,ii::NTuple{N,Int} where N)
    out=Vector{Union{typeof.(X)...}}(undef,length(ii))
    for (iout,i) in enumerate(ii)
        out[iout] = X[i]
    end
    return(out)
end
# function aA_(x::Array) x[[1,5]] end
# function aT_(x::Array) x[(1,5)] end
#
# @btime for _ in 1:1e3 aA_([1,3,5,65,6]) end
# @btime for _ in 1:1e3 aT_([1,3,5,65,6]) end

struct ExtendableVector{T}
    x::Vector{T}
end

str="""
$INLINE function getindex(X::ExtendableVector{Int},i::Int)
    $INBOUNDstart
    if length(X.x) < i
        append!(X.x,zeros(Int,i - length(X.x)))
    end
    return(X.x[i])
    $INBOUNDend
end"""
eval(Meta.parse(str))
str="""
$INLINE function setindex!(X::ExtendableVector{Int}, val::Int, i::Int)
    $INBOUNDstart
    if length(X.x) < i
        append!(X.x,fill(NaN,i - length(X.x)))
    end
    X.x[i] = val
    $INBOUNDend
end"""
eval(Meta.parse(str))

str="""
$INLINE function getindex(X::ExtendableVector{Symbol},i::Int)
    $INBOUNDstart
    if length(X.x) < i
        append!(X.x,fill(Symbol(),i - length(X.x)))
    end
    return(X.x[i])
    $INBOUNDend
end"""
eval(Meta.parse(str))
str="""
$INLINE function setindex!(X::ExtendableVector{Symbol}, val::Symbol, i::Int)
    $INBOUNDstart
    if length(X.x) < i
        append!(X.x,fill(NaN,i - length(X.x)))
    end
    X.x[i] = val
    $INBOUNDend
end"""
eval(Meta.parse(str))


str="""
$INLINE function getindex(X::ExtendableVector{Float},i::Int)
    $INBOUNDstart
    if length(X.x) < i
        append!(X.x,fill(NaN,i - length(X.x)))
    end
    return(X.x[i])
    $INBOUNDend
end"""
eval(Meta.parse(str))
str="""
$INLINE function setindex!(X::ExtendableVector{Float}, val::Float, i::Int)
    $INBOUNDstart
    if length(X.x) < i
        append!(X.x,fill(NaN,i - length(X.x)))
    end
    X.x[i] = val
    $INBOUNDend
end"""
eval(Meta.parse(str))

str="""
$INLINE function getindex(X::ExtendableVector{ExtendableVector{T}},i::Int) where T
    $INBOUNDstart
    if length(X.x) < i
        append!(X.x,[ExtendableVector{T}([]) for _ in 1:(i - length(X.x)) ])
    end
    return(X.x[i])
    $INBOUNDend
end"""
eval(Meta.parse(str))

length(X::ExtendableVector) = length(X.x)
str="""
$INLINE function setindex!(X::ExtendableVector{Int}, val::Int, i::Int)
    $INBOUNDstart
    if length(X.x) < i
        append!(X.x,zeros(Int,i - length(X.x)))
    end
    X.x[i] = val
    $INBOUNDend
end"""
eval(Meta.parse(str))

struct BiDirVector{T}
    neg::ExtendableVector{T}
    pos::ExtendableVector{T}
end
BiDirVector() = BiDirVector(ExtendableVector{Int}([]),ExtendableVector{Int}([]))

str="""
$INLINE function getindex(X::BiDirVector,i::Int)
    $INBOUNDstart
    if i>0
        return(X.pos[i])
    else
        return(X.neg[-i])
    end
    $INBOUNDend
end"""
eval(Meta.parse(str))

str="""
$INLINE function setindex!(X::BiDirVector, val::Int, i::Int)
    $INBOUNDstart
    if i>0
        X.pos[i] = val
    else
        X.neg[-i]= val
    end
    $INBOUNDend
end"""
eval(Meta.parse(str))

length(X::BiDirVector) = length(X.neg) + length(X.pos)
length(X::BiDirVector,sign::Val{-1}) = length(X.neg)
length(X::BiDirVector,sign::Val{1 }) = length(X.pos)

function NaNrm(X::Vector{T}) where T <: Number
    X[.!isequal.(X, NaN)]
end

function RecursiveArrayToTuple(a,i=0)
    A = deepcopy(a)
    a = Array{Any,1}()
    for aa in A push!(a,aa) end
    for i in eachindex(a)
        if typeof(a[i])<:Vector
            a[i] = RecursiveArrayToTuple(a[i])
        end
    end
    a = tuple(a...)
    return a
end

str = """
$INLINE function RandMultivarHyperGeometric(N::Int,Pop::Vector{Int})
    $FASTMATHstart
    $INBOUNDstart
    sample = Vector{Int}(undef,length(Pop))
    CumPopDiff = cumsum(Pop) .- Pop
    for i in length(Pop):-1:2 #         # successes in Pop  # failures in Pop      sample size
        sample[i] = rand(Hypergeometric(  Pop[i]           ,CumPopDiff[i]      ,   N    ))
        N -= sample[i]
    end
    sample[1] = N
    return(sample)
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

str = """
$INLINE function RandMultivarHyperGeometric!(N::Int,Pop::Vector{Int},sample::Vector{Int})
    $FASTMATHstart
    $INBOUNDstart
    CumPopDiff = cumsum(Pop) .- Pop
    for i in length(Pop):-1:2 #         # successes in Pop  # failures in Pop      sample size
        sample[i] = rand(Hypergeometric(  Pop[i]           ,CumPopDiff[i]      ,   N    ))
        N -= sample[i]
    end
    sample[1] = N
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))


str = """
function RandMultivarHyperGeometric_!(sample::Vector{Int},CumPopDiff::Vector{Int},N::Int,Pop::Vector{Int})
    $FASTMATHstart
    $INBOUNDstart
    cumsum!(CumPopDiff,Pop)
    CumPopDiff .-= Pop
    for i in length(Pop):-1:2 #         # successes in Pop  # failures in Pop      sample size
        sample[i] = rand(Hypergeometric(  Pop[i]           ,CumPopDiff[i]      ,   N    ))
        N -= sample[i]
    end
    sample[1] = N
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

str = """
$INLINE function RecursiveRandMultivarHyperGeometric(N::Int,N__::Vector{Vector{Int}})
    $FASTMATHstart
    $INBOUNDstart
    $PRINT if !all(sum.(N__) .== N) error("Some vectors of N__ do not sum to N") end
    out::Vector{Vector{Int}} = [[i] for i in eachindex(N__[1]) if N__[1][i] > 0]
    N_::Vector{Int} = Vector{Int}(filter(n -> n > 0,N__[1]))
    Sample::Vector{Int}       = Vector{Int}(undef,length(N_))
    CumPopDiff::Vector{Int}   = Vector{Int}(undef,length(N_))
    out2::Vector{Vector{Int}} = Vector{Vector{Int}}()
    N_2::Vector{Int}          = Vector{Int}()
    # RecursiveRandMultivarHyperGeometric!(view(N__,2:length(N__)),out,N_,Sample, CumPopDiff, out2, N_2)
    RecursiveRandMultivarHyperGeometric!(N__[2:end]             ,out,N_,Sample, CumPopDiff, out2, N_2)
    return(N_,out)
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

str = """
$INLINE function RecursiveRandMultivarHyperGeometric!(N__::Vector{Vector{Int}}, out::Vector{Vector{Int}}, N_::Vector{Int}, Sample::Vector{Int}, CumPopDiff::Vector{Int}, out2::Vector{Vector{Int}}, N_2::Vector{Int}  )
# function RecursiveRandMultivarHyperGeometric!(N__::SubArray{Array{Int64,1},1,Array{Array{Int64,1},1}} , out::Vector{Vector{Int}}, N_::Vector{Int}, Sample::Vector{Int}, CumPopDiff::Vector{Int}, out2::Vector{Vector{Int}}, N_2::Vector{Int}  )
    $FASTMATHstart
    $INBOUNDstart
    empty!(out2)
    empty!(N_2)
    resize!(Sample    , length(N_))
    resize!(CumPopDiff, length(N_))
    # Ndone::Vector{Int} = fill(0,length(N_))
    for i in 1:length(N__[1])
        if N__[1][i] > 0
            RandMultivarHyperGeometric_!(Sample,CumPopDiff,N__[1][i],N_)
            for (ii,n) in enumerate(Sample)
                if n > 0
                    # push!(out2,vcat(out[ii],i))
                    push!(out2,[X for X in out[ii]])
                    push!(out2[end],i)
                    push!(N_2,n)
                    # Ndone[ii] += n
                    N_[ii] -= n
    end ; end ; end ; end
    # out = deepcopy(out2)
    # N_  = deepcopy(N_2)
    empty!(out) ; append!(out,out2)
    empty!(N_ ) ; append!(N_,N_2)
    if length(N__) > 1
        # RecursiveRandMultivarHyperGeometric!(view(N__,2:length(N__)),out,N_,Sample, CumPopDiff, out2, N_2)
        RecursiveRandMultivarHyperGeometric!(N__[2:end],out,N_,Sample, CumPopDiff, out2, N_2)
    end
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

str="""
$INLINE function RecursiveRandMultinom(N::Int,P::Vector{Vector{Float}})
    $FASTMATHstart
    $INBOUNDstart
    out::Vector{Vector{Int}} = Vector{Vector{Int}}()
    N_::Vector{Int} = Vector{Int}()
    for (val,n) in enumerate(rand(Multinomial(N,P[1])))
        RecursiveRandMultinom!(n,P[2:end],out,N_,[val])
    end
    return(N_,out)
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

str="""
$INLINE function RecursiveRandMultinom!(N::Int, P::Vector{Vector{Float}}                             , out::Vector{Vector{Int}}, N_::Vector{Int}, valS::Vector{Int})
# function RecursiveRandMultinom!(N::Int, P::SubArray{Array{Int64,1},1,Array{Array{Float,1},1}}, out::Vector{Vector{Int}}, N_::Vector{Int}, valS::Vector{Int})
    $FASTMATHstart
    $INBOUNDstart
    if length(P) > 1
        for (val,n) in enumerate(rand(Multinomial(N,P[1])))
            if n>0
                valS2::Vector{Int} = [x for x in valS]
                push!(valS2, val)
                # RecursiveRandMultinom!(n,view(P,2:length(P)),out,N_,valS2)
                RecursiveRandMultinom!(n,P[2:end],out,N_,valS2)
            end
        end
    else
        for (val,n) in enumerate(rand(Multinomial(N,P[1])))
            if n>0
                # push!(out,vcat(valS, val))
                push!(out,[X for X in valS])
                push!(out[end],val)
                push!(N_,n)
            end
        end
    end
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

##################################
##### Set the data structure #####
##################################
if GfG const if_GfG_1_else_0 = 1 else const if_GfG_1_else_0 = 0 end
if GfG addOneIfGfG_str = " + 1" else addOneIfGfG_str = "" end


const NallelesPerTrait_infection_success_AllSp = (NallelesPerTrait + DistBetweenHostSp) * NhostSp
temp = [] ; x = 1 ; for sp in 1:NhostSp global xend = x+(NallelesPerTrait-1) ; push!(temp,x:xend) ; global x += NallelesPerTrait+DistBetweenHostSp end ; const AlleleRangePerSp_infection_success = tuple(temp...)

NallelesPerTrait_EachTrait = tuple([if (trait <= N_traits_infection_success) NallelesPerTrait_infection_success_AllSp else NallelesPerTrait end for trait in 1:Ntraits]...)

if (((UInt128(NallelesPerTrait_infection_success_AllSp)^(Ntraits*2)) == 0)  ||  ((UInt128(NallelesPerTrait_infection_success_AllSp)^(Ntraits*2))) == 0)
        error("Choose a lower value for one of these parameters:\nN_traits_infection_success, N_traits_immunity, NallelesPerTrait, DistBetweenHostSp")
end

const MaxIntValues = ((UInt128(2)^63-1), Int128(2)^127-1, UInt128(4)^127-1)
const IntTypes     = (Int64, Int128, UInt128)
const Int_IDh = IntTypes[minimum(filter(i -> ((UInt128(NallelesPerTrait_infection_success_AllSp + if_GfG_1_else_0)^(Ntraits*2))) < UInt128(MaxIntValues[i]), 1:3))]
const Int_IDp = IntTypes[minimum(filter(i -> ((UInt128(NallelesPerTrait_infection_success_AllSp + if_GfG_1_else_0)^(Ntraits  ))) < UInt128(MaxIntValues[i]), 1:3))]

# MaxDist
const MaxDist_infection = NallelesPerTrait_infection_success_AllSp
const MaxDist_immunity  = NallelesPerTrait

#### Circular distance     infection
str = """
function CircularDist_infection(x::Int, y::Int)
    if x>y
        min(x-y   ,$NallelesPerTrait_infection_success_AllSp-x + y  )
    elseif x<y
        min(y-x   ,$NallelesPerTrait_infection_success_AllSp-y + x  )
    else
        0
    end
end"""
eval(Meta.parse(str))

CircularDist_ = []
for x in 1:NallelesPerTrait_infection_success_AllSp
    push!(CircularDist_,[])
    for y in 1:NallelesPerTrait_infection_success_AllSp
        push!(CircularDist_[x],CircularDist_infection(x,y))
    end
end
const Circular_Dist_infection = RecursiveArrayToTuple(CircularDist_)

str="""
$INLINE function CircularDist_infection(x::Int, y::Int)
    $INBOUNDstart
    $(ifelse(GfG,"
    if (x === 0) | (y === 0)
        return $MaxDist_infection
    else
    ",""))
        Circular_Dist_infection[x][y]
    $(ifelse(GfG,"end",""))
    $INBOUNDend
end"""
eval(Meta.parse(str))


#### Circular distance     immunity
str = """
function CircularDist_immunity(x::Int, y::Int)
    if x>y
        min(x-y           ,$NallelesPerTrait -x + y  )
    elseif x<y
        min(y-x           ,$NallelesPerTrait -y + x  )
    else
        0
    end
end"""
eval(Meta.parse(str))

CircularDist_ = []
for x in 1:NallelesPerTrait
    push!(CircularDist_,[])
    for y in 1:NallelesPerTrait
        push!(CircularDist_[x],CircularDist_immunity(x,y))
    end
end
const Circular_Dist_immunity = RecursiveArrayToTuple(CircularDist_)

str="""
$INLINE function CircularDist_immunity(x::Int, y::Int)
    $INBOUNDstart
    $(ifelse(GfG,"
    if (x === 0) | (y === 0)
        return $MaxDist_immunity
    else
    ",""))
        return Circular_Dist_immunity[x][y]
    $(ifelse(GfG,"end",""))
    $INBOUNDend
end"""
eval(Meta.parse(str))

CircularDist_ = nothing

# str="""
# $INLINE function GetOverallDist_infection_success(Dist_::Vector{Int})
#     $INBOUNDstart
#     insertsort!(Dist_)
#     Dist_[$Min_Ntraits_matching_for_infection_success]
#     $INBOUNDend
# end"""
# eval(Meta.parse(str))
#
# str="""
# $INLINE function GetOverallDist_recoveryinnateimmu(Dist_::Vector{Int})
#     $INBOUNDstart
#     insertsort!(Dist_)
#     Dist_[$Min_Ntraits_matching_for_recovery_innate_immunity]
#     $INBOUNDend
# end"""
# eval(Meta.parse(str))
GetOverallDist_infection_success = GetOverallDist_recoveryinnateimmu = minimum

if PhenotypicStochasticity
    global const PhenRobustnessLevels = (1e-15):(1-1e-14)/(NallelesPerTrait_infection_success_AllSp-1):(1-1e-15)
    if length(PhenRobustnessLevels) != NallelesPerTrait_infection_success_AllSp error("length(PhenRobustnessLevels) != NallelesPerTrait_infection_success_AllSp") end
    if (maximum(PhenRobustnessLevels)+0.001) < 1 error("When setting PhenRobustnessLevels of phenotypic noise : (maximum(PhenRobustnessLevels)+0.001) < 1") end
    if (minimum(PhenRobustnessLevels)-0.001) > 0 error("When setting PhenRobustnessLevels of phenotypic noise : (minimum(PhenRobustnessLevels)-0.001) > 0") end
    #
    function GetSymetricalGeom(p::Float,N::Int)
       P = map(k -> (1-p)^(k)*p , 1:N)
       P = P ./ sum(P)
       P[1] = P[1]*2 # not changing is not duplicated
       P = [reverse(P[2:end])...,P...]
       P = P ./ 2
       P = P ./ sum(P)
       return(P)
    end
    const SymetricalProbilitiesGeomDistr_Parasit = RecursiveArrayToTuple([GetSymetricalGeom(p,NallelesPerTrait)                for p in PhenRobustnessLevels])
    const SymetricalProbilitiesGeomDistr_Host    = RecursiveArrayToTuple([GetSymetricalGeom(p,NallelesPerTrait)                for p in PhenRobustnessLevels])
end

for typeLoop in [Vector{Int},NTuple{Ntraits*2,Int_IDp}]
    strLoop = """
        $INLINE function GetHostID(genotype::$typeLoop)
        $FASTMATHstart
        $INBOUNDstart
#         x = 0
#         for i in 2:length(genotype)
#             x += ((NallelesPerTrait_infection_success_AllSp + if_GfG_1_else_0)^(i-1))*(genotype[i] )
#         end
#         x += genotype[1]
#         return(x)
    return(Int_IDh(genotype[1]+ $(
        join( ["$(Int_IDh(NallelesPerTrait_infection_success_AllSp+if_GfG_1_else_0)^(i-1))*(genotype[$i] )" for i in 2:(Ntraits*2)], "+" )
    )))
    $INBOUNDend
    $FASTMATHend
    end"""
eval(Meta.parse(strLoop));
end

for typeLoop in [Vector{Int_IDp},NTuple{Ntraits,Int_IDp}]
    strLoop = """
        $INLINE function GetParasitID(genotype::$typeLoop)
        $FASTMATHstart
        $INBOUNDstart
#       x = 0
#       for i in 2:length(genotype)
#           x += ((NallelesPerTrait_infection_success_AllSp + if_GfG_1_else_0)^(i-1))*(genotype[i] )
#       end
#       x += genotype[1]
#       return(x)
        return(Int_IDp(genotype[1]+ $(
            join( ["$(Int_IDp(NallelesPerTrait_infection_success_AllSp+if_GfG_1_else_0)^(i-1))*(genotype[$i] )" for i in 2:Ntraits], "+" )
        )))
        $INBOUNDend
        $FASTMATHend
        end"""
    eval(Meta.parse(strLoop));
end

# offspring will always be non immunised
str = """
$INLINE function GetCategoryID(IDh::Int_IDh, IDhPhen::Int_IDh,  IDpGenPhen_HostShiftHistory::Symbol)
return(Symbol(:IDcat_, IDh,:_ , IDhPhen,:_, IDpGenPhen_HostShiftHistory))
end"""
eval(Meta.parse(str));

if SimulateImmunisation
    str = """
        $INLINE function GetCategoryID(IDh::Int_IDh, IDhPhen::Int_IDh,  IDpGenPhen_HostShiftHistory::Symbol,  IDi::Symbol)
            return(Symbol(:IDcat_, IDh,:_ , IDhPhen,:_, IDpGenPhen_HostShiftHistory,:_, IDi))
        end"""
    eval(Meta.parse(str));
end

str = """
    $INLINE function GetPgenXphenID_HostShiftHistory(IDp::Int_IDp,  IDpPhen::Int_IDp, HostShiftHistory::Symbol )
        return(Symbol(:IDparasitCat_, IDp,:_, IDpPhen,:_HostShiftHistory_,HostShiftHistory))
    end"""
eval(Meta.parse(str));

str = """
struct HostsGenotypeS
    IDhS::Vector{Int_IDh}
    IDhSphen::Vector{Int_IDh}
    posiInHostsGenotypeS_List::Vector{Int} # POSITION OF THE PHENOTYPE ! not of the genotype
    posiInParasitGenotypeS_List::Vector{Int} # POSITION OF THE PHENOTYPE ! not of the genotype
    IDpS::Vector{Int_IDp}
    IDpSphen::Vector{Int_IDp}
    IDpGenPhen_HostShiftHistoryS::Vector{Symbol}
    IDp_HostShiftHistoryS::Vector{Symbol}
    $(ifelse(SimulateImmunisation
    ,"IDiS::Vector{Symbol} # POSI of the GENOTYPE that have infected the individuals"
    ,""))
    IDcategorieS::Vector{Symbol}
    $(ifelse(GfG
    ,"rate_of_gametes_prodS::Vector{Float}
      rate_of_parasites_emissionS::Vector{Float}"
    ,""))
    virulenceS::Vector{Float}
    Precoveryinnateimmu::Vector{Float}
    $(ifelse(SimulateImmunisation,"Precoveryacquiredimmu::Vector{Float}",""))
    N::Ref{Int}
    N_::Vector{Int}
    dN_::Vector{Int}
    AllelesPool::NTuple{Ntraits,Vector{Float}}
    Ngametes::Vector{Int}
    NparasitesEmitted::Vector{Int}
    nMut::Vector{Int} # Number of mutants for the current allele and trait that decrease [1] and that increase [2]
    MutatedAlleles::NTuple{Ntraits,Vector{Int}} # mutation
    AllelesChanges::Vector{Int} # mutation Vector of size (NallelesPerTrait-1)
    Nrecovery::Vector{Int}
    NmigrEachPop::Vector{Int}
    NmigrEachGenotype::Vector{Int}
end"""
eval(Meta.parse(str));

if SimulateImmunisation
    if GfG
        const ListFieldsOneValPerCat = (:IDhS, :IDhSphen, :IDpS, :IDpSphen, :IDpGenPhen_HostShiftHistoryS, :IDp_HostShiftHistoryS, :IDiS, :IDcategorieS, :rate_of_gametes_prodS, :rate_of_parasites_emissionS, :virulenceS, :Precoveryinnateimmu, :Precoveryacquiredimmu, :N_, :dN_, :posiInHostsGenotypeS_List, :posiInParasitGenotypeS_List)
    else
        const ListFieldsOneValPerCat = (:IDhS, :IDhSphen, :IDpS, :IDpSphen, :IDpGenPhen_HostShiftHistoryS, :IDp_HostShiftHistoryS, :IDiS, :IDcategorieS                                                      , :virulenceS, :Precoveryinnateimmu, :Precoveryacquiredimmu, :N_, :dN_, :posiInHostsGenotypeS_List, :posiInParasitGenotypeS_List)
    end
else
    if GfG
        const ListFieldsOneValPerCat = (:IDhS, :IDhSphen, :IDpS, :IDpSphen, :IDpGenPhen_HostShiftHistoryS, :IDp_HostShiftHistoryS,        :IDcategorieS, :rate_of_gametes_prodS, :rate_of_parasites_emissionS, :virulenceS, :Precoveryinnateimmu,                         :N_, :dN_, :posiInHostsGenotypeS_List, :posiInParasitGenotypeS_List)
    else
        const ListFieldsOneValPerCat = (:IDhS, :IDhSphen, :IDpS, :IDpSphen, :IDpGenPhen_HostShiftHistoryS, :IDp_HostShiftHistoryS,        :IDcategorieS                                                      , :virulenceS, :Precoveryinnateimmu,                         :N_, :dN_, :posiInHostsGenotypeS_List, :posiInParasitGenotypeS_List)
    end
end

str = """
function println(X::HostsGenotypeS; Full = false, Return = false)
    println()
    println()
    A = DataFrame(IDhS=X.IDhS, IDhSphen=X.IDhSphen, IDpS=X.IDpS, IDpSphen=X.IDpSphen, IDp_HostShiftHistoryS=X.IDp_HostShiftHistoryS
    $(ifelse(SimulateImmunisation
    ,", IDiS=X.IDiS"
    ,""))
    , IDcategorieS=X.IDcategorieS
    $(ifelse(GfG
    ,",rate_of_gametes_prodS=X.rate_of_gametes_prodS, rate_of_parasites_emissionS=X.rate_of_parasites_emissionS"
    ,""))
    , virulenceS=X.virulenceS, Precoveryinnateimmu=X.Precoveryinnateimmu, Precoveryacquiredimmu=X.Precoveryacquiredimmu, N_=X.N_, dN_=X.dN_, posiInHostsGenotypeS_List=X.posiInHostsGenotypeS_List, posiInParasitGenotypeS_List=X.posiInParasitGenotypeS_List ,copycols = false)
    if Full         println(A)    else        show(A) ; println()    end
    if Return return(A) end
end"""
eval(Meta.parse(str));

str = """
struct ParasitesPop
    IDpS::Vector{Int_IDp}
    IDpSphen::Vector{Int_IDp}
    IDpGenPhen_HostShiftHistoryS::Vector{Symbol}
    IDp_HostShiftHistoryS::Vector{Symbol}
    posiInParasitGenotypeS_List::Vector{Int}
    ParasitFreq::Vector{Float16}
    N::Ref{Int}
    N_::Vector{Int}
    Nmutating::Vector{Int}
    NmutatingEachTrait::Vector{Int}
    NmuSize::NTuple{Ntraits,Vector{Int}}
    $(ifelse(GfG
    ,"NmutEachAllelesAppearingTrait::NTuple{Ntraits,Vector{Int}}
      rate_of_parasites_emissionS::Vector{Float}"
    ,""))
    NewGenotype::Vector{Int}
    NewId::Ref{Int}
end"""
eval(Meta.parse(str));

str = """
function println(X::ParasitesPop; Full = false, Return = false)
    println()
    println()
    A = DataFrame(IDp_HostShiftHistoryS=X.IDp_HostShiftHistoryS, IDpS=X.IDpS, IDpSphen = X.IDpSphen ,N_ = X.N_ ,posiInParasitGenotypeS_List = X.posiInParasitGenotypeS_List, ParasitFreq = X.ParasitFreq $(ifelse(GfG,", rate_of_parasites_emissionS = X.rate_of_parasites_emissionS","")) ,copycols = false)
    if Full         println(A)      else        show(A)  ; println()  end
    if Return return(A) end
end"""
eval(Meta.parse(str));

struct HostSexS
    HostSexS::NTuple{2,HostsGenotypeS} # one entry for each host sex [always (1,2) for male, female]
    IDsexS::NTuple{2,Int}
    N::Ref{Int}
    N_::Vector{Int}
    K::Ref{Float} # this needs to be updatable to simulate changes in population size
    density_dependence::Ref{Float} # this needs to be updatable to simulate changes in population size
    Sociality::Float
    Ht::Float
    Ht_Sociality::Float
    Ht_one_Sociality_K::Ref{Float} # this needs to be updatable to simulate changes in population size
    ParasitesPop::ParasitesPop
    # Temporary storages for calculus
    NcontactsInPop::Ref{Int} # For the ongoing pop, number of contact between individuals for males and females
    NcontactsCrossSp::Vector{Int}
    P_NcontactsInPop_1_2_3_ect_parasites::Vector{Float16}
    ParasitFreq_x_Pinfection1contact::Vector{Float16}
    Pinfection::Vector{Float}
    virulenceScompar::Vector{Float}
    P_HigherVirulenceS::Vector{Float16}
    P_SameVirulenceS::Vector{Float16}
    Posi_HigherOrEqualVirulenceS_comparToVirulenceInHost::Vector{Int}
    NsuccessfullInfectionEachParasit::Vector{Int}
    RefPhen::Vector{Int}
    HtCrossSp::Ref{Float}
    SpInteract::Ref{Float}
end
##################################
str = """
struct ParasitesGenotype # assuming no noise
    IDp::Int_IDp
    Alleles::NTuple{Ntraits,Int}
    traits_infection_success::NTuple{N_traits_infection_success,Int}
    traits_immunity::NTuple{N_traits_immunity,Int}
    $(ifelse(GfG,"rate_of_parasites_emission::Float",""))
end"""
eval(Meta.parse(str))

if PhenotypicStochasticity
    str = """ ### whatever PhenotypicStochasticity is when we use RandParasitesGenotype, the  Phenotype == Genotype (easier to program).
	function RandParasitesGenotype()
		Alleles = [if (trait <= N_traits_infection_success) sample(1:NallelesPerTrait) else sample(1:NallelesPerTrait) end for trait in 1:NtraitsPhenOnly
        ] $(ifelse(GfG,".* rand(Binomial(1,0.5),$(NtraitsPhenOnly)","")))
		NoiseAlleles = rand(1:NallelesPerTrait, NtraitsPhenOnly)
        $(ifelse(GfG,"NoiseAlleles[Alleles .=== 0] .= 0",""))
		Genotype = tuple(Alleles...,NoiseAlleles...)
		ParasitesGenotype(
    		GetParasitID(Genotype)
    		,Genotype # Alleles::NTuple{N_traits_infection_success + N_traits_immunity,Int}
    		,tuple(Alleles[1:N_traits_infection_success]...) # traits_infection_success::NTuple{N_traits_infection_success,Int}
    		,tuple(Alleles[(N_traits_infection_success+1):end]...) # traits_immunity::NTuple{N_traits_immunity,Int}
            $(ifelse(GfG,",GetRate_of_parasites_emission(Genotype)","")) # rate_of_parasites_emission
		)
	end"""
	eval(Meta.parse(str))
else
    str = """
        function RandParasitesGenotype()
            Alleles = (sample(1:NallelesPerTrait,NtraitsPhenOnly) $(ifelse(GfG,".* rand(Binomial(1,0.5),$(NtraitsPhenOnly))","")))
    		Genotype = tuple(Alleles...)
            ParasitesGenotype(
                GetParasitID(Genotype)
                ,Genotype # Alleles::Array{N_traits_infection_success + N_traits_immunity,Int}
                ,tuple(Alleles[1:N_traits_infection_success]...) # traits_infection_success::Array{N_traits_infection_success,Int}
                ,tuple(Alleles[(N_traits_infection_success+1):end]...) # traits_immunity::Array{N_traits_immunity,Int}
                $(ifelse(GfG,",GetRate_of_parasites_emission(Genotype)","")) # rate_of_parasites_emission
            )
        end"""
    eval(Meta.parse(str))
end

str = """
struct HostsGenotype
    IDh::Int_IDh
    Alleles::NTuple{Ntraits*2, Int}
    traits_infection_success::NTuple{N_traits_infection_success*2,Int}
    traits_immunity::NTuple{N_traits_immunity*2,Int}
    traits_infection_success_Symb::Symbol
    traits_immunity_Symb::Symbol
    $(ifelse(GfG,"rate_of_gametes_prod::Float",""))
end"""
eval(Meta.parse(str))

if PhenotypicStochasticity
	str = """
	function RandHostsGenotype(sp::Int)
        Alleles      = [$(join([ if trait <= N_traits_infection_success "[sample(AlleleRangePerSp_infection_success[sp]) for _ in 1:2]..." else "[sample(1:NallelesPerTrait, 1) for _ in 1:2]..." end for trait in 1:NtraitsPhenOnly],", "))
        ] $(ifelse(GfG,".* rand(Binomial(1,0.5),$(NtraitsPhenOnly*2))",""))
        for traitPosi in 2:2:$(Ntraits*2)
            if Alleles[traitPosi-1] > Alleles[traitPosi]
                Alleles[(traitPosi-1):traitPosi] .= Alleles[traitPosi:-1:(traitPosi-1)]
            end
        end
		NoiseAlleles =                                                                                                          [$(join(["[sample(1:NallelesPerTrait, 1) for _ in 1:2]..."     for trait in 1:NtraitsPhenOnly],", "))]
        $(ifelse(GfG,"NoiseAlleles[Alleles .=== 0] .= 0","")))
		Genotype = (Alleles..., NoiseAlleles...)
		HostsGenotype(
    		GetHostID(Genotype) # IDh
    		,Genotype # Alleles::NTuple{NtraitsPhenOnly*2, Int}
            ,tuple(             Alleles[1:(2*N_traits_infection_success                    )]...) # traits_infection_success::NTuple{N_traits_infection_success,Int}
    		,tuple(             Alleles[(2*N_traits_infection_success+1):(2*NtraitsPhenOnly)]...) # traits_immunity::NTuple{N_traits_immunity,Int}
            ,Symbol(Symbol.("_",Alleles[1:(2*N_traits_infection_success                    )])...)  # traits_infection_success_Symb
    		,Symbol(Symbol.("_",Alleles[(2*N_traits_infection_success+1):(2*NtraitsPhenOnly)])...)  # traits_immunity_Symb
            $(ifelse(GfG,",GetRate_of_gametes_prod(Genotype)","")) # rate_of_gametes_prod
        )
	end"""
	eval(Meta.parse(str))
else
    str = """
        function RandHostsGenotype(sp::Int)
            Alleles      = [$(join([ if trait <= N_traits_infection_success "[sample(AlleleRangePerSp_infection_success[sp]) for _ in 1:2]..." else "[sample(1:NallelesPerTrait) for _ in 1:2]..." end for trait in 1:NtraitsPhenOnly],", "))
            ] $(ifelse(GfG,".* rand(Binomial(1,0.5),$(NtraitsPhenOnly*2))",""))
            for traitPosi in 2:2:$(Ntraits*2)
                if Alleles[traitPosi-1] > Alleles[traitPosi]
                    Alleles[(traitPosi-1):traitPosi] .= Alleles[traitPosi:-1:(traitPosi-1)]
                end
            end
            Genotype = (Alleles...,)
            HostsGenotype(
                GetHostID(Genotype) # IDh
                ,Genotype # Alleles::Array{NtraitsPhenOnly, Array{2,Int}}
                ,tuple(             Alleles[1:(2*N_traits_infection_success                    )]...) # traits_infection_success::NTuple{N_traits_infection_success,Int}
        		,tuple(             Alleles[(2*N_traits_infection_success+1):(2*NtraitsPhenOnly)]...) # traits_immunity::NTuple{N_traits_immunity,Int}
                ,Symbol(Symbol.("_",Alleles[1:(2*N_traits_infection_success                    )])...)  # traits_infection_success_Symb
        		,Symbol(Symbol.("_",Alleles[(2*N_traits_infection_success+1):(2*NtraitsPhenOnly)])...)  # traits_immunity_Symb
                $(ifelse(GfG,",GetRate_of_gametes_prod(Genotype)","")) # rate_of_gametes_prod
            )
        end"""
    eval(Meta.parse(str))
end

Int() = 0
Float16() = NaN16
Float32() = NaN32
Float64() = NaN64
struct ExtendableDict{K<:Value,V}
    x::Dict{K,V}
end
#
function setindex!(X::ExtendableDict{K,V},k::K,v::V) where K<:Value where V
    X.x[k] = v
end
#
function getindex(X::ExtendableDict{K,V},k::K) where K<:Value where V
    if !haskey(X.x, k)
        X.x[k] =  V()
    end
    return(X.x[k])
end
#
ExtendableDict{Int,V}() where V = ExtendableDict(Dict{Int,V}())

keys(ExtendableDict) = keys(ExtendableDict.x)

str = """
struct HostsGenotypeS_List
    IDhS::Vector{Int_IDh}
    GenotypeS::Vector{HostsGenotype}
    HostSexOffspringMap::$(join(["ExtendableDict{Int," for traitPosi in 1:(Ntraits*2)]))Int$(join(["}" for traitPosi in 1:(Ntraits*2)]))  # POSITION OF GENOTYPES NOT ID
end"""
eval(Meta.parse(str))

struct ParasitesGenotypeS_List
    IDpS::Vector{Int_IDp}
    GenotypeS::Vector{ParasitesGenotype}
end
##################################
struct MetaPop
    PopS::Vector{HostSexS} # one entry for each host population
    IDpopS::Vector{Int}
    IDsp::Int
    N::Ref{Int}
    N_::Vector{Int}
    MigrRate::Ref{Float}
    PmigrSettlementEachPop::Vector{Float}
end


if SimulateImmunisation
    struct Immunisations
        History::NTuple{n, Int} where n
        Map::Dict{Int,Symbol} # For each immunisationsType, if it recovers from parasit genotype which is at position 1, 2, 3 ..., which immunisationsType this gives ?
    end

    str = """
    $INLINE function getindex(X::Dict{Symbol,Immunisations},immun::Symbol,i::Int)
        $INBOUNDstart
        x = X[immun]
        if !haskey(x.Map, i)
            out::Vector{Int} = collect(x.History)
            # push!(out,i)
            # sort!(out)
            at = findfirst(out .> i)
            if isnothing(at)
                push!(out,i)
            else
                insert!(out, at, i)
            end
            outTuple::NTuple{n,Int_IDp} where n = tuple(out...)
            outSymbol::Symbol = Symbol("_",join(outTuple,"°"))
            x.Map[i] = outSymbol
            if !haskey(X, outSymbol)
                X[outSymbol] = Immunisations(outTuple, Dict{Int,Symbol}())
            end
            return(outSymbol)
        else
            return(x.Map[i])
        end
        $INBOUNDend
    end"""
    eval(Meta.parse(str))
end

str = """
struct HPinteractions
    TraitsValues_infection_success::Dict{Symbol, NTuple{ N_traits_immunity*2  , Int}}
    TraitsValues_immunity::Dict{  Symbol, NTuple{ N_traits_infection_success*2, Int}}
    Pinfection1contact::Dict{Symbol, Vector{Float}}        #          for each set of values of traits_infection_success_Symb, what is Pinfection1contact     on each parasite ?
    Pinfection1contactFloat16::Dict{Symbol, Vector{Float}} # Float16( for each set of values of traits_infection_success_Symb, what is Pinfection1contact     on each parasite ? )
    virulenceS::Dict{Symbol, Vector{Float}}                #          for each set of values of traits_infection_success_Symb, what is virulenceS             on each parasite ?
    Precoveryinnateimmu::Dict{Symbol, Vector{Float}}       #          for each set of values of traits_immunity_Symb         , what is P_recovery_innate_immu on each parasite ?
    $(ifelse(SimulateImmunisation
    ,"# for each value of each trait_immunity of the immunising parasite, for each value of each trait_immunity of the infecting parasite, what is Precoveryacquiredimmu ?
      Precoveryacquiredimmu::$(join(["ExtendableDict{Int," for traitPosi in 1:(N_traits_immunity*2)]))Float$(join(["}" for traitPosi in 1:(N_traits_immunity*2)]))
      immunisationS::Dict{Symbol,Immunisations}
      TempDistImmunisingInfecting::Vector{Int}"
    ,""))
end"""
eval(Meta.parse(str))

struct MetaCommunity
    MetaPopS::NTuple{NhostSp,MetaPop} # one entry for each host species
    IDspS::NTuple{NhostSp,Int}
    N::Ref{Int}
    N_::Vector{Int}
    hostsGenotypeS_List::HostsGenotypeS_List
    parasitesGenotypeS_List::ParasitesGenotypeS_List
    #
    HPinteractions::HPinteractions
    ParasitMutatingMap::NTuple{Ntraits,Vector{ExtendableDict{Int,Int}}} # POSITION OF GENOTYPES NOT ID
    Storage::NamedTuple
end


str = #### GetNewHostGenotype!  Here the phenotypic stochastisity don't matter since the it has already played to chose "IDh"
"""
$INLINE function GetNewHostGenotype!(X::MetaCommunity, genotype::NTuple{$(Ntraits*2), Int}, IDh::Int_IDh)
    $FASTMATHstart
    $INBOUNDstart
    $PRINT println("GetNewHostGenotype!")
    $PRINT if !all([genotype[traitPosi] <= genotype[traitPosi+1] for traitPosi in 1:2:$(Ntraits*2)])
    $PRINT      error("In GetNewHostGenotype!, some alleles are not ordered
    $PRINT Genotype to be constructed was \$(genotype)")
    $PRINT end
    $PRINT println(genotype)
    h_traits_infection_success     = tuple($(join(["genotype[$traitPosi]" for traitPosi in 1:(N_traits_infection_success*2)],", ")))
    h_traits_immunity              = tuple($(join(["genotype[$traitPosi]" for traitPosi in (N_traits_infection_success*2+1):(NtraitsPhenOnly*2)],", ")))
    h_traits_infection_success_Symb = Symbol(Symbol.("_",h_traits_infection_success)...)
    h_traits_immunity_Symb          = Symbol(Symbol.("_",h_traits_immunity         )...)
    if !haskey(X.HPinteractions.Pinfection1contact, h_traits_infection_success_Symb)
        X.HPinteractions.TraitsValues_infection_success[h_traits_infection_success_Symb] = h_traits_infection_success
        X.HPinteractions.Pinfection1contact[            h_traits_infection_success_Symb] = [GetPinfection1contact(  GetOverallDist_infection_success( [$(join(["min(CircularDist_infection(X.HPinteractions.TraitsValues_infection_success[h_traits_infection_success_Symb][$i    ], p.traits_infection_success[$(Int(i/2))])
                                                                                                                                                                  , CircularDist_infection(X.HPinteractions.TraitsValues_infection_success[h_traits_infection_success_Symb][$(i-1)], p.traits_infection_success[$(Int(i/2))]))" for i in 2:2:(2*N_traits_infection_success)]," , ")) ] ))       for p in X.parasitesGenotypeS_List.GenotypeS]
        X.HPinteractions.Pinfection1contactFloat16[     h_traits_infection_success_Symb] = Float16.(X.HPinteractions.Pinfection1contact[h_traits_infection_success_Symb])   # Pinfection1contactFloat16
        X.HPinteractions.virulenceS[                    h_traits_infection_success_Symb] = GetVirulence.(X.HPinteractions.Pinfection1contact[h_traits_infection_success_Symb]) # virulenceS
    end
    if !haskey(X.HPinteractions.Precoveryinnateimmu, h_traits_immunity_Symb)
        X.HPinteractions.TraitsValues_immunity[h_traits_immunity_Symb] = h_traits_immunity
        X.HPinteractions.Precoveryinnateimmu[  h_traits_immunity_Symb] = [GetPrecoveryinnateimmu( GetOverallDist_recoveryinnateimmu([ $(join(["min(CircularDist_immunity(X.HPinteractions.TraitsValues_immunity[h_traits_immunity_Symb][$(i-1)], p.traits_immunity[$(Int(i/2))])
                                                                                                                                                 , CircularDist_immunity(X.HPinteractions.TraitsValues_immunity[h_traits_immunity_Symb][$i    ], p.traits_immunity[$(Int(i/2))]))" for i in 2:2:(2*N_traits_immunity)]," , ")) ] ))   for p in X.parasitesGenotypeS_List.GenotypeS] # P_recovery_innate_immu
    end
    push!(X.hostsGenotypeS_List.GenotypeS,
                       HostsGenotype(IDh # IDh
                                    ,genotype # Alleles::NTuple{Ntraits*2,Int}
                                    ,h_traits_infection_success # traits_infection_success::NTuple{N_traits_infection_success,Int}
                                    ,h_traits_immunity          # traits_immunity::NTuple{N_traits_immunity,Int}
                                    ,h_traits_infection_success_Symb
                                    ,h_traits_immunity_Symb
                      $(ifelse(GfG,",GetRate_of_gametes_prod(genotype)","")) # rate_of_gametes_prod
                                    ))
    push!(X.hostsGenotypeS_List.IDhS,IDh)
    X.hostsGenotypeS_List.HostSexOffspringMap$(join(["[genotype[$traitPosi]$addOneIfGfG_str]" for traitPosi in 1:(Ntraits*2-1)])).x[genotype[$(Ntraits*2)]$addOneIfGfG_str] = length(X.hostsGenotypeS_List.IDhS)
    $PRINT println("GetNewHostGenotype! Done")
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

str = """
$INLINE function GetNewParasGenotype!(HPinteractions::HPinteractions, parasitesGenotypeS_List::ParasitesGenotypeS_List, hostsGenotypeS_List::HostsGenotypeS_List, ParasitMutatingMap::NTuple{Ntraits,Vector{ExtendableDict{Int,Int}}}, genotypes::NTuple{$Ntraits,Int}, IDp::Int_IDp)
    $FASTMATHstart
    $INBOUNDstart
    $PRINT println("GetNewParasGenotype!")
    p = ParasitesGenotype(
		IDp
		,genotypes # Alleles::NTuple{N_traits_infection_success + N_traits_immunity,Int}
		,tuple(genotypes[1:$N_traits_infection_success]...) # traits_infection_success::NTuple{N_traits_infection_success,Int}
		,tuple(genotypes[$(N_traits_infection_success+1):$(N_traits_infection_success+N_traits_immunity)]...) # traits_immunity::NTuple{N_traits_immunity,Int}
        $(ifelse(GfG,",GetRate_of_parasites_emission(genotypes)","")) # rate_of_parasites_emission
		)
	push!(parasitesGenotypeS_List.IDpS     , IDp)
    $PRINT println("p = ",p)
	push!(parasitesGenotypeS_List.GenotypeS, p )
	# For each trait, create the new genotype in ParasitMutatingMap
	$SIMD for trait in 1:$Ntraits
		push!(ParasitMutatingMap[trait],ExtendableDict{Int,Int}())
	end
    # For each Host traits_infection_success value, update   Pinfection1contact, Pinfection1contactFloat16, virulenceS
    for k in keys(HPinteractions.Pinfection1contact)
        $PRINT println("Pinfection1contact   k = ",k)
        $PRINT println("HPinteractions.TraitsValues_infection_success[k] = ", HPinteractions.TraitsValues_infection_success[k])
        $PRINT println("p.traits_infection_success = ", p.traits_infection_success)
        $PRINT println("GetOverallDist_infection_success = ",
        $PRINT [ $(join(["min(CircularDist_infection(HPinteractions.TraitsValues_infection_success[k][$(i-1)], p.traits_infection_success[$(Int(i/2))])
        $PRINT              , CircularDist_infection(HPinteractions.TraitsValues_infection_success[k][$i    ], p.traits_infection_success[$(Int(i/2))]))" for i in 2:2:(2*N_traits_infection_success)]," , ")) ] )
        push!(HPinteractions.Pinfection1contact[k]         , GetPinfection1contact( GetOverallDist_infection_success(  [ $(join(["min(CircularDist_infection(HPinteractions.TraitsValues_infection_success[k][$(i-1)], p.traits_infection_success[$(Int(i/2))])
                                                                                                                                    , CircularDist_infection(HPinteractions.TraitsValues_infection_success[k][$i    ], p.traits_infection_success[$(Int(i/2))]))" for i in 2:2:(2*N_traits_infection_success)]," , ")) ] )))
        $PRINT println("HPinteractions.Pinfection1contact[k][end] = ", HPinteractions.Pinfection1contact[k][end])
        push!(HPinteractions.Pinfection1contactFloat16[k]  , Float16(HPinteractions.Pinfection1contact[k][end]))
        $PRINT println("HPinteractions.Pinfection1contactFloat16[k][end] = ", HPinteractions.Pinfection1contactFloat16[k][end])
        push!(HPinteractions.virulenceS[k]                 , GetVirulence(HPinteractions.Pinfection1contact[k][end]) )
        $PRINT println("HPinteractions.virulenceS[k][end] = ", HPinteractions.virulenceS[k][end])
    end
    # For each Host traits_immunity          value, update   Precoveryinnateimmu
    for k in keys(HPinteractions.Precoveryinnateimmu)
        $PRINT println("Precoveryinnateimmu    k = ",k)
        $PRINT println("HPinteractions.TraitsValues_immunity[k] = ", HPinteractions.TraitsValues_immunity[k])
        $PRINT println("p.traits_immunity = ", p.traits_immunity)
        push!(HPinteractions.Precoveryinnateimmu[k]        , GetPrecoveryinnateimmu( GetOverallDist_recoveryinnateimmu([ $(join(["min(CircularDist_immunity( HPinteractions.TraitsValues_immunity[k][$(i-1)], p.traits_immunity[$(Int(i/2))])
                                                                                                                                    , CircularDist_immunity( HPinteractions.TraitsValues_immunity[k][    $i], p.traits_immunity[$(Int(i/2))]))" for i in 2:2:(2*N_traits_immunity)]," , ")) ] )))
    end
    $PRINT println("GetNewParasGenotype! Done")
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))
#
#
# [ min(CircularDist_infection(metaCommunity.HPinteractions.TraitsValues_infection_success[k][1], p.traits_infection_success[1])
#     , CircularDist_infection(metaCommunity.HPinteractions.TraitsValues_infection_success[k][2], p.traits_infection_success[1])) , min(CircularDist_infection(metaCommunity.HPinteractions.TraitsValues_infection_success[k][3], p.traits_infection_success[2])
#     , CircularDist_infection(metaCommunity.HPinteractions.TraitsValues_infection_success[k][4], p.traits_infection_success[2])) , min(CircularDist_infection(metaCommunity.HPinteractions.TraitsValues_infection_success[k][5], p.traits_infection_success[3])
#     , CircularDist_infection(metaCommunity.HPinteractions.TraitsValues_infection_success[k][6], p.traits_infection_success[3])) ]


## Cost of having or not having some traits for hosts and parasites

# Equivalent R funtions :
# fN = function(x,y_Min,y_Max,x_Min,x_Max,a){
#   if(a>0){ y_Max-((y_Max-y_Min)*(1-(x_Max-x)/(x_Max-x_Min )))/(1+a/(x_Max-x_Min )*((x_Max-x)))
#     }else{ y_Min+((y_Max-y_Min)*(1-(x-x_Min)/(x_Max-x_Min )))/(1-a/(x_Max-x_Min )*((x-x_Min)))
# }}
#
# fP = function(x,y_Min,y_Max,x_Min,x_Max,a){ # Positive relationship
#   if(a>0){ y_Max-((y_Max-y_Min)*(1-(x-x_Min)/(x_Max-x_Min )))/(1+a/(x_Max-x_Min )*((x-x_Min)))
#     }else{ y_Min+((y_Max-y_Min)*(1-(x_Max-x)/(x_Max-x_Min )))/(1-a/(x_Max-x_Min )*((x_Max-x)))
# }}


if GfG
    if MaxHostCostOfMissing_traits_infection_success > Maxrate_of_gametes_prod error("MaxHostCostOfMissing_traits_infection_success > Maxrate_of_gametes_prod") end
    if MaxHostCostOfHaving_traits_immunity           > Maxrate_of_gametes_prod error("MaxHostCostOfHaving_traits_immunity > Maxrate_of_gametes_prod") end

    y_Min_inf = MaxHostCostOfMissing_traits_infection_success / sqrt(Maxrate_of_gametes_prod)
    y_Min_imm = MaxHostCostOfHaving_traits_immunity           / sqrt(Maxrate_of_gametes_prod)

    str = """
    $INLINE function GetRate_of_gametes_prod(genotype::NTuple{Ntraits*2,Int})
        $FASTMATHstart
        $INBOUNDstart
        NtraitsInfectionSuccessPresent = $(join(["(genotype[$traitPosi] !== 0)" for traitPosi in 1:(N_traits_infection_success*2)                      ]," + "))
        NtraitsImmunityPresent         = $(join(["(genotype[$traitPosi] !== 0)" for traitPosi in ((N_traits_infection_success*2)+1):(NtraitsPhenOnly*2)]," + "))
        #        y_Min+((y_Max-y_Min)*(1-(x-x_Min)/(x_Max-x_Min )))/(1+a*((x-x_Min)/(x_Max-x_Min )))
        #              ((y_Max-y_Min)*(1-(x_Max-x)/(x_Max-x_Min )))/(1-a*((x_Max-x)/(x_Max-x_Min )))
        #        y_Max-((y_Max-y_Min)*(1-(x_Max-x)/(x_Max-x_Min )))/(1-a*((x_Max-x)/(x_Max-x_Min )))
        # y_Max-(y_Min+((y_Max-y_Min)*(1-(x-x_Min)/(x_Max-x_Min )))/(1+a*((x-x_Min)/(x_Max-x_Min ))))

        return (
            (   $(ifelse(TolParasCostOfMissing_traits_immunity >= 0
                    # With x_Min = 0 (no traits)
                    #                y_Max         -((            y_Max           -  y_Min   )*(1-           (x-x_Min)          /       (x_Max-x_Min)           ))/(1+                    a                  /       (x_Max-x_Min )          *(           (x-x_Min)          ))
                    ,"$(sqrt(Maxrate_of_gametes_prod))-(($(sqrt(Maxrate_of_gametes_prod)-y_Min_inf))*(1-NtraitsInfectionSuccessPresent/$(N_traits_infection_success*2)))/(1+$(TolParasCostOfMissing_traits_immunity/(N_traits_infection_success*2))*(NtraitsInfectionSuccessPresent))"
                    #     y_Min +((              y_Max          -   y_Min  )*(1-(            x_Max          -             x                )/        (x_Max-x_Min)          ))/(1-                       a                       /      (x_Max-x_Min)            *((          x_Max                -              x               )))
                    ,"($y_Min_inf+(($(sqrt(Maxrate_of_gametes_prod)-y_Min_inf))*(1-($N_traits_infection_success-NtraitsInfectionSuccessPresent)/$(N_traits_infection_success*2)))/(1-$(TolHostCostOfMissing_traits_infection_success/(N_traits_infection_success*2))*(($(N_traits_infection_success*2)-NtraitsInfectionSuccessPresent)))"
                ))
            ) * (
                $(ifelse(TolHostCostOfHaving_traits_immunity >= 0
                    #    y_Min  +((           y_Max            -   y_Min  )*(1-    (x-x_Min)         /    (x_Max-x_Min )    ))/(1+                   a                 /  (x_Max-x_Min )  *(     (x-x_Min)        ))
                    ,"$y_Min_imm+(($(sqrt(Maxrate_of_gametes_prod)-y_Min_imm))*(1-NtraitsImmunityPresent/$(N_traits_immunity*2)))/(1+$(TolHostCostOfHaving_traits_immunity/N_traits_immunity)*(NtraitsImmunityPresent))"
                    #                y_Max         -((            y_Max           -   y_Min  )*(1-(          x_Max       -          x           )/      (x_Max-x_Min )  ))/(1-                   a                 /  (x_Max-x_Min )      *((         x_Max        -          x           )))
                    ,"$(sqrt(Maxrate_of_gametes_prod))-(($(sqrt(Maxrate_of_gametes_prod)-y_Min_imm))*(1-($(N_traits_immunity*2)-NtraitsImmunityPresent)/$(N_traits_immunity*2)))/(1-$(TolHostCostOfHaving_traits_immunity/(N_traits_immunity*2))*(($(N_traits_immunity*2)-NtraitsImmunityPresent)))"
                ))
            )
        )
        $INBOUNDend
        $FASTMATHend
    end"""
    eval(Meta.parse(str))
end


if GfG
    if MaxParasCostOfMissing_traits_immunity         > MaxRate_of_parasites_emission error("MaxParasCostOfMissing_traits_immunity > MaxRate_of_parasites_emission") end
    if MaxParasCostOfHaving_traits_infection_success > MaxRate_of_parasites_emission error("MaxParasCostOfHaving_traits_infection_success > MaxRate_of_parasites_emission") end

    y_Min_imm = MaxParasCostOfMissing_traits_immunity         / sqrt(MaxRate_of_parasites_emission)
    y_Min_inf = MaxParasCostOfHaving_traits_infection_success / sqrt(MaxRate_of_parasites_emission)

    str = """
    $INLINE function GetRate_of_parasites_emission(genotype::NTuple{Ntraits,Int})
        $FASTMATHstart
        $INBOUNDstart
        NtraitsInfectionSuccessPresent = $(join(["(genotype[$traitPosi] !== 0)" for traitPosi in 1:N_traits_infection_success                  ]," + "))
        NtraitsImmunityPresent         = $(join(["(genotype[$traitPosi] !== 0)" for traitPosi in (N_traits_infection_success+1):NtraitsPhenOnly]," + "))
        return ( # positive relationship with NtraitsImmunityPresent
            (   $(ifelse(TolParasCostOfMissing_traits_immunity >= 0
                    # With x_Min = 0 (no traits)
                    #                y_Max              -((                y_Max            -  y_Min   )*(1-      (x-x_Min)       /  (x_Max-x_Min)  ))/(1+                   a                   /  (x_Max-x_Min )  *(     (x-x_Min)        ))
                    ,"$(sqrt(MaxRate_of_parasites_emission))-(($(sqrt(MaxRate_of_parasites_emission)-y_Min_imm))*(1-NtraitsImmunityPresent/N_traits_immunity))/(1+$(TolParasCostOfMissing_traits_immunity/N_traits_immunity)*(NtraitsImmunityPresent))"
                    #    y_Min  +((                y_Max            -   y_Min  )*(1-(      x_Max       -             x        )/   (x_Max-x_Min)  ))/(1-                     a                 /  (x_Max-x_Min)   *(      x_Max       -         x            ))
                    ,"$y_Min_imm+(($(sqrt(MaxRate_of_parasites_emission)-y_Min_imm))*(1-($N_traits_immunity-NtraitsImmunityPresent)/$N_traits_immunity))/(1-$(TolParasCostOfMissing_traits_immunity/N_traits_immunity)*($N_traits_immunity-NtraitsImmunityPresent))"
                ))
            ) * ( # negative relationship with NtraitsInfectionSuccessPresent
                $(ifelse(TolParasCostOfHaving_traits_infection_success >= 0
                    #    y_Min  +((              y_Max              -  y_Min   )*(1-    (x-x_Min)                 /     (x_Max-x_Min )       ))/(1+                    a                          /      (x_Max-x_Min )       *(         (x-x_Min)            ))
                    ,"$y_Min_inf+(($(sqrt(MaxRate_of_parasites_emission)-y_Min_inf))*(1-NtraitsInfectionSuccessPresent/N_traits_infection_success))/(1+$(TolParasCostOfHaving_traits_infection_success/N_traits_infection_success)*(NtraitsInfectionSuccessPresent))"
                    #                 y_Max             -((             y_Max               -   y_Min  )*(1-(           x_Max           -               x              )/     ( x_Max-x_Min )       ))/(1-                        a                      /  (x_Max-x_Min )           *(         x_Max             -               x              ))
                    ,"$(sqrt(MaxRate_of_parasites_emission))-(($(sqrt(MaxRate_of_parasites_emission)-y_Min_inf))*(1-($N_traits_infection_success-NtraitsInfectionSuccessPresent)/$N_traits_infection_success))/(1-$(TolParasCostOfHaving_traits_infection_success/N_traits_infection_success)*($N_traits_infection_success-NtraitsInfectionSuccessPresent))"
                ))
            )
        )
        $INBOUNDend
        $FASTMATHend
    end"""
    eval(Meta.parse(str))
end

### SetHostParasitInteractionParameters!
str = """
$INLINE function GetPinfection1contact(Dist::Union{Float,Int})
    $FASTMATHstart
        if (iszero(Dist)) $Max_Pinfection1contact else
            $(ifelse(Tol_infection >= 0
            # y_Max                  -(       (y_Max-y_Min)*(1-(x_Max-x)/(x_Max-x_Min )))/(1+a/(x_Max-x_Min )*((x_Max-x)))
            ,"$Max_Pinfection1contact-($(Max_Pinfection1contact/MaxDist_infection)*Dist)/(1+$(Tol_infection/MaxDist_infection)*($MaxDist_infection-Dist))"
            # y_Min+(    (y_Max-y_Min)      *(1-(x-x_Min)/(x_Max-x_Min )    ))/(1   -a            /    (x_Max-x_Min )*((x-x_Min)))
            ,"      ($Max_Pinfection1contact*(1-Dist     /$MaxDist_infection))/(1+$(-Tol_infection/MaxDist_infection)*Dist)"
            ))
        end
    $FASTMATHend
    end"""
eval(Meta.parse(str))

str = """
$INLINE function GetVirulence(Pinfection1contact::Float)
    $FASTMATHstart
        if (iszero(Pinfection1contact)) $Min_virul else
            $(ifelse(Tol_virul >= 0
            # y_Max     -((y_Max-y_Min)         *(1-(x-x_Min)                  /(x_Max-x_Min )         ))/(1+a          /(x_Max-x_Min )         *(x-x_Min))
            ,"$Max_virul-($(Max_virul-Min_virul)×(1-Pinfection1contact         /$Max_Pinfection1contact))/(1+$(Tol_virul/Max_Pinfection1contact)× Pinfection1contact)"
            # y_Min     +((y_Max-y_Min)         *(1-(x_Max -x                 )/(x_Max-x_Min )         ))/(1   -a        /(x_Max-x_Min )         *(x_Max                  -x))
            ,"$Min_virul+($(Max_virul-Min_virul)×(          Pinfection1contact /$Max_Pinfection1contact))/(1+$(-Tol_virul/Max_Pinfection1contact)×($Max_Pinfection1contact-Pinfection1contact))"
            ))
        end
    $FASTMATHend
    end"""
eval(Meta.parse(str))

str = """
$INLINE function GetPrecoveryinnateimmu(Dist::Union{Float,Int})
    $FASTMATHstart
    if (iszero(Dist)) $Max_Precovery_innate_immu else
        $(ifelse(Tol_InnateImmu >= 0
        ,"$Max_Precovery_innate_immu*(1-(Dist/$MaxDist_immunity)/(1+$Tol_InnateImmu×(1-Dist/$MaxDist_immunity) ))"
        ,"$Max_Precovery_innate_immu*($MaxDist_immunity/Dist-1)/($MaxDist_immunity/Dist - $Tol_InnateImmu)"
        ))
    end
    $FASTMATHend
end"""
eval(Meta.parse(str))

if SimulateImmunisation
    str = """
    $INLINE function GetPrecoveryacquiredimmu_fromMeanDist(Dist::Float)
        $FASTMATHstart
        $INBOUNDstart
        # fN = function(x,y_Min,y_Max,x_Min,x_Max,a){ # Negative relationship
        #   if(a>0){ y_Min+((y_Max-y_Min)*(1-(x-x_Min)/(x_Max-x_Min )))/(1+a/(x_Max-x_Min )*((x-x_Min)))
        #     }else{ y_Max-((y_Max-y_Min)*(1-(x_Max-x)/(x_Max-x_Min )))/(1-a/(x_Max-x_Min )*((x_Max-x)))
        # }}
        return (
            $(ifelse(Tol_CrossImmun >= 0
            # y_Max                        -((y_Max-y_Min)                *(1-(x_Max-x)               / (x_Max-x_Min )  ))/(1+  a             /(x_Max-x_Min )   *(( x_Max           -x)))
            ,"$Max_Precovery_acquiered_immu-($Max_Precovery_acquiered_immu*(1-($MaxDist_immunity-Dist)/$MaxDist_immunity))/(1+$(Tol_CrossImmun/MaxDist_immunity)*(($MaxDist_immunity-Dist)))"
            # y_Min+((y_Max-y_Min)                *(1-(x-x_Min)/(  x_Max-x_Min )  ))/(1   -a             /  (x_Max-x_Min ) *(x-x_Min))
            ,"      ($Max_Precovery_acquiered_immu*(1- Dist    / $MaxDist_immunity))/(1+$(-Tol_CrossImmun/MaxDist_immunity)*Dist     )"
            ))
        )
        $INBOUNDend
        $FASTMATHend
    end"""
    eval(Meta.parse(str))

    str = """
      $INLINE function GetPrecoveryacquiredimmu!(Precoveryacquiredimmu::$(join(["ExtendableDict{Int," for traitPosi in 1:(N_traits_immunity*2)]))Float$(join(["}" for traitPosi in 1:(N_traits_immunity*2)]))
                                              ,posiImmunisingParas::Int ,  posiInfectingParas::Int
                                              ,parasitesGenotypeS_List::ParasitesGenotypeS_List
                                              ,TempDistImmunisingInfecting::Vector{Int})
        $FASTMATHstart
        $INBOUNDstart
        if posiImmunisingParas === 0
            return(  0.0  )
        else
            ImmunisingParas_traits_immunity::NTuple{$N_traits_immunity,Int} = parasitesGenotypeS_List.GenotypeS[posiImmunisingParas].traits_immunity
            InfectingParas_traits_immunity::NTuple{$N_traits_immunity,Int}  = parasitesGenotypeS_List.GenotypeS[posiInfectingParas ].traits_immunity                                                                                                                     # -1 -->> !! => get the last ExtendableDict
            Precoveryacquiredimmu_ = Precoveryacquiredimmu$(  join(["[ImmunisingParas_traits_immunity[$traitPosi]$addOneIfGfG_str]" for traitPosi in 1:N_traits_immunity])    )$(    join(["[InfectingParas_traits_immunity[ $traitPosi]$addOneIfGfG_str]" for traitPosi in 1:(N_traits_immunity-1)])  )
            posiLastAllele = InfectingParas_traits_immunity[end]$addOneIfGfG_str
            if isequal(Precoveryacquiredimmu_[posiLastAllele] , NaN)
                map!(i -> CircularDist_immunity(ImmunisingParas_traits_immunity[i], InfectingParas_traits_immunity[i]), TempDistImmunisingInfecting, 1:$N_traits_immunity)
                Precoveryacquiredimmu_.x[posiLastAllele] = GetPrecoveryacquiredimmu_fromMeanDist(  mean(TempDistImmunisingInfecting)  )
            end
            return(  Precoveryacquiredimmu_.x[posiLastAllele]  )
        end
        $INBOUNDend
        $FASTMATHend
        end"""
    eval(Meta.parse(str))
end

str = """
function SetHostParasitInteractionParameters!(X::MetaCommunity)
    $FASTMATHstart
    $INBOUNDstart
    for k in keys(X.HPinteractions.Pinfection1contact)
        empty!(X.HPinteractions.Pinfection1contact[k])
        empty!(X.HPinteractions.Pinfection1contactFloat16[k])
        empty!(X.HPinteractions.virulenceS[k])
    end
    for k in keys(X.HPinteractions.Precoveryinnateimmu)
        empty!(X.HPinteractions.Precoveryinnateimmu[k])
    end
    for p in X.parasitesGenotypeS_List.GenotypeS
        for k in (g.traits_infection_success_Symb for g in X.hostsGenotypeS_List.GenotypeS)
            push!(X.HPinteractions.Pinfection1contact[k]         , GetPinfection1contact( GetOverallDist_infection_success(  [ $(join(["min(CircularDist_infection(X.HPinteractions.TraitsValues_infection_success[k][$(i-1)], p.traits_infection_success[$(Int(i/2))])
                                                                                                                                          , CircularDist_infection(X.HPinteractions.TraitsValues_infection_success[k][$i    ], p.traits_infection_success[$(Int(i/2))]))" for i in 2:2:(2*N_traits_infection_success)]," , ")) ] )))
            push!(X.HPinteractions.Pinfection1contactFloat16[k]  , Float16(X.HPinteractions.Pinfection1contact[k][end]))
            push!(X.HPinteractions.virulenceS[k]                 , GetVirulence(X.HPinteractions.Pinfection1contact[k][end]) )
        end
        # For each Host traits_immunity          value, update   Precoveryinnateimmu
        for k in (g.traits_immunity_Symb for g in X.hostsGenotypeS_List.GenotypeS)
            push!(X.HPinteractions.Precoveryinnateimmu[k]        , GetPrecoveryinnateimmu( GetOverallDist_recoveryinnateimmu([ $(join(["min(CircularDist_immunity( X.HPinteractions.TraitsValues_immunity[k][$(i-1)], p.traits_immunity[$(Int(i/2))])
                                                                                                                                          , CircularDist_immunity( X.HPinteractions.TraitsValues_immunity[k][    $i], p.traits_immunity[$(Int(i/2))]))" for i in 2:2:(2*N_traits_immunity)]," , ")) ] )))
        end
    end
    $(ifelse(SimulateImmunisation  # Symbol("_",join([0],","))
    ,"X.HPinteractions.immunisationS[            :_0          ] = Immunisations((0,), Dict{Int,Symbol}())"
    ,""))
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

#### HostParasitesPhenotypes_INI! fill HostsGenotypeS accounting for genotype stochastisity
# for SP in 1:NhostSp
    # function HostParasitesPhenotypes_INI!(MetaComm::MetaCommunity, hostsGenotypeS::HostsGenotypeS, hostGenotype::NTuple{$(Ntraits*2), Int}, parasitGenotype::NTuple{$Ntraits,Int}, N::Int, sp::Val{$SP})
    str = """
        function HostParasitesPhenotypes_INI!(MetaComm::MetaCommunity, hostsGenotypeS::HostsGenotypeS, hostGenotype::NTuple{$(Ntraits*2), Int}, parasitGenotype::NTuple{$Ntraits,Int}, N::Int, hSp::Int )
            $INBOUNDstart
            $PRINT println("HostParasitesPhenotypes_INI! with parasitGenotype")
            HostID = GetHostID(      hostGenotype)
            ParaID = GetParasitID(parasitGenotype)
            # Check the existence of the genotypes
            if !any(MetaComm.hostsGenotypeS_List.IDhS .== HostID)
                GetNewHostGenotype!(MetaComm, hostGenotype, HostID)
            end
            if !any(MetaComm.parasitesGenotypeS_List.IDpS .== ParaID)
                $PRINT println("parasitGenotype = ", parasitGenotype)
                $PRINT println("ParaID = ", ParaID)
                $PRINT println("HostParasitesPhenotypes_INI! : GetNewParasGenotype!(Genotype)")
                GetNewParasGenotype!(MetaComm.HPinteractions, MetaComm.parasitesGenotypeS_List, MetaComm.hostsGenotypeS_List ,MetaComm.ParasitMutatingMap, parasitGenotype, ParaID)
                $PRINT println("HostParasitesPhenotypes_INI! : GetNewParasGenotype!(Genotype) Done")
            end
            if PhenotypicStochasticity
        		hostNoiseAlleles    = hostGenotype[  $(NtraitsPhenOnly*2+1):$(Ntraits*2)]
        		parasitNoiseAlleles = parasitGenotype[$(NtraitsPhenOnly+1):$Ntraits]
                #                 HostTrait1chr1,HostTrait1chr2,...,ParasTrait1,ParasTrait2,...
        		N_::Vector{Int} , PhenS::Vector{Vector{Int}}                                    = RecursiveRandMultinom(N, vcat(SymetricalProbilitiesGeomDistr_Host[hostNoiseAlleles],SymetricalProbilitiesGeomDistr_Parasit[parasitNoiseAlleles]))
    		else
                N_ , PhenS = [N], [[$(join(["$(NallelesPerTrait),$(NallelesPerTrait)" for trait in 1:NtraitsPhenOnly],", ")*"   ,     "*join(["$NallelesPerTrait" for trait in 1:NtraitsPhenOnly],", ") )]]
    		end
            RefPhen = vcat(hostGenotype[1:$(NtraitsPhenOnly*2)], parasitGenotype[1:$NtraitsPhenOnly])
            RefPhen .-= $NallelesPerTrait
            #
            RefPhen[$(NtraitsPhenOnly*2+1):end] .-= $(NallelesPerTrait) #  .- NallelesPerTrait{in parasites}
            HostsIDphenS = Ref{Int_IDh}()
            ParasIDphenS = Ref{Int_IDp}()
            for i in eachindex(PhenS)
                PhenS[i] .+= RefPhen # .- Npossible alleles {in hosts} || .- NallelesPerTrait{in parasites}
                # Host
                for chr in 0:1
                    for trait in 1:$NtraitsPhenOnly
                        traitChr = trait*2-chr
                        if PhenS[i][traitChr] < 1
                            PhenS[i][traitChr] += $(NallelesPerTrait_EachTrait)[trait] # circular trait /!| we are conscidering the phenotype !
                        elseif PhenS[i][traitChr] > $(NallelesPerTrait_EachTrait)[trait]
                            PhenS[i][traitChr] -= $(NallelesPerTrait_EachTrait)[trait] # circular trait /!| we are conscidering the phenotype !
                        end
                    end
                end
                for trait in 1:2:$(NtraitsPhenOnly*2) # order chromosoms
                    if PhenS[i][trait] > PhenS[i][trait+1]
                        PhenS[i][trait:(trait+1)] = PhenS[i][(trait+1):-1:trait]
                    end
                end
                # Parasite
                for traitPosiParas in $(NtraitsPhenOnly*2+1):$(NtraitsPhenOnly*2+NtraitsPhenOnly)
                    if PhenS[i][traitPosiParas] < 1
                        PhenS[i][traitPosiParas] += $(NallelesPerTrait_EachTrait)[traitPosiParas-$(NtraitsPhenOnly*2)]
                    elseif PhenS[i][traitPosiParas] > $(NallelesPerTrait_EachTrait)[traitPosiParas-$(NtraitsPhenOnly*2)]
                        PhenS[i][traitPosiParas] -= $(NallelesPerTrait_EachTrait)[traitPosiParas-$(NtraitsPhenOnly*2)]
                    end
                end
                #
                $(ifelse(PhenotypicStochasticity
                ,"PhenS[i] = vcat(PhenS[i][1:$(NtraitsPhenOnly*2)],  hostNoiseAlleles,  PhenS[i][$(NtraitsPhenOnly*2+1):end],  parasitNoiseAlleles)"
                ,"PhenS[i] = vcat(PhenS[i][1:$(NtraitsPhenOnly*2)],                     PhenS[i][$(NtraitsPhenOnly*2+1):end]                      )"
                ))
                $PRINT println(); println(); println(); println("hostGenotype = ", hostGenotype  ,"
                $PRINT  ;  parasitGenotype=",parasitGenotype,"
                $PRINT  ;  PhenS=",PhenS); println(); println(); println()
                HostsIDphenS[] = GetHostID(   PhenS[i][1:$(Ntraits*2)])
                ParasIDphenS[] = GetParasitID(PhenS[i][$(Ntraits*2+1):end])
                # Is this combinaison of host and parasit genotype and phenotype already existing in the local population ?
                IDcategory = GetCategoryID(HostID, HostsIDphenS[], GetPgenXphenID_HostShiftHistory(ParaID, ParasIDphenS[], Symbol(:IDp°,  ParaID,  :__Sp°,  hSp ,  :°°Time°,  0.0,  :°°IDh°,  HostsIDphenS[])))
                posiLocal = findall(hostsGenotypeS.IDcategorieS .== IDcategory)
                $PRINT if length(posiLocal)>1 error("the choice between categories that matches host and parasit genotype and phenotype has not yet been implemented") end
                if length(posiLocal)==0 # it does not exist
                    push!(hostsGenotypeS.IDhS                        ,HostID)
                    push!(hostsGenotypeS.IDhSphen                    ,HostsIDphenS[])
                    push!(hostsGenotypeS.IDpS                        ,ParaID)
                    push!(hostsGenotypeS.IDpSphen                    ,ParasIDphenS[])
                    push!(hostsGenotypeS.IDpGenPhen_HostShiftHistoryS,GetPgenXphenID_HostShiftHistory(ParaID, ParasIDphenS[], Symbol(:IDp°,  ParaID,  :__Sp°,  hSp ,  :°°Time°,  0.0,  :°°IDh°,  HostsIDphenS[])))
                    push!(hostsGenotypeS.IDp_HostShiftHistoryS       ,                                                        Symbol(:IDp°,  ParaID,  :__Sp°,  hSp ,  :°°Time°,  0.0,  :°°IDh°,  HostsIDphenS[]))
                    $(ifelse(SimulateImmunisation
                    ,"push!(hostsGenotypeS.IDiS    ,:_0)"
                    ,""    ))
                    push!(hostsGenotypeS.IDcategorieS, IDcategory)
                    $(ifelse(GfG
                    ,"push!(hostsGenotypeS.rate_of_gametes_prodS      ,GetRate_of_gametes_prod(      tuple(PhenS[i][1:$(Ntraits*2)]...)))
                      push!(hostsGenotypeS.rate_of_parasites_emissionS,GetRate_of_parasites_emission(tuple(PhenS[i][$(Ntraits*2+1):end]...)))"
                    ,""))
                    #
                    h_traits_infection_success = tuple($(join(["PhenS[i][$traitPosi]" for traitPosi in 1:(N_traits_infection_success*2)],", ")))
                    p_traits_infection_success = PhenS[i][$(Ntraits*2+1):$(Ntraits*2+N_traits_infection_success)]
                    push!(hostsGenotypeS.virulenceS, GetVirulence(
                                                                  GetPinfection1contact(  GetOverallDist_infection_success([
                                                                        $(join(["min(CircularDist_infection(h_traits_infection_success[$(i-1)], p_traits_infection_success[$(Int(i/2))]) , CircularDist_infection(h_traits_infection_success[$i], p_traits_infection_success[$(Int(i/2))]))" for i in 2:2:(2*N_traits_immunity)]," , ")) ])
                                                                )))
                    #
                    h_traits_immunity = tuple($(join(["PhenS[i][$traitPosi]" for traitPosi in (N_traits_infection_success*2+1):(NtraitsPhenOnly*2)],", ")))
                    p_traits_immunity = PhenS[i][$(Ntraits*2+N_traits_infection_success+1):$(Ntraits*2+N_traits_infection_success+N_traits_immunity)]

                    push!(hostsGenotypeS.Precoveryinnateimmu, GetPrecoveryinnateimmu( GetOverallDist_recoveryinnateimmu([  $(join(["min(CircularDist_immunity(h_traits_immunity[$(i-1)], p_traits_immunity[$(Int(i/2))]) , CircularDist_immunity(h_traits_immunity[$i], p_traits_immunity[$(Int(i/2))]))" for i in 2:2:(2*N_traits_immunity)]," , ")) ])  ))
                    $(ifelse(SimulateImmunisation
                    ,"push!(hostsGenotypeS.Precoveryacquiredimmu, 0.0  )" # GetPrecoveryacquiredimmu!(MetaComm.HPinteractions.Precoveryacquiredimmu, posiParasGeno, 0, MetaComm.parasitesGenotypeS_List, MetaComm.HPinteractions.TempDistImmunisingInfecting)
                    ,""))
                    #
                    push!(hostsGenotypeS.N_      ,N_[i])
                    push!(hostsGenotypeS.dN_     ,0    )
                    hostsGenotypeS.N[] += N_[i]
                    #
                    # Is this host *phenotype* already existing in the MetaCommunity ?
                    posiHostPhen = findfirst(xx -> xx === HostsIDphenS[], MetaComm.hostsGenotypeS_List.IDhS)
                    if isnothing(posiHostPhen) # it does not exist
                        GetNewHostGenotype!(MetaComm, tuple(PhenS[i][1:$(Ntraits*2)]...), HostsIDphenS[])
                        push!(hostsGenotypeS.posiInHostsGenotypeS_List,length(MetaComm.hostsGenotypeS_List.IDhS))
                    else
                        push!(hostsGenotypeS.posiInHostsGenotypeS_List,posiHostPhen)
                    end
                    # Is this parasit *phenotype* already existing in the MetaCommunity ?
                    posiParasPhen = findfirst(xx -> xx === ParasIDphenS[], MetaComm.parasitesGenotypeS_List.IDpS)
                    if isnothing(posiParasPhen) # it does not exist
                        $PRINT println("HostParasitesPhenotypes_INI! : GetNewParasGenotype!(Phenotype)")
                        GetNewParasGenotype!(MetaComm.HPinteractions, MetaComm.parasitesGenotypeS_List, MetaComm.hostsGenotypeS_List, MetaComm.ParasitMutatingMap, tuple(PhenS[i][$(Ntraits*2+1):$(Ntraits*3)]...), ParasIDphenS[])
                        $PRINT println("HostParasitesPhenotypes_INI! : GetNewParasGenotype!(Phenotype) Done")
                        push!(hostsGenotypeS.posiInParasitGenotypeS_List,length(MetaComm.parasitesGenotypeS_List.IDpS))
                    else
                        push!(hostsGenotypeS.posiInParasitGenotypeS_List,posiParasPhen)
                    end
                else
                    hostsGenotypeS.N_[posiLocal[1]] += N_[i]
                    hostsGenotypeS.N[] += N_[i]
                end
            end
        $INBOUNDend
        end"""
    eval(Meta.parse(str))
# end


# for SP in 1:NhostSp # hostGenotype::Nothing => uninfected hosts
      # function HostParasitesPhenotypes_INI!(MetaComm::MetaCommunity, hostsGenotypeS::HostsGenotypeS, hostGenotype::NTuple{$(Ntraits*2), Int}, parasitGenotype::Nothing, N::Int, sp::Val{$SP})
      str = """
        function HostParasitesPhenotypes_INI!(MetaComm::MetaCommunity, hostsGenotypeS::HostsGenotypeS, hostGenotype::NTuple{$(Ntraits*2), Int}, parasitGenotype::Nothing, N::Int, hSp::Int              )
            $INBOUNDstart
            HostID = GetHostID(   hostGenotype )
            # Check the existence of the genotypes
            if !any(MetaComm.hostsGenotypeS_List.IDhS .== HostID)
                GetNewHostGenotype!(MetaComm, hostGenotype, HostID)
            end
            $PRINT println("HostParasitesPhenotypes_INI! : GetNewHostGenotype!(MetaComm, hostGenotype, HostID) DONE")
            if PhenotypicStochasticity
                hostNoiseAlleles    = hostGenotype[  $(NtraitsPhenOnly*2+1):$(Ntraits*2)]
                #
                #                 HostTrait1chr1, HostTrait1chr2, ...
                N_::Vector{Int} , PhenS::Vector{Vector{Int}}                                    = RecursiveRandMultinom(N, SymetricalProbilitiesGeomDistr_Host[hostNoiseAlleles])
                $PRINT println("HostParasitesPhenotypes_INI! : RecursiveRandMultinom DONE")
            else
                N_ , PhenS = [N], [[$(join(["$(NallelesPerTrait),$(NallelesPerTrait)" for trait in 1:NtraitsPhenOnly],", ") )]]
            end
            RefPhen = Vector{Int}(undef,$(NtraitsPhenOnly*2))  ;  $(join(["RefPhen[$i] = hostGenotype[$i]" for i in 1:(NtraitsPhenOnly*2)], "; ")) # RefPhen = hostGenotype[1:$(NtraitsPhenOnly*2)]
            RefPhen .-= $NallelesPerTrait
            #
            $PRINT println("HostParasitesPhenotypes_INI! : RefPhen DONE")
            $PRINT if length(hostsGenotypeS.IDhS) != length(hostsGenotypeS.IDpS) global hostsGenotypeS__ = hostsGenotypeS end
            HostsIDphenS = Ref{Int_IDh}()
            $PRINT println("HostParasitesPhenotypes_INI! : PosiGenotypeInHostsGenotypeS DONE")
            $PRINT println("HostParasitesPhenotypes_INI! : Entering in eachindex(PhenS)")
            for i in eachindex(PhenS)
                $PRINT println("HostParasitesPhenotypes_INI! : Computing actual phenotype")
                PhenS[i] .+= RefPhen # .- Npossible alleles {in hosts} || .- NallelesPerTrait{in parasites}
                # Host
                $PRINT println("HostParasitesPhenotypes_INI! : checking AllelesRangePerSp")
                for chr in 0:1
                    for trait in 1:$NtraitsPhenOnly
                        traitChr = trait*2-chr
                        if PhenS[i][traitChr] < 1
                            PhenS[i][traitChr] += $(NallelesPerTrait_EachTrait)[trait] # circular trait /!| we are conscidering the phenotype !
                        elseif PhenS[i][traitChr] > $(NallelesPerTrait_EachTrait)[trait] # circular trait /!| we are conscidering the phenotype !
                            PhenS[i][traitChr] -= $(NallelesPerTrait_EachTrait)[trait] # circular trait /!| we are conscidering the phenotype !
                        end
                    end
                end
                for trait in 1:2:$(NtraitsPhenOnly*2) # order chromosoms
                    if PhenS[i][trait] > PhenS[i][trait+1]
                        PhenS[i][trait:(trait+1)] = PhenS[i][(trait+1):-1:trait]
                    end
                end
                #
                $(ifelse(PhenotypicStochasticity
                ,"append!(PhenS[i],hostNoiseAlleles)"
                ,""
                ))
                HostsIDphenS[] = GetHostID( PhenS[i][1:$(Ntraits*2)] )
                # Is this combinaison of host and parasit genotype and phenotype already existing in the local population ?
                posiLocal = findall((hostsGenotypeS.IDhSphen .== HostsIDphenS[]) .& (hostsGenotypeS.IDhS .== HostID) .& (hostsGenotypeS.IDpS .== $(Int_IDp(0))))
                $PRINT println("HostParasitesPhenotypes_INI! : posiLocal DONE")
                $PRINT if length(posiLocal)>1 error("the choice between categories that matches host and parasit genotype and phenotype has not yet been implemented") end
                if length(posiLocal) === 0 # it does not exist
                    $PRINT println("HostParasitesPhenotypes_INI! : Locally creating the host category")
                    push!(hostsGenotypeS.IDhS    ,HostID)
                    push!(hostsGenotypeS.IDhSphen,HostsIDphenS[])
                    push!(hostsGenotypeS.IDpS    ,Int_IDp(0))
                    push!(hostsGenotypeS.IDpSphen,Int_IDp(0))
                    push!(hostsGenotypeS.IDpGenPhen_HostShiftHistoryS,:IDparasitCat_0_0)
                    push!(hostsGenotypeS.IDp_HostShiftHistoryS       ,:IDparasitCat_0_0)
                    $(ifelse(SimulateImmunisation
                    ,"push!(hostsGenotypeS.IDiS    ,:_0)"
                    ,""    ))
                    push!(hostsGenotypeS.IDcategorieS, GetCategoryID(HostID, HostsIDphenS[], :IDparasitCat_0_0))
                    $(ifelse(GfG
                    ,"push!(hostsGenotypeS.rate_of_gametes_prodS      ,GetRate_of_gametes_prod(hostGenotype))
                      push!(hostsGenotypeS.rate_of_parasites_emissionS,0.0)"
                    ,""))
                    push!(hostsGenotypeS.virulenceS ,0.0)
                    push!(hostsGenotypeS.Precoveryinnateimmu, 0.0  )
                    $(ifelse(SimulateImmunisation
                    ,"push!(hostsGenotypeS.Precoveryacquiredimmu, 0.0  )" # GetPrecoveryacquiredimmu!(MetaComm.HPinteractions.Precoveryacquiredimmu, posiParasGeno, 0, MetaComm.parasitesGenotypeS_List, MetaComm.HPinteractions.TempDistImmunisingInfecting)
                    ,""))
                    push!(hostsGenotypeS.posiInParasitGenotypeS_List, Int_IDp(0))
                    push!(hostsGenotypeS.N_      ,N_[i])
                    push!(hostsGenotypeS.dN_     ,0    )
                    hostsGenotypeS.N[] += N_[i]
                    #
                    # Is this host *phenotype* already existing in the MetaCommunity ?
                    posiHostPhen = findfirst(xx -> xx === HostsIDphenS[],MetaComm.hostsGenotypeS_List.IDhS)
                    if isnothing(posiHostPhen) # it does not exist
                        $PRINT println("HostParasitesPhenotypes_INI! : Globally creating the host category")
                        GetNewHostGenotype!(MetaComm, tuple(PhenS[i][1:$(Ntraits*2)]...), HostsIDphenS[])
                        push!(hostsGenotypeS.posiInHostsGenotypeS_List,length(MetaComm.hostsGenotypeS_List.IDhS))
                    else
                        push!(hostsGenotypeS.posiInHostsGenotypeS_List,posiHostPhen)
                    end
                else
                    $PRINT println("HostParasitesPhenotypes_INI! : Updating hostsGenotypeS.N_")
                    hostsGenotypeS.N_[posiLocal[1]] += N_[i]
                    $PRINT println("HostParasitesPhenotypes_INI! : Updating hostsGenotypeS.N[]")
                    hostsGenotypeS.N[] += N_[i]
                    $PRINT println("HostParasitesPhenotypes_INI! : hostsGenotypeS.N_ Updated")
                end
            end
        $INBOUNDend
        end"""
    eval(Meta.parse(str))
# end


if PhenotypicStochasticity
	# for SP in 1:NhostSp
    # function GetOffspring!(MetaComm::MetaCommunity, hostsGenotypeS::HostsGenotypeS, hostGenotype::NTuple{$Ntraits,Int}, N::Int, sp::Val{$SP})
	str = """
$INLINE function GetOffspring!(MetaComm::MetaCommunity, hostsGenotypeS::HostsGenotypeS, hostGenotype::NTuple{$Ntraits,Int}, N::Int)
        $FASTMATHstart
        $INBOUNDstart
		$PRINT println("GetOffspring!")
		HostID = GetHostID(hostGenotype)
		# Check the existence of the genotypes
		# if !any(MetaComm.hostsGenotypeS_List.IDhS .== HostID)
        if MetaComm.hostsGenotypeS_List.HostSexOffspringMap$(join(["[hostGenotype[$traitPosi]$addOneIfGfG_str]" for traitPosi in 1:(Ntraits*2)])) === 0
		    GetNewHostGenotype!(MetaComm, hostGenotype, HostID)
		end
		$PRINT println("GetOffspring! : GetNewHostGenotype!(MetaComm, hostGenotype, HostID) DONE")
        # $(join(["hostNoiseAlleles[\$i] = hostGenotype[\$ii]" for (i,ii) in zip(1:length((NtraitsPhenOnly*2+1):(Ntraits*2)), (NtraitsPhenOnly*2+1):(Ntraits*2))], "; ")))
		hostNoiseAlleles    = hostGenotype[$(NtraitsPhenOnly*2+1):$(Ntraits*2)]
        #
        #                 HostTrait1chr1,HostTrait1chr2,...,ParasTrait1,ParasTrait2,...
		N_::Vector{Int} , PhenS::Vector{Vector{Int}}                                     = RecursiveRandMultinom(N, SymetricalProbilitiesGeomDistr_Host[hostNoiseAlleles])
		$PRINT println("GetOffspring! : RecursiveRandMultinom DONE")
        RefPhen = Vector{Int}(undef,$(NtraitsPhenOnly*2))  ;  $(join(["RefPhen[$i] = hostGenotype[$i]" for i in 1:(NtraitsPhenOnly*2)], "; ")) # RefPhen = hostGenotype[1:$(NtraitsPhenOnly*2)]
        RefPhen .-= $NallelesPerTrait
		$PRINT println("GetOffspring! : RefPhen DONE")
		$PRINT if length(hostsGenotypeS.IDhS) != length(hostsGenotypeS.IDpS) global hostsGenotypeS__ = hostsGenotypeS end
		HostsIDphenS = Ref{Int_IDh}()
		$PRINT println("GetOffspring! : PosiGenotypeInHostsGenotypeS DONE")
		$PRINT println("GetOffspring! : Entering in eachindex(PhenS)")
		for i in eachindex(PhenS)
			$PRINT println("GetOffspring! : Computing actual phenotype")
			PhenS[i] .+= RefPhen # .- Npossible alleles {in hosts} || .- NallelesPerTrait{in parasites}
			# Host
			$SIMD for chr in 0:1
			    $SIMD for trait in 1:$NtraitsPhenOnly
				traitChr = trait*2-chr
                    if PhenS[i][traitChr] < 1
                        PhenS[i][traitChr] += $(NallelesPerTrait_EachTrait)[trait] # circular trait /!| we are conscidering the phenotype !
                    elseif PhenS[i][traitChr] > $(NallelesPerTrait_EachTrait)[trait] # circular trait /!| we are conscidering the phenotype !
                        PhenS[i][traitChr] -= $(NallelesPerTrait_EachTrait)[trait] # circular trait /!| we are conscidering the phenotype !
                    end
			    end
			end
			$SIMD for trait in 1:2:$(NtraitsPhenOnly*2) # order chromosomes
			    if PhenS[i][trait] > PhenS[i][trait+1]
                    PhenS[i][trait:(trait+1)] = PhenS[i][(trait+1):-1:trait]
			    end
			end
			$PRINT println("GetOffspring! : AllelesRangePerSp checked")
			#
			append!(PhenS[i],hostNoiseAlleles)
			HostsIDphenS[] = GetHostID(   PhenS[i][1:$(Ntraits*2)])
			# Is this combinaison of host and parasit genotype and phenotype already existing in the local population ?
			CategoryID = GetCategoryID(HostID, HostsIDphenS[], :IDparasitCat_0_0)
    		posiLocal = findfirst(xx -> xx === CategoryID, hostsGenotypeS.IDcategorieS)
			$PRINT println("GetOffspring! : posiLocal DONE")
			$PRINT if length(posiLocal)>1 error("the choice between categories that matches host and parasit genotype and phenotype has not been implemented") end
			if isnothing(posiLocal) # it does not exist
			    $PRINT println("GetOffspring! : Locally creating the host category")
			    push!(hostsGenotypeS.IDhS                        ,HostID)
			    push!(hostsGenotypeS.IDhSphen                    ,HostsIDphenS[])
			    push!(hostsGenotypeS.IDpS                        ,Int_IDp(0))
                push!(hostsGenotypeS.IDpSphen                    ,Int_IDp(0))
                push!(hostsGenotypeS.IDpGenPhen_HostShiftHistoryS,:IDparasitCat_0_0)
                push!(hostsGenotypeS.IDp_HostShiftHistoryS       ,:IDparasitCat_0_0)
                $(ifelse(SimulateImmunisation
                    ,"push!(hostsGenotypeS.IDiS                        ,:_0)"
                    ,""    ))
			    push!(hostsGenotypeS.posiInParasitGenotypeS_List ,0)
			    push!(hostsGenotypeS.IDcategorieS                ,CategoryID)
                $(ifelse(GfG
                ,"push!(hostsGenotypeS.rate_of_gametes_prodS      ,GetRate_of_gametes_prod(PhenS[i][1:$(Ntraits*2)]))
                  push!(hostsGenotypeS.rate_of_parasites_emissionS,0.0                     )"
                ,""))
			    push!(hostsGenotypeS.virulenceS                  ,0.0  )
			    push!(hostsGenotypeS.Precoveryinnateimmu         ,0.0  )
                $(ifelse(SimulateImmunisation
                ,"push!(hostsGenotypeS.Precoveryacquiredimmu, 0.0  )" # GetPrecoveryacquiredimmu!(MetaComm.HPinteractions.Precoveryacquiredimmu, posiParasGeno, 0, MetaComm.parasitesGenotypeS_List, MetaComm.HPinteractions.TempDistImmunisingInfecting)
                ,""))
                push!(hostsGenotypeS.N_                          ,0    )
                push!(hostsGenotypeS.dN_                         ,N_[i])
			    #
			    # Is this host *phenotype* already existing in the *MetaCommunity* ?
			    # posiHostPhen = findfirst(xx -> xx === HostsIDphenS[], MetaComm.hostsGenotypeS_List.IDhS)
                posiHostPhen = MetaComm.hostsGenotypeS_List.HostSexOffspringMap$(join(["[PhenS[i][$traitPosi]$addOneIfGfG_str]" for traitPosi in 1:(Ntraits*2)]))
			    # if isnothing(posiHostPhen) # it does not exist
                if posiHostPhen === 0
                    $PRINT println("GetOffspring! : Globally creating the host category")
                    GetNewHostGenotype!(MetaComm, tuple(PhenS[i][1:$(Ntraits*2)]...), HostsIDphenS[])
                    push!(hostsGenotypeS.posiInHostsGenotypeS_List,length(MetaComm.hostsGenotypeS_List.IDhS))
			    else
                    push!(hostsGenotypeS.posiInHostsGenotypeS_List,posiHostPhen)
			    end
			else
			    $PRINT println("GetOffspring! : Updating hostsGenotypeS.N_")
			    hostsGenotypeS.dN_[posiLocal] += N_[i]
			    $PRINT println("GetOffspring! : Updating hostsGenotypeS.N[]")
			    $PRINT println("GetOffspring! : hostsGenotypeS.N_ Updated")
			end
	    end
    $INBOUNDend
    $FASTMATHend
	end"""
	eval(Meta.parse(str))
	# end
else
    # for SP in 1:NhostSp
    # function GetOffspring!(MetaComm::MetaCommunity, hostsGenotypeS::HostsGenotypeS, hostGenotype::NTuple{$(Ntraits*2),Int}, N::Int, sp::Val{$SP})
    str = """
    $INLINE function GetOffspring!(MetaComm::MetaCommunity, hostsGenotypeS::HostsGenotypeS, hostGenotype::NTuple{$(Ntraits*2),Int}, N::Int)
        $FASTMATHstart
        $INBOUNDstart
            $PRINT println("GetOffspring!")
            HostID = GetHostID(hostGenotype)
            CategoryID = GetCategoryID(HostID, HostID, :IDparasitCat_0_0)
            # Is this combinaison of host and parasit genotype and phenotype already existing in the local population ?
            posiLocal = findfirst(xx -> xx === CategoryID, hostsGenotypeS.IDcategorieS)
            $PRINT println("GetOffspring! : posiLocal DONE")
            if isnothing(posiLocal) # it does not exist
                $PRINT println("GetOffspring! : Locally creating the host category")
                push!(hostsGenotypeS.IDhS                        ,HostID)
                push!(hostsGenotypeS.IDhSphen                    ,HostID)
                push!(hostsGenotypeS.IDpS                        ,Int_IDp(0))
                push!(hostsGenotypeS.IDpSphen                    ,Int_IDp(0))
                push!(hostsGenotypeS.IDpGenPhen_HostShiftHistoryS,:IDparasitCat_0_0)
                push!(hostsGenotypeS.IDp_HostShiftHistoryS       ,:IDparasitCat_0_0)
                $(ifelse(SimulateImmunisation
                    ,"push!(hostsGenotypeS.IDiS                        ,:_0)"
                    ,""    ))
                push!(hostsGenotypeS.posiInParasitGenotypeS_List ,0)
                push!(hostsGenotypeS.IDcategorieS                ,CategoryID)
                $(ifelse(GfG
                ,"push!(hostsGenotypeS.rate_of_gametes_prodS      ,GetRate_of_gametes_prod(hostGenotype))
                  push!(hostsGenotypeS.rate_of_parasites_emissionS,0.0         )"
                ,""))
                push!(hostsGenotypeS.virulenceS                  ,0.0  )
                push!(hostsGenotypeS.Precoveryinnateimmu         ,0.0  )
                $(ifelse(SimulateImmunisation
                ,"push!(hostsGenotypeS.Precoveryacquiredimmu, 0.0  )" # GetPrecoveryacquiredimmu!(MetaComm.HPinteractions.Precoveryacquiredimmu, posiParasGeno, 0, MetaComm.parasitesGenotypeS_List, MetaComm.HPinteractions.TempDistImmunisingInfecting)
                ,""))
                push!(hostsGenotypeS.dN_                         ,N   )
                push!(hostsGenotypeS.N_                          ,0   )
                #
                # Check the existence of the genotypes == phenotype
                # posiHost = findfirst(xx -> xx === HostID, MetaComm.hostsGenotypeS_List.IDhS)
                posiHost = MetaComm.hostsGenotypeS_List.HostSexOffspringMap$(join(["[hostGenotype[$traitPosi]$addOneIfGfG_str]" for traitPosi in 1:(Ntraits*2)]))
                # if isnothing(posiHost)
                if posiHost === 0
                    GetNewHostGenotype!(MetaComm, hostGenotype, HostID)
                    push!(hostsGenotypeS.posiInHostsGenotypeS_List, length(MetaComm.hostsGenotypeS_List.IDhS))
                else
                    push!(hostsGenotypeS.posiInHostsGenotypeS_List, posiHost)
                end
            else
                $PRINT println("GetOffspring! : Updating hostsGenotypeS.N_")
                hostsGenotypeS.dN_[posiLocal] += N
                $PRINT println("GetOffspring! : Updating hostsGenotypeS.N[]")
                $PRINT println("GetOffspring! : hostsGenotypeS.N_ Updated")
            end
            $INBOUNDend
            $FASTMATHend
        end"""
        eval(Meta.parse(str))
end

function SetParasitMutatingMap!(X::MetaCommunity)
    for trait in 1:Ntraits
        for g in eachindex(X.parasitesGenotypeS_List.IDpS)
            push!(X.ParasitMutatingMap[trait],ExtendableDict{Int,Int}() )
        end
    end
end

str = """
function Check_ParasitMutatingMap(X::MetaCommunity,IDpS::Union{Nothing,Vector{Int_IDp}}=nothing)
    $PRINT println("In Check_ParasitMutatingMap")
    if isnothing(IDpS)
        IDpS = X.parasitesGenotypeS_List.IDpS
    end
    W = filter(i -> in(X.parasitesGenotypeS_List.IDpS[i],IDpS), eachindex(X.parasitesGenotypeS_List.IDpS))
    for g in W
        global CurrentGenotype = X.parasitesGenotypeS_List.GenotypeS[g].Alleles
        for trait in 1:Ntraits
            for mut in keys(X.ParasitMutatingMap[trait][g])
                if !iszero( X.ParasitMutatingMap[trait][g][mut])
                    NewGenotypeShouldBe = collect(CurrentGenotype)
                    NewGenotypeShouldBe[trait] += mut
                    NewGenotypeIs = collect(X.parasitesGenotypeS_List.GenotypeS[X.ParasitMutatingMap[trait][g][mut]].Alleles)
                    if NewGenotypeShouldBe != NewGenotypeIs
                        global Xerror = CurrentGenotype
                        global MetaComm_error = X
                        error("Check_ParasitMutatingMap found this error:
In ParasitMutatingMap, at trait=\$trait and posi=\$g, got genotype \$(X.ParasitMutatingMap[trait][g][mut])
which is                        \$(NewGenotypeIs)
while CurrentGenotype is        \$(CurrentGenotype)
therfore it should actually be  \$(NewGenotypeShouldBe)"
                        )
                    end
                end
            end
        end
    end
    println("Check_ParasitMutatingMap found no errors !")
end"""
eval(Meta.parse(str))

#### MutateParasites!
# Probilities Geometric distribution
if MuEffectParasites * NallelesPerTrait_infection_success_AllSp < 1.0 error("MuEffectParasites which is $MuEffectParasites must be above 1/NallelesPerTrait_infection_success_AllSp which is $(1/NallelesPerTrait_infection_success_AllSp)") end
if MuEffectParasites * NallelesPerTrait                         < 1.0 error("MuEffectParasites which is $MuEffectParasites must be above 1/NallelesPerTrait                         which is $(1/NallelesPerTrait                        )") end

function GetP(Nalleles::Int)
    p = 1/(MuEffectParasites*Nalleles)
    p = map(k -> (1-p)^(k-1)*p , 1:Nalleles)
    return(p/sum(p))
end
p = tuple(map(Nalleles -> GetP(Nalleles-1), NallelesPerTrait_EachTrait )...)

if !GfG global MuParasit_trait_presence = NaN  ; global MuHost_trait_presence = fill(NaN,NhostSp) end
str = """
$(ifelse(GfG
,ifelse(PhenotypicStochasticity
    ,"$INLINE function MutateParasites!(HPinteractions::HPinteractions, parasitesGenotypeS_List::ParasitesGenotypeS_List, hostsGenotypeS_List::HostsGenotypeS_List, ParasitMutatingMap::NTuple{Ntraits,Vector{ExtendableDict{Int,Int}}}, N_IDp_HostShiftHistoryS::Vector{Int}, IDpS::Vector{Int}, IDp_HostShiftHistoryS::Vector{Symbol}, posiInParasitGenotypeS_List::Vector{Int} ,Nmutating::Vector{Int}, NmutatingEachTrait::Vector{Int}, NmuSize::NTuple{$Ntraits,Vector{Int}}, NmutEachAllelesAppearingTrait::NTuple{$Ntraits,Vector{Int}},                                             NewGenotype::Vector{Int}, NewId::Ref{Int}, dt::Float)"
    ,"$INLINE function MutateParasites!(HPinteractions::HPinteractions, parasitesGenotypeS_List::ParasitesGenotypeS_List, hostsGenotypeS_List::HostsGenotypeS_List, ParasitMutatingMap::NTuple{Ntraits,Vector{ExtendableDict{Int,Int}}}, N_IDp_HostShiftHistoryS::Vector{Int}, IDpS::Vector{Int}, IDp_HostShiftHistoryS::Vector{Symbol}, posiInParasitGenotypeS_List::Vector{Int} ,Nmutating::Vector{Int}, NmutatingEachTrait::Vector{Int}, NmuSize::NTuple{$Ntraits,Vector{Int}}, NmutEachAllelesAppearingTrait::NTuple{$Ntraits,Vector{Int}}, rate_of_parasites_emissionS::Vector{Float}, NewGenotype::Vector{Int}, NewId::Ref{Int}, dt::Float)")
,    "$INLINE function MutateParasites!(HPinteractions::HPinteractions, parasitesGenotypeS_List::ParasitesGenotypeS_List, hostsGenotypeS_List::HostsGenotypeS_List, ParasitMutatingMap::NTuple{Ntraits,Vector{ExtendableDict{Int,Int}}}, N_IDp_HostShiftHistoryS::Vector{Int}, IDpS::Vector{Int}, IDp_HostShiftHistoryS::Vector{Symbol}, posiInParasitGenotypeS_List::Vector{Int} ,Nmutating::Vector{Int}, NmutatingEachTrait::Vector{Int}, NmuSize::NTuple{$Ntraits,Vector{Int}},                                                                                                          NewGenotype::Vector{Int}, NewId::Ref{Int}, dt::Float)"
))
    $FASTMATHstart
    $INBOUNDstart
    $PRINT println("MutateParasites!")
    for idxMutatingIDp in eachindex(IDpS)
        PosiMutatingIDp::Int =   posiInParasitGenotypeS_List[idxMutatingIDp]
        $PRINT println("PosiMutatingIDp = ",PosiMutatingIDp)
        genotype::NTuple{$Ntraits,Int} = parasitesGenotypeS_List.GenotypeS[PosiMutatingIDp].Alleles
        $PRINT println("genotype = ",genotype)
        $(ifelse(GfG
        ,"Ntrait_Present::Int = sum(genotype .!== 0)
          $PRINT println(:Ntrait_Present_,Ntrait_Present)
          Mu_trait_alleles::Float = $MuParasit_trait_alleles * Ntrait_Present
          $PRINT println(:Mu_trait_alleles_,Mu_trait_alleles)
          $PRINT $(ifelse(GfG,"println((Mu_trait_alleles + $(MuParasit_trait_presence*Ntraits))*dt)",""))
          Nmutating[1] = rand(Binomial(N_IDp_HostShiftHistoryS[idxMutatingIDp],(Mu_trait_alleles + $(MuParasit_trait_presence*Ntraits))*dt)) # total number of mutants"
        ,"Nmutating[1] = rand(Binomial(N_IDp_HostShiftHistoryS[idxMutatingIDp],$(MuParasit_trait_alleles)*dt)) # total number of mutants"
        ))
        $PRINT println("Nmutating[1] = ",Nmutating[1])
        if Nmutating[1] > 0
            $(ifelse(GfG
            ,"Nmutating[3] = rand(Binomial(Nmutating[1], $(MuParasit_trait_presence*Ntraits) / (Mu_trait_alleles + $(MuParasit_trait_presence*Ntraits)))) # number that swap the presence-absence of a trait
              Nmutating[1] -= Nmutating[3]"
            ,""
            ))
            $PRINT println("Nmutating[1] = ",Nmutating[1])
            $PRINT $(ifelse(GfG,"println(\"Nmutating[3] = \",Nmutating[3])",""))
              Nmutating[2] = rand(Binomial(Nmutating[1])) # number that increase for a trait
              Nmutating[1] -= Nmutating[2] # number that decrease for a trait
              $PRINT println("Nmutating[1] = ",Nmutating[1])
              $PRINT println("Nmutating[2] = ",Nmutating[2])
              $PRINT $(ifelse(GfG,"println(\"Nmutating[3] = \",Nmutating[3])",""))
            for dir in eachindex(Nmutating) # 1 = decrease ;  2 = increase ;  3 = swap presence of the trait
                $PRINT println("dir = ",dir)
                # NmutatingEachTrait : every trait can mutate several times -> sampling with resampling
                $(ifelse(GfG
                ,"if dir === 3
                    rand!(Multinomial( Nmutating[dir], $(fill(1/Ntraits, Ntraits))                                               ), NmutatingEachTrait) # every trait can have its present that is swapp
                else
                    rand!(Multinomial( Nmutating[dir], [if (genotype[t] === 0) 0.0 else 1/Ntrait_Present end for t in 1:$Ntraits]), NmutatingEachTrait) # only traits that are present can mutate
                end"
                ,"  rand!(Multinomial( Nmutating[dir], $(fill(1/Ntraits, Ntraits))                                               ), NmutatingEachTrait) # every trait can mutate several times -> sampling with resampling"
                ))
                for mutatingTrait in 1:$(Ntraits)
                    $PRINT println("mutatingTrait = ",mutatingTrait)
                    # Which mutation for the mutatingTrait ?
                    if NmutatingEachTrait[mutatingTrait] > 0
                        if dir === 3 # Choose the alleles of the traits that is appearing
                            if genotype[mutatingTrait] === 0 # the trait is appearing
                                rand!(Multinomial(NmutatingEachTrait[mutatingTrait],fill(1/$NallelesPerTrait_EachTrait[mutatingTrait], $NallelesPerTrait_EachTrait[mutatingTrait])), NmutEachAllelesAppearingTrait[mutatingTrait])
                            else                             # the trait is disappearing (only one type of mutation is happening)
                                fill!(NmutEachAllelesAppearingTrait[mutatingTrait],0)
                                NmutEachAllelesAppearingTrait[mutatingTrait][genotype[mutatingTrait]] = NmutatingEachTrait[mutatingTrait]
                            end
                        else
                            rand!(Multinomial(NmutatingEachTrait[mutatingTrait],$p[mutatingTrait]), NmuSize[mutatingTrait]) # in a geometric distribution
                            $PRINT println("NmuSize[mutatingTrait] = ",NmuSize[mutatingTrait])
                        end
                        for muSize in (  dir===3  ?  (1:$NallelesPerTrait_EachTrait[mutatingTrait])  :  (1:($NallelesPerTrait_EachTrait[mutatingTrait]-1))  ) # dir = 3  ==>>  Either the current allele is 0 ( => size of the mutation = value of the new allele) or the current allele is [[ index of NmutEachAllelesAppearingTrait that is not zero ]] (=> size of the mutation = value of the current allele)
                            $PRINT println("muSize = ",muSize)
                            Nmut::Int = (dir===3  ?  NmutEachAllelesAppearingTrait[mutatingTrait][muSize]  :  NmuSize[mutatingTrait][muSize])
                            if Nmut > 0  # if there are some mutants with that muSize
                                if     dir === 1 muSizeSign = -muSize
                                elseif dir === 2 muSizeSign = muSize
                                else # dir === 3 => if the trait disappear    muSizeSign = -muSize   else (the trait is appearing)  muSizeSign = muSize
                                    if genotype[mutatingTrait] === 0 # the trait is appearing
                                        muSizeSign =  muSize
                                    else
                                        muSizeSign = -genotype[mutatingTrait]
                                    end
                                end
                                $PRINT println("genotype[mutatingTrait] = ",genotype[mutatingTrait])
                                $PRINT println("muSizeSign = ",muSizeSign)
                                $PRINT println("dir = ",dir)
                                if dir!==3
                                    if  (genotype[mutatingTrait] + muSizeSign) < 1 # Circular mutation
                                        muSizeSign += $NallelesPerTrait_EachTrait[mutatingTrait]
                                        muSize = abs(muSizeSign)
                                    elseif (genotype[mutatingTrait] + muSizeSign) > $NallelesPerTrait_EachTrait[mutatingTrait] # Circular mutation
                                        muSizeSign -= $NallelesPerTrait_EachTrait[mutatingTrait]
                                        muSize = abs(muSizeSign)
                                    end
                                end
                                if muSize !== 0
                                    $PRINT println("-> genotype[mutatingTrait] = ",genotype[mutatingTrait])
                                    $PRINT println("-> muSizeSign = ",muSizeSign)
                                    posiGlobalNewId = ParasitMutatingMap[mutatingTrait][PosiMutatingIDp][ muSizeSign]
                                    $PRINT println("posiGlobalNewId = ",posiGlobalNewId)
                                    if iszero(posiGlobalNewId) # we need to update the ParasitMutatingMap
                                        $PRINT println("genotype = ",genotype)
                                        NewGenotype .= genotype ; NewGenotype[mutatingTrait] += muSizeSign
                                        $PRINT println("NewGenotype = ",NewGenotype)
                                        $PRINT if any(NewGenotype .< 0) error("any(NewGenotype .< 0)") end
                                        NewId[] = GetParasitID(NewGenotype)
                                        $PRINT println("NewId[] = ",NewId[])
                                        # Is the new genotype absent in the metaCommunity ?
                                        posiGlobalNewId = findfirst(xx -> xx === NewId[], parasitesGenotypeS_List.IDpS)
                                        $PRINT println("posiGlobalNewId = ",posiGlobalNewId)
                                        if isnothing(posiGlobalNewId) # the new genotype is absent
                                            $PRINT println("Getting GetNewParasGenotype!")
                                            GetNewParasGenotype!(HPinteractions, parasitesGenotypeS_List, hostsGenotypeS_List, ParasitMutatingMap, tuple(NewGenotype...), NewId[])
                                            posiGlobalNewId = length(parasitesGenotypeS_List.IDpS)
                                        end
                                        # update the ParasitMutatingMap
                                        ParasitMutatingMap[mutatingTrait][PosiMutatingIDp].x[ muSizeSign] = posiGlobalNewId # the mutation that is happening
                                        ParasitMutatingMap[mutatingTrait][posiGlobalNewId].x[-muSizeSign] = PosiMutatingIDp # the revert mutation
                                        $PRINT println("mutatingTrait = ",mutatingTrait )
                                        $PRINT println("muSizeSign = "   ,muSizeSign    )
                                        for futureMutatingTrait in 1:$Ntraits
                                            $PRINT println("futureMutatingTrait = "   ,futureMutatingTrait    )
                                            for futureChange in keys(ParasitMutatingMap[futureMutatingTrait][PosiMutatingIDp])
                                                $PRINT println("futureChange = ",futureChange )
                                                if !iszero(ParasitMutatingMap[futureMutatingTrait][PosiMutatingIDp][futureChange])
                                                    $PRINT if iszero(ParasitMutatingMap[futureMutatingTrait][PosiMutatingIDp][futureChange])
                                                    $PRINT     error("iszero(ParasitMutatingMap[futureMutatingTrait][PosiMutatingIDp][futureChange])")
                                                    $PRINT end
                                                    IfFuturChangeHappenInFormerGenot = ParasitMutatingMap[futureMutatingTrait][PosiMutatingIDp][futureChange]
                                                    $PRINT println("IfFuturChangeHappenInFormerGenot = ",IfFuturChangeHappenInFormerGenot )
                                                    If_IfFuturChangeHappenInFormerGenot_mutatingTrait_mut = ParasitMutatingMap[mutatingTrait][IfFuturChangeHappenInFormerGenot][muSizeSign]
                                                    $PRINT println("If_IfFuturChangeHappenInFormerGenot_mutatingTrait_mut = ",If_IfFuturChangeHappenInFormerGenot_mutatingTrait_mut )
                                                    if !iszero(If_IfFuturChangeHappenInFormerGenot_mutatingTrait_mut)
                                                        ParasitMutatingMap[futureMutatingTrait][posiGlobalNewId].x[futureChange] = If_IfFuturChangeHappenInFormerGenot_mutatingTrait_mut
                                                        $PRINT println("Check_ParasitMutatingMap futureChange mutatingTrait=\$mutatingTrait, futureChange=\$futureChange")
                                                    end
                                                end
                                            end
                                        end
                                    else
                                        NewId[] = parasitesGenotypeS_List.IDpS[posiGlobalNewId]
                                    end
                                    posiNewId = findfirst(xx -> xx === NewId[], IDpS)
                                    if isnothing(posiNewId) # the new genotype is absent from the ParasitesPop
                                         push!(IDpS                       , NewId[])
                                         $PRINT println("Symbol(:IDp_,NewId[], __ .*split(string(IDp_HostShiftHistoryS[idxMutatingIDp]), __ )[2])  =  ")
                                         $PRINT println(Symbol(:IDp_,NewId[],"__".*split(string(IDp_HostShiftHistoryS[idxMutatingIDp]),"__")[2]))
                                         push!(IDp_HostShiftHistoryS      , Symbol(:IDp_,NewId[],"__".*split(string(IDp_HostShiftHistoryS[idxMutatingIDp]),"__")[2])   )
                                         $PRINT println("muSize = ", muSize)
                                         push!(N_IDp_HostShiftHistoryS    , Nmut )
                                         $PRINT println("posiGlobalNewId = ", posiGlobalNewId)
                                         push!(posiInParasitGenotypeS_List, posiGlobalNewId)
                                         $PRINT println("  parasitesGenotypeS_List.GenotypeS[posiGlobalNewId].Alleles    = ", parasitesGenotypeS_List.GenotypeS[posiGlobalNewId].Alleles   )
                                         $PRINT $(ifelse(GfG&(!PhenotypicStochasticity),"println(\"GetRate_of_parasites_emission(   parasitesGenotypeS_List.GenotypeS[posiGlobalNewId].Alleles   ) = \", GetRate_of_parasites_emission(   parasitesGenotypeS_List.GenotypeS[posiGlobalNewId].Alleles   ))",""))
$(ifelse(GfG&(!PhenotypicStochasticity),"push!(rate_of_parasites_emissionS , GetRate_of_parasites_emission(   parasitesGenotypeS_List.GenotypeS[posiGlobalNewId].Alleles   ))",""))
                                    else # the new genotype is present
                                        N_IDp_HostShiftHistoryS[posiNewId] += Nmut
                                    end
                                    N_IDp_HostShiftHistoryS[idxMutatingIDp] -= Nmut # these mutations actually happen
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    $PRINT println("MutateParasites! Done")
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

#### SetNParasites!  This also implements the parasit PhenotypicStochasticity
str = """
function SetNParasites!(X::MetaCommunity, mutate::Bool = false, dt::Float = 0.0)
    $FASTMATHstart
    $PRINT println("SetNParasites: ",MetaCommunity)
    for (sp,x) in enumerate(getfield(X,:$(fieldname(MetaCommunity,1))))
        SetNParasites_MetaPop!(x, X.HPinteractions, X.parasitesGenotypeS_List, X.hostsGenotypeS_List, X.ParasitMutatingMap, mutate, dt)
        X.N_[sp] = x.N[]
    end
    X.N[] = $(join(["X.N_[$i]" for i in 1:NhostSp], " + "))
    $PRINT println("After SetNParasites! Check_ParasitMutatingMap(X)")
    $PRINT Check_ParasitMutatingMap(X)
    $FASTMATHend
end"""
eval(Meta.parse(str))

str = """
$INLINE function SetNParasites_MetaPop!(X::MetaPop, HPinteractions::HPinteractions, parasitesGenotypeS_List::ParasitesGenotypeS_List, hostsGenotypeS_List::HostsGenotypeS_List, ParasitMutatingMap::NTuple{Ntraits,Vector{ExtendableDict{Int,Int}}}, mutate::Bool, dt::Float)
    $FASTMATHstart
    $PRINT println("SetNParasites: ",MetaPop)
    for (pop,x) in enumerate(getfield(X,:$(fieldname(MetaPop,1))))
        SetNParasites_HostSexS!(x, HPinteractions, parasitesGenotypeS_List, hostsGenotypeS_List, ParasitMutatingMap, mutate, dt)
        X.N_[pop] = x.N[]
    end
    X.N[] = $(join(["X.N_[$i]" for i in 1:NPops], " + "))    $FASTMATHend
end"""
eval(Meta.parse(str))


if PhenotypicStochasticity
	str = """
$INLINE function SetNParasites_HostSexS!(X::HostSexS, HPinteractions::HPinteractions, parasitesGenotypeS_List::ParasitesGenotypeS_List, hostsGenotypeS_List::HostsGenotypeS_List, ParasitMutatingMap::NTuple{Ntraits,Vector{ExtendableDict{Int,Int}}}, mutate::Bool, dt::Float)
    $FASTMATHstart
    $INBOUNDstart
        for sex in 1:2
            X.N_[sex] = X.HostSexS[sex].N[] = sum(X.HostSexS[sex].N_)
        end
        X.N[] = X.N_[1] + X.N_[2]
        if X.N[] > 0
    	    $PRINT println("SetNParasites: ",HostSexS)
    	    empty!(X.ParasitesPop.IDpS)
    	    empty!(X.ParasitesPop.IDpSphen)
            empty!(X.ParasitesPop.IDpGenPhen_HostShiftHistoryS)
            empty!(X.ParasitesPop.IDp_HostShiftHistoryS)
            empty!(X.ParasitesPop.N_)
            empty!(X.ParasitesPop.posiInParasitGenotypeS_List)
            $(ifelse(GfG,"empty!(X.ParasitesPop.rate_of_parasites_emissionS)",""))
            $PRINT println("IDp_HostShiftHistoryS = ",IDp_HostShiftHistoryS)
            #
            IDp_HostShiftHistoryS       = Vector{Symbol }()
            N_IDp_HostShiftHistoryS     = Vector{Int    }()
            IDpS                        = Vector{Int_IDp}()
            posiInParasitGenotypeS_List = Vector{Int    }()
            #
            for sex in 1:2
                posiParasites = findall(X.HostSexS[sex].IDpS .!== Int_IDp(0))
                resize!(X.HostSexS[sex].NparasitesEmitted, length(posiParasites))
                $(ifelse(GfG
                ,"X.HostSexS[sex].NparasitesEmitted .= rand.(Poisson.(X.HostSexS[sex].N_[posiParasites] .* X.HostSexS[sex].rate_of_parasites_emissionS[posiParasites]))"
                ,"X.HostSexS[sex].NparasitesEmitted .= rand.(Poisson.(X.HostSexS[sex].N_[posiParasites] .* $MaxRate_of_parasites_emission))"
                ))
                for (ii,i) in enumerate(posiParasites)
                    posiInParasitesPop = findfirst(xx -> xx === X.HostSexS[sex].IDp_HostShiftHistoryS[i],  IDp_HostShiftHistoryS)
                    if isnothing(posiInParasitesPop)
                        push!(IDp_HostShiftHistoryS      , X.HostSexS[sex].IDp_HostShiftHistoryS[      i])
                        push!(N_IDp_HostShiftHistoryS    , X.HostSexS[sex].NparasitesEmitted[         ii])
                        push!(posiInParasitGenotypeS_List, X.HostSexS[sex].posiInParasitGenotypeS_List[i])
                        push!(IDpS                       , X.HostSexS[sex].IDpS[                       i])
                    else
                        N_IDp_HostShiftHistoryS[posiInParasitesPop] += X.HostSexS[sex].NparasitesEmitted[ii])
                    end
                end
            end
            #
            if length(IDp_HostShiftHistoryS) !== 0 # if there are some parasites in the current population
                if mutate
                    $PRINT println("Entring: MutateParasites!")
                    $(ifelse(GfG
                     # where MutateParasites! uptates N_IDp_HostShiftHistoryS AND, _only in case of creation of new genotypes_, IDpS as well as IDp_HostShiftHistoryS
                    ,"MutateParasites!(HPinteractions, parasitesGenotypeS_List, hostsGenotypeS_List, ParasitMutatingMap, N_IDp_HostShiftHistoryS, IDpS, IDp_HostShiftHistoryS, posiInParasitGenotypeS_List ,X.ParasitesPop.Nmutating, X.ParasitesPop.NmutatingEachTrait, X.ParasitesPop.NmuSize, X.ParasitesPop.NmutEachAllelesAppearingTrait                                          , X.ParasitesPop.NewGenotype, X.ParasitesPop.NewId, dt::Float)"
                    ,"MutateParasites!(HPinteractions, parasitesGenotypeS_List, hostsGenotypeS_List, ParasitMutatingMap, N_IDp_HostShiftHistoryS, IDpS, IDp_HostShiftHistoryS, posiInParasitGenotypeS_List ,X.ParasitesPop.Nmutating, X.ParasitesPop.NmutatingEachTrait, X.ParasitesPop.NmuSize,                                                                                         X.ParasitesPop.NewGenotype, X.ParasitesPop.NewId, dt::Float)"
                    ))
                    # For some very rare parasites genotype, all individuals might have mutated (e.g. if X.ParasitesPop.N_[...] === 1)
                    Del = findall(N_IDp_HostShiftHistoryS .=== 0)
                    if length(Del) > 0
                        deleteat!(IDpS  ,Del)
                        deleteat!(N_IDp_HostShiftHistoryS    ,Del)
                        deleteat!(IDp_HostShiftHistoryS      ,Del)
                        deleteat!(posiInParasitGenotypeS_List,Del)
                    end
            	end
                # implement the PhenotypicStochasticity
            	for (i,IDp) in enumerate(IDpS)
            	    Genot = parasitesGenotypeS_List.GenotypeS[posiInParasitGenotypeS_List[i]].Alleles
                    $(ifelse(GfG
                    ,"MissingTraits = findall(Genot .== 0)"
                    ,""))
                    $(join(["X.RefPhen[$i] = Genot[$i]" for i in 1:NtraitsPhenOnly], "; ")) # X.RefPhen = Genot[1:NtraitsPhenOnly]
            	    N_::Vector{Int} , PhenS::Vector{Vector{Int}} = RecursiveRandMultinom(N_IDp_HostShiftHistoryS[i], SymetricalProbilitiesGeomDistr_Parasit[Genot[$(NtraitsPhenOnly+1):end]])
            	    X.RefPhen .-= $(NallelesPerTrait) # .- Npossible alleles {in hosts} || .- NallelesPerTrait{in parasites}
            	    ParasIDphen = Ref{Int_IDp}()
            	    ParasIDpGenPhen_HostShiftHistoryS = Ref{Symbol}()
            	    for i in eachindex(PhenS)
                        $(ifelse(GfG
                        ,"PhenS[i][MissingTraits] .= 0"
                        ,""))
                        PhenS[i] .+= RefPhen # .- Npossible alleles {in hosts} || .- NallelesPerTrait{in parasites}
                		for traitPosi in 1:$(NtraitsPhenOnly) ## circular traits
                		    if PhenS[i][traitPosi] < 1
                		        PhenS[i][traitPosi]  += $NallelesPerTrait_EachTrait[traitPosi]
                		    end
                		    if PhenS[i][traitPosi] > $NallelesPerTrait_EachTrait[traitPosi]
                		        PhenS[i][traitPosi]  -= $NallelesPerTrait_EachTrait[traitPosi]
                		    end
                		end
                		#
                		PhenS[i] = vcat(PhenS[i],  Genot[$(NtraitsPhenOnly+1):end])
                		ParasIDphen[] = GetParasitID(PhenS[i])
                		ParasIDpGenPhen_HostShiftHistoryS[] = GetPgenXphenID_HostShiftHistory(IDp, ParasIDphen[], IDp_HostShiftHistoryS[i])
                		# This combinaison parasit genotype and phenotype always absent in ParasitesPop
                		push!(X.ParasitesPop.IDpS      , IDp)
                		push!(X.ParasitesPop.IDpSphen  , ParasIDphen[])
                		push!(X.ParasitesPop.IDpGenPhen_HostShiftHistoryS, ParasIDpGenPhen_HostShiftHistoryS[])
                        push!(X.ParasitesPop.IDp_HostShiftHistoryS, IDp_HostShiftHistoryS[i])
                		push!(X.ParasitesPop.N_        , N_[i])
          $(ifelse(GfG,"push!(X.ParasitesPop.rate_of_parasites_emissionS),GetRate_of_parasites_emission(PhenS[i]))",""))
               		# Is this parasit *phenotype* already existing in the MetaCommunity ?
                		posiParasPhen = findfirst(xx -> xx === ParasIDphen[],parasitesGenotypeS_List.IDpS)
                		if isnothing(posiParasPhen) # it does not exist
                		    GetNewParasGenotype!(HPinteractions, parasitesGenotypeS_List, hostsGenotypeS_List, ParasitMutatingMap, tuple(PhenS[i]...), ParasIDphen[])
                		    push!(X.ParasitesPop.posiInParasitGenotypeS_List, length(parasitesGenotypeS_List.IDpS))
                		else
                		    push!(X.ParasitesPop.posiInParasitGenotypeS_List, posiParasPhen)
                		end
            	    end
            	end
            end
        	X.ParasitesPop.N[] = sum(X.ParasitesPop.N_)
            resize!(X.ParasitesPop.ParasitFreq, length(X.ParasitesPop.N_) )
            X.ParasitesPop.ParasitFreq .= X.ParasitesPop.N_ ./ X.ParasitesPop.N[] # RELATIVE frequency of each parasit
        end
        $PRINT println("SetNParasites DONE: ",HostSexS)
        $INBOUNDend
        $FASTMATHend
	end"""
	eval(Meta.parse(str))
else
    str = """
    $INLINE function SetNParasites_HostSexS!(X::HostSexS, HPinteractions::HPinteractions, parasitesGenotypeS_List::ParasitesGenotypeS_List, hostsGenotypeS_List::HostsGenotypeS_List, ParasitMutatingMap::NTuple{Ntraits,Vector{ExtendableDict{Int,Int}}}, mutate::Bool, dt::Float)
    $FASTMATHstart
    $INBOUNDstart
    for sex in 1:2
        X.N_[sex] = X.HostSexS[sex].N[] = sum(X.HostSexS[sex].N_)
    end
    X.N[] = X.N_[1] + X.N_[2]
    if X.N[] > 0
        $PRINT println("SetNParasites: ",HostSexS)
        empty!(X.ParasitesPop.IDp_HostShiftHistoryS)
        empty!(X.ParasitesPop.IDpSphen)
        #
        empty!(X.ParasitesPop.IDpGenPhen_HostShiftHistoryS)
        empty!(X.ParasitesPop.N_                          )
        empty!(X.ParasitesPop.IDpS                        )
        empty!(X.ParasitesPop.posiInParasitGenotypeS_List )
        $(ifelse(GfG
      ,"empty!(X.ParasitesPop.rate_of_parasites_emissionS)",""))
        #
        for sex in 1:2
            posiParasites = findall(X.HostSexS[sex].IDpS .!== Int_IDp(0))
            resize!(X.HostSexS[sex].NparasitesEmitted, length(posiParasites))
            $(ifelse(GfG
            ,"X.HostSexS[sex].NparasitesEmitted .= rand.(Poisson.(X.HostSexS[sex].N_[posiParasites] .* X.HostSexS[sex].rate_of_parasites_emissionS[posiParasites]))"
            ,"X.HostSexS[sex].NparasitesEmitted .= rand.(Poisson.(X.HostSexS[sex].N_[posiParasites] .* $MaxRate_of_parasites_emission))"
            ))
            for (ii,i) in enumerate(posiParasites)
                posiInParasitesPop = findfirst(xx -> xx === X.HostSexS[sex].IDp_HostShiftHistoryS[i],  X.ParasitesPop.IDp_HostShiftHistoryS)
                if isnothing(posiInParasitesPop)
                    push!(X.ParasitesPop.IDp_HostShiftHistoryS      , X.HostSexS[sex].IDp_HostShiftHistoryS[      i])
                    push!(X.ParasitesPop.N_                         , X.HostSexS[sex].NparasitesEmitted[         ii])
                    push!(X.ParasitesPop.posiInParasitGenotypeS_List, X.HostSexS[sex].posiInParasitGenotypeS_List[i])
                    push!(X.ParasitesPop.IDpS                       , X.HostSexS[sex].IDpS[                       i])
      $(ifelse(GfG,"push!(X.ParasitesPop.rate_of_parasites_emissionS, X.HostSexS[sex].rate_of_parasites_emissionS[i])",""))
                else
                    X.ParasitesPop.N_[posiInParasitesPop] += X.HostSexS[sex].NparasitesEmitted[ii]
                end
            end
        end
        #
        if length(X.ParasitesPop.IDp_HostShiftHistoryS) !== 0 # if there is some parasites
            if mutate
                $PRINT println("Entring: MutateParasites!")
                $(ifelse(GfG
                # where MutateParasites! uptates N_IDp_HostShiftHistoryS AND, _only in case of creation of new genotypes_, IDpS as well as IDp_HostShiftHistoryS
                ,"MutateParasites!(HPinteractions, parasitesGenotypeS_List, hostsGenotypeS_List, ParasitMutatingMap, X.ParasitesPop.N_, X.ParasitesPop.IDpS, X.ParasitesPop.IDp_HostShiftHistoryS, X.ParasitesPop.posiInParasitGenotypeS_List ,X.ParasitesPop.Nmutating, X.ParasitesPop.NmutatingEachTrait, X.ParasitesPop.NmuSize, X.ParasitesPop.NmutEachAllelesAppearingTrait, X.ParasitesPop.rate_of_parasites_emissionS, X.ParasitesPop.NewGenotype, X.ParasitesPop.NewId, dt::Float)"
                ,"MutateParasites!(HPinteractions, parasitesGenotypeS_List, hostsGenotypeS_List, ParasitMutatingMap, X.ParasitesPop.N_, X.ParasitesPop.IDpS, X.ParasitesPop.IDp_HostShiftHistoryS, X.ParasitesPop.posiInParasitGenotypeS_List ,X.ParasitesPop.Nmutating, X.ParasitesPop.NmutatingEachTrait, X.ParasitesPop.NmuSize,                                                                                           X.ParasitesPop.NewGenotype, X.ParasitesPop.NewId, dt::Float)"
                ))
                $PRINT println("MutateParasites! DONE")
                # For some very rare parasites genotype, all in individual might have mutated
                Del = findall(X.ParasitesPop.N_ .=== 0)
                if length(Del) > 0
                    deleteat!(X.ParasitesPop.IDpS                        ,Del)
                    deleteat!(X.ParasitesPop.N_                          ,Del)
                    deleteat!(X.ParasitesPop.IDp_HostShiftHistoryS       ,Del)
                    deleteat!(X.ParasitesPop.posiInParasitGenotypeS_List ,Del)
      $(ifelse(GfG,"deleteat!(X.ParasitesPop.rate_of_parasites_emissionS ,Del)",""))
                end
        	end
            append!(X.ParasitesPop.IDpSphen                    , X.ParasitesPop.IDpS)
            append!(X.ParasitesPop.IDpGenPhen_HostShiftHistoryS, GetPgenXphenID_HostShiftHistory.(X.ParasitesPop.IDpS, X.ParasitesPop.IDpSphen, X.ParasitesPop.IDp_HostShiftHistoryS))
        end
        X.ParasitesPop.N[] = sum(X.ParasitesPop.N_)
        resize!(X.ParasitesPop.ParasitFreq, length(X.ParasitesPop.N_) )
        X.ParasitesPop.ParasitFreq .= X.ParasitesPop.N_ ./ X.ParasitesPop.N[] # RELATIVE frequency of each parasit
    end
    $PRINT println("In   SetNParasites_HostSexS!")
    $PRINT println(X.ParasitesPop)
    $INBOUNDend
    $FASTMATHend
    end"""
    eval(Meta.parse(str))
end

##################
##### EVOLVE #####
##################
# UpdateHostAllelesPool_and_Mutate!
# Probilities Geometric distribution
if MuEffectHost*NallelesPerTrait < 1.0 error("MuEffectHost which is $MuEffectHost must be above 1/NallelesPerTrait which is $(1/NallelesPerTrait)") end
p = 1/(MuEffectHost * NallelesPerTrait)
p = map(k -> (1-p)^(k-1)*p , 1:(NallelesPerTrait-1))
p=p/sum(p)

for MUTATE in [false,true]
    strLoop = """
    $(ifelse(MUTATE
    ,"function UpdateHostAllelesPool_and_Mutate!(X::MetaCommunity)"
    ,"function UpdateHostAllelesPool!(           X::MetaCommunity)"
    ))
    $FASTMATHstart
    $INBOUNDstart
        for (sp,x) in enumerate(X.MetaPopS)
            $(ifelse(MUTATE
            ,"UpdateHostAllelesPool_and_Mutate__MetaPop!(x,X.hostsGenotypeS_List, sp)"
            ,"UpdateHostAllelesPool__MetaPop!(           x,X.hostsGenotypeS_List    )"
            ))
        end
        $PRINT println("UpdateHostAllelesPool_and_Mutate! Done")
    $INBOUNDend
    $FASTMATHend
    end"""
    eval(Meta.parse(strLoop))
end
####
for MUTATE in [false,true]
    typeLoop = MetaPop
    strLoop = """$INLINE $(ifelse(MUTATE
    ,"function UpdateHostAllelesPool_and_Mutate__MetaPop!(X::$typeLoop,hostsGenotypeS_List::HostsGenotypeS_List, sp::Int)"
    ,"function UpdateHostAllelesPool__MetaPop!(           X::$typeLoop,hostsGenotypeS_List::HostsGenotypeS_List         )"
    ))
    $FASTMATHstart
    $INBOUNDstart
        $PRINT println("UpdateHostAllelesPool_and_Mutate!", $typeLoop)
        for x in getfield(X,:$(fieldname(typeLoop,1)))
            if x.N[] > 0
                $(ifelse(MUTATE
                ,"UpdateHostAllelesPool_and_Mutate__HostSexS!(x, hostsGenotypeS_List, sp)"
                ,"UpdateHostAllelesPool__HostSexS!(           x, hostsGenotypeS_List    )"
                ))
            end
        end
    $INBOUNDend
    $FASTMATHend
    end"""
    eval(Meta.parse(strLoop))
end
for MUTATE in [false,true]
    typeLoop = HostSexS
    strLoop = """$INLINE $(ifelse(MUTATE
    ,"function UpdateHostAllelesPool_and_Mutate__HostSexS!(X::$typeLoop,hostsGenotypeS_List::HostsGenotypeS_List, sp::Int)"
    ,"function UpdateHostAllelesPool__HostSexS!(           X::$typeLoop,hostsGenotypeS_List::HostsGenotypeS_List         )"))
    $FASTMATHstart
    $INBOUNDstart
        $PRINT println("UpdateHostAllelesPool_and_Mutate!", $typeLoop)
        for x in getfield(X,:$(fieldname(typeLoop,1)))
            $(ifelse(MUTATE
            ,"UpdateHostAllelesPool_and_Mutate__HostsGenotypeS!(x, hostsGenotypeS_List, Val(sp))"
            ,"UpdateHostAllelesPool__HostsGenotypeS!(           x, hostsGenotypeS_List         )"
            ))
        end
    $INBOUNDend
    $FASTMATHend
    end"""
    eval(Meta.parse(strLoop))
end
####
for MUTATE in [false,true]
    for sp in 1:NhostSp
        strLoop = """
        $INLINE $(ifelse(MUTATE
        ,"function UpdateHostAllelesPool_and_Mutate__HostsGenotypeS!(X::HostsGenotypeS, hostsGenotypeS_List::HostsGenotypeS_List, sp::Val{$sp})"
        ,"function UpdateHostAllelesPool__HostsGenotypeS!(           X::HostsGenotypeS, hostsGenotypeS_List::HostsGenotypeS_List              )"
        ))
        $FASTMATHstart
        $INBOUNDstart
            $PRINT println("UpdateHostAllelesPool_and_Mutate! HostsGenotypeS")
                # Reset AllelesPool
                $SIMD for trait in 1:$Ntraits
                    fill!(X.AllelesPool[trait],0.0) # for each trait for each allele number of individual with this allele
                end
                $PRINT println("X.N_ = ",X.N_)
                resize!(X.Ngametes,length(X.N_))
                $PRINT if GfG println("X.rate_of_gametes_prodS = ",X.rate_of_gametes_prodS) end
                $(ifelse(GfG
                    ,"X.Ngametes .= rand.(Poisson.(X.N_ .* X.rate_of_gametes_prodS))"
                    ,"X.Ngametes .= rand.(Poisson.(X.N_ .* $Maxrate_of_gametes_prod))"
                ))
                $PRINT println("X.Ngametes = ",X.Ngametes)
                $SIMD for i in findall(X.Ngametes .!== 0)
                    genotype = hostsGenotypeS_List.GenotypeS[  X.posiInHostsGenotypeS_List[i]  ].Alleles
                    $SIMD for traitPosi in 2:2:$(Ntraits*2)
                        X.AllelesPool[div(traitPosi,2)][genotype[traitPosi-1]$addOneIfGfG_str] += X.Ngametes[i] # Chromosome 1
                        X.AllelesPool[div(traitPosi,2)][genotype[traitPosi  ]$addOneIfGfG_str] += X.Ngametes[i] # Chromosome 2
                    end
                end
                $PRINT println(\"before Mutate   X.AllelesPool = \", X.AllelesPool)
                $(ifelse(MUTATE
                ,"
                # Mutate
                $PRINT println(\"Mutating the host\")
                $(ifelse(GfG
                ,"AlleleRange = $([0,collect(AlleleRangePerSp_infection_success[sp])...])"
                ,"AlleleRange = $(AlleleRangePerSp_infection_success[sp])"))
                $PRINT println(\"AlleleRange = \",AlleleRange)
                for mutatingTrait in 1:$Ntraits
                    if mutatingTrait === $(N_traits_infection_success+1) # change the AlleleRange once the traits_infection_success are done
                        $(ifelse(GfG
                        ,"AlleleRange = 0:$NallelesPerTrait"
                        ,"AlleleRange = 1:$NallelesPerTrait"))
                    end
                    $PRINT println(\"AlleleRange = \",AlleleRange)
                    fill!(X.MutatedAlleles[mutatingTrait],0)
                    for allele in AlleleRange
                        $PRINT println(\"allele = \",allele)
                        if X.AllelesPool[mutatingTrait][allele$addOneIfGfG_str] > 0.0
                            # Here we mutate a proportion of the allele pool whatever the value of dt because if dt only a small proportion of the allele pool will be used by HaveSex!, and this portion will be function of dt
                            $(ifelse(GfG
                            ,"if allele === 0 # the trait can only appear
                                  X.nMut[1] = X.nMut[2] = 0
                                  X.nMut[3] = rand(Binomial(Int(X.AllelesPool[mutatingTrait][allele$addOneIfGfG_str]), $(MuHost_trait_presence[sp]))) # total number of mutants === number of mutants for which the trait appear
                              else
                                  X.nMut[1]  = rand(Binomial(Int(X.AllelesPool[mutatingTrait][allele$addOneIfGfG_str]), $(MuHost_trait_alleles[sp] + MuHost_trait_presence[sp]) )) # total number of mutants
                                  X.nMut[3]  = rand(Binomial(X.nMut[1] , $(  MuHost_trait_presence[sp] / (MuHost_trait_alleles[sp] + MuHost_trait_presence[sp])))) # number of mutants for which the trait disappear
                                  X.nMut[1] -= X.nMut[3] # number of mutants that have the mutatingTrait present and that mutate it
                                  X.nMut[2]  = rand(Binomial(X.nMut[1])) # number of mutants that decrease
                                  X.nMut[1] -= X.nMut[2] # number of mutants that increase
                              end"
                            ,"X.nMut[1]  = rand(Binomial(Int(X.AllelesPool[mutatingTrait][allele$addOneIfGfG_str]), $(MuHost_trait_alleles[sp]))) # total number of mutants
                              X.nMut[2]  = rand(Binomial(X.nMut[1])) # number of mutants that decrease
                              X.nMut[1] -= X.nMut[2] # number of mutants that increase
                              "
                            ))
                            $PRINT println(\"X.nMut = \",X.nMut)
                            $PRINT TotalAdded = 0
                            $(ifelse(GfG, "if allele !==0", ""))
                                # Decrease
                                $PRINT println(\"Decrease\")
                                $PRINT println(\"rand!(Multinomial(X.nMut[1], $p), X.AllelesChanges)\")
                                rand!(Multinomial(X.nMut[1], $p), X.AllelesChanges) # Size of mutations in a geometric distribution
                                for musize in 1:$(NallelesPerTrait-1)
                                    if X.AllelesChanges[musize] > 0
                                        $PRINT println(\"allele = \",allele)
                                        $PRINT println(\"musize = \",musize)
                                        $PRINT println(\"mutatingTrait = \",mutatingTrait)
                                        if mutatingTrait <= $N_traits_infection_success
                                            if (allele-musize) >= $(AlleleRangePerSp_infection_success[sp][1]) # just apply the mutation
                                                X.MutatedAlleles[mutatingTrait][allele - musize $addOneIfGfG_str] += X.AllelesChanges[musize]
                                                $PRINT TotalAdded += X.AllelesChanges[musize]
                                            else # don't apply the mutation (reflexive boundaries)
                                                X.nMut[1] -= X.AllelesChanges[musize]
                                            end
                                        else
                                            if (allele-musize) >= 1 # just apply the mutation
                                                X.MutatedAlleles[mutatingTrait][allele - musize $addOneIfGfG_str] += X.AllelesChanges[musize]
                                                $PRINT TotalAdded += X.AllelesChanges[musize]
                                            else # circular trait
                                                X.MutatedAlleles[mutatingTrait][allele - musize + $(NallelesPerTrait + if_GfG_1_else_0) ] += X.AllelesChanges[musize]
                                                $PRINT TotalAdded += X.AllelesChanges[musize]
                                            end
                                        end
                                    end
                                end
                                # Increase
                                $PRINT println(\"Increase\")
                                $PRINT println(\"rand!(Multinomial(X.nMut[2], $p), X.AllelesChanges)\")
                                rand!(Multinomial(X.nMut[2], $p), X.AllelesChanges) # Size of mutations in a geometric distribution
                                $PRINT println(\"X.AllelesChanges = \", X.AllelesChanges)
                                for musize in 1:$(NallelesPerTrait-1)
                                    if X.AllelesChanges[musize] > 0
                                        $PRINT println(\"allele = \",allele)
                                        $PRINT println(\"musize = \",musize)
                                        $PRINT println(\"mutatingTrait = \",mutatingTrait)
                                        if mutatingTrait <= $N_traits_infection_success
                                            if (allele+musize) <= $(AlleleRangePerSp_infection_success[sp][end]) # just apply the mutation
                                                X.MutatedAlleles[mutatingTrait][allele + musize $addOneIfGfG_str] += X.AllelesChanges[musize]
                                                $PRINT TotalAdded += X.AllelesChanges[musize]
                                            else # don't apply the mutation (reflexive boundaries)
                                                X.nMut[2] -= X.AllelesChanges[musize]
                                            end
                                        else
                                            if (allele+musize) <= $(NallelesPerTrait) # just apply the mutation
                                                X.MutatedAlleles[mutatingTrait][allele + musize $addOneIfGfG_str] += X.AllelesChanges[musize]
                                                $PRINT TotalAdded += X.AllelesChanges[musize]
                                            else # circular trait              #                -                    -  = +
                                                X.MutatedAlleles[mutatingTrait][allele + musize - $(NallelesPerTrait - if_GfG_1_else_0) ] += X.AllelesChanges[musize]
                                                $PRINT TotalAdded += X.AllelesChanges[musize]
                                            end
                                        end
                                    end
                                end
                            $(ifelse(GfG, "end", ""))
                            $(ifelse(GfG
                            ,"# Swapping the trait presence
                            if allele === 0 # the trait appear
                                X.MutatedAlleles[mutatingTrait][AlleleRange[2:end]] .+= rand(Multinomial(X.nMut[3],( fill(1/(length(AlleleRange)-1), length(AlleleRange)-1 ))))
                                $PRINT TotalAdded += X.nMut[3]
                            else # the trait disappear
                                X.MutatedAlleles[mutatingTrait][1] += X.nMut[3]
                                $PRINT TotalAdded += X.nMut[3]
                            end"
                            ,""
                            ))
                            X.MutatedAlleles[mutatingTrait][allele $addOneIfGfG_str] -= (X.nMut[1] + X.nMut[2]$(ifelse(GfG," + X.nMut[3]",""))) # remove the mutants of the current allele
                            $PRINT println(\"TotalAdded =\", TotalAdded)
                            $PRINT println(\"sum(X.nMut) =\", sum(X.nMut))
                            $PRINT if sum(X.MutatedAlleles[mutatingTrait]) !==0 println(\"mutatingTrait = \", mutatingTrait) ; println(\"allele = \", allele) ; global Xerror = X ; error(\"In UpdateHostAllelesPool_and_Mutate!(X::HostsGenotypeS,hostsGenotypeS_List::HostsGenotypeS_List)  :   sum(X.MutatedAlleles[mutatingTrait]) !== 0 \") end
                        end
                    end
                    X.AllelesPool[mutatingTrait] .+= X.MutatedAlleles[mutatingTrait] # where X.MutatedAlleles have positive and negative values that sum to zero.
                end
                "
                ,""
                ))
                # Convert to frequencies
                $PRINT if length(unique(sum.(X.AllelesPool)))>1 global Xerror = X ; error("In UpdateHostAllelesPool_and_Mutate!(X::HostsGenotypeS,hostsGenotypeS_List::HostsGenotypeS_List)  :   length(unique(sum.(X.AllelesPool)))>0") end
                $SIMD for trait in 1:$(Ntraits)
                    $PRINT println(\"before ./= sum    X.AllelesPool[trait] = \", X.AllelesPool[trait])
                    Tot = sum(X.AllelesPool[trait])
                    if (Tot > 0) X.AllelesPool[trait] ./= Tot end
                    $PRINT println(\"after ./= sum    X.AllelesPool[trait] = \", X.AllelesPool[trait])
                end
            $INBOUNDend
            $FASTMATHend
            end"""
        eval(Meta.parse(strLoop))
    end
end
####
# print(str)


##################
# HaveSex!
type = MetaCommunity
str = """
function HaveSex!(X::$type,dt::Float)
    $FASTMATHstart
    $INBOUNDstart
    $PRINT println("HaveSex: ",$type)
    $PRINT println("UpdateHostAllelesPool_and_Mutate")
    UpdateHostAllelesPool_and_Mutate!(X)
    for (sp, x) in enumerate(getfield(X,:$(fieldname(type,1))))
        HaveSex_MetaPop!(x,X,dt,sp)
    end
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))
##
type = MetaPop
str = """
$INLINE function HaveSex_MetaPop!(X::$type,MetaComm::MetaCommunity,dt::Float,sp::Int)
    $FASTMATHstart
    $INBOUNDstart
    $PRINT println("HaveSex!: MetaPop")
    for (pop,x) in enumerate(getfield(X,:$(fieldname(type,1))))
        if all(x.N_ .> 0) # if there is both males and females in the population
            HaveSex_HostSexS!(x,MetaComm,dt,pop)
        end
    end
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))
##
str = """
$INLINE function HaveSex_HostSexS!(X::HostSexS, MetaComm::MetaCommunity, dt::Float, pop::Int)
    $FASTMATHstart
    $INBOUNDstart
    $PRINT println("HaveSex!: HostSexS")
    E_Noff::Float = X.N[] * ($b0-X.density_dependence[] × X.N[]) * dt / 2 # expected number of sons and daughters by individuals, thus divide by two
    if E_Noff > 0.0
        Noff::Vector{Int} = rand(Poisson( E_Noff ),2)
        $PRINT println("Noff", Noff)
        for sex in 1:2 # sons (1) then daughters (2)
            if Noff[sex] !== 0
                $PRINT println("sex ", sex)
                N_geno, Genotypes = RecursiveRandMultinom(Noff[sex],
                    [$(join(["X.HostSexS[1].AllelesPool[$trait],X.HostSexS[2].AllelesPool[$trait]" for trait in 1:Ntraits],",  "))]
                )
                $PRINT println("N_geno =",N_geno)
                $PRINT println("Genotypes =",Genotypes)
                # Add offspring
                $PRINT println("Add offspring")
                for i in eachindex(N_geno)
                    $(ifelse(GfG
                    ,"Genotypes[i] .-= 1"
                    ,""
                    ))
                    $SIMD for posi in (1:2:$(2*Ntraits)) # order alleles
                        if Genotypes[i][posi] > Genotypes[i][posi+1]
                            Genotypes[i][posi:(posi+1)] = Genotypes[i][(posi+1):-1:posi]
                        end
                    end
                    $PRINT println("GetOffspring! in HaveSex!")
                    GetOffspring!(MetaComm, X.HostSexS[sex], tuple(Genotypes[i]...), N_geno[i])
                end
            end
        end
    end
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

#################
# Die!
str = """
function Die!(X::MetaCommunity,dt::Float)
    $FASTMATHstart
    $INBOUNDstart
    $PRINT println("Die! ",typeof(X))
    for (sp, x) in enumerate(getfield(X, :$(fieldname(MetaCommunity,1))))
        Die_MetaPop!(x,dt,sp)
    end
    $PRINT println("Die! ",typeof(X),"Done")
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

str = """
$INLINE function Die_MetaPop!(X::MetaPop, dt::Float, sp::Int)
    $FASTMATHstart
    $INBOUNDstart
    $PRINT println("Die! ",typeof(X))
    for (pop,x) in enumerate(getfield(X,:$(fieldname(MetaPop,1))))
        if x.N[] > 0
            Die_HostSexS!(x,dt,pop)
        end
    end
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

str = """
$INLINE function Die_HostSexS!(X::HostSexS, dt::Float, pop::Int)
    $FASTMATHstart
    $INBOUNDstart
    $PRINT println("Die! ",typeof(X))
    NdeadTot::Ref{Int} = Ref(0)
    for (sex, x) in enumerate(getfield(X,:$(fieldname(HostSexS,1))))
        Die_HostsGenotypeS!(x, X.density_dependence[],X.N[] , sex, dt)
    end
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

str = """
$INLINE function Die_HostsGenotypeS!(X::HostsGenotypeS, density_dependence::Float, NindPop::Int, sex::Int, dt::Float)
    $FASTMATHstart
    $INBOUNDstart
    $PRINT println("Die! ",typeof(X))
    p = min.(((X.virulenceS .+ ($d0 + density_dependence × NindPop)).*dt),1.0) # for each category, for each individual, probability to die  (d0+D×N+virulence)×dt
    $PRINT println("p ",p)
    $PRINT println("X.N_ =",X.N_)
    Ndead_ = rand.(Binomial.(X.N_,p))
    $PRINT println(2)
    X.dN_ .-= Ndead_
    $PRINT println(3)
    $PRINT println("Die! ",typeof(X),"Done")
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

##################
# InfectHosts

# str = """
# function GetPinfection!(Pinfection1contact::Float,Ncontacts::Int,N::Int)
#     sum([((1-(1-Pinfection1contact)^n)   *
#         pdf(Binomial(Ncontacts  ,1/N ),n)
#         # Float( binomial(big(Ncontacts),big(n))×(1/N)^n×(1-1/N )^(Ncontacts-n))
#         )
#     for n in 1:min($NcoinfectionsMax,Ncontacts)])
# end"""
# eval(Meta.parse(str))

# Get all Float16 values from 0 to 1
try # only create  P_0_1_16f  if it does not already exists
    P_0_1_16f[1]
catch
    global P_0_1_16f = Vector{Float16}()
    global x = 0.0
    while x < 1.1
        push!(P_0_1_16f,Float16(x))
        global x += 0.003/(2^16)
    end
    P_0_1_16f = unique(P_0_1_16f)
    # P_0_1_16f = (P_0_1_16f...,) ;
end

function GetPinfection_pdfBinomialPart___ParasitFreq_x_Pinfection1contact__p(ParasitFreq_x_Pinfection1contact__p::Float16, Max_NcontactsMax::Int)
    out::Vector{Float} = Vector{Float}(undef,Max_NcontactsMax)
    for n in 1:Max_NcontactsMax
        out[n] = ( 1 - pdf(Binomial(n, Float(ParasitFreq_x_Pinfection1contact__p)),0))
    end
    return(out)
end
function GetPinfection_pdfBinomialPart___P_HigherVirulence(P_HigherVirulence::Float16, Max_NcontactsMax::Int)
    out::Vector{Float} = Vector{Float}(undef,Max_NcontactsMax)
    for n in 1:Max_NcontactsMax
        out[n] = (pdf(Binomial((n-1), Float(P_HigherVirulence) ),0 ))
    end
    return(out)
end
# Float( binomial(big((n-1)),big(nn))×(P_SameVirulence)^nn×(1-P_SameVirulence )^((n-1)-nn))
function GetPinfection_pdfBinomialPart___P_SameVirulence(P_SameVirulence::Float16, Max_NcontactsMax::Int)
    temp::Vector{Float} = Vector{Float}(undef,Max_NcontactsMax)
    out::Vector{Float} = Vector{Float}(undef,Max_NcontactsMax)
    for n in 1:Max_NcontactsMax
        map!(nn -> pdf(Binomial((n-1), Float(P_SameVirulence) ), nn) / (nn+1) ,temp,1:n)
        out[n] = sum(temp[1:n])
    end
    return(out)
end

function GetAllComb_GetPinfection_pdfBinomialPart(N::Int,P::Vector{Float16})
    # ParasitFreq_x_Pinfection1contact__p, P_SameVirulence, P_HigherVirulence = Dict{Float16, NTuple{N,Float}}(), Dict{Float16, NTuple{N,Float}}(), Dict{Float16, NTuple{N,Float}}()
    ParasitFreq_x_Pinfection1contact__p::Dict{Float16, NTuple{N,Float64}}, P_SameVirulence::Dict{Float16, NTuple{N,Float64}}, P_HigherVirulence::Dict{Float16, NTuple{N,Float64}}   =   Dict{Float16, NTuple{N,Float64}}(), Dict{Float16, NTuple{N,Float64}}(), Dict{Float16, NTuple{N,Float64}}()
    # ParasitFreq_x_Pinfection1contact__p::Dict{Float16, Vector{Float}}, P_SameVirulence::Dict{Float16, Vector{Float}}, P_HigherVirulence::Dict{Float16, Vector{Float}}   =   Dict{Float16, Vector{Float}}(), Dict{Float16, Vector{Float}}(), Dict{Float16, Vector{Float}}()
    @async @inbounds Threads.@threads for p in P
        ParasitFreq_x_Pinfection1contact__p[p] = tuple(GetPinfection_pdfBinomialPart___ParasitFreq_x_Pinfection1contact__p(min(p,one(Float16)),N)...) ;
        P_SameVirulence[p]                     = tuple(GetPinfection_pdfBinomialPart___P_SameVirulence(min(p,one(Float16)),N)...) ;
        P_HigherVirulence[p]                   = tuple(GetPinfection_pdfBinomialPart___P_HigherVirulence(min(p,one(Float16)),N)...) ;
        # ParasitFreq_x_Pinfection1contact__p[p] = GetPinfection_pdfBinomialPart___ParasitFreq_x_Pinfection1contact__p(min(p,one(Float16)),N)
        # P_SameVirulence[p]                     = GetPinfection_pdfBinomialPart___P_SameVirulence(min(p,one(Float16)),N)
        # P_HigherVirulence[p]                   = GetPinfection_pdfBinomialPart___P_HigherVirulence(min(p,one(Float16)),N)
    end
    return( ParasitFreq_x_Pinfection1contact__p, P_SameVirulence, P_HigherVirulence )
end

temp = load("./p_Float16_.jld")
ParasitFreq_x_Pinfection1contact__p_Float16_ = temp["ParasitFreq_x_Pinfection1contact__p_Float16_"] ;
P_SameVirulence_Float16_ = temp["P_SameVirulence_Float16_"] ;
P_HigherVirulence_Float16_ = temp["P_HigherVirulence_Float16_"] ;

temp = try # only create  ParasitFreq_x_Pinfection1contact__p_Float16 ECT  if they do not already exists
    ParasitFreq_x_Pinfection1contact__p_Float16_[zero(Float16)]
catch end
if isnothing(temp) || (round(maximum(Ht)*10,digits=0) > length(ParasitFreq_x_Pinfection1contact__p_Float16_[zero(Float16)]))
    ParasitFreq_x_Pinfection1contact__p_Float16_, P_SameVirulence_Float16_, P_HigherVirulence_Float16_ = GetAllComb_GetPinfection_pdfBinomialPart(Int(maximum(Ht)*10),P_0_1_16f) ;
    save("/DISQUE/0_Grenoble/Biodiv_zoonoses/p_Float16_.jld", "ParasitFreq_x_Pinfection1contact__p_Float16_", ParasitFreq_x_Pinfection1contact__p_Float16_, "P_SameVirulence_Float16_", P_SameVirulence_Float16_, "P_HigherVirulence_Float16_", P_HigherVirulence_Float16_)
end

const ParasitFreq_x_Pinfection1contact__p_Float16 = ParasitFreq_x_Pinfection1contact__p_Float16_ ;
const P_SameVirulence_Float16 = P_SameVirulence_Float16_ ;
const P_HigherVirulence_Float16 = P_HigherVirulence_Float16_ ;


str = """
$INLINE function GetPinfection!(ParasitFreq_x_Pinfection1contact__p::Float16, P_SameVirulence::Float16, P_HigherVirulence::Float16, P_NcontactsInPop_1_2_3_ect_parasites::Vector{Float16}, NcontactsMax::Int)
    $FASTMATHstart
    $INBOUNDstart
    return(
    sum(view(P_NcontactsInPop_1_2_3_ect_parasites,1:NcontactsMax) .* view(ParasitFreq_x_Pinfection1contact__p_Float16[ParasitFreq_x_Pinfection1contact__p],1:NcontactsMax) .* ( view(P_SameVirulence_Float16[P_HigherVirulence],1:NcontactsMax)     .+ view(P_SameVirulence_Float16[P_SameVirulence],1:NcontactsMax) ))
    # sum(P_NcontactsInPop_1_2_3_ect_parasites[1:NcontactsMax] .* ParasitFreq_x_Pinfection1contact__p_Float16[ParasitFreq_x_Pinfection1contact__p][1:NcontactsMax] .* ( P_SameVirulence_Float16[P_HigherVirulence][1:NcontactsMax]     .+ P_SameVirulence_Float16[P_SameVirulence][1:NcontactsMax] ))
    )
    # # # temp = Vector{Float}(undef,NcontactsMax)
    # # Tot = 0.0
    # # for n in 1:NcontactsMax
    # #     # map!(nn -> pdf(Binomial((n-1), P_SameVirulence ), nn) / (nn+1) ,temp,1:n)
    # #     # Tot += P_NcontactsInPop_1_2_3_ect_parasites[n] * ( 1 - pdf(Binomial(n, ParasitFreq_x_Pinfection1contact__p),0))                      * ((pdf(Binomial((n-1), P_HigherVirulence ),0 )) + sum(temp[1:n])                                      )
    # #       Tot += P_NcontactsInPop_1_2_3_ect_parasites[n] * ParasitFreq_x_Pinfection1contact__p_Float16[ParasitFreq_x_Pinfection1contact__p][n] * ( P_SameVirulence_Float16[P_HigherVirulence][n]      + P_SameVirulence_Float16[P_SameVirulence][n] )
    # # end
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))


# @btime GetPinfection_!($Float16(0.564), $Float16(0.564), $Float16(0.564), $[Float16(0.564),Float16(0.564),Float16(0.564),Float16(0.564),Float16(0.564),Float16(0.564)], 6);
# @btime GetPinfection!( $Float16(0.564), $Float16(0.564), $Float16(0.564), $[Float16(0.564),Float16(0.564),Float16(0.564),Float16(0.564),Float16(0.564),Float16(0.564)], 6);


# for NcontactsMax in 1:Int(maximum(Ht)*10) # Max_NcontactsMax
#     str = """
#     $NOINLINE function GetPinfection!(ParasitFreq_x_Pinfection1contact__p::Float16, P_SameVirulence::Float16, P_HigherVirulence::Float16, P_NcontactsInPop_1_2_3_ect_parasites::Vector{Float16}, NcontactsMax::Val{$NcontactsMax})
#         $INBOUNDstart
#         return(
#             $(join(["P_NcontactsInPop_1_2_3_ect_parasites[$Ncontacts] * ParasitFreq_x_Pinfection1contact__p_Float16[ParasitFreq_x_Pinfection1contact__p][$Ncontacts] * ( P_SameVirulence_Float16[P_HigherVirulence][$Ncontacts] + P_SameVirulence_Float16[P_SameVirulence][$Ncontacts] )" for Ncontacts in 1:NcontactsMax], " + "))
#         )
#         # (
#         # sum(P_NcontactsInPop_1_2_3_ect_parasites[1:NcontactsMax] .* ParasitFreq_x_Pinfection1contact__p_Float16[ParasitFreq_x_Pinfection1contact__p][1:NcontactsMax] .* ( P_SameVirulence_Float16[P_HigherVirulence][1:NcontactsMax]     .+ P_SameVirulence_Float16[P_SameVirulence][1:NcontactsMax] ))
#         # )
#         # # # temp = Vector{Float}(undef,NcontactsMax)
#         # # Tot = 0.0
#         # # for n in 1:NcontactsMax
#         # #     # map!(nn -> pdf(Binomial((n-1), P_SameVirulence ), nn) / (nn+1) ,temp,1:n)
#         # #     # Tot += P_NcontactsInPop_1_2_3_ect_parasites[n] * ( 1 - pdf(Binomial(n, ParasitFreq_x_Pinfection1contact__p),0))                      * ((pdf(Binomial((n-1), P_HigherVirulence ),0 )) + sum(temp[1:n])                                      )
#         # #       Tot += P_NcontactsInPop_1_2_3_ect_parasites[n] * ParasitFreq_x_Pinfection1contact__p_Float16[ParasitFreq_x_Pinfection1contact__p][n] * ( P_SameVirulence_Float16[P_HigherVirulence][n]      + P_SameVirulence_Float16[P_SameVirulence][n] )
#         # # end
#         $INBOUNDend
#     end"""
#     eval(Meta.parse(str))
# end


str = """
function InfectHosts!(X::MetaCommunity,dt::Float)
    $FASTMATHstart
    $INBOUNDstart
    SetNParasites!(X,true,dt) # true means that mutation happen
    $PRINT println("InfectHosts : ", MetaCommunity)
    for (spRecipient, x) in enumerate(getfield(X,:$(fieldname(MetaCommunity,1))))
        InfectHosts_MetaPop!(x, X , X.hostsGenotypeS_List, $(ifelse(SimulateImmunisation,"X.parasitesGenotypeS_List, ",""))spRecipient, dt)
    end
    $PRINT println("SetNParasites")
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

str = """
$INLINE function InfectHosts_MetaPop!(X::MetaPop, MetaComm::MetaCommunity , hostsGenotypeS_List::HostsGenotypeS_List, $(ifelse(SimulateImmunisation,"parasitesGenotypeS_List::ParasitesGenotypeS_List, ",""))spRecipient::Int, dt::Float)
    $FASTMATHstart
    $INBOUNDstart
    $PRINT println("InfectHosts : ", MetaPop)
    for (pop, x) in enumerate(getfield(X,:$(fieldname(MetaPop,1))))
        $PRINT if x.ParasitesPop.N[] !== sum(x.ParasitesPop.N_) global xx_ = x ; global XX_ = X ; error("x.ParasitesPop.N[] !== sum(x.ParasitesPop.N_)") end
        $PRINT println("x.N[] = ",x.N[])
        if x.N[] > 0
            if x.ParasitesPop.N[] > 0
                $PRINT println("InfectHosts : ", HostSexS)
                InfectHosts_HostSexS!(x, MetaComm.HPinteractions, hostsGenotypeS_List, $(ifelse(SimulateImmunisation,"parasitesGenotypeS_List, ",""))              x.ParasitesPop,                                                                          dt)
            end
            for spDonor in 1:$NhostSp
                if spDonor !== spRecipient
                    DonorPop::HostSexS = MetaComm.MetaPopS[spDonor].PopS[pop]
                    if DonorPop.ParasitesPop.N[] > 0
                        $PRINT println("InfectHostsShift_HostSexS")
                        InfectHostsShift_HostSexS!(x, MetaComm.HPinteractions, hostsGenotypeS_List, $(ifelse(SimulateImmunisation,"parasitesGenotypeS_List, ","")) DonorPop.ParasitesPop, DonorPop.N[], DonorPop.K[], spRecipient, spDonor, MetaComm.Storage.TimeEvolved[], dt)
                    end
                end
            end
        end
    end
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

for HostShift in [true,false]
    strLoop = """
    $NOINLINE $(ifelse(HostShift
        ,"function InfectHostsShift_HostSexS!(X::HostSexS, HPinteractions::HPinteractions, hostsGenotypeS_List::HostsGenotypeS_List, $(ifelse(SimulateImmunisation,"parasitesGenotypeS_List::ParasitesGenotypeS_List, ",""))parasitesPop::ParasitesPop, Nd::Int, Kd::Float, spRecipient::Int, spDonor::Int, Time::Float, dt::Float)"
        ,"function InfectHosts_HostSexS!(     X::HostSexS, HPinteractions::HPinteractions, hostsGenotypeS_List::HostsGenotypeS_List, $(ifelse(SimulateImmunisation,"parasitesGenotypeS_List::ParasitesGenotypeS_List, ",""))parasitesPop::ParasitesPop,                                                                  dt::Float)"  ))
        $FASTMATHstart
        $INBOUNDstart
        $PRINT println("X = ")
        $PRINT println.(X.HostSexS)
        $PRINT println("X.N[] = ",X.N[])
        $PRINT println("parasitesPop.N[] = ",parasitesPop.N[])
        $PRINT println("parasitesPop = ")
        $PRINT println("typeof(parasitesPop) = ", typeof(parasitesPop))
        $PRINT println(parasitesPop)
        NparasCategories = length(parasitesPop.ParasitFreq)
        $PRINT println("NparasCategories = ",NparasCategories)
        #
        $(ifelse(HostShift
        ,"
        $PRINT println(\"Nd = \",Nd)
        $PRINT println(\"(parasitesPop.N[]/Nd) * X.N[] * X.HtCrossSp[] *(X.SpInteract[] + Nd/Kd*X.N[]/X.K[]*(1-X.SpInteract[])) = \",(parasitesPop.N[]/Nd) * X.N[] * X.HtCrossSp[] *(X.SpInteract[] + Nd/Kd*X.N[]/X.K[]*(1-X.SpInteract[])))
        X.NcontactsInPop[] = # number of contacts during dt whatever the sex of individuals
        rand(Poisson( (parasitesPop.N[]/Nd) * X.N[] * X.HtCrossSp[] *(X.SpInteract[] + Nd/Kd*X.N[]/X.K[]*(1-X.SpInteract[]))
        *dt
        )) # For each host genotype total number of contacts with a parasit it recieves
        "
        ,"X.NcontactsInPop[] = # number of contacts during dt whatever the sex of individuals
        rand(Poisson( ((X.Ht_Sociality + X.Ht_one_Sociality_K[] * X.N[]) # H_T×Sociality + (N×H_T×(1-Sociality))/K
        *dt * parasitesPop.N[]
        ))) # For each host genotype total number of contacts with a parasit it recieves
        "))
        if X.NcontactsInPop[] !== 0
            $(ifelse(HostShift
            ,"X.NcontactsCrossSp[spDonor] += X.NcontactsInPop[]"
            ,""))
            $PRINT println("Poisson X.NcontactsInPop = ",X.NcontactsInPop)
            # Some hosts will undergo several contacts
            resize!(X.ParasitFreq_x_Pinfection1contact , NparasCategories)
            resize!(X.virulenceScompar                 , NparasCategories)
            resize!(X.Pinfection                       , NparasCategories)
            resize!(X.NsuccessfullInfectionEachParasit , NparasCategories)
            resize!(X.P_HigherVirulenceS               , NparasCategories)
            resize!(X.P_SameVirulenceS                 , NparasCategories)
            #
            empty!(X.P_NcontactsInPop_1_2_3_ect_parasites)
            one_Nhostes::Float = 1/X.N[]
            $PRINT println("one_Nhostes = ", one_Nhostes)
            # Approximation of the probability that a host undergoes 1, 2, 3 ,ect , contacts with parasites (X.P_NcontactsInPop_1_2_3_ect_parasites >~ thresholdApproxBinomial)
            n = 0 ; Tot = pdf(Binomial(X.NcontactsInPop[]  ,one_Nhostes ),0) # Skip the probability that n = 0
            while (Tot < $thresholdApproxBinomial) || (length(X.P_NcontactsInPop_1_2_3_ect_parasites) < 2)
                n += 1
                temp::Float = pdf(Binomial(X.NcontactsInPop[]  , one_Nhostes),n)
                Tot += temp
                push!(X.P_NcontactsInPop_1_2_3_ect_parasites, Float16(temp))
            end
            if (Tot > one(Float16)) X.P_NcontactsInPop_1_2_3_ect_parasites ./= sum(X.P_NcontactsInPop_1_2_3_ect_parasites) end
            $PRINT if length(X.P_NcontactsInPop_1_2_3_ect_parasites) === 0 global X_ = X ; println("X.NcontactsInPop = ",X.NcontactsInPop) ; error("length(X.P_NcontactsInPop_1_2_3_ect_parasites) === 0") end
            #
            NcontactsMax = length(X.P_NcontactsInPop_1_2_3_ect_parasites)
            $PRINT println("NcontactsMax = ",NcontactsMax)
            for sex in 1:2
                x = X.HostSexS[sex]
                # println("x.N[] = ",x.N[])
                # Start = x.N[]
                for h in eachindex(x.N_) # For each host genotype
                    $PRINT println("sex = ",sex," h = ",h)
                    $PRINT println("length(X.ParasitFreq_x_Pinfection1contact) = ",length(X.ParasitFreq_x_Pinfection1contact))
                    $PRINT println("parasitesPop.ParasitFreq = ",parasitesPop.ParasitFreq)
                    X.ParasitFreq_x_Pinfection1contact .= parasitesPop.ParasitFreq .* HPinteractions.Pinfection1contact[     hostsGenotypeS_List.GenotypeS[x.posiInHostsGenotypeS_List[h]].traits_infection_success_Symb    ][    parasitesPop.posiInParasitGenotypeS_List    ]
                    $PRINT println("X.ParasitFreq_x_Pinfection1contact = ",X.ParasitFreq_x_Pinfection1contact)
                    # order the virulence according to parasitesPop
                    $PRINT println("length(X.virulenceScompar) = ",length(X.virulenceScompar))
                    X.virulenceScompar                 .=                             HPinteractions.virulenceS[             hostsGenotypeS_List.GenotypeS[x.posiInHostsGenotypeS_List[h]].traits_infection_success_Symb    ][    parasitesPop.posiInParasitGenotypeS_List    ]
                    $PRINT println("X.virulenceScompar = ",X.virulenceScompar)
                    $PRINT println("1")
                    #
                    # Handle multiple infections ############################
                    P_SameVirulence = P_HigherVirulence = zero(Float16) # ... than the current virulence induced by the current parasit of the host (if any)
                    empty!(X.Posi_HigherOrEqualVirulenceS_comparToVirulenceInHost)
                    fill!(X.P_HigherVirulenceS, Float16(0.0)) # ... than the other parasite that will coinfect the same host individual
                    fill!(X.P_SameVirulenceS  , Float16(0.0)) # ...   as the other parasite that will coinfect the same host individual
                    for p in eachindex(X.virulenceScompar) # focal parasit
                        for p2 in eachindex(X.virulenceScompar) # eventual co-infecting parasites
                            if X.virulenceScompar[p2] === X.virulenceScompar[p]
                                X.P_SameVirulenceS[p] += X.ParasitFreq_x_Pinfection1contact[p2]
                            elseif X.virulenceScompar[p2] > X.virulenceScompar[p]
                                X.P_HigherVirulenceS[p] += X.ParasitFreq_x_Pinfection1contact[p2]
                            end
                        end
                        if X.virulenceScompar[p] === x.virulenceS[h]
                            P_SameVirulence += X.ParasitFreq_x_Pinfection1contact[p]
                            push!(X.Posi_HigherOrEqualVirulenceS_comparToVirulenceInHost,p)
                        elseif X.virulenceScompar[p] > x.virulenceS[h]
                            P_HigherVirulence += X.ParasitFreq_x_Pinfection1contact[p]
                            push!(X.Posi_HigherOrEqualVirulenceS_comparToVirulenceInHost,p)
                        # else # we are in the case : X.virulenceScompar[p] < x.virulenceS[h]
                        #     push!(X.Posi_HigherOrEqualVirulenceS_comparToVirulenceInHost,p)
                        end
                    end
                    $PRINT println("X.P_HigherVirulenceS = ",X.P_HigherVirulenceS)
                    $PRINT println("X.P_SameVirulenceS = ",X.P_SameVirulenceS)
                    $PRINT println("X.P_HigherVirulenceS = ",X.P_HigherVirulenceS)
                    $PRINT println("NcontactsMax = ",NcontactsMax)
                    $PRINT println("X.ParasitFreq_x_Pinfection1contact = ",X.ParasitFreq_x_Pinfection1contact)
                    $PRINT println("X.P_NcontactsInPop_1_2_3_ect_parasites = ",X.P_NcontactsInPop_1_2_3_ect_parasites)#
                    # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                    Pinfected_for_each_host_indiv = sum( X.P_NcontactsInPop_1_2_3_ect_parasites .* [(1-(1-(   P_HigherVirulence   +   P_SameVirulence/2 ))^n) for n in 1:NcontactsMax] ) # For individuals of genotype 'h', probability for each of them to became infected
                    #                                    Proba of 1,2,3 or N contacts           x  |  |  |Proba that an infecting parasit is successful|| ||
                    #                                    Proba of 1,2,3 or N contacts           x  |  |Proba that an infecting parasit fails ---------->| ||
                    #                                    Proba of 1,2,3 or N contacts           x  |  |Proba that the n parasites that have a contact fail>||
                    #                                    Proba of 1,2,3 or N contacts           x  |Proba that at least one parasit succeed -------------->|
                    # _____________________________________________________________________________________________________________________________________|
                    $PRINT if Pinfected_for_each_host_indiv > 1.0 ; global x_ = x ; global X_ = X ; println("X.P_NcontactsInPop_1_2_3_ect_parasites = ",X.P_NcontactsInPop_1_2_3_ect_parasites) ; println("P_HigherVirulence = ",P_HigherVirulence) ; println("P_SameVirulence = ",P_SameVirulence) ; println("X.ParasitFreq_x_Pinfection1contact = ",X.ParasitFreq_x_Pinfection1contact) ; println("NcontactsMax = ",NcontactsMax) ; error("Pinfected_for_each_host_indiv > 1.0"); end
                    N_hosts_infected_current_genotype = rand(Binomial(x.N_[h], Float(min(Pinfected_for_each_host_indiv, 1.0))))
                    # $PRINT println("N_hosts_infected_current_genotype = ",N_hosts_infected_current_genotype)
                    $PRINT println("X.P_NcontactsInPop_1_2_3_ect_parasites = ",X.P_NcontactsInPop_1_2_3_ect_parasites)
                    if N_hosts_infected_current_genotype > 0
                        fill!(X.Pinfection, 0.0) # for each parasite genotype, what is the probability of infecting the current host geneotype ?
                        for p in X.Posi_HigherOrEqualVirulenceS_comparToVirulenceInHost # only conscider cases where the infecting parasit has a higher virulence
                            X.Pinfection[p] = GetPinfection!(X.ParasitFreq_x_Pinfection1contact[p], X.P_SameVirulenceS[p], X.P_HigherVirulenceS[p], X.P_NcontactsInPop_1_2_3_ect_parasites, NcontactsMax)
                        end
                        $PRINT println("X.Pinfection = ",X.Pinfection)
                        if sum(X.Pinfection) !== 0.0
                            $PRINT println("X.Pinfection = ",X.Pinfection)
                            $PRINT println("X.Pinfection = ",X.Pinfection)
                            X.Pinfection ./= sum(X.Pinfection)
                            rand!(Multinomial(N_hosts_infected_current_genotype, X.Pinfection), X.NsuccessfullInfectionEachParasit)
                            $PRINT println("X.NsuccessfullInfectionEachParasit= ",X.NsuccessfullInfectionEachParasit)
                            $PRINT println("sum(X.NsuccessfullInfectionEachParasit) = ",sum(X.NsuccessfullInfectionEachParasit))
                            #
                            # Apply the infections ######################
                            x.N_[h] -= N_hosts_infected_current_genotype
                            for p in 1:NparasCategories
                                if X.NsuccessfullInfectionEachParasit[p] > 0
                                    # Is the new host categogy existing ?
                                    $(ifelse(SimulateImmunisation
                                        ,"idCat = GetCategoryID(x.IDhS[h], x.IDhSphen[h], parasitesPop.IDpGenPhen_HostShiftHistoryS[p], x.IDiS[h])"
                                        ,"idCat = GetCategoryID(x.IDhS[h], x.IDhSphen[h], parasitesPop.IDpGenPhen_HostShiftHistoryS[p]           )"  ))
                                    $PRINT println("idCat = ",idCat)
                                    $(ifelse(HostShift
                                        ,"localPosi = nothing"
                                        ,"localPosi = findfirst(xx -> xx === idCat, x.IDcategorieS)"))
                                    $PRINT println("localPosi = ",localPosi)
                                    if isnothing(localPosi)
                                        push!(x.IDhS                        , x.IDhS[                                   h])
                                        push!(x.IDhSphen                    , x.IDhSphen[                               h])
                                        push!(x.posiInHostsGenotypeS_List   , x.posiInHostsGenotypeS_List[              h])
                                        push!(x.posiInParasitGenotypeS_List , parasitesPop.posiInParasitGenotypeS_List[ p])
                                        push!(x.IDpS                        , parasitesPop.IDpS[                        p])
                                        push!(x.IDpSphen                    , parasitesPop.IDpSphen[                    p])
                                        push!(x.IDpGenPhen_HostShiftHistoryS, parasitesPop.IDpGenPhen_HostShiftHistoryS[p])
                                        push!(x.IDp_HostShiftHistoryS       , parasitesPop.IDp_HostShiftHistoryS[       p])
                                        $(ifelse(SimulateImmunisation
                                            ,"push!(x.IDiS                  , x.IDiS[                                   h])"
                                            ,""  ))
                                        push!(x.IDcategorieS                , idCat                                       )
                                        $(ifelse(GfG
                                        ,"push!(x.rate_of_gametes_prodS     , x.rate_of_gametes_prodS[                  h])
                                          push!(x.rate_of_parasites_emissionS, parasitesPop.rate_of_parasites_emissionS[p])"
                                        ,""))
                                        push!(x.virulenceS                  , HPinteractions.virulenceS[             hostsGenotypeS_List.GenotypeS[   x.posiInHostsGenotypeS_List[h]   ].traits_infection_success_Symb     ][    parasitesPop.posiInParasitGenotypeS_List[p]   ])
                                        push!(x.Precoveryinnateimmu         , HPinteractions.Precoveryinnateimmu[    hostsGenotypeS_List.GenotypeS[   x.posiInHostsGenotypeS_List[h]   ].traits_immunity_Symb              ][    parasitesPop.posiInParasitGenotypeS_List[p]   ])
                                        $(ifelse(SimulateImmunisation
                                        ,"push!(x.Precoveryacquiredimmu,
                                            maximum( map(posiImmunisingParas -> GetPrecoveryacquiredimmu!(HPinteractions.Precoveryacquiredimmu                # Precoveryacquiredimmu
                                                                             , posiImmunisingParas # posiImmunisingParas
                                                                             , parasitesPop.posiInParasitGenotypeS_List[p]                  # posiInfectingParas
                                                                             , parasitesGenotypeS_List
                                                                             , HPinteractions.TempDistImmunisingInfecting
                                                                             ), HPinteractions.immunisationS[  x.IDiS[h]  ].History ))
                                         )"
                                        ,""))
                                        push!(x.N_                          , X.NsuccessfullInfectionEachParasit[p] )
                                        push!(x.dN_                         , 0 )
                                    else
                                        x.N_[localPosi] += X.NsuccessfullInfectionEachParasit[p]
                                    end
                                    $(ifelse(HostShift
                                    ,"add::Symbol = Symbol(Symbol(\">IDp°\"),  parasitesPop.IDpS[p],  :__Sp°,  spRecipient ,  :°°Time°,  round(Time,digits=3),  :°°IDh°,  x.IDhSphen[h])
                                      x.IDpGenPhen_HostShiftHistoryS[end]=Symbol(x.IDpGenPhen_HostShiftHistoryS[end], add )
                                      x.IDp_HostShiftHistoryS[end]       =Symbol(x.IDp_HostShiftHistoryS[       end], add )"
                                    ,""
                                    ))
                                end
                            end
                        end
                    end
                end
            end
        end
        $INBOUNDend
        $FASTMATHend
    end"""
    eval(Meta.parse(strLoop))
end

str="""
function RecoverInnateimmu!(X::MetaCommunity,dt::Float)
    $FASTMATHstart
    $INBOUNDstart
    $PRINT println("RecoverInnateimmu!")
    for x in getfield(X,:$(fieldname(MetaCommunity,1)))
        RecoverInnateimmu_MetaPop!(x,$(ifelse(SimulateImmunisation,"X.HPinteractions.immunisationS,",""))dt)
    end
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

str="""
function RecoverInnateimmu_MetaPop!(X::MetaPop,$(ifelse(SimulateImmunisation,"immunisationS::Dict{Symbol,Immunisations},",""))dt::Float)
    $FASTMATHstart
    $INBOUNDstart
    $PRINT println("RecoverInnateimmu_MetaPop!")
    for x in getfield(X,:$(fieldname(MetaPop,1)))
        if x.N[] > 0
            RecoverInnateimmu_HostSexS!(x,$(ifelse(SimulateImmunisation,"immunisationS,",""))dt)
        end
    end
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

str="""
function RecoverInnateimmu_HostSexS!(X::HostSexS,$(ifelse(SimulateImmunisation,"immunisationS::Dict{Symbol,Immunisations},",""))dt::Float)
    $FASTMATHstart
    $INBOUNDstart
    $PRINT println("RecoverInnateimmu_HostSexS!")
    for x in getfield(X,:$(fieldname(HostSexS,1)))
        if x.N[] > 0
            Recover_InnateORAcquired_immu_HostsGenotypeS!(x,$(ifelse(SimulateImmunisation,"immunisationS,",""))x.Precoveryinnateimmu,dt)
        end
    end
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

str="""
function Recover_InnateORAcquired_immu_HostsGenotypeS!(X::HostsGenotypeS,$(ifelse(SimulateImmunisation,"immunisationS::Dict{Symbol,Immunisations},",""))Precovery::Vector{Float},dt::Float)
    $FASTMATHstart
    $INBOUNDstart
        $PRINT println("Recover_InnateORAcquired_immu_HostsGenotypeS!")
        resize!(X.Nrecovery,length(X.N_))
        $PRINT println("Precovery = ", Precovery)
        $PRINT println("1.0-(1.0-Precovery)^dt = ", 1.0.-(1.0.-Precovery).^dt)
        ## checking : function Nreco!(N::Int,p::Float,dt::Float,End::Float);    Nreco = 0.0;    for t in 0:dt:End;        recov = rand(Binomial(N, 1.0-(1.0-p)^dt ));        N -= recov;        Nreco += recov;        recov;    end;    return Nreco;end                  ;            dt = 0.05 ;  p = 0.9 ; mean([Nreco!(100,p,dt,1.5) for rep in 1:100000])
        map!(i -> rand(Binomial(X.N_[i], 1.0-(1.0-Precovery[i])^dt)), X.Nrecovery, eachindex(X.N_))
        $PRINT println("X.Nrecovery = ",X.Nrecovery)
        for h in findall(X.Nrecovery .!== 0)
            X.N_[h] -= X.Nrecovery[h]
            $(ifelse(SimulateImmunisation
                    ,"newIDi = immunisationS[ X.IDiS[h], X.posiInParasitGenotypeS_List[h]]
                      # println(\"X.IDiS[h] = \", X.IDiS[h], \" X.posiInParasitGenotypeS_List[h] = \", X.posiInParasitGenotypeS_List[h], \" newIDi = \",newIDi)
                      NewCategoryID = GetCategoryID(X.IDhS[h], X.IDhSphen[h], :IDparasitCat_0_0, newIDi )"
                    ,"NewCategoryID = GetCategoryID(X.IDhS[h], X.IDhSphen[h], :IDparasitCat_0_0         )"        ))
            localPosi = findfirst(xx -> xx ===NewCategoryID, X.IDcategorieS)
            $PRINT println("localPosi = ",localPosi)
            if isnothing(localPosi)
                $PRINT println("X.IDhS[h] = ",X.IDhS[h])
                push!(X.IDhS                         , X.IDhS[h]                      )
                $PRINT println("X.IDhSphen[h] = ",X.IDhSphen[h])
                push!(X.IDhSphen                     , X.IDhSphen[h]                  )
                $PRINT println("X.posiInHostsGenotypeS_List[h] = ",X.posiInHostsGenotypeS_List[h])
                push!(X.posiInHostsGenotypeS_List    , X.posiInHostsGenotypeS_List[h] )
                push!(X.posiInParasitGenotypeS_List  , 0                              )
                push!(X.IDpS                         , Int_IDp(0)                     )
                push!(X.IDpSphen                     , Int_IDp(0)                     )
                push!(X.IDpGenPhen_HostShiftHistoryS , :IDparasitCat_0_0              )
                push!(X.IDp_HostShiftHistoryS        , :IDparasitCat_0_0              )
                $(ifelse(SimulateImmunisation
                        ,"push!(X.IDiS, newIDi        )"
                        ,""
                ))
                push!(X.IDcategorieS                 , NewCategoryID                  )
                $(ifelse(GfG
                ,"$PRINT println(\"X.rate_of_gametes_prodS= \",X.rate_of_gametes_prodS)
                  push!(X.rate_of_gametes_prodS     , X.rate_of_gametes_prodS[h]      )
                  push!(X.rate_of_parasites_emissionS, 0.0                            )"
                ,""))
                push!(X.virulenceS                   , 0.0                            )
                push!(X.Precoveryinnateimmu          , 0.0                            )
                $(ifelse(SimulateImmunisation
                ,"push!(X.Precoveryacquiredimmu      , 0.0                            )"
                ,""))
                $PRINT println("X.Nrecovery[h] = ",X.Nrecovery[h] )
                push!(X.N_                           , X.Nrecovery[h]                 )
                push!(X.dN_                          , 0                              )
            else
                $PRINT println("X.Nrecovery[h] = ",X.Nrecovery[h] )
                X.N_[localPosi] += X.Nrecovery[h]
            end
        end
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

## RecoverAcquiredimmu!  ########
if SimulateImmunisation
    str="""
    function RecoverAcquiredimmu!(X::MetaCommunity,dt::Float)
        $FASTMATHstart
        $INBOUNDstart
        for x in getfield(X,:$(fieldname(MetaCommunity,1)))
            RecoverAcquiredimmu_MetaPop!(x,X.HPinteractions.immunisationS,dt)
        end
        $INBOUNDend
        $FASTMATHend
    end"""
    eval(Meta.parse(str))

    str="""
    function RecoverAcquiredimmu_MetaPop!(X::MetaPop,immunisationS::Dict{Symbol,Immunisations},dt::Float)
        $FASTMATHstart
        $INBOUNDstart
        for x in getfield(X,:$(fieldname(MetaPop,1)))
            if x.N[] > 0
                RecoverAcquiredimmu_HostSexS!(x,immunisationS,dt)
            end
        end
        $INBOUNDend
        $FASTMATHend
    end"""
    eval(Meta.parse(str))

    str="""
    function RecoverAcquiredimmu_HostSexS!(X::HostSexS,immunisationS::Dict{Symbol,Immunisations},dt::Float)
        $FASTMATHstart
        $INBOUNDstart
        for x in getfield(X,:$(fieldname(HostSexS,1)))
            if x.N[] > 0
                Recover_InnateORAcquired_immu_HostsGenotypeS!(x, immunisationS, x.Precoveryacquiredimmu, dt)
            end
        end
        $INBOUNDend
        $FASTMATHend
    end"""
    eval(Meta.parse(str))
end


# Migrate Hosts Nmigr ~ Binomial ; population of settement ~ Hypergeometric(K/N)
str = """
function MigrateHosts!(X::MetaCommunity,dt::Float)
    for (sp,x) in enumerate(X.MetaPopS)
        $PRINT println("MigrateHosts!  sp = ",sp)
        MigrateHosts_MetaPop!(x, dt)
        X.N_[sp] = x.N[]
    end
    X.N[] = $(join(["X.N_[$i]" for i in 1:NhostSp], " + "))
    $PRINT println("X.N[] = ",X.N[])
    $PRINT println("MigrateHosts!  Done")
end"""
eval(Meta.parse(str))

str = """
function MigrateHosts_MetaPop!(X::MetaPop,dt::Float)
    $FASTMATHstart
    $INBOUNDstart
    map!(pop -> if (pop.K[] === 0.0) 0.0 else pop.K[] / pop.N[] end, X.PmigrSettlementEachPop, X.PopS)
    IsInf = X.PmigrSettlementEachPop .=== Inf
    if any(IsInf) X.PmigrSettlementEachPop[IsInf] .= 1 ; X.PmigrSettlementEachPop[.!IsInf] .= 0 end
    X.PmigrSettlementEachPop ./= ($(join(    ["X.PmigrSettlementEachPop[$pop]" for pop in 1:NPops], "+"    )))
    $PRINT println("X.PmigrSettlementEachPop = ",X.PmigrSettlementEachPop)
    for (currentPop,x) in enumerate(X.PopS)
        if x.N[] > 0
            MigrateHosts_HostSexS!(x, X.MigrRate[], X.PopS, X.PmigrSettlementEachPop, currentPop ,dt)
        end
        X.N_[currentPop] = x.N[]
        $PRINT println("currentPop = ",currentPop)
        $PRINT println("X.N_ = ",X.N_)
    end
    $PRINT println("X.N[] = ",X.N[])
    X.N[] = $(join(["X.N_[$i]" for i in 1:NPops], " + "))
    $PRINT println("X.N[] = ",X.N[])
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

function MigrateHosts_HostSexS!(X::HostSexS, MigrRate::Float, PopS::Vector{HostSexS}, PmigrSettlementEachPop::Vector{Float}, currentPop::Int,  dt::Float)
    for sex in 1:2
        MigrateHosts_HostsGenotypeS!(X.HostSexS[sex], MigrRate, PopS, PmigrSettlementEachPop, currentPop, sex, dt)
    end
end

str= """
function MigrateHosts_HostsGenotypeS!(X::HostsGenotypeS,MigrRate::Float,PopS::Vector{HostSexS}, PmigrSettlementEachPop::Vector{Float},currentPop::Int,sex::Int,dt::Float)
    $FASTMATHstart
    $INBOUNDstart
    resize!(X.NmigrEachGenotype,length(X.N_))
    if MigrRate*dt > 1.0 MigrRate = 1.0 else MigrRate *= dt end
    map!(n -> rand(Binomial(n,MigrRate)), X.NmigrEachGenotype, X.N_)
    $PRINT println("X.NmigrEachGenotype = ", X.NmigrEachGenotype)
    $PRINT println("findall(X.N_ !== 0) = ", findall(X.N_ !== 0))
    for g in findall(X.N_ !== 0)
        $PRINT println("g = ", g)
        rand!(Multinomial(X.NmigrEachGenotype[g], PmigrSettlementEachPop),X.NmigrEachPop)
        $PRINT println("X.NmigrEachPop = ", X.NmigrEachPop)
        $PRINT println("sum(X.NmigrEachPop[1:end .!== currentPop]) = ", sum(X.NmigrEachPop[1:end .!== currentPop]))
        if sum(X.NmigrEachPop[1:end .!== currentPop]) > 0
            # Apply migrations
            $PRINT println("X.NmigrEachGenotype[g] = ", X.NmigrEachGenotype[g])
            X.NmigrEachGenotype[g] -= X.NmigrEachPop[currentPop]
            for popSettlement in 1:$NPops
                $PRINT println("popSettlement = ",popSettlement)
                if (currentPop != popSettlement) & (X.NmigrEachPop[popSettlement] > 0) # Ignore migrations in your own population
                    $PRINT println("X.NmigrEachPop[popSettlement] = ",X.NmigrEachPop[popSettlement])
                    $PRINT println("X.N_ = ",X.N_)
                    # is the genotype existing ?
                    posiPopSettlement = findfirst(xx -> xx === X.IDcategorieS[g], PopS[popSettlement].HostSexS[sex].IDcategorieS)
                    if isnothing(posiPopSettlement)
                        push!(PopS[popSettlement].HostSexS[sex].IDhS                        , X.IDhS[g])
                        push!(PopS[popSettlement].HostSexS[sex].IDhSphen                    , X.IDhSphen[g])
                        push!(PopS[popSettlement].HostSexS[sex].IDpS                        , X.IDpS[g])
                        push!(PopS[popSettlement].HostSexS[sex].IDpSphen                    , X.IDpSphen[g])
                        push!(PopS[popSettlement].HostSexS[sex].IDpGenPhen_HostShiftHistoryS, X.IDpGenPhen_HostShiftHistoryS[g])
                        push!(PopS[popSettlement].HostSexS[sex].IDp_HostShiftHistoryS       , X.IDp_HostShiftHistoryS[g])
                        $(ifelse(SimulateImmunisation
                        ,"push!(PopS[popSettlement].HostSexS[sex].IDiS                      , X.IDiS[g])"
                        ,""
                        ))
                        push!(PopS[popSettlement].HostSexS[sex].IDcategorieS                , X.IDcategorieS[g])
                        $(ifelse(GfG
                        ,"push!(PopS[popSettlement].HostSexS[sex].rate_of_gametes_prodS     , X.rate_of_gametes_prodS[g]       )
                          push!(PopS[popSettlement].HostSexS[sex].rate_of_parasites_emissionS, X.rate_of_parasites_emissionS[g]  )"
                        ,""))
                        push!(PopS[popSettlement].HostSexS[sex].virulenceS                  , X.virulenceS[g])
                        push!(PopS[popSettlement].HostSexS[sex].Precoveryinnateimmu         , X.Precoveryinnateimmu[g])
                        $(ifelse(SimulateImmunisation
                        ,"push!(PopS[popSettlement].HostSexS[sex].Precoveryacquiredimmu     , X.Precoveryacquiredimmu[g]  )"
                        ,""))
                        push!(PopS[popSettlement].HostSexS[sex].posiInHostsGenotypeS_List   , X.posiInHostsGenotypeS_List[g])
                        push!(PopS[popSettlement].HostSexS[sex].posiInParasitGenotypeS_List , X.posiInParasitGenotypeS_List[g])

                        push!(PopS[popSettlement].HostSexS[sex].N_ , X.NmigrEachPop[popSettlement]  )
                        push!(PopS[popSettlement].HostSexS[sex].dN_, 0  )
                    else
                        PopS[popSettlement].HostSexS[sex].N_[posiPopSettlement] += X.NmigrEachPop[popSettlement]
                    end
                end
            end
        end
        X.N_ .-= X.NmigrEachGenotype
        X.N[] = sum(X.N_)
    end
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

##################
# CheckIfStop__apply_dN__RecordState!
str= """
function GetMean_trait_infection_success_H_UnInfected(X::MetaCommunity,pop::Int,sp::Int,trait::Int = 1)
    $FASTMATHstart
    $INBOUNDstart
    W = [X.MetaPopS[sp].PopS[pop].HostSexS[sex].IDpS .=== 0 for sex in 1:2]
    unInf_ = sum([ sum(   [(X.hostsGenotypeS_List.GenotypeS[posi].traits_infection_success[trait*2-1] + X.hostsGenotypeS_List.GenotypeS[posi].traits_infection_success[trait*2])/2 for posi in X.MetaPopS[sp].PopS[pop].HostSexS[sex].posiInHostsGenotypeS_List[W[sex]]]     .* X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_[W[sex]]) for sex in 1:2]) / (sum(X.MetaPopS[sp].PopS[pop].HostSexS[1].N_[W[1]])  +  sum(X.MetaPopS[sp].PopS[pop].HostSexS[2].N_[W[2]]))
    W = [X.MetaPopS[sp].PopS[pop].HostSexS[sex].IDpS .!== 0 for sex in 1:2]
      Inf_ = sum([ sum(   [(X.hostsGenotypeS_List.GenotypeS[posi].traits_infection_success[trait*2-1] + X.hostsGenotypeS_List.GenotypeS[posi].traits_infection_success[trait*2])/2 for posi in X.MetaPopS[sp].PopS[pop].HostSexS[sex].posiInHostsGenotypeS_List[W[sex]]]     .* X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_[W[sex]]) for sex in 1:2]) / (sum(X.MetaPopS[sp].PopS[pop].HostSexS[1].N_[W[1]])  +  sum(X.MetaPopS[sp].PopS[pop].HostSexS[2].N_[W[2]]))
    return [unInf_, Inf_]
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

str= """
function GetSD_trait_infection_success_H_UnInfected(X::MetaCommunity,pop::Int,sp::Int,trait::Int = 1)
    $FASTMATHstart
    $INBOUNDstart
    # Mean
    W = [X.MetaPopS[sp].PopS[pop].HostSexS[sex].IDpS .=== 0 for sex in 1:2]
    unInf_ =      sum([ sum(   [(X.hostsGenotypeS_List.GenotypeS[posi].traits_infection_success[trait*2-1] + X.hostsGenotypeS_List.GenotypeS[posi].traits_infection_success[trait*2])/2 for posi in X.MetaPopS[sp].PopS[pop].HostSexS[sex].posiInHostsGenotypeS_List[W[sex]]]                 .* X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_[W[sex]]) for sex in 1:2]) / (sum(X.MetaPopS[sp].PopS[pop].HostSexS[1].N_[W[1]])  +  sum(X.MetaPopS[sp].PopS[pop].HostSexS[2].N_[W[2]]))
    W = [X.MetaPopS[sp].PopS[pop].HostSexS[sex].IDpS .!== 0 for sex in 1:2]
      Inf_ =      sum([ sum(   [(X.hostsGenotypeS_List.GenotypeS[posi].traits_infection_success[trait*2-1] + X.hostsGenotypeS_List.GenotypeS[posi].traits_infection_success[trait*2])/2 for posi in X.MetaPopS[sp].PopS[pop].HostSexS[sex].posiInHostsGenotypeS_List[W[sex]]]                 .* X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_[W[sex]]) for sex in 1:2]) / (sum(X.MetaPopS[sp].PopS[pop].HostSexS[1].N_[W[1]])  +  sum(X.MetaPopS[sp].PopS[pop].HostSexS[2].N_[W[2]]))
    # Sd
    W = [X.MetaPopS[sp].PopS[pop].HostSexS[sex].IDpS .=== 0 for sex in 1:2]
    unInf_ = sqrt(sum([ sum(  ([(X.hostsGenotypeS_List.GenotypeS[posi].traits_infection_success[trait*2-1] + X.hostsGenotypeS_List.GenotypeS[posi].traits_infection_success[trait*2])/2 for posi in X.MetaPopS[sp].PopS[pop].HostSexS[sex].posiInHostsGenotypeS_List[W[sex]]] .- unInf_)^2    .* X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_[W[sex]]) for sex in 1:2]) / (sum(X.MetaPopS[sp].PopS[pop].HostSexS[1].N_[W[1]])  +  sum(X.MetaPopS[sp].PopS[pop].HostSexS[2].N_[W[2]])))
    W = [X.MetaPopS[sp].PopS[pop].HostSexS[sex].IDpS .!== 0 for sex in 1:2]
      Inf_ = sqrt(sum([ sum(  ([(X.hostsGenotypeS_List.GenotypeS[posi].traits_infection_success[trait*2-1] + X.hostsGenotypeS_List.GenotypeS[posi].traits_infection_success[trait*2])/2 for posi in X.MetaPopS[sp].PopS[pop].HostSexS[sex].posiInHostsGenotypeS_List[W[sex]]] .-   Inf_)^2    .* X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_[W[sex]]) for sex in 1:2]) / (sum(X.MetaPopS[sp].PopS[pop].HostSexS[1].N_[W[1]])  +  sum(X.MetaPopS[sp].PopS[pop].HostSexS[2].N_[W[2]])))
    return [unInf_, Inf_]
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

str= """
function GetMean_trait_immunity_H_UnInfected(X::MetaCommunity,pop::Int,sp::Int,trait::Int = 1)
    $FASTMATHstart
    $INBOUNDstart
    W = [X.MetaPopS[sp].PopS[pop].HostSexS[sex].IDpS .=== 0 for sex in 1:2]
    unInf_ = sum([ sum(   [(X.hostsGenotypeS_List.GenotypeS[posi].traits_immunity[trait*2-1] + X.hostsGenotypeS_List.GenotypeS[posi].traits_immunity[trait*2])/2 for posi in X.MetaPopS[sp].PopS[pop].HostSexS[sex].posiInHostsGenotypeS_List[W[sex]]]     .* X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_[W[sex]]) for sex in 1:2]) / (sum(X.MetaPopS[sp].PopS[pop].HostSexS[1].N_[W[1]])  +  sum(X.MetaPopS[sp].PopS[pop].HostSexS[2].N_[W[2]]))
    W = [X.MetaPopS[sp].PopS[pop].HostSexS[sex].IDpS .!== 0 for sex in 1:2]
      Inf_ = sum([ sum(   [(X.hostsGenotypeS_List.GenotypeS[posi].traits_immunity[trait*2-1] + X.hostsGenotypeS_List.GenotypeS[posi].traits_immunity[trait*2])/2 for posi in X.MetaPopS[sp].PopS[pop].HostSexS[sex].posiInHostsGenotypeS_List[W[sex]]]     .* X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_[W[sex]]) for sex in 1:2]) / (sum(X.MetaPopS[sp].PopS[pop].HostSexS[1].N_[W[1]])  +  sum(X.MetaPopS[sp].PopS[pop].HostSexS[2].N_[W[2]]))
    return [unInf_, Inf_]
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

str= """
function GetSd_trait_immunity_H_UnInfected(X::MetaCommunity,pop::Int,sp::Int,trait::Int = 1)
    $FASTMATHstart
    $INBOUNDstart
    # Mean
    W = [X.MetaPopS[sp].PopS[pop].HostSexS[sex].IDpS .=== 0 for sex in 1:2]
    unInf_ =      sum([ sum(   [(X.hostsGenotypeS_List.GenotypeS[posi].traits_immunity[trait*2-1] + X.hostsGenotypeS_List.GenotypeS[posi].traits_immunity[trait*2])/2 for posi in X.MetaPopS[sp].PopS[pop].HostSexS[sex].posiInHostsGenotypeS_List[W[sex]]]                  .* X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_[W[sex]]) for sex in 1:2]) / (sum(X.MetaPopS[sp].PopS[pop].HostSexS[1].N_[W[1]])  +  sum(X.MetaPopS[sp].PopS[pop].HostSexS[2].N_[W[2]]))
    W = [X.MetaPopS[sp].PopS[pop].HostSexS[sex].IDpS .!== 0 for sex in 1:2]
      Inf_ =      sum([ sum(   [(X.hostsGenotypeS_List.GenotypeS[posi].traits_immunity[trait*2-1] + X.hostsGenotypeS_List.GenotypeS[posi].traits_immunity[trait*2])/2 for posi in X.MetaPopS[sp].PopS[pop].HostSexS[sex].posiInHostsGenotypeS_List[W[sex]]]                  .* X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_[W[sex]]) for sex in 1:2]) / (sum(X.MetaPopS[sp].PopS[pop].HostSexS[1].N_[W[1]])  +  sum(X.MetaPopS[sp].PopS[pop].HostSexS[2].N_[W[2]]))
    # Sd
    W = [X.MetaPopS[sp].PopS[pop].HostSexS[sex].IDpS .=== 0 for sex in 1:2]
    unInf_ = sqrt(sum([ sum(  ([(X.hostsGenotypeS_List.GenotypeS[posi].traits_immunity[trait*2-1] + X.hostsGenotypeS_List.GenotypeS[posi].traits_immunity[trait*2])/2 for posi in X.MetaPopS[sp].PopS[pop].HostSexS[sex].posiInHostsGenotypeS_List[W[sex]]] .- unInf_)^2     .* X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_[W[sex]]) for sex in 1:2]) / (sum(X.MetaPopS[sp].PopS[pop].HostSexS[1].N_[W[1]])  +  sum(X.MetaPopS[sp].PopS[pop].HostSexS[2].N_[W[2]])))
    W = [X.MetaPopS[sp].PopS[pop].HostSexS[sex].IDpS .!== 0 for sex in 1:2]
      Inf_ = sqrt(sum([ sum(  ([(X.hostsGenotypeS_List.GenotypeS[posi].traits_immunity[trait*2-1] + X.hostsGenotypeS_List.GenotypeS[posi].traits_immunity[trait*2])/2 for posi in X.MetaPopS[sp].PopS[pop].HostSexS[sex].posiInHostsGenotypeS_List[W[sex]]] .-   Inf_)^2     .* X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_[W[sex]]) for sex in 1:2]) / (sum(X.MetaPopS[sp].PopS[pop].HostSexS[1].N_[W[1]])  +  sum(X.MetaPopS[sp].PopS[pop].HostSexS[2].N_[W[2]])))
    return [unInf_, Inf_]
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

str= """
function GetMean_trait_infection_success_H(X::MetaCommunity,pop::Int,sp::Int,trait::Int = 1)
    $FASTMATHstart
    $INBOUNDstart
    sum([ sum(   [(X.hostsGenotypeS_List.GenotypeS[posi].traits_infection_success[trait*2-1] + X.hostsGenotypeS_List.GenotypeS[posi].traits_infection_success[trait*2])/2 for posi in X.MetaPopS[sp].PopS[pop].HostSexS[sex].posiInHostsGenotypeS_List]     .* X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_) for sex in 1:2]) / X.MetaPopS[sp].PopS[pop].N[]
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))
str= """
function GetSd_trait_infection_success_H(X::MetaCommunity,pop::Int,sp::Int,trait::Int = 1)
    $FASTMATHstart
    $INBOUNDstart
    M =  sum([ sum(   [(X.hostsGenotypeS_List.GenotypeS[posi].traits_infection_success[trait*2-1] + X.hostsGenotypeS_List.GenotypeS[posi].traits_infection_success[trait*2])/2 for posi in X.MetaPopS[sp].PopS[pop].HostSexS[sex].posiInHostsGenotypeS_List]         .* X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_) for sex in 1:2]) / X.MetaPopS[sp].PopS[pop].N[]
    sqrt(sum([ sum(  ([(X.hostsGenotypeS_List.GenotypeS[posi].traits_infection_success[trait*2-1] + X.hostsGenotypeS_List.GenotypeS[posi].traits_infection_success[trait*2])/2 for posi in X.MetaPopS[sp].PopS[pop].HostSexS[sex].posiInHostsGenotypeS_List] .- M)^2 .* X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_) for sex in 1:2]) / X.MetaPopS[sp].PopS[pop].N[])
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

str= """
function GetMean_trait_immunity_H(X::MetaCommunity,pop::Int,sp::Int,trait::Int = 1)
    $FASTMATHstart
    $INBOUNDstart
    sum([ sum(   [(X.hostsGenotypeS_List.GenotypeS[posi].traits_immunity[         trait*2-1] + X.hostsGenotypeS_List.GenotypeS[posi].traits_immunity[         trait*2])/2 for posi in X.MetaPopS[sp].PopS[pop].HostSexS[sex].posiInHostsGenotypeS_List]     .* X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_) for sex in 1:2]) / X.MetaPopS[sp].PopS[pop].N[]
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))
str= """
function GetSd_trait_immunity_H(X::MetaCommunity,pop::Int,sp::Int,trait::Int = 1)
    $FASTMATHstart
    $INBOUNDstart
    M =  sum([ sum(   [(X.hostsGenotypeS_List.GenotypeS[posi].traits_immunity[         trait*2-1] + X.hostsGenotypeS_List.GenotypeS[posi].traits_immunity[         trait*2])/2 for posi in X.MetaPopS[sp].PopS[pop].HostSexS[sex].posiInHostsGenotypeS_List]         .* X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_) for sex in 1:2]) / X.MetaPopS[sp].PopS[pop].N[]
    sqrt(sum([ sum(  ([(X.hostsGenotypeS_List.GenotypeS[posi].traits_immunity[         trait*2-1] + X.hostsGenotypeS_List.GenotypeS[posi].traits_immunity[         trait*2])/2 for posi in X.MetaPopS[sp].PopS[pop].HostSexS[sex].posiInHostsGenotypeS_List] .- M)^2 .* X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_) for sex in 1:2]) / X.MetaPopS[sp].PopS[pop].N[])
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

str= """
function GetMean_trait_infection_success_P(X::MetaCommunity,pop::Int,sp::Int,SetNParasites::Bool,trait::Int = 1)
    $FASTMATHstart
    $INBOUNDstart
    if SetNParasites SetNParasites!(X, false) end
    sum( [X.parasitesGenotypeS_List.GenotypeS[posi].traits_infection_success[trait] for posi in X.MetaPopS[sp].PopS[pop].ParasitesPop.posiInParasitGenotypeS_List] .* X.MetaPopS[sp].PopS[pop].ParasitesPop.N_) / X.MetaPopS[sp].PopS[pop].ParasitesPop.N[]
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))
str= """
function GetSd_trait_infection_success_P(X::MetaCommunity,pop::Int,sp::Int,SetNParasites::Bool,trait::Int = 1)
    $FASTMATHstart
    $INBOUNDstart
    if SetNParasites SetNParasites!(X, false) end
    M =  sum( [X.parasitesGenotypeS_List.GenotypeS[posi].traits_infection_success[trait] for posi in X.MetaPopS[sp].PopS[pop].ParasitesPop.posiInParasitGenotypeS_List]         .* X.MetaPopS[sp].PopS[pop].ParasitesPop.N_) / X.MetaPopS[sp].PopS[pop].ParasitesPop.N[]
    sqrt(sum(([X.parasitesGenotypeS_List.GenotypeS[posi].traits_infection_success[trait] for posi in X.MetaPopS[sp].PopS[pop].ParasitesPop.posiInParasitGenotypeS_List] .- M)^2 .* X.MetaPopS[sp].PopS[pop].ParasitesPop.N_) / X.MetaPopS[sp].PopS[pop].ParasitesPop.N[])
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

str= """
function GetMean_trait_immunity_P(X::MetaCommunity,pop::Int,sp::Int,SetNParasites::Bool,trait::Int = 1)
    $FASTMATHstart
    $INBOUNDstart
    if SetNParasites SetNParasites!(X, false) end
    sum( [X.parasitesGenotypeS_List.GenotypeS[posi].traits_immunity[trait] for posi in X.MetaPopS[sp].PopS[pop].ParasitesPop.posiInParasitGenotypeS_List] .* X.MetaPopS[sp].PopS[pop].ParasitesPop.N_) / X.MetaPopS[sp].PopS[pop].ParasitesPop.N[]
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))
str= """
function GetSd_trait_immunity_P(X::MetaCommunity,pop::Int,sp::Int,SetNParasites::Bool,trait::Int = 1)
    $FASTMATHstart
    $INBOUNDstart
    if SetNParasites SetNParasites!(X, false) end
    M =  sum( [X.parasitesGenotypeS_List.GenotypeS[posi].traits_immunity[trait] for posi in X.MetaPopS[sp].PopS[pop].ParasitesPop.posiInParasitGenotypeS_List]         .* X.MetaPopS[sp].PopS[pop].ParasitesPop.N_) / X.MetaPopS[sp].PopS[pop].ParasitesPop.N[]
    sqrt(sum(([X.parasitesGenotypeS_List.GenotypeS[posi].traits_immunity[trait] for posi in X.MetaPopS[sp].PopS[pop].ParasitesPop.posiInParasitGenotypeS_List] .- M)^2 .* X.MetaPopS[sp].PopS[pop].ParasitesPop.N_) / X.MetaPopS[sp].PopS[pop].ParasitesPop.N[])
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

str= """
function GetFreq_trait_H(X::MetaCommunity, pop::Int, sp::Int, allele::Int, trait::Int)
    $FASTMATHstart
    $INBOUNDstart
    # (X.MetaPopS[sp].PopS[pop].HostSexS[1].AllelesPool[trait][allele$addOneIfGfG_str] * X.MetaPopS[sp].PopS[pop].HostSexS[1].N[]  +  X.MetaPopS[sp].PopS[pop].HostSexS[2].AllelesPool[trait][allele$addOneIfGfG_str] * X.MetaPopS[sp].PopS[pop].HostSexS[2].N[]) / X.MetaPopS[sp].PopS[pop].N[]
    # problem witht the strategy above is that AllelesPool is a random sample of the genotypes
    
    
    
    # nallele = 0
    # for sex in 1:2
    #         for posiOnChr in 0:1
    #                 for posi in unique(X.MetaPopS[sp].PopS[pop].HostSexS[sex].posiInHostsGenotypeS_List)
    #                         
    #                 end
    #                 
    #                 
    #             nallele += 
    #             sum( 
    #                 X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_[
    #                     [X.hostsGenotypeS_List.GenotypeS[posi].Alleles[trait*2-posiOnChr]
    #             for posi in X.MetaPopS[sp].PopS[pop].HostSexS[sex].posiInHostsGenotypeS_List] .== allele ]
    #         )
    #         end
    # end
            
    nallele = 0
    for sex in 1:2
        for posiOnChr in 0:1
            for posi in unique(X.MetaPopS[sp].PopS[pop].HostSexS[sex].posiInHostsGenotypeS_List)
                if X.hostsGenotypeS_List.GenotypeS[posi].Alleles[trait*2-posiOnChr] == allele
                    nallele += 
                        sum(
                            X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_[
                                X.MetaPopS[sp].PopS[pop].HostSexS[sex].posiInHostsGenotypeS_List .== posi
                                ])
                end
            end
        end
    end
            
    return(nallele / (X.MetaPopS[sp].PopS[pop].N[]*2)) # 2 chromosomes
    
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

str= """
function GetFreq_trait_P(X::MetaCommunity, pop::Int, sp::Int, allele::Int, trait::Int)
    $FASTMATHstart
    $INBOUNDstart
    sum(X.MetaPopS[sp].PopS[pop].ParasitesPop.N_[[X.parasitesGenotypeS_List.GenotypeS[posi].Alleles[trait] for posi in X.MetaPopS[sp].PopS[pop].ParasitesPop.posiInParasitGenotypeS_List] .== allele]) / X.MetaPopS[sp].PopS[pop].ParasitesPop.N[]
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))


str= """
function GetFitn_trait_H(X::MetaCommunity, pop::Int, sp::Int, allele::Int, trait::Int)
    $FASTMATHstart
    $INBOUNDstart
    $PRINT println("GetFitn_trait_H")
    $PRINT println("pop = ",pop)
    $PRINT println("sp = ",sp)
    $PRINT println("allele = ",allele)
    $PRINT println("trait = ",trait)
    # # is the allele existing in this host pop
    if X.MetaPopS[sp].PopS[pop].HostSexS[1].AllelesPool[trait][allele$addOneIfGfG_str] + X.MetaPopS[sp].PopS[pop].HostSexS[2].AllelesPool[trait][allele$addOneIfGfG_str] === 0.0
        return(NaN)
    else
        if trait <= $N_traits_infection_success
            return(
            1 - sum([GetPinfection1contact( CircularDist_infection(allele, ParasitAllele)) for ParasitAllele in (1-$if_GfG_1_else_0):NallelesPerTrait_infection_success_AllSp]  .*  [x[end] for x in X.Storage.Recorded.ParasTraits_infection_success_freq[sp][pop][trait                            ]] )
            )
        else
            return(
                sum([GetPrecoveryinnateimmu(CircularDist_immunity( allele, ParasitAllele)) for ParasitAllele in (1-$if_GfG_1_else_0):NallelesPerTrait                        ]  .*  [x[end] for x in X.Storage.Recorded.ParasTraits_immunity_freq[         sp][pop][trait-$N_traits_infection_success]] )
            )
        end
    end
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

str= """
function GetFitn_trait_P(X::MetaCommunity, pop::Int, sp::Int, allele::Int, trait::Int)
    $FASTMATHstart
    $INBOUNDstart
    $PRINT println("GetFitn_trait_P")
    $PRINT println("pop = ",pop)
    $PRINT println("sp = ",sp)
    $PRINT println("allele = ",allele)
    $PRINT println("trait = ",trait)
    if trait <= $N_traits_infection_success
        # is the allele existing in this parasit pop
        if any([X.parasitesGenotypeS_List.GenotypeS[posi].traits_infection_success[trait                            ] for posi in X.MetaPopS[sp].PopS[pop].ParasitesPop.posiInParasitGenotypeS_List] .=== allele)
            return(
                sum([GetPinfection1contact( CircularDist_infection(allele, HostAllele)) for HostAllele in $(ifelse(GfG,"[0,$AlleleRangePerSp_infection_success[sp]...]","$AlleleRangePerSp_infection_success[sp]"))                      ]  .*  [x[end] for x in X.Storage.Recorded.HostsTraits_infection_success_freq[sp][pop][trait                            ]][$AlleleRangePerSp_infection_success[sp]] )
            )
        else
            return(NaN)
        end
     else
        # is the allele existing in this parasit pop
        if any([X.parasitesGenotypeS_List.GenotypeS[posi].traits_immunity[         trait-$N_traits_infection_success] for posi in X.MetaPopS[sp].PopS[pop].ParasitesPop.posiInParasitGenotypeS_List] .=== allele)
            return(
              1 - sum([GetPrecoveryinnateimmu(CircularDist_immunity( allele, HostAllele)) for HostAllele in (1-$if_GfG_1_else_0):NallelesPerTrait]  .*  [x[end] for x in X.Storage.Recorded.HostsTraits_immunity_freq[         sp][pop][trait-$N_traits_infection_success]] )
            )
        else
            return(NaN)
        end
    end
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))


str= """
function GetParasitesHe(X::MetaCommunity,sp::Int,pop::Int,trait::Int)
    if length( X.MetaPopS[sp].PopS[pop].ParasitesPop.IDpS) > 0
        Freq::Vector{Union{Int,Float}} = [ sum([if (g.Alleles[trait] === allele) sum(X.MetaPopS[sp].PopS[pop].ParasitesPop.N_[X.MetaPopS[sp].PopS[pop].ParasitesPop.IDpS .=== g.IDp]) else 0 end for g in X.parasitesGenotypeS_List.GenotypeS]) for allele in (1-$if_GfG_1_else_0):$NallelesPerTrait_EachTrait[trait]]
        Freq ./= sum(Freq)
        return 1-sum(Freq.^2)
    else
        return 0.0
    end
end"""
eval(Meta.parse(str))


for Dt in [1,2]
    for ParasitesMustNotDisapear in [ [Nothing], [Nothing,Vector{NTuple{2,Int}}] ][ Dt ]
        strLoop = """
        function CheckIfStop__apply_dN__RecordState!(X::MetaCommunity $(["",", dt::Float, RecordEvery::Float, SpPop_WhereParasitesMustNotDisapear::$ParasitesMustNotDisapear"][Dt]))
            $FASTMATHstart
            $INBOUNDstart
            $(["Record = true", "Record = X.Storage.TimeEvolved[] >= X.Storage.NextTimeRecord[] ; X.Storage.TimeEvolved[] += dt
               if Record println(X.Storage.TimeEvolved[]) ; X.Storage.NextTimeRecord[] += RecordEvery end"][Dt])
            $(["push!(X.Storage.Recorded.Time, 0.0)","if Record push!(X.Storage.Recorded.Time, X.Storage.TimeEvolved[] + dt) end"][Dt])
            # ExtinctionParasites = Vector{Int}()
            if Record
                SetNParasites!(X)
                $PRINT for sp in 1:length(X.MetaPopS) for pop in 1:length(X.MetaPopS[sp].PopS) println("1_") ; println(X.MetaPopS[sp].PopS[pop].ParasitesPop) end end
                UpdateHostAllelesPool!(X)
                push!(X.Storage.Recorded.RealWordTime, Dates.now() )
                Equilibrium = Vector{Bool}()
                $(ifelse(ALLOWSforLARGE_INACCURATE_dt
                ,""
                ,"tooLargeDt = false"))
            end
            for sp in 1:$NhostSp
                $PRINT println("sp=",sp)
                for pop in 1:$NPops
                    $PRINT println("pop=",pop)
                    for sex in 1:2
                        $PRINT println("sex=",sex)
                        # X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_ .+= X.MetaPopS[sp].PopS[pop].HostSexS[sex].dN_
                        $(ifelse(ALLOWSforLARGE_INACCURATE_dt
                        ,""
                        ,"  if any(X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_ .< 0)  tooLargeDt = true end
                          # if abs(sum(X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_[   findall(X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_ .< 0)  ])) > (sum(abs.(X.MetaPopS[sp].PopS[pop].HostSexS[sex].dN_[   findall(X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_ .> 0)  ]))/10) tooLargeDt = true end"))
                        # Delete empty categories
                        Del = findall(X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_ .<= 0)
                        if length(Del) > 0
                            $(join(["deleteat!(X.MetaPopS[sp].PopS[pop].HostSexS[sex]."*string(field)*",Del) \n  " for field in ListFieldsOneValPerCat]) )
                        end
                        X.MetaPopS[sp].PopS[pop].N_[sex] = X.MetaPopS[sp].PopS[pop].HostSexS[sex].N[] = sum(X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_)
                        # fill!(X.MetaPopS[sp].PopS[pop].HostSexS[sex].dN_,0)
                    end
                    X.MetaPopS[sp].N_[pop] = X.MetaPopS[sp].PopS[pop].N[] = sum(X.MetaPopS[sp].PopS[pop].N_)
                    $PRINT println("X.MetaPopS[sp].N_[pop] = ",X.MetaPopS[sp].N_[pop])
                    #
                    $(ifelse(ParasitesMustNotDisapear == Nothing
                    ,""
                    ,"if any(SpPop_WhereParasitesMustNotDisapear .=== [(sp,pop)]) && X.MetaPopS[sp].PopS[pop].ParasitesPop.N[] === 0
                        X.Storage.Extinction[] = true
                      end"))
                    # Demography
                    $PRINT println("Demography")
                    # push!(ExtinctionParasites, X.MetaPopS[sp].PopS[pop].ParasitesPop.N[])
                    if Record
                        push!(X.Storage.Recorded.Parasites[sp][pop]  , sum(X.MetaPopS[sp].PopS[pop].HostSexS[1].N_[X.MetaPopS[sp].PopS[pop].HostSexS[1].IDpS.!==0]) + sum(X.MetaPopS[sp].PopS[pop].HostSexS[2].N_[X.MetaPopS[sp].PopS[pop].HostSexS[2].IDpS.!==0]))
                        push!(X.Storage.Recorded.Hosts[sp][pop]     , X.MetaPopS[sp].PopS[pop].N[])
                        # Evolution
                        for trait_infection_success in 1:N_traits_infection_success
                            $PRINT println("trait_infection_success = ",trait_infection_success)
                            for allele in $(1-if_GfG_1_else_0):$NallelesPerTrait_infection_success_AllSp
                                $PRINT println("allele = ",allele)
                                push!(X.Storage.Recorded.HostsTraits_infection_success_freq[sp][pop][trait_infection_success][allele$addOneIfGfG_str], GetFreq_trait_H(X, pop, sp, allele, trait_infection_success))
                                push!(X.Storage.Recorded.ParasTraits_infection_success_freq[sp][pop][trait_infection_success][allele$addOneIfGfG_str], GetFreq_trait_P(X, pop, sp, allele, trait_infection_success))
                            end
                            for allele in $(1-if_GfG_1_else_0):$NallelesPerTrait_infection_success_AllSp  # GetFitn_... uses the output of GetFreq_...
                                push!(X.Storage.Recorded.HostsTraits_infection_success_fitn[sp][pop][trait_infection_success][allele$addOneIfGfG_str], GetFitn_trait_H(X, pop, sp, allele, trait_infection_success))
                                push!(X.Storage.Recorded.ParasTraits_infection_success_fitn[sp][pop][trait_infection_success][allele$addOneIfGfG_str], GetFitn_trait_P(X, pop, sp, allele, trait_infection_success))
                            end
                        end
                        $PRINT println("trait_infection_success = Done")
                        for trait_immunity in 1:N_traits_immunity
                            $PRINT println("trait_immunity = ",trait_immunity)
                            for allele in $(1-if_GfG_1_else_0):$NallelesPerTrait
                                $PRINT println("allele = ",allele)
                                push!(X.Storage.Recorded.HostsTraits_immunity_freq[sp][pop][trait_immunity][allele$addOneIfGfG_str], GetFreq_trait_H(X, pop, sp, allele, $N_traits_infection_success+trait_immunity))
                                push!(X.Storage.Recorded.ParasTraits_immunity_freq[sp][pop][trait_immunity][allele$addOneIfGfG_str], GetFreq_trait_P(X, pop, sp, allele, $N_traits_infection_success+trait_immunity))
                            end
                            for allele in $(1-if_GfG_1_else_0):$NallelesPerTrait  # GetFitn_... uses the output of GetFreq_...
                                push!(X.Storage.Recorded.HostsTraits_immunity_fitn[sp][pop][trait_immunity][allele$addOneIfGfG_str], GetFitn_trait_H(X, pop, sp, allele, $N_traits_infection_success+trait_immunity))
                                push!(X.Storage.Recorded.ParasTraits_immunity_fitn[sp][pop][trait_immunity][allele$addOneIfGfG_str], GetFitn_trait_P(X, pop, sp, allele, $N_traits_infection_success+trait_immunity))
                            end
                        end
                        $PRINT println("Pinfectionsuccess = Done")
                        $PRINT println("sp = ",sp)
                        $PRINT println("pop = ",pop)
                        push!(X.Storage.Recorded.Pinfectionsuccess[  sp][pop] , mean([GetPinfection1contact(  GetOverallDist_infection_success( [ $(join(["CircularDist_infection(  sample($(1-if_GfG_1_else_0):$NallelesPerTrait_infection_success_AllSp ,Weights([x[end] for x in X.Storage.Recorded.HostsTraits_infection_success_freq[sp][pop][$trait]])) , sample($(1-if_GfG_1_else_0):$NallelesPerTrait_infection_success_AllSp, Weights([x[end] for x in X.Storage.Recorded.ParasTraits_infection_success_freq[sp][pop][$trait]]))  )" for trait in 1:N_traits_infection_success]," , "))])) for rep in 1:100]) )
                        $PRINT println("trait_immunity = Done")
                        W = [findall(X.MetaPopS[sp].PopS[pop].HostSexS[sex].IDpS .!== 0) for sex in 1:2]
                        push!(X.Storage.Recorded.Precoveryinnateimmu[sp][pop]  , mean([GetPrecoveryinnateimmu( GetOverallDist_recoveryinnateimmu([ $(join(["CircularDist_immunity(   sample($(1-if_GfG_1_else_0):$NallelesPerTrait                        ,Weights([x[end] for x in X.Storage.Recorded.HostsTraits_immunity_freq[         sp][pop][$trait]])) , sample($(1-if_GfG_1_else_0):$NallelesPerTrait                        , Weights([x[end] for x in X.Storage.Recorded.ParasTraits_immunity_freq[         sp][pop][$trait]]))  )" for trait in 1:N_traits_immunity         ]," , "))])) for rep in 1:100]) )
                        $PRINT println("Precoveryinnateimmu = Done")
                        if X.MetaPopS[sp].PopS[pop].ParasitesPop.N[] === 0
                            push!(X.Storage.Recorded.Precoveryacquiredimmu[sp][pop], NaN)
                        else
                            push!(X.Storage.Recorded.Precoveryacquiredimmu[sp][pop], mean([
                            maximum( map(posiImmunisingParas -> GetPrecoveryacquiredimmu!(X.HPinteractions.Precoveryacquiredimmu                # Precoveryacquiredimmu
                                                             , posiImmunisingParas # posiImmunisingParas
                                                             , sample(X.MetaPopS[sp].PopS[pop].ParasitesPop.posiInParasitGenotypeS_List, Weights(X.MetaPopS[sp].PopS[pop].ParasitesPop.N_))  # posiInfectingParas
                                                             , X.parasitesGenotypeS_List
                                                             , X.HPinteractions.TempDistImmunisingInfecting
                                                            ), X.HPinteractions.immunisationS[  sample([X.MetaPopS[sp].PopS[pop].HostSexS[1].IDiS..., X.MetaPopS[sp].PopS[pop].HostSexS[2].IDiS...], Weights([X.MetaPopS[sp].PopS[pop].HostSexS[1].N_..., X.MetaPopS[sp].PopS[pop].HostSexS[2].N_...]))  ].History ))
                                                              for rep in 1:100])
                            )
                        end
                        $PRINT println("Precoveryacquiredimmu = Done")
                        #
                        push!(X.Storage.Recorded.HostsHe[    sp][pop] , sum([ 1-sum((sum([X.MetaPopS[sp].PopS[pop].HostSexS[sex].AllelesPool[trait] for sex in 1:2])./2).^2) for trait in 1:$Ntraits])/$Ntraits )
                        push!(X.Storage.Recorded.ParasitesHe[sp][pop] , sum([ GetParasitesHe(X,sp,pop,trait)                                                                 for trait in 1:$Ntraits])/$Ntraits )
                        #
                        L = length(X.Storage.Recorded.HostsHe[   sp][pop])
                        push!(Equilibrium, (L >= 100) && abs(FastLinearModel(X.Storage.Recorded.HostsHe[    sp][pop][Int(round(L*0.5)):end])) < $(0.0001*NhostSp*NPops) )
                        push!(Equilibrium, (L >= 100) && abs(FastLinearModel(X.Storage.Recorded.ParasitesHe[sp][pop][Int(round(L*0.5)):end])) < $(0.0001*NhostSp*NPops) )
                    end
                end
                X.N_[sp] = X.MetaPopS[sp].N[] = sum(X.MetaPopS[sp].N_)
            end
            X.N[] = sum(X.N_)
            #
            if all(X.N_ .=== 0) # | all(ExtinctionParasites .=== 0)
                X.Storage.Extinction[] = true
            end
            if Record
                X.Storage.EquilibriumAllPop[] = all(Equilibrium) | X.Storage.EquilibriumAllPop[]
                $(ifelse(ALLOWSforLARGE_INACCURATE_dt
                ,""
                ,"if tooLargeDt error(\"any(X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_ .< 0)\")
                              # error(\"To keep computations accurate, you must use a smaller dt value.\\nCurrent dt value was \$dt\")
                    end"))
            end
            $INBOUNDend
            $FASTMATHend
        end"""
        eval(Meta.parse(strLoop))
    end
end

str = """
function Implement_dN!(X::MetaCommunity)
    $FASTMATHstart
    $INBOUNDstart
    # Implement results : N = N + dN
    for sp in 1:$NhostSp
        $PRINT println("sp=",sp)
        for pop in 1:$NPops
            $PRINT println("pop=",pop)
            for sex in 1:2
                X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_ .+= X.MetaPopS[sp].PopS[pop].HostSexS[sex].dN_
                # Delete empty categories
                Del = findall(X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_ .<= 0)
                if length(Del) > 0
                    $(join(["deleteat!(X.MetaPopS[sp].PopS[pop].HostSexS[sex]."*string(field)*",Del) \n  " for field in ListFieldsOneValPerCat]) )
                end
                X.MetaPopS[sp].PopS[pop].N_[sex] = X.MetaPopS[sp].PopS[pop].HostSexS[sex].N[] = sum(X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_)
                fill!(X.MetaPopS[sp].PopS[pop].HostSexS[sex].dN_,0)
            end
            X.MetaPopS[sp].N_[pop] = X.MetaPopS[sp].PopS[pop].N[] = sum(X.MetaPopS[sp].PopS[pop].N_)
        end
        X.N_[sp] = X.MetaPopS[sp].N[] = sum(X.MetaPopS[sp].N_)
    end
    X.N[] = sum(X.N_)
    $INBOUNDend
    $FASTMATHend
end"""
eval(Meta.parse(str))

function Demography!(X::MetaCommunity, dt::Float)
    # reproduction and death need to hapen simultaneously because otherwise, when dt is large, the population size depend on dt.
    # This is because population size is always too high after reptoduction and too low after death.
    # For the other actions, we need to update directly the number of individuals to avoid negative number of individuals !
    HaveSex!(X,dt)
    Die!(X,dt)
    Implement_dN!(X)
end

##################
Actions      = (( InfectHosts! ,  Demography! ,  RecoverInnateimmu! ,  RecoverAcquiredimmu! ,  MigrateHosts! ))
# Evolve!
str = """
function Evolve!(X::MetaCommunity, dt::Float, RecordEvery::Float, SpPop_WhereParasitesMustNotDisapear::Union{Nothing,Vector{NTuple{2,Int}}}
    # , Actions::NTuple{$(length(Actions)),Function} = Actions
    , ActionOrder::Vector{Int} = collect(1:$(length(Actions))))
    shuffle!(ActionOrder)
    for act in ActionOrder
        $PRINT println("act = ",act)
        SetNParasites!(X, false, dt)
        # Np = sum([sum([sum(X.MetaPopS[sp].PopS[pop].ParasitesPop.N_ ) for pop in 1:$NPops]) for sp in 1:$NhostSp])
#        Actions[act](X,dt)
        if act==1
            InfectHosts!(X,dt)
        elseif act==2
            Demography!(X,dt)
        elseif act==3
            RecoverInnateimmu!(X,dt)
        elseif act==4
            RecoverAcquiredimmu!(X,dt)
        else act==5
            MigrateHosts!(X,dt)
        end
    end
    $PRINT println("starting CheckIfStop__apply_dN__RecordState!")
    CheckIfStop__apply_dN__RecordState!(X,dt,RecordEvery,SpPop_WhereParasitesMustNotDisapear)
    # $PRINT println()
    # $PRINT println("dN_INI = ",    sum([sum([sum([sum(X.MetaPopS[sp].PopS[pop].HostSexS[sex].dN_) for sex in 1:2]) for pop in 1:$NPops]) for sp in 1:$NhostSp])*dt    )
    # $PRINT println(" N_INI = ",    sum([sum([sum([sum(X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_ ) for sex in 1:2]) for pop in 1:$NPops]) for sp in 1:$NhostSp])*dt    )
    # InfectHosts!(X,dt)
    # $PRINT println("dN_POSTInfectHosts = ",    sum([sum([sum([sum(X.MetaPopS[sp].PopS[pop].HostSexS[sex].dN_) for sex in 1:2]) for pop in 1:$NPops]) for sp in 1:$NhostSp])    )
    # $PRINT println(" N_POSTInfectHosts = ",    sum([sum([sum([sum(X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_ ) for sex in 1:2]) for pop in 1:$NPops]) for sp in 1:$NhostSp])    )
    # RecoverInnateimmu!(X,dt)
    # $PRINT println("dN_POSTHaveSex = ",    sum([sum([sum([sum(X.MetaPopS[sp].PopS[pop].HostSexS[sex].dN_) for sex in 1:2]) for pop in 1:$NPops]) for sp in 1:$NhostSp])    )
    # $PRINT println(" N_POSTHaveSex = ",    sum([sum([sum([sum(X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_ ) for sex in 1:2]) for pop in 1:$NPops]) for sp in 1:$NhostSp])    )
    # $PRINT println("dN_POSTDie = ",    sum([sum([sum([sum(X.MetaPopS[sp].PopS[pop].HostSexS[sex].dN_) for sex in 1:2]) for pop in 1:$NPops]) for sp in 1:$NhostSp])    )
    # $PRINT println(" N_POSTDie = ",    sum([sum([sum([sum(X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_ ) for sex in 1:2]) for pop in 1:$NPops]) for sp in 1:$NhostSp])    )
    # $PRINT println("dN_POSTRecoverInnateimmu = ",    sum([sum([sum([sum(X.MetaPopS[sp].PopS[pop].HostSexS[sex].dN_) for sex in 1:2]) for pop in 1:$NPops]) for sp in 1:$NhostSp])    )
    # $PRINT println(" N_POSTRecoverInnateimmu = ",    sum([sum([sum([sum(X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_ ) for sex in 1:2]) for pop in 1:$NPops]) for sp in 1:$NhostSp])    )
    # MigrateHosts!(X,dt)
    # $PRINT println("dN_POSTMigrateHosts = ",    sum([sum([sum([sum(X.MetaPopS[sp].PopS[pop].HostSexS[sex].dN_) for sex in 1:2]) for pop in 1:$NPops]) for sp in 1:$NhostSp])    )
    # $PRINT println(" N_POSTMigrateHosts = ",    sum([sum([sum([sum(X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_ ) for sex in 1:2]) for pop in 1:$NPops]) for sp in 1:$NhostSp])    )
    # $PRINT println("CheckIfStop__apply_dN__RecordState")
    # $PRINT println("dN_POST_apply_dN = ",    sum([sum([sum([sum(X.MetaPopS[sp].PopS[pop].HostSexS[sex].dN_) for sex in 1:2]) for pop in 1:$NPops]) for sp in 1:$NhostSp])    )
    # $PRINT println(" N_POST_apply_dN = ",    sum([sum([sum([sum(X.MetaPopS[sp].PopS[pop].HostSexS[sex].N_ ) for sex in 1:2]) for pop in 1:$NPops]) for sp in 1:$NhostSp])    )
end"""
eval(Meta.parse(str))

str = """
function EvolveOver!(X::MetaCommunity;dt::Float,Duration::Float,RecordEvery::Float,SpPop_WhereParasitesMustNotDisapear::Union{Nothing,Vector{NTuple{2,Int}},Symbol})
    if SpPop_WhereParasitesMustNotDisapear == :All
        SpPop_WhereParasitesMustNotDisapear = $([(sp, pop) for sp in 1:NhostSp for pop in 1:NPops])
    end
    Time = 0.0
    while Time < Duration
        if X.Storage.Extinction[] break end
        Time += dt
        Evolve!(X,dt,RecordEvery,SpPop_WhereParasitesMustNotDisapear)
    end
end"""
eval(Meta.parse(str))

str = """
function EvolveOver_Or_UntilEquilibrium!(X::MetaCommunity;dt::Float,MaxDuration::Float,RecordEvery::Float,SpPop_WhereParasitesMustNotDisapear::Union{Nothing,Vector{NTuple{2,Int}},Symbol})
    if SpPop_WhereParasitesMustNotDisapear == :All
        SpPop_WhereParasitesMustNotDisapear = $([(sp, pop) for sp in 1:NhostSp for pop in 1:NPops])
    end
    Time = 0.0
    while (Time < MaxDuration) & !X.Storage.EquilibriumAllPop[]
        if X.Storage.Extinction[] break end
        Time += dt
        Evolve!(X,dt,RecordEvery,SpPop_WhereParasitesMustNotDisapear)
    end
end"""
eval(Meta.parse(str))


function PlotParameters(;legendfontsize=7)
    default(legendfontsize = legendfontsize # , titlefont = (20, "times"), guidefont = (18, :darkgreen), tickfont = (12, :orange), guide = "x", framestyle = :zerolines, yminorgrid = true
    , legend = :outerbottom)
    #### SUPR
    # plotS = []
    # TraitsNames = flat([["trait infection success "*string(t) for t in 1:N_traits_infection_success], ["trait immunity "*string(t) for t in 1:N_traits_immunity]])
    # for t in 1:Ntraits
    #     push!(plotS, heatmap(DistanceBetweenAlleles[t],title = TraitsNames[t]))
    # end
    # while length(plotS) < round(sqrt(Ntraits),RoundUp)^2
    #     push!(plotS, plot())
    # end
    # plot(plotS..., layout = Plots.GridLayout(Int(round(sqrt(Ntraits),RoundUp)), Int(round(sqrt(Ntraits),RoundUp))),yaxis=:none,size = (Int(round(sqrt(Ntraits),RoundUp)), Int(round(sqrt(Ntraits),RoundUp))).*400,guidefontsize=8)

    TextSize = 10
    plotS = []
    Kmax = maximum(maximum.(K))
    N = 1:0.1:Kmax
    # Generation time (time to get two offspring)
    Y = map(n-> 2/((b0-(b0-d0)/2×n/Kmax) ),N)
    push!(plotS, plot(N,Y,xlab="Population size",ylab="Expected generation time (time to get two offspring)",lab="",annotations=(0,maximum(Y),Plots.text("Assuming in all pop and sp K=$Kmax", :left,TextSize))))
    #
    if GfG
        # rate_of_gametes_prod - N traits infection success
        push!(plotS, plot(xlab="N traits infection success",ylab="Rate of gametes production"))
        for x2 in unique([0, 2, N_traits_immunity, N_traits_immunity*2-2, N_traits_immunity*2])
            g = [tuple(vcat(fill(1,x1),fill(0,N_traits_infection_success*2-x1), fill(1,x2),fill(0,N_traits_immunity*2-x2))...) for x1 in 0:(N_traits_infection_success*2)]
            Y = GetRate_of_gametes_prod.(g)
            plotS[end] = plot!(0:0.5:N_traits_infection_success,Y,lab="N traits immunity presents = $(x2/2)")
        end
        plotS[end] = scatter!([0                ], [MaxHostCostOfMissing_traits_infection_success], lab = "MaxHostCostOfMissing_traits_infection_success")
        plotS[end] = scatter!([N_traits_immunity], [Maxrate_of_gametes_prod                         ], lab = "Maxrate_of_gametes_prod")
        # plot(plotS[end])
        # rate_of_gametes_prod - N traits immunity
        push!(plotS, plot(xlab="N traits immunity",ylab="Rate of gametes production"  ))
        for x1 in unique([0, 2, N_traits_infection_success, N_traits_infection_success*2-2, N_traits_infection_success*2])
            g = [tuple(vcat(fill(1,x1),fill(0,N_traits_immunity*2-x1), fill(1,x2),fill(0,N_traits_immunity*2-x2))...) for x2 in 0:(N_traits_immunity*2)]
            Y = GetRate_of_gametes_prod.(g)
            plotS[end] = plot!(0:0.5:N_traits_immunity,Y,lab="N traits infection success = $(x1/2)")
        end
        plotS[end] = scatter!([N_traits_infection_success], [MaxHostCostOfHaving_traits_immunity], lab = "MaxHostCostOfHaving_traits_immunity")
        plotS[end] = scatter!([0                         ], [Maxrate_of_gametes_prod               ], lab = "Maxrate_of_gametes_prod")
        # plot(plotS[end])
        # rate_of_parasites_emission - N traits infection success
        push!(plotS, plot(xlab="N traits infection success",ylab="Rate of parasites emission"))
        for x2 in unique([0, 1, Int(round(N_traits_immunity/2)), N_traits_immunity-1, N_traits_immunity])
            g = [tuple(vcat(fill(1,x1),fill(0,N_traits_infection_success-x1), fill(1,x2),fill(0,N_traits_immunity-x2))...) for x1 in 0:N_traits_infection_success]
            Y = GetRate_of_parasites_emission.(g)
            plotS[end] = plot!(0:1:N_traits_infection_success,Y,lab="N traits immunity presents = $x2)")
        end
        plotS[end] = scatter!([N_traits_infection_success], [MaxParasCostOfHaving_traits_infection_success], lab = "MaxParasCostOfHaving_traits_infection_success")
        plotS[end] = scatter!([0                         ], [MaxRate_of_parasites_emission                    ], lab = "MaxRate_of_parasites_emission")
        # plot(plotS[end])
        # rate_of_parasites_emission - N traits immunity
        push!(plotS, plot(xlab="N traits immunity",ylab="Rate of parasites emission"  ))
        for x1 in unique([0, 1, Int(round(N_traits_immunity/2)), N_traits_immunity-1, N_traits_immunity])
            g = [tuple(vcat(fill(1,x1),fill(0,N_traits_immunity-x1), fill(1,x2),fill(0,N_traits_immunity-x2))...) for x2 in 0:N_traits_immunity]
            Y = GetRate_of_parasites_emission.(g)
            plotS[end] = plot!(0:1:N_traits_immunity,Y,lab="N traits infection success = $x1")
        end
        plotS[end] = scatter!([0                ], [MaxParasCostOfMissing_traits_immunity], lab = "MaxParasCostOfMissing_traits_immunity")
        plotS[end] = scatter!([N_traits_immunity], [MaxRate_of_parasites_emission            ], lab = "MaxRate_of_parasites_emission")
        # plot(plotS[end])
    end
    # Number of host-parasit contacts induced by each infected individuals during one time unit
    # P(infection|1 contact)
    push!(plotS, plot())
    X = 0:(MaxDist_infection/100):MaxDist_infection
    for sp in 1:NhostSp
        plotS[end] = plot!(X,GetPinfection1contact.(X),xlab="Distance between\nan host and a parasite for traits_infection_success",ylab="P(successful infection | 1 contact)",lab="sp$sp",line=flat([[:solid, :dash, :dot, :dashdot] for rep in 1:NhostSp])[sp])
    end
    # plot(plotS[end])
    # Virulence
    X = 0:0.01:Max_Pinfection1contact
    push!(plotS, plot(X,GetVirulence.(X),xlab="P(successful infection | 1 contact)",ylab="Virulence",lab="",line=:solid))
    # plot(plotS[end])
    # P(innate immu.->recovery)
    push!(plotS, plot())
    X = 0:(MaxDist_immunity/100):MaxDist_immunity
    for sp in 1:NhostSp
        plotS[end] = plot!(X,GetPrecoveryinnateimmu.(X),xlab="Distance between\nan host and a parasite for traits_immunity",ylab="P(innate immu.⇒recovery)",lab="sp$sp",line=flat([[:solid, :dash, :dot, :dashdot] for rep in 1:NhostSp])[sp])
    end
    # plot(plotS[end])
#     # P(recovery_Pg | acquired immu._Ph)
    if SimulateImmunisation
        X = 0:0.01:MaxDist_immunity
        push!(plotS, plot(X, GetPrecoveryacquiredimmu_fromMeanDist.(X), xlab="Distance between an infecting parasite and the most\nsimilar immunising parasit the host has undergone",ylab="P(acquired immu.⇒recovery)",lab=""))
    else
        push!(plotS, plot())
    end
    if GfG
        plot(plotS..., layout = Plots.GridLayout(3, 3),yaxis=:none,size = (3, 3).*350,guidefontsize=9)
    else
        push!(plotS, plot())
        plot(plotS..., layout = Plots.GridLayout(2, 3),yaxis=:none,size = (3, 2).*350,guidefontsize=9)
    end
end



# Demography
str = """
function Plot_Demography(X::MetaCommunity;Time::UnitRange{Int64}=1:1, ylog::Bool,scaleX::Float=1.0,scaleY::Float=1.0)
    if Time === 1:1 Time = X.Storage.Recorded.Time[1]:X.Storage.Recorded.Time[end] end
    W = (X.Storage.Recorded.Time .>= Time[1]) .& (X.Storage.Recorded.Time .<= Time[end])
    plotS = []
    for sp in eachindex(X.MetaPopS)
        for pop in eachindex(X.MetaPopS[sp].PopS)
            if any(X.Storage.Recorded.Hosts[   sp][pop][W] .> 0) # & any(X.Storage.Recorded.Parasites[   sp][pop][W] .> 0)
                push!(plotS, plot( X.Storage.Recorded.Time[W], X.Storage.Recorded.Hosts[   sp][pop][W],color = :blue, title = "pop = \$pop; sp = \$sp",lab = "",legendfont = font(8),ylim=(1,maximum(X.Storage.Recorded.Hosts[   sp][pop][W]))))
                plotS[end] = plot!(X.Storage.Recorded.Time[W], X.Storage.Recorded.Parasites[sp][pop][W],color = :red ,lab = "")
            else
                push!(plotS, plot([1,10,100],[1,10,100],title = "pop = \$pop; sp = \$sp",color = :white ,lab = ""))
            end
        end
    end
    if ylog
        plot(plotS..., layout = Plots.GridLayout($NhostSp, $NPops),yaxis=:log ,size = ($NPops*scaleX, $NhostSp*scaleY).*400)
    else
        plot(plotS..., layout = Plots.GridLayout($NhostSp, $NPops),yaxis=:none,size = ($NPops*scaleX, $NhostSp*scaleY).*400)
    end
end"""
eval(Meta.parse(str))
# Ngenotypes
str = """
function Plot_He(X::MetaCommunity;Time::UnitRange{Int64}=1:1, ylog::Bool)
    if Time === 1:1 Time = X.Storage.Recorded.Time[1]:X.Storage.Recorded.Time[end] end
    W = (X.Storage.Recorded.Time .>= Time[1]) .& (X.Storage.Recorded.Time .<= Time[end])
    plotS = []
    for sp in eachindex(X.MetaPopS)
        $PRINT println("sp =", sp)
        for pop in eachindex(X.MetaPopS[sp].PopS)
            $PRINT println("pop =", pop)
            if any(X.Storage.Recorded.Hosts[   sp][pop][W] .> 0)# & any(X.Storage.Recorded.Parasites[ sp][pop][W] .> 0)
                push!(plotS, plot( X.Storage.Recorded.Time[W], X.Storage.Recorded.HostsHe[ sp][pop][W],color = :blue, title = "pop = \$pop; sp = \$sp",lab = "",legendfont = font(8)))
                if sum(X.Storage.Recorded.ParasitesHe[sp][pop][W] .> 0) > 2
                    Y = X.Storage.Recorded.ParasitesHe[sp][pop][W]
                    Y[Y .=== 0] .= 1
                    plotS[end] = plot!(X.Storage.Recorded.Time[W], Y ,color = :red ,lab = "")
                    # plotS[end] = scatter!([X.Storage.Recorded.ParasitesHe[sp][pop].TimeMaxReached[]], [X.Storage.Recorded.ParasitesHe[sp][pop].MaxReached[]],color = :black ,lab = "")
                end
            else
                push!(plotS, plot([1,10,100],[1,10,100],title = "pop = \$pop; sp = \$sp",color = :white ,lab = ""))
            end
        end
    end
    if ylog
        plot(plotS..., layout = Plots.GridLayout($NhostSp, $NPops),yaxis=:log ,size = ($NPops, $NhostSp).*400)
    else
        plot(plotS..., layout = Plots.GridLayout($NhostSp, $NPops),yaxis=:none,size = ($NPops, $NhostSp).*400)
    end
end"""
eval(Meta.parse(str))

# Evolution Freq alleles traits_infection_success
str = """
function Plot_traits_infection_success_FREQ_ALLELES(X::MetaCommunity;Time::UnitRange{Int64}=1:1,trait::Int=1,scaleX::Float=1.0,scaleY::Float=1.0)
    if Time === 1:1 Time = X.Storage.Recorded.Time[1]:X.Storage.Recorded.Time[end] end
    W = (X.Storage.Recorded.Time .>= Time[1]) .& (X.Storage.Recorded.Time .<= Time[end])
    Xaxis = X.Storage.Recorded.Time[W]
    Xaxis_ = vcat(Xaxis, reverse(Xaxis))
    colors = range(HSV(0,1,1), stop=HSV(360,1,1), length=NallelesPerTrait_infection_success_AllSp + 1)[2:end]
    # colors = palette(:hsv,NallelesPerTrait_infection_success_AllSp + 1)[2:end]
    $(ifelse(GfG, "pushfirst!(colors,RGB())", ""))
    # \$(ifelse(NallelesPerTrait_infection_success_AllSp > 10
    # ,"colors = repeat(collect(palette(:tab20).colors.colors.data), $NallelesPerTrait_infection_success_AllSp)[1:$NallelesPerTrait_infection_success_AllSp]"
    # ,"colors = repeat(collect(palette(:tab10).colors.colors.data), $NallelesPerTrait_infection_success_AllSp)[1:$NallelesPerTrait_infection_success_AllSp]"
    # ))
    plotS = []
    for sp in eachindex(X.MetaPopS)
        for pop in eachindex(X.MetaPopS[sp].PopS)
            if any(X.Storage.Recorded.Hosts[   sp][pop][W] .> 0) # & any(X.Storage.Recorded.Parasites[   sp][pop][W] .> 0)
                # for trait in 1:$N_traits_infection_success
                    push!(plotS, plot(ylims = (0,1),  title = "Host traits_infection_success  pop = \$pop; sp = \$sp",titlefont = font(9) ))
                    Yaxis = fill(0.0,sum(W))
                    for allele in $(1-if_GfG_1_else_0):$NallelesPerTrait_infection_success_AllSp
                        YaxisUP = Yaxis .+ X.Storage.Recorded.HostsTraits_infection_success_freq[sp][pop][trait][allele$addOneIfGfG_str][W]
                        plotS[end] = plot!(Xaxis_, vcat(Yaxis, reverse(YaxisUP)) , seriestype=:shape,opacity=.9,fillcolor=colors[allele$addOneIfGfG_str])
                        Yaxis = YaxisUP
                    end
                # end
            else
                push!(plotS, plot([1,10,100],[1,10,100],title = "pop = \$pop; sp = \$sp",color = :white ,lab = ""))
            end
        end
    end
    for sp in eachindex(X.MetaPopS)
        for pop in eachindex(X.MetaPopS[sp].PopS)
            if any(X.Storage.Recorded.Hosts[   sp][pop][W] .> 0) # & any(X.Storage.Recorded.Parasites[   sp][pop][W] .> 0)
            # for trait in 1:$N_traits_infection_success
                push!(plotS, plot(ylims = (0,1),titlefont = font(9)))
                Yaxis = fill(0.0,sum(W))
                Lags = collect(-(sum(W)-3):(sum(W)-3))
                for allele in $(1-if_GfG_1_else_0):$NallelesPerTrait_infection_success_AllSp
                    YaxisUP = Yaxis .+ X.Storage.Recorded.ParasTraits_infection_success_freq[sp][pop][trait][allele$addOneIfGfG_str][W]
                    plotS[end] = plot!(Xaxis_, vcat(Yaxis, reverse(YaxisUP)) , seriestype=:shape,opacity=.9,fillcolor=colors[allele$addOneIfGfG_str])
                    Yaxis = YaxisUP
                end
                #
                Null = []
                title!(plotS[end], "Parasit traits_infection_success  pop = \$pop; sp = \$sp" )
            # end
            else
                push!(plotS, plot([1,10,100],[1,10,100],title = "pop = \$pop; sp = \$sp",color = :white ,lab = ""))
            end
        end
    end
    plot(plotS..., layout = Plots.GridLayout($NPops*2, $NhostSp),yaxis=:none,size = ($(NPops*NhostSp)*scaleX, 1*scaleY).*400, legend = false)
end"""
eval(Meta.parse(str))

# Evolution Freq alleles traits_immunity
str = """
function Plot_traits_immunity_FREQ_ALLELES(X::MetaCommunity;Time::UnitRange{Int64}=1:1,trait::Int=1,scaleX::Float=1.0,scaleY::Float=1.0)
    if Time == 1:1 Time = X.Storage.Recorded.Time[1]:X.Storage.Recorded.Time[end] end
    W = (X.Storage.Recorded.Time .>= Time[1]) .& (X.Storage.Recorded.Time .<= Time[end])
    Xaxis = X.Storage.Recorded.Time[W]
    Xaxis_ = vcat(Xaxis, reverse(Xaxis))
    colors = range(HSV(0,1,1), stop=HSV(360,1,1), length=NallelesPerTrait_infection_success_AllSp + 1)[2:end]
    # colors = palette(:hsv,NallelesPerTrait_infection_success_AllSp + 1)[2:end]
    $(ifelse(GfG, "pushfirst!(colors,RGB())", ""))
    # \$(ifelse(NallelesPerTrait>10
    # ,"colors = repeat(collect(palette(:tab20).colors.colors.data), $NallelesPerTrait)[1:$NallelesPerTrait]"
    # ,"colors = repeat(collect(palette(:tab10).colors.colors.data), $NallelesPerTrait)[1:$NallelesPerTrait]"
    # ))
    plotS = []
    for sp in eachindex(X.MetaPopS)
        for pop in eachindex(X.MetaPopS[sp].PopS)
            if any(X.Storage.Recorded.Hosts[   sp][pop][W] .> 0) # & any(X.Storage.Recorded.Parasites[   sp][pop][W] .> 0)
                # for trait in 1:$(N_traits_immunity)
                    push!(plotS, plot(ylims = (0,1),  title = "Host traits_immunity    pop = \$pop; sp = \$sp",titlefont = font(9) ))
                    Yaxis = fill(0.0,sum(W))
                    for allele in $(1-if_GfG_1_else_0):$NallelesPerTrait
                        YaxisUP = Yaxis .+ X.Storage.Recorded.HostsTraits_immunity_freq[sp][pop][trait][allele$addOneIfGfG_str][W]
                        plotS[end] = plot!(Xaxis_, vcat(Yaxis, reverse(YaxisUP)) , seriestype=:shape,opacity=.9,fillcolor=colors[allele$addOneIfGfG_str])
                        Yaxis = YaxisUP
                    end
                # end
            else
                push!(plotS, plot([1,10,100],[1,10,100],title = "pop = \$pop; sp = \$sp",color = :white ,lab = ""))
            end
        end
    end
    for sp in eachindex(X.MetaPopS)
        for pop in eachindex(X.MetaPopS[sp].PopS)
            if any(X.Storage.Recorded.Hosts[   sp][pop][W] .> 0) # & any(X.Storage.Recorded.Parasites[   sp][pop][W] .> 0)
                # for trait in 1:$(N_traits_immunity)
                    push!(plotS, plot(ylims = (0,1),  title = "Parasit traits_immunity    pop = \$pop; sp = \$sp",titlefont = font(9) ))
                    Yaxis = fill(0.0,sum(W))
                    Lags = collect(-(sum(W)-3):(sum(W)-3))
                    for allele in $(1-if_GfG_1_else_0):$NallelesPerTrait
                        YaxisUP = Yaxis .+ X.Storage.Recorded.ParasTraits_immunity_freq[sp][pop][trait][allele$addOneIfGfG_str][W]
                        plotS[end] = plot!(Xaxis_, vcat(Yaxis, reverse(YaxisUP)) , seriestype=:shape,opacity=.9,fillcolor=colors[allele$addOneIfGfG_str])
                        Yaxis = YaxisUP
                    end
                    #
                    Null = []
                    title!(plotS[end], "Parasit traits_immunity  pop = \$pop; sp = \$sp" )
                # end
            else
                push!(plotS, plot([1,10,100],[1,10,100],title = "pop = \$pop; sp = \$sp",color = :white ,lab = ""))
            end
        end
    end
    plot(plotS..., layout = Plots.GridLayout($NPops*2, $NhostSp),yaxis=:none,size = ($(NPops*NhostSp)*scaleX, 2*scaleY).*400, legend = false)
end"""
eval(Meta.parse(str))



function BiDirLog(x::Vector)
    W = x.>0 ; Min = minimum(x[W]) ; x[W] .= log.( x[W]) ; x[W] .+= -minimum(x[W]) + Min
    W = x.<0 ; Min = maximum(x[W]) ; x[W] .= log.(-x[W]) ; x[W] .+= -minimum(x[W]) + Min ; x[W] .= .-x[W]
    return(x)
end
function covNaN(x::Vector{T},y::Vector{T}) where T<:Number
    W = ((.!isequal.(x , NaN)) .& (.!isequal.(y , NaN)))
    cov(x[W],y[W])
end
# Covariation fitness - evolution
str = """
function Plot_Covariation_fitn_evol_traits_infection_success(X::MetaCommunity;Time::UnitRange{Int64}=1:1,trait::Int=1)
    if Time === 1:1 Time = X.Storage.Recorded.Time[1]:X.Storage.Recorded.Time[end] end
    W = (X.Storage.Recorded.Time .>= Time[1]) .& (X.Storage.Recorded.Time .<= Time[end])
    plotS = []
    for sp in eachindex(X.MetaPopS)
        println("sp = ",sp)
        for pop in eachindex(X.MetaPopS[sp].PopS)
            println("pop = ",pop)
            if any(X.Storage.Recorded.Hosts[   sp][pop][W] .> 0) # & any(X.Storage.Recorded.Parasites[   sp][pop][W] .> 0)
                push!(plotS, plot(X.Storage.Recorded.Time[W][[1,end]], [0,0], color = :grey, title = "Cov(fitness, evolution)    pop = \$pop; sp = \$sp",legendfont = font(8) ,lab = "" ))
                # for trait in 1:$N_traits_infection_success
                    Yh = [  covNaN(
                               [X.Storage.Recorded.HostsTraits_infection_success_fitn[sp][pop][trait][allele][t]                                                                                      for allele in 1:$(NallelesPerTrait_infection_success_AllSp+if_GfG_1_else_0)]
                              ,[X.Storage.Recorded.HostsTraits_infection_success_freq[sp][pop][trait][allele][t+1] - X.Storage.Recorded.HostsTraits_infection_success_freq[sp][pop][trait][allele][t] for allele in 1:$(NallelesPerTrait_infection_success_AllSp+if_GfG_1_else_0)]
                              )
                        for t in findall(W)[1:end-1]]
                    Yp = [  covNaN(
                               [X.Storage.Recorded.ParasTraits_infection_success_fitn[sp][pop][trait][allele][t]                                                                                      for allele in 1:$(NallelesPerTrait_infection_success_AllSp+if_GfG_1_else_0)]
                              ,[X.Storage.Recorded.ParasTraits_infection_success_freq[sp][pop][trait][allele][t+1] - X.Storage.Recorded.ParasTraits_infection_success_freq[sp][pop][trait][allele][t] for allele in 1:$(NallelesPerTrait_infection_success_AllSp+if_GfG_1_else_0)]
                              )
                        for t in findall(W)[1:end-1]]
                    plotS[end] = plot!(X.Storage.Recorded.Time[W][1:end-1], Yh, color = :blue, line = :solid, lab = "H")
                    plotS[end] = plot!(X.Storage.Recorded.Time[W][1:end-1], Yp, color = :red , line = :solid, lab = "P")
                    Wbx = ((.!isequal.(Yh, NaN)) .& (.!isequal.(Yp, NaN)))
                    push!(plotS, boxplot(vcat(fill("Hosts",sum(Wbx)),fill("Parasites",sum(Wbx))), (BiDirLog(vcat(Yh[Wbx],Yp[Wbx]))), lab = "", title = "±log(±Cov(fitness, evolution))",whisker_width=0.5))
                    Plots.abline!(0,0, lab = "", line=:dash)
                # end
            else
                push!(plotS, plot([1,10,100],[1,10,100],title = "Cov(fitness, evolution)    pop = \$pop; sp = \$sp",color = :white ,lab = ""))
                push!(plotS, plot([1,10,100],[1,10,100],title = "Cov(fitness, evolution)    pop = \$pop; sp = \$sp",color = :white ,lab = ""))
            end
        end
    end
    plot(plotS..., layout = Plots.GridLayout($NhostSp, $NPops*2),yaxis=:none,size = (max($NPops, $NPops*maximum(Time)/4e3)*2, $NhostSp).*400)
end"""
eval(Meta.parse(str))

str = """
function Plot_Covariation_fitn_evol_traits_immunity(X::MetaCommunity;Time::UnitRange{Int64}=1:1,trait::Int=1)
    if Time === 1:1 Time = X.Storage.Recorded.Time[1]:X.Storage.Recorded.Time[end] end
    W = (X.Storage.Recorded.Time .>= Time[1]) .& (X.Storage.Recorded.Time .<= Time[end])
    plotS = []
    for sp in eachindex(X.MetaPopS)
        for pop in eachindex(X.MetaPopS[sp].PopS)
            if any(X.Storage.Recorded.Hosts[   sp][pop][W] .> 0) # & any(X.Storage.Recorded.Parasites[   sp][pop][W] .> 0)
                push!(plotS, plot(X.Storage.Recorded.Time[W][[1,end]], [0,0], color = :grey, title = "Cov(fitness, evolution)    pop = \$pop; sp = \$sp",legendfont = font(8) ,lab = "" ))
                # for trait in 1:$N_traits_immunity
                    Yh = [  covNaN(
                               [X.Storage.Recorded.HostsTraits_immunity_fitn[sp][pop][trait][allele][t]                                                                                      for allele in 1:$(NallelesPerTrait+if_GfG_1_else_0)]
                              ,[X.Storage.Recorded.HostsTraits_immunity_freq[sp][pop][trait][allele][t+1] - X.Storage.Recorded.HostsTraits_immunity_freq[sp][pop][trait][allele][t] for allele in 1:$(NallelesPerTrait+if_GfG_1_else_0)]
                              )
                        for t in findall(W)[1:end-1]]
                    Yp = [  covNaN(
                               [X.Storage.Recorded.ParasTraits_immunity_fitn[sp][pop][trait][allele][t]                                                                                      for allele in 1:$(NallelesPerTrait+if_GfG_1_else_0)]
                              ,[X.Storage.Recorded.ParasTraits_immunity_freq[sp][pop][trait][allele][t+1] - X.Storage.Recorded.ParasTraits_immunity_freq[sp][pop][trait][allele][t] for allele in 1:$(NallelesPerTrait+if_GfG_1_else_0)]
                              )
                        for t in findall(W)[1:end-1]]
                    plotS[end] = plot!(X.Storage.Recorded.Time[W][1:end-1], Yh, color = :blue, line = :solid, lab = "H")
                    plotS[end] = plot!(X.Storage.Recorded.Time[W][1:end-1], Yp, color = :red , line = :solid, lab = "P")
                    Wbx = ((.!isequal.(Yh, NaN)) .& (.!isequal.(Yp, NaN)))
                    push!(plotS, boxplot(vcat(fill("Hosts",sum(Wbx)),fill("Parasites",sum(Wbx))), (BiDirLog(vcat(Yh[Wbx],Yp[Wbx]))), lab = "", title = "±log(±Cov(fitness, evolution))",whisker_width=0.5))
                    Plots.abline!(0,0, lab = "", line=:dash)
                # end
            else
                push!(plotS, plot([1,10,100],[1,10,100],title = "Cov(fitness, evolution)    pop = \$pop; sp = \$sp",color = :white ,lab = ""))
                push!(plotS, plot([1,10,100],[1,10,100],title = "Cov(fitness, evolution)    pop = \$pop; sp = \$sp",color = :white ,lab = ""))
            end
        end
    end
    plot(plotS..., layout = Plots.GridLayout($NhostSp, $NPops*2),yaxis=:none,size = (max($NPops, $NPops*maximum(Time)/4e3)*2, $NhostSp).*400)
end"""
eval(Meta.parse(str))

# Evolution traits_infection_success_MEAN
str = """
function Plot_traits_infection_success_MEAN(X::MetaCommunity;Time::UnitRange{Int64}=1:1)
    if Time === 1:1 Time = X.Storage.Recorded.Time[1]:X.Storage.Recorded.Time[end] end
    W = (X.Storage.Recorded.Time .>= Time[1]) .& (X.Storage.Recorded.Time .<= Time[end])
    plotS = []
    for sp in eachindex(X.MetaPopS)
        for pop in eachindex(X.MetaPopS[sp].PopS)
            if any(X.Storage.Recorded.Hosts[   sp][pop][W] .> 0) # & any(X.Storage.Recorded.Parasites[   sp][pop][W] .> 0)
                push!(plotS, plot(ylims = (0.95,$(NallelesPerTrait_infection_success_AllSp+0.05)),  title = "traits_infection_success    pop = \$pop; sp = \$sp",legendfont = font(8) ))
                for trait in 1:$N_traits_infection_success
                    plotS[end] = plot!(vcat(X.Storage.Recorded.Time[W], reverse(X.Storage.Recorded.Time[W]))
                    ,vcat(X.Storage.Recorded.HostsTraits_infection_success_mean[sp][pop][trait][W] .+ X.Storage.Recorded.HostsTraits_infection_success_sd[sp][pop][trait][W], reverse(X.Storage.Recorded.HostsTraits_infection_success_mean[sp][pop][trait][W] .- X.Storage.Recorded.HostsTraits_infection_success_sd[sp][pop][trait][W]))
                    , seriestype=:shape,opacity=.08,fillcolor=:blue)  # plot polygons
                    plotS[end] = plot!(vcat(X.Storage.Recorded.Time[W], reverse(X.Storage.Recorded.Time[W]))
                    ,vcat(X.Storage.Recorded.ParasTraits_infection_success_mean[sp][pop][trait][W] .+ X.Storage.Recorded.ParasTraits_infection_success_sd[sp][pop][trait][W], reverse(X.Storage.Recorded.ParasTraits_infection_success_mean[sp][pop][trait][W] .- X.Storage.Recorded.ParasTraits_infection_success_sd[sp][pop][trait][W]))
                    , seriestype=:shape,opacity=.08,fillcolor=:red)  # plot polygons
                    plotS[end] = plot!(X.Storage.Recorded.Time[W], X.Storage.Recorded.HostsTraits_infection_success_mean[sp][pop][trait][W],color = :blue,line = [:solid,:dash,:dot][trait],lab = "H")
                    plotS[end] = plot!(X.Storage.Recorded.Time[W], X.Storage.Recorded.ParasTraits_infection_success_mean[sp][pop][trait][W],color = :red ,line = [:solid,:dash,:dot][trait],lab = "P")
                    ############ SUPR ############
                    # plotS[end] = plot!(X.Storage.Recorded.Time[findall(W)[[1,end]]], AllelesRangePerSp[sp][trait][[1,1]]  ,color = :grey ,line = :dash,lab = "")
                    # plotS[end] = plot!(X.Storage.Recorded.Time[findall(W)[[1,end]]], AllelesRangePerSp[sp][trait][[end,end]],color = :grey ,line = :dash,lab = "")
                end
            else
                push!(plotS, plot([1,10,100],[1,10,100],title = "pop = \$pop; sp = \$sp",color = :white ,lab = ""))
            end
        end
    end
    plot(plotS..., layout = Plots.GridLayout($NhostSp, $NPops),yaxis=:none,size = (max($NPops, $NPops*maximum(Time)/4e3)*2, $NhostSp).*400)
end"""
eval(Meta.parse(str))
# Evolution infection_success MEAN
str = """
function Plot_infection_success_MEAN(X::MetaCommunity;Time::UnitRange{Int64}=1:1)
    if Time == 1:1 Time = X.Storage.Recorded.Time[1]:X.Storage.Recorded.Time[end] end
    W = (X.Storage.Recorded.Time .>= Time[1]) .& (X.Storage.Recorded.Time .<= Time[end])
    plotS = []
    for sp in eachindex(X.MetaPopS)
        for pop in eachindex(X.MetaPopS[sp].PopS)
            if any(X.Storage.Recorded.Hosts[   sp][pop][W] .> 0) # & any(X.Storage.Recorded.Parasites[   sp][pop][W] .> 0)
                push!(plotS, plot(ylims = (0,1)))
                plotS[end] = plot!(X.Storage.Recorded.Time[W], X.Storage.Recorded.Pinfectionsuccess[sp][pop][W],color = :blue, line = :solid, lab = "",  title = "P(infection success) pop = \$pop; sp = \$sp",legendfont = font(8))
            else
                push!(plotS, plot([1,10,100],[1,10,100],title = "pop = \$pop; sp = \$sp",color = :white ,lab = ""))
            end
        end
    end
    plot(plotS..., layout = Plots.GridLayout($NhostSp, $NPops),yaxis=:none,size = ($NPops, $NhostSp).*400)
end"""
eval(Meta.parse(str))
# Evolution traits_infection_success SD
str = """
function Plot_traits_infection_success_SD(X::MetaCommunity;Time::UnitRange{Int64}=1:1)
    if Time == 1:1 Time = X.Storage.Recorded.Time[1]:X.Storage.Recorded.Time[end] end
    W = (X.Storage.Recorded.Time .>= Time[1]) .& (X.Storage.Recorded.Time .<= Time[end])
    plotS = []
    Yrange = []
    for sp in eachindex(X.MetaPopS)
        for pop in eachindex(X.MetaPopS[sp].PopS)
            if any(X.Storage.Recorded.Hosts[   sp][pop][W] .> 0) # & any(X.Storage.Recorded.Parasites[   sp][pop][W] .> 0)
                push!(plotS, plot())
                for trait in 1:$N_traits_infection_success
                    Yrange = [minimum([Yrange..., X.Storage.Recorded.HostsTraits_infection_success_sd[sp][pop][trait][W]..., X.Storage.Recorded.ParasTraits_infection_success_sd[sp][pop][trait][W]...]), maximum([Yrange..., X.Storage.Recorded.HostsTraits_infection_success_sd[sp][pop][trait][W]..., X.Storage.Recorded.ParasTraits_infection_success_sd[sp][pop][trait][W]...])]
                    plotS[end] = plot!(X.Storage.Recorded.Time[W], X.Storage.Recorded.HostsTraits_infection_success_sd[sp][pop][trait][W],color = :blue,line = [:solid,:dash,:dot][trait],lab = "",  title = "pop = \$pop; sp = \$sp",legendfont = font(8))
                    plotS[end] = plot!(X.Storage.Recorded.Time[W], X.Storage.Recorded.ParasTraits_infection_success_sd[sp][pop][trait][W],color = :red ,line = [:solid,:dash,:dot][trait],lab = "")
                end
            else
                push!(plotS, plot([1,10,100],[1,10,100],title = "pop = \$pop; sp = \$sp",color = :white ,lab = ""))
            end
        end
    end
    plot(plotS..., layout = Plots.GridLayout($NhostSp, $NPops),yaxis=:none,ylim=Yrange,size = ($NPops, $NhostSp).*400)
end"""
eval(Meta.parse(str))
# Evolution traits_immunity MEAN
str = """
function Plot_traits_immunity_MEAN(X::MetaCommunity;Time::UnitRange{Int64}=1:1)
    if Time == 1:1 Time = X.Storage.Recorded.Time[1]:X.Storage.Recorded.Time[end] end
    W = (X.Storage.Recorded.Time .>= Time[1]) .& (X.Storage.Recorded.Time .<= Time[end])
    plotS = []
    for sp in eachindex(X.MetaPopS)
        for pop in eachindex(X.MetaPopS[sp].PopS)
            if any(X.Storage.Recorded.Hosts[   sp][pop][W] .> 0) # & any(X.Storage.Recorded.Parasites[   sp][pop][W] .> 0)
                push!(plotS, plot(ylims = (0.95,$(NallelesPerTrait+0.05)),  title = "traits_immunity     pop = \$pop; sp = \$sp",legendfont = font(8) ))
                for trait in 1:$N_traits_immunity
                    plotS[end] = plot!(vcat(X.Storage.Recorded.Time[W], reverse(X.Storage.Recorded.Time[W]))
                    ,vcat(X.Storage.Recorded.HostsTraits_immunity_mean[sp][pop][trait][W] .+ X.Storage.Recorded.HostsTraits_immunity_sd[sp][pop][trait][W], reverse(X.Storage.Recorded.HostsTraits_immunity_mean[sp][pop][trait][W] .- X.Storage.Recorded.HostsTraits_immunity_sd[sp][pop][trait][W]))
                    , seriestype=:shape,opacity=.08,fillcolor=:blue)  # plot polygons
                    plotS[end] = plot!(vcat(X.Storage.Recorded.Time[W], reverse(X.Storage.Recorded.Time[W]))
                    ,vcat(X.Storage.Recorded.ParasTraits_immunity_mean[sp][pop][trait][W] .+ X.Storage.Recorded.ParasTraits_immunity_sd[sp][pop][trait][W], reverse(X.Storage.Recorded.ParasTraits_immunity_mean[sp][pop][trait][W] .- X.Storage.Recorded.ParasTraits_immunity_sd[sp][pop][trait][W]))
                    , seriestype=:shape,opacity=.08,fillcolor=:red)  # plot polygons
                    plotS[end] = plot!(X.Storage.Recorded.Time[W], X.Storage.Recorded.HostsTraits_immunity_mean[sp][pop][trait][W],color = :blue,line = [:solid,:dash,:dot][trait],lab = "")
                    plotS[end] = plot!(X.Storage.Recorded.Time[W], X.Storage.Recorded.ParasTraits_immunity_mean[sp][pop][trait][W],color = :red ,line = [:solid,:dash,:dot][trait],lab = "")
                    ############ SUPR ############
                    # plotS[end] = plot!(X.Storage.Recorded.Time[findall(W)[[1,end]]], AllelesRangePerSp[sp][trait+$N_traits_infection_success][[1,1]]  ,color = :grey ,line = :dash,lab = "")
                    # plotS[end] = plot!(X.Storage.Recorded.Time[findall(W)[[1,end]]], AllelesRangePerSp[sp][trait+$N_traits_infection_success][[end,end]],color = :grey ,line = :dash,lab = "")
                end
            else
                push!(plotS, plot([1,10,100],[1,10,100],title = "pop = \$pop; sp = \$sp",color = :white ,lab = ""))
            end
        end
    end
    plot(plotS..., layout = Plots.GridLayout($NhostSp, $NPops),yaxis=:none,size = ($NPops, $NhostSp).*400)
end"""
eval(Meta.parse(str))
# Evolution P_recovery_innate_immu MEAN
str = """
function Plot_P_recovery_innate_immu(X::MetaCommunity;Time::UnitRange{Int64}=1:1)
    if Time == 1:1 Time = X.Storage.Recorded.Time[1]:X.Storage.Recorded.Time[end] end
    W = (X.Storage.Recorded.Time .>= Time[1]) .& (X.Storage.Recorded.Time .<= Time[end])
    plotS = []
    for sp in eachindex(X.MetaPopS)
        for pop in eachindex(X.MetaPopS[sp].PopS)
            if any(X.Storage.Recorded.Hosts[   sp][pop][W] .> 0) # & any(X.Storage.Recorded.Parasites[   sp][pop][W] .> 0)
                push!(plotS, plot(ylims = (0,1)))
                plotS[end] = plot!(X.Storage.Recorded.Time[W], X.Storage.Recorded.Precoveryinnateimmu[sp][pop][W],color = :blue, line = :solid, lab = "",  title = "P(recovery innate immu) pop = \$pop; sp = \$sp",legendfont = font(8))
            else
                push!(plotS, plot([1,10,100],[1,10,100],title = "pop = \$pop; sp = \$sp",color = :white ,lab = ""))
            end
        end
    end
    plot(plotS..., layout = Plots.GridLayout($NhostSp, $NPops),yaxis=:none,size = ($NPops, $NhostSp).*400)
end"""
eval(Meta.parse(str))
# Evolution P_recovery_acquired_immu MEAN
str = """
function Plot_P_recovery_acquired_immu(X::MetaCommunity;Time::UnitRange{Int64}=1:1)
    if Time == 1:1 Time = X.Storage.Recorded.Time[1]:X.Storage.Recorded.Time[end] end
    W = (X.Storage.Recorded.Time .>= Time[1]) .& (X.Storage.Recorded.Time .<= Time[end])
    plotS = []
    for sp in eachindex(X.MetaPopS)
        for pop in eachindex(X.MetaPopS[sp].PopS)
            if any(X.Storage.Recorded.Hosts[   sp][pop][W] .> 0) # & any(X.Storage.Recorded.Parasites[   sp][pop][W] .> 0)
                push!(plotS, plot(ylims = (0,1)))
                plotS[end] = plot!(X.Storage.Recorded.Time[W], X.Storage.Recorded.Precoveryacquiredimmu[sp][pop][W],color = :blue, line = :solid, lab = "",  title = "P(recovery acquired immu) pop = \$pop; sp = \$sp",legendfont = font(8))
            else
                push!(plotS, plot([1,10,100],[1,10,100],title = "pop = \$pop; sp = \$sp",color = :white ,lab = ""))
            end
        end
    end
    plot(plotS..., layout = Plots.GridLayout($NhostSp, $NPops),yaxis=:none,size = ($NPops, $NhostSp).*400)
end"""
eval(Meta.parse(str))
# Evolution N_traits_immunity SD
str = """
function Plot_traits_immunity_SD(X::MetaCommunity;Time::UnitRange{Int64}=1:1,scaleX::Float,scaleY::Float)
    if Time == 1:1 Time = X.Storage.Recorded.Time[1]:X.Storage.Recorded.Time[end] end
    W = (X.Storage.Recorded.Time .>= Time[1]) .& (X.Storage.Recorded.Time .<= Time[end])
    plotS = []
    Yrange = []
    for sp in eachindex(X.MetaPopS)
        for pop in eachindex(X.MetaPopS[sp].PopS)
            if any(X.Storage.Recorded.Hosts[   sp][pop][W] .> 0) # & any(X.Storage.Recorded.Parasites[   sp][pop][W] .> 0)
                push!(plotS, plot())
                for trait in 1:$N_traits_immunity
                    Yrange = [minimum([Yrange..., NaNrm(X.Storage.Recorded.HostsTraits_immunity_sd[sp][pop][trait][W])..., NaNrm(X.Storage.Recorded.ParasTraits_immunity_sd[sp][pop][trait][W])...]), maximum([Yrange..., NaNrm(X.Storage.Recorded.HostsTraits_immunity_sd[sp][pop][trait][W])..., NaNrm(X.Storage.Recorded.ParasTraits_immunity_sd[sp][pop][trait][W])...])]
                    plotS[end] = plot!(X.Storage.Recorded.Time[W], X.Storage.Recorded.HostsTraits_immunity_sd[sp][pop][trait][W],color = :blue,line = [:solid,:dash,:dot][trait],lab = "",  title = "pop = \$pop; sp = \$sp",legendfont = font(8))
                    plotS[end] = plot!(X.Storage.Recorded.Time[W], X.Storage.Recorded.ParasTraits_immunity_sd[sp][pop][trait][W],color = :red ,line = [:solid,:dash,:dot][trait],lab = "")
                end
            else
                push!(plotS, plot([1,10,100],[1,10,100],title = "pop = \$pop; sp = \$sp",color = :white ,lab = ""))
            end
        end
    end
    plot(plotS..., layout = Plots.GridLayout($NhostSp, $NPops),yaxis=:none,ylim=Yrange,size = ($NPops*scaleX, $NhostSp*scaleY).*400)
end"""
eval(Meta.parse(str))
