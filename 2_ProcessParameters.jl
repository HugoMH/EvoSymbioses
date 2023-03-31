
#####################################################
############### Processing parameters ###############
#####################################################
const NhostSp = Int(NhostSp)
const NallelesPerTrait = Int(NallelesPerTrait)
const N_traits_infection_success = Int(N_traits_infection_success)
const N_traits_immunity = Int(N_traits_immunity)
const DistBetweenHostSp = Int(DistBetweenHostSp)

const NtraitsPhenOnly = N_traits_infection_success + N_traits_immunity

if PhenotypicStochasticity
    if OnePhenNoiseTraitPerTrait
        global const Ntraits = NtraitsPhenOnly * 2
    else
        global const Ntraits = NtraitsPhenOnly + 1
    end
else
    global const Ntraits = NtraitsPhenOnly
end

const N_gen_X_Phen_X_H_X_P = Ntraits*3  *2 # Number of characteristics of each category (exept immunisations) : each trait for each host chromosome + each trait for parasit + the same for the phenotypes (=> *2)


include("./3_functions.jl");
# Evolution

# Host population dynamics
const density_dependence = [[(b0-d0)/(2×K[sp][pop]) for pop in 1:NPops] for sp in 1:NhostSp] # carrying capacity [for each population, one value per species]
const dt_Max = 1/((d0+ maximum(maximum.([d.*k for (d, k) in zip(density_dependence , K)])) ) )

##################################
##### Set the data structure #####
##################################

hostGenotype_sp_pop    = [[RandHostsGenotype(   sp) for pop in 1:NPops] for sp in 1:NhostSp]
parasitGenotype_sp_pop = [[RandParasitesGenotype( ) for pop in 1:NPops] for sp in 1:NhostSp]
##########################
str = """
hostsGenotypes_pop_sp = [[HostsGenotypeS(
 Vector{Int_IDh}() # IDhS::Vector{Int_IDh}
,Vector{Int_IDh}() # IDhSphen::Vector{Int_IDh}
,Vector{Int}()     # posiInHostsGenotypeS_List::Vector{Int}
,Vector{Int}()     # posiInParasitGenotypeS_List::Vector{Int}
,Vector{Int_IDp}() # IDpS::Vector{Int_IDp}
,Vector{Int_IDp}() # IDpSphen::Vector{Int_IDp}
,Vector{Symbol}()  # IDpGenPhen_HostShiftHistoryS::Vector{Symbol}
,Vector{Symbol}()  # IDp_HostShiftHistoryS::Vector{Symbol}
$(ifelse(SimulateImmunisation
,",Vector{Symbol}()" # IDiS::Vector{Int}
,""))
,Vector{Symbol}() # IDcategorieS
$(ifelse(GfG
,",Vector{Float}() # rate_of_gametes_prod
  ,Vector{Float}() # rate_of_parasites_emission"
,""
))
,Vector{Float}()  # virulenceS
,Vector{Float}()  # Precoveryinnateimmu
$(ifelse(SimulateImmunisation
,",Vector{Float}()" # Precoveryacquiredimmu
,""))
,Ref{Int}()       # N
,Vector{Int}()    # N_
,Vector{Int}()    # dN_
$(ifelse(GfG # AllelesPool
    ,",tuple([if (trait <= N_traits_infection_success) Vector{Float}(undef,NallelesPerTrait_infection_success_AllSp + 1) else Vector{Float}(undef,NallelesPerTrait + 1) end for trait in 1:Ntraits]...)" # include the zero allele
    ,",tuple([if (trait <= N_traits_infection_success) Vector{Float}(undef,NallelesPerTrait_infection_success_AllSp    ) else Vector{Float}(undef,NallelesPerTrait    ) end for trait in 1:Ntraits]...)"))
,Vector{Int}()    # Ngametes
,Vector{Int}()    # NparasitesEmitted
$(ifelse(GfG # nMut: Number of mutants for the current allele and trait that decrease [1] and that increase [2] or for which the trait presence is swapped [3]
    ,",zeros(Int,3)"
    ,",zeros(Int,2)"))
$(ifelse(GfG # MutatedAlleles
,",tuple([if (trait <= N_traits_infection_success) zeros(Int,NallelesPerTrait_infection_success_AllSp + 1) else zeros(Int,NallelesPerTrait + 1) end for trait in 1:Ntraits]...)"
,",tuple([if (trait <= N_traits_infection_success) zeros(Int,NallelesPerTrait_infection_success_AllSp    ) else zeros(Int,NallelesPerTrait    ) end for trait in 1:Ntraits]...)"))
,Vector{Int}(undef,NallelesPerTrait-1) # AllelesChanges
,Vector{Int}() # Nrecovery
,Vector{Int}(undef,NPops) # NmigrEachPop
,Vector{Int}() # NmigrEachGenotype
) for pop in 1:NPops
] for sp in 1:NhostSp
]"""
eval(Meta.parse(str))

##########################
str = """
parasitesPop_pop_sp = [[ParasitesPop(
Vector{Int_IDp}() # parasitesGenotypes
,Vector{Int_IDp}() # parasitesPhenotypes
,Vector{Symbol}()  # IDpGenPhen_HostShiftHistoryS
,Vector{Symbol}()  # IDp_HostShiftHistoryS
,Vector{Int}()     # posiInParasitGenotypeS_List
,Vector{Float16}() # ParasitFreq : RELATIVE frequency of each parasit
,Ref(0) # N
,[0] # N_
$(ifelse(GfG, ",Vector{Int}(undef,3)", ",Vector{Int}(undef,2)" )) # Nmutating  [decrease, increase, swap the presence-absence of a trait]
,Vector{Int}(undef,Ntraits) # NmutatingEachTrait
,tuple([Vector{Int}(undef,(Nalleles-1)) for Nalleles in NallelesPerTrait_EachTrait]...) # NmuSize
$(ifelse(GfG
,",tuple([Vector{Int}(undef,(Nalleles  )) for Nalleles in NallelesPerTrait_EachTrait]...) # NmutEachAllelesAppearingTrait
  ,Vector{Int}() #rate_of_parasites_emissionS"
,""))
,Vector{Int}(undef,Ntraits) # NewGenotype
,Ref(0) # NewId
) for pop in 1:NPops
] for sp in 1:NhostSp
]"""
eval(Meta.parse(str))
##########################
hostSexS_pop_sp = [[HostSexS(
tuple(deepcopy(hostsGenotypes_pop_sp[sp][pop]),deepcopy(hostsGenotypes_pop_sp[sp][pop])) # HostsGenotypeS
,tuple(1,2) # 1, 2 = males, females
,Ref(0) # N
,[0,0]# N_
,Ref(Float(K[sp][pop])) # K
,Ref(Float(density_dependence[sp][pop])) # density_dependence
,Float(Sociality[sp]) # Sociality
,Float(Ht[sp]) # Ht
,Float(Ht[sp]*(Sociality[sp])) # Ht_Sociality
,Ref(Float(Ht[sp]*(1-Sociality[sp])/K[sp][pop])) # Ht_one_Sociality_K
,parasitesPop_pop_sp[sp][pop] # ParasitesPop
,Ref{Int}() # NcontactsInPop For the ongoing pop, number of contact between individuals
,fill(0,NhostSp) # NcontactsCrossSp
,Vector{Float}() # P_NcontactsInPop_1_2_3_ect_parasites
,Vector{Float}() # ParasitFreq_x_Pinfection1contact
,Vector{Float}() # Pinfection
,Vector{Float}() # virulenceScompar
,Vector{Float16}() # P_HigherVirulenceS
,Vector{Float16}() # P_SameVirulenceS
,Vector{Int}() # Posi_LowerVirulenceS_comparToVirulenceInHost
,Vector{Int}() # NsuccessfullInfectionEachParasit
,Vector{Int}(undef,NtraitsPhenOnly) # RefPhen
,Ref(HtCrossSp) # HtCrossSp
,Ref(SpInteract) # SpInteract
) for pop in 1:NPops
] for sp in 1:NhostSp
]##########################

if length(MigrRate) != NhostSp
    error("length of MigrRate must be the same as NhostSp")
end

metaPop = [MetaPop(
[deepcopy(hostSexS_pop_sp[sp][pop])    for pop in 1:NPops] # PopS
,[pop                                   for pop in 1:NPops] # IDpopS
,sp # IDsp
,Ref(0) # N
,zeros(Int,NPops) # N_
,Ref(MigrRate[sp]) # MigrRate
,Vector{Float}(undef,NPops) # PmigrSettlementEachPop
) for sp in 1:NhostSp
]##########################
str = """
hostsGenotypeS_List = HostsGenotypeS_List(
Vector{Int}()
,Vector{HostsGenotype}()
,$(join(["ExtendableDict(Dict(0 => " for traitPosi in 1:((Ntraits*2)-1)])) ExtendableDict(Dict{Int,Int}( $(join(["))" for traitPosi in 1:(Ntraits*2)])) # HostSexOffspringMap
)"""
eval(Meta.parse(str))
##########################
parasitesGenotypeS_List = ParasitesGenotypeS_List(
Vector{Int}()
,Vector{ParasitesGenotype}()
)##########################
str = """
hpinteractions = HPinteractions(
     Dict{Symbol, NTuple{ N_traits_immunity         , Int}}() # TraitsValues_infection_success
    ,Dict{Symbol, NTuple{ N_traits_infection_success, Int}}() # TraitsValues_immunity
    ,Dict{Symbol, Vector{Float}}() # Pinfection1contact                  for each ParasitesGenotype, what is its Pinfection1contact on this HostsGenotype ?
    ,Dict{Symbol, Vector{Float}}() # Pinfection1contactFloat16  Float16( for each ParasitesGenotype, what is its Pinfection1contact on this HostsGenotype ? )
    ,Dict{Symbol, Vector{Float}}() # virulenceS                          for each ParasitesGenotype, what is its virulence on this HostsGenotype ?
    ,Dict{Symbol, Vector{Float}}() # Precoveryinnateimmu                 for each ParasitesGenotype, what is its P_recovery_innate_immu on this HostsGenotype ?
    $(ifelse(SimulateImmunisation
    ,",$(join(["ExtendableDict(Dict(0 => " for traitPosi in 1:((N_traits_immunity*2)-1)])) ExtendableDict(Dict{Int,Float}( $(join(["))" for traitPosi in 1:(N_traits_immunity*2)])) # Precoveryacquiredimmu
      ,Dict{Symbol,Immunisations}() # immunisationS : For each immunisationsType, if it recovers from parasit genotype 1, 2, 3 ..., this which immunisationsType this gives ?
      ,Vector{Int}(undef,N_traits_immunity) # TempDistImmunisingInfecting"
    ,""))
)"""
eval(Meta.parse(str))
##########################
str = """
metaCommunity = MetaCommunity(
    tuple([deepcopy(metaPop[sp]) for sp in 1:NhostSp]...) # MetaPopS
    ,tuple([sp                   for sp in 1:NhostSp]...) # IDspS
    ,Ref(0) # N
    ,zeros(Int,NhostSp) # N_
    ,deepcopy(hostsGenotypeS_List)# HostsGenotypeS_List
    ,deepcopy(parasitesGenotypeS_List) # ParasitesGenotypeS_List
    #
    ,deepcopy(hpinteractions)
    #
    ,ntuple(trait -> Vector{ExtendableDict{Int,Int}}() , Ntraits) # ParasitMutatingMap  [trait][current_genotype][sens_of_the_mutation_1=decrease_2=increase] #  POSITION OF GENOTYPES NOT ID
    ,(# Storage
        TimeEvolved = Ref(0.0)
        ,NextTimeRecord = Ref(0.0)
        ,Recorded = (           RealWordTime                       =     Vector{DateTime}()
                               ,Time                               =     Vector{  Float }()
                               ,Parasites                          = [[  Vector{Int}(     )                                                                                                                                           for pop in 1:NPops] for sp in 1:NhostSp]
                               ,Hosts                              = [[  Vector{Int}(     )                                                                                                                                           for pop in 1:NPops] for sp in 1:NhostSp]
                               ,Precoveryinnateimmu                = [[  Vector{Float}(   )                                                                                                                                           for pop in 1:NPops] for sp in 1:NhostSp]
$(ifelse(SimulateImmunisation,",Precoveryacquiredimmu              = [[  Vector{Float}(   )                                                                                                                                           for pop in 1:NPops] for sp in 1:NhostSp]",""))
                               ,Pinfectionsuccess                  = [[  Vector{Float}(   )                                                                                                                                           for pop in 1:NPops] for sp in 1:NhostSp]
                               ,HostsTraits_infection_success_freq = [[[[Vector{Float}() for allele in 1:(NallelesPerTrait_EachTrait[trait]                            + if_GfG_1_else_0)] for trait in 1:N_traits_infection_success] for pop in 1:NPops] for sp in 1:NhostSp]
                               ,ParasTraits_infection_success_freq = [[[[Vector{Float}() for allele in 1:(NallelesPerTrait_EachTrait[trait]                            + if_GfG_1_else_0)] for trait in 1:N_traits_infection_success] for pop in 1:NPops] for sp in 1:NhostSp]
                               ,HostsTraits_immunity_freq          = [[[[Vector{Float}() for allele in 1:(NallelesPerTrait_EachTrait[trait+N_traits_infection_success] + if_GfG_1_else_0)] for trait in 1:N_traits_immunity         ] for pop in 1:NPops] for sp in 1:NhostSp]
                               ,ParasTraits_immunity_freq          = [[[[Vector{Float}() for allele in 1:(NallelesPerTrait_EachTrait[trait+N_traits_infection_success] + if_GfG_1_else_0)] for trait in 1:N_traits_immunity         ] for pop in 1:NPops] for sp in 1:NhostSp]
                               ,HostsTraits_infection_success_fitn = [[[[Vector{Float}() for allele in 1:(NallelesPerTrait_EachTrait[trait]                            + if_GfG_1_else_0)] for trait in 1:N_traits_infection_success] for pop in 1:NPops] for sp in 1:NhostSp]
                               ,ParasTraits_infection_success_fitn = [[[[Vector{Float}() for allele in 1:(NallelesPerTrait_EachTrait[trait]                            + if_GfG_1_else_0)] for trait in 1:N_traits_infection_success] for pop in 1:NPops] for sp in 1:NhostSp]
                               ,HostsTraits_immunity_fitn          = [[[[Vector{Float}() for allele in 1:(NallelesPerTrait_EachTrait[trait+N_traits_infection_success] + if_GfG_1_else_0)] for trait in 1:N_traits_immunity         ] for pop in 1:NPops] for sp in 1:NhostSp]
                               ,ParasTraits_immunity_fitn          = [[[[Vector{Float}() for allele in 1:(NallelesPerTrait_EachTrait[trait+N_traits_infection_success] + if_GfG_1_else_0)] for trait in 1:N_traits_immunity         ] for pop in 1:NPops] for sp in 1:NhostSp]
                               ,HostsHe                            = [[  Vector{Float}()                                                                                                                                              for pop in 1:NPops] for sp in 1:NhostSp]
                               ,ParasitesHe                        = [[  Vector{Float}()                                                                                                                                              for pop in 1:NPops] for sp in 1:NhostSp]
        )
    ,Extinction        = Ref(false)
    ,EquilibriumAllPop = Ref(false)
    )
)"""
eval(Meta.parse(str))
##########################

## reset these objects to ensure that functions are not using them
metaPop = hostsGenotypeS_List = parasitesGenotypeS_List = parasitesGenotypeS_List = nothing

for sp in 1:NhostSp
    for pop in 1:NPops
        for sex in 1:2
            println(sp,pop,sex)
            HostParasitesPhenotypes_INI!(metaCommunity # uninfecteds
            ,metaCommunity.MetaPopS[sp].PopS[pop].HostSexS[sex] # hostsGenotypeS
            ,hostGenotype_sp_pop[sp][pop].Alleles # hostGenotype
            ,nothing # parasitGenotype
            ,Int(round(N0[sp][pop]/2)) # N
            ,sp
            )
            HostParasitesPhenotypes_INI!(metaCommunity # infecteds
            ,metaCommunity.MetaPopS[sp].PopS[pop].HostSexS[sex] # hostsGenotypeS
            ,hostGenotype_sp_pop[sp][pop].Alleles # hostGenotype
            ,parasitGenotype_sp_pop[sp][pop].Alleles # parasitGenotype
            ,Int(round(N1[sp][pop]/2)) # N
            ,sp
            )
            println(metaCommunity.MetaPopS[sp].PopS[pop].HostSexS[sex])
        end
    end
end

SetHostParasitInteractionParameters!(metaCommunity) # Pinfection1contact, virulenceS, P_recovery_innate_immu
SetNParasites!(metaCommunity)

SetParasitMutatingMap!(metaCommunity)
CheckIfStop__apply_dN__RecordState!(metaCommunity)



# create the output file
if !isnothing(RecordedStatistics)
    open(RecordedStatistics, "w") do f
        write(f, join(NamesOfRecordedStatistics,"\t"))
    end
end
