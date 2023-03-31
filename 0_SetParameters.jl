# To validate this model
# Red queen dynamique :
# 	- Effet de la taille de la pop
# 	- Effet contradictoire sur le polymorphisme (effet seuil de la taille de la pop, en grande taille la taille de pop augmente; en petite taille, elle diminue)
# 	- Effet des interactions GfG, MA et leurs inverses
# 	  Géographie (mosaic theory of coevolution : augmente la divergence entre pop - traveaux de Sylvain Gandon)
# To test with this model
# - Effect of population size with continous matching alleles modèle
# - Effect of population size depending on sociality

# /DISQUE/FAC/Julia/julia-1.4.2-linux-x86_64/julia-1.4.2/bin/julia

# Set working directory
cd("/DISQUE/0_Grenoble/Biodiv_zoonoses/EvoSymbioses/")
# Set working directory
include("./1_LoadPkgs_importFunctions.jl");

#############################
# Compilation setting
# PRINT = ""  # Debbuging (printing MANY information during the computions to help identifying the origin of an error)
PRINT = "#" # No debbuging

# SIMD = "" # never use simd
SIMD = "@simd" # use simd for some loops (very few)

# FASTMATH = ""
FASTMATH = "@fastmath"

# NO_INLINE = ""
NO_INLINE = "@inline / @noinline"

# INBOUND = ""
INBOUND = "@inbounds"

# ALLOWSforLARGE_INACCURATE_dt = true
ALLOWSforLARGE_INACCURATE_dt = false
#############################

##################################################
############### Setting parameters ###############
##################################################
const NhostSp = 1
const NPops = 1

# Host population dynamics
const b0 = 5   # basal birth (the probability to have an offspring during one time step is  b0 * 2)
const d0 = 0.05 # death rates
const K = [[Int(1e2) for pop in 1:NPops] for sp in 1:NhostSp] # carrying capacity [for each species, one value per population]
# const K = [[Int(5e3), Int(0)] ,   [Int(1e3), Int(1e3)]] # carrying capacity [for each species, one value per population]
# K[1][1] *= 5
# K[2][1] *= 2

const thresholdApproxBinomial = 0.99

const PrevalenceIni = 0.75 # Start with a monomorphic parasit population with a prevalence of PrevalenceIni.
const N0 = [rand.(Poisson.(Int.(round.(K_./2 .* (1-PrevalenceIni))))) for K_ in K] # [for each sp[for each pop initial number of uninfected individuals]]
const N1 = [rand.(Poisson.(Int.(round.(K_./2 .* (  PrevalenceIni))))) for K_ in K] # [for each sp[for each pop initial number of   infected individuals]]
# N0[2][1] , N1[2][1] = 0, 0 # Initialy there is no species 2 in the population 1 (and vice versa but that is already done by the setting of the K object)
# N0 .+= N1 ; N1 -= N1

##  Transmission
const Ht = fill(50.0,NhostSp) # One value per host species. The average number of contacts between the parasites of one infected host and some other hosts when N=K
const HtCrossSp = 0.0 # The average number of contacts a host individual of the species SP_r (r for recipient) undergoes with the species SP_d (d for donor) when for both species N=K or whatever the values of N if 〖Sp〗_interact=1
const SpInteract = 0.5 # The proportion of contact induced by interspecies interactions that are independent of the species densities.
const Sociality = fill(0.75,NhostSp) # One value per host species. The proportion of these contacts that are independant of the host density (i.e. induced by social interactions)
const N_traits_infection_success = 1 # Number of trait which matching between the host and the parasites will determin P(infection|1 contact) noted Pinfection1contact

const Tol_infection = -10.0 # Shape of the relationship between Pinfection1contact and the distance between the host and parasit for the traits_infection_success
const Max_Pinfection1contact = 1.0
# const Min_Pinfection1contact = 0.0 by construction, it is set to zero.

##  Virulence
const Tol_virul = -500.0 # Shape of the relationship between Pinfection1contact and the virulence (the increase in the probability of death because of the infection)
const Max_virul = 0.1 # [0,1] Maximal virulence (this will be added to the death rate, the probability to die)
const Min_virul = 0.1 # [0,1] Minimal virulence (this will be added to the death rate, the probability to die)

## Host recovery
# Innate immunity
const N_traits_immunity = 1 # Number of trait which
# i)  matching between the host and the parasites will determine the probability to recover thank to innate immune system P(innate immu.->recovery) noted P_recovery_innate_immu
# ii) matching between the parasites genotypes Pg and Ph will determine the probability to recover thank to cross immuniny P(recovery_Pg|acquired immu._Ph) which will be stored in the tuple of tuples P_recovery_acquired_immu
const Tol_InnateImmu = -10.0 # Shape of the relationship between P_recovery_innate_immu and the minimum of the distances between the host and parasit for the traits_immunity
const Max_Precovery_innate_immu = 0.95

# Adaptive immune system and Cross-immunity
const SimulateImmunisation = true # Should the Adaptive immunity be simulated ?
const Tol_CrossImmun = 0.0 # Shape of the relationship between P(recovery_Pg|acquired immu._Ph) and the distance between the parasites genotypes Pg and Ph for the traits_immunity
const Max_Precovery_acquiered_immu = 0.25

### Evolution
const Maxrate_of_gametes_prod = 1.0       # Expected number of gamètes   produced during one time step by one host              which fitness                    is maximal. If GfG is false, then fitness is always maximal.
const MaxRate_of_parasites_emission = 1.0 # Expected number of parasites produced during one time step by one infected host for which the fitness of the parasit is maximal. If GfG is false, then fitness is always maximal.

const GfG = false# should an (inverse) Gene for Gene dynamic be modeled by allowing traits to appear and disappear ?
if GfG
    const MaxHostCostOfMissing_traits_infection_success = 0.0 # Value of rate_of_gametes_prod       if NtraitsInfectionSuccessPresent = 0                          when NtraitsImmunityPresent         = 0                            (NtraitsImmunityPresent         = 0                            maximize the value of rate_of_gametes_prod )
    const MaxHostCostOfHaving_traits_immunity = 2.0           # Value of rate_of_gametes_prod       if NtraitsImmunityPresent         = N_traits_immunity*2        when NtraitsInfectionSuccessPresent = N_traits_infection_success*2 (NtraitsInfectionSuccessPresent = N_traits_infection_success*2 maximize the value of rate_of_gametes_prod )
    const MaxParasCostOfHaving_traits_infection_success = 1.0 # Value of rate_of_parasites_emission if NtraitsInfectionSuccessPresent = N_traits_infection_success when NtraitsImmunityPresent         = N_traits_immunity            (NtraitsImmunityPresent         = N_traits_immunity            maximize the value of rate_of_parasites_emission)
    const MaxParasCostOfMissing_traits_immunity = 2.0         # Value of rate_of_parasites_emission if NtraitsImmunityPresent         = 0                          when NtraitsInfectionSuccessPresent = 0                            (NtraitsInfectionSuccessPresent = 0                            maximize the value of rate_of_parasites_emission)

    const TolHostCostOfMissing_traits_infection_success = 0.0 # Shape of the relationship between NtraitsInfectionSuccessPresent and Maxrate_of_gametes_prod
    const TolParasCostOfMissing_traits_immunity = 0.0         # Shape of the relationship between NtraitsImmunityPresent         and rate_of_gametes_prod
    const TolHostCostOfHaving_traits_immunity = 0.0           # Shape of the relationship between NtraitsImmunityPresent         and Maxrate_of_gametes_prod
    const TolParasCostOfHaving_traits_infection_success = 0.0 # Shape of the relationship between NtraitsInfectionSuccessPresent and rate_of_gametes_prod
end

## Evolution
const NallelesPerTrait = 6
const DistBetweenHostSp = 0 # number of alleles that separate each host species for the traits_infection_success
const PhenotypicStochasticity = false # Do we want to modèle PhenotypicStochasticity ?
# The following line only matters if PhenotypicStochasticity is set to true
const OnePhenNoiseTraitPerTrait = true # Should each trait have its own Phenotypic Noise ?

# Mutation rate for each trait
const MuHost_trait_alleles = [0.02 for sp in 1:NhostSp] # one mutation rate per species
const MuParasit_trait_alleles = 0.05 # mutation rate
const MuEffectHost = const MuEffectParasites = 0.5 # Average proportion of the range of alleles values covered by only one mutation
if GfG
    const MuHost_trait_presence = [0.001 for sp in 1:NhostSp] # one mutation rate per species
    const MuParasit_trait_presence = 0.001 # mutation rate
end
# For each species, probability that an individual migrate (for simplicity we assum that individuals can migrate to there own population - with a probability that depend on carrying capacities)
# No population structure corresponds to MigrRate = 1.0
MigrRate = [0.0 for sp in 1:NhostSp]


## Simulation parameters
# dt : time increment
dt = 0.01 # Float
# Maximal duration of the simulation to reach the equilibrium (if the equilibrium has not been reached after that time, the simulation stops any way)
# Whether the equilibrium has been reached or not is recorded in MetaCommunity.Storage.EquilibriumAllPop[]
const MaxTimeToEquilibrium = 20000.0 # Float
# Duration of the simulation once equilibrium as been reached (same unit as dt)
const TimeAfterEquilibrium = 1000.0 # Float
# AfterEquilibrium, at which frequency the state of the population should be recorded ?
const RecordEvery = dt*10 # 5.0 # Float
# Should only the dynamic after Equilibrium be kept ?
OnlyRecordDynamicAfterEquilibrium = true # true

# Number of repeatition of each simulation
Nsimulations = 1
# Number of repeatition of each simulation
DiscardeSimulationIfEquilibriumNoReached = true

## RecordStatistics : function that will be applied at the end of each simulation to record some statistics
RecordedStatistics = "./out.txt" # name of the file where the statistics will be written
NamesOfRecordedStatistics = flat([[["Hosts$(sp)He_immunity_pop$(pop)","Hosts$(sp)He_infection_success_pop$(pop)","Parasites_in_H$(sp)_He_immunity_pop$(pop)","Parasites_in_H$(sp)_He_infection_success_pop$(pop)"] for sp in 1:NhostSp] for pop in 1:NPops])
# set these variables to 'nothing' if you don't want to record some statistics
function RecordStatistics(X)
    function GetP_He(trait::Int,X,pop::Int,sp::Int)
        alleles = unique([g.Alleles[trait] for g in X.parasitesGenotypeS_List.GenotypeS])
        if length( X.MetaPopS[sp].PopS[pop].ParasitesPop.IDpS ) > 0
            Freq::Vector{Union{Int,Float}} = [ sum([if (g.Alleles[trait] === allele) sum(X.MetaPopS[sp].PopS[pop].ParasitesPop.N_[X.MetaPopS[sp].PopS[pop].ParasitesPop.IDpS .=== g.IDp]) else 0 end for g in X.parasitesGenotypeS_List.GenotypeS]) for allele in alleles ]
            Freq ./= sum(Freq)
            return(1-sum(Freq.^2))
        else
            return 0.0
        end
    end
    out = Vector{Float}()
    for pop in 1:NPops
        for sp in 1:NhostSp
            # "Hosts$(sp)He_immunity"
            push!(out, mean([ 1-sum((sum([X.MetaPopS[sp].PopS[pop].HostSexS[sex].AllelesPool[trait] for sex in 1:2])./2).^2) for trait in 1:N_traits_immunity]) )
            # "Hosts$(sp)He_infection_success"
            push!(out, mean([ 1-sum((sum([X.MetaPopS[sp].PopS[pop].HostSexS[sex].AllelesPool[trait] for sex in 1:2])./2).^2) for trait in N_traits_immunity:(N_traits_immunity+N_traits_infection_success)]))
            # "Parasites_in_H$(sp)_He_immunity"
            push!(out, mean([ GetP_He(trait,X,pop,sp) for trait in 1:N_traits_immunity]))
            # "Parasites_in_H$(sp)_He_infection_success"
            push!(out, mean([ GetP_He(trait,X,pop,sp) for trait in N_traits_immunity:(N_traits_immunity+N_traits_infection_success)]))
        end
    end
    # The output of this function must be a Vector which values will be appended to the file "RecordedStatistics".
    # They will be separated by tabulations ("\t").
    return(out)
end

include("./2_ProcessParameters.jl");

PlotParameters(legendfontsize=7)

# Run the simulations
include("./4_RunSimulations.jl");

# Explore the simulation (assuming Nsimulations == 1, otherwise you are only exploring the last simulation)
X.Storage.EquilibriumAllPop[]

X.Storage.Extinction


Plot_Demography(X;ylog = false, Time=2:1000,scaleX=5.0,scaleY=0.75)
Plot_traits_infection_success_FREQ_ALLELES(X,trait=1,scaleX=5.0,scaleY=0.75)
Plot_traits_immunity_FREQ_ALLELES(X,trait=1,scaleX=5.0,scaleY=0.75)


plot( X.Storage.Recorded.HostsTraits_infection_success_freq[1][1][1][1][2:end],trait=1)
plot!(X.Storage.Recorded.ParasTraits_infection_success_freq[1][1][1][1][2:end],trait=1)

plot( X.Storage.Recorded.HostsTraits_immunity_freq[1][1][1][1][2:end],trait=1)
plot!(X.Storage.Recorded.ParasTraits_immunity_freq[1][1][1][1][2:end],trait=1)

### store Alleles frequecies
File = "./Alleles_frequecies.csv" ,
open(File, truncate =true) do f,
    write(f, join(,
            ["Nhosts",
             "FreqImmH".* string.(1:NallelesPerTrait)...,
             "FreqInfH".* string.(1:NallelesPerTrait)...,
             "Nparasits",
             "FreqImmP".* string.(1:NallelesPerTrait)...,
             "FreqInfP".* string.(1:NallelesPerTrait)...,
            ],
            ,";"))
    
     for i in 1:length(X.Storage.Recorded.Hosts[1][1])    ,
        write(f, *join( [ X.Storage.Recorded.Hosts[1][1][i],
                            ,[X.Storage.Recorded.HostsTraits_immunity_freq[1         ][1][1][allele][i] for allele in 1:NallelesPerTrait]...,
                            ,[X.Storage.Recorded.HostsTraits_infection_success_freq[1][1][1][allele][i] for allele in 1:NallelesPerTrait]...,
                            , X.Storage.Recorded.Parasites[1][1][i],
                            ,[X.Storage.Recorded.ParasTraits_immunity_freq[         1][1][1][allele][i] for allele in 1:NallelesPerTrait]...,
                            ,[X.Storage.Recorded.ParasTraits_infection_success_freq[1][1][1][allele][i] for allele in 1:NallelesPerTrait]...,
                ]
            ,";" ) )
     end
end


Plot_He(X;ylog = false)
Plot_Demography(X;ylog = false, Time=2:1000)
Plot_traits_infection_success_FREQ_ALLELES(X,trait=1)
Plot_traits_infection_success_FREQ_ALLELES(X,trait=2)
Plot_traits_infection_success_FREQ_ALLELES(X,trait=3)

Plot_infection_success_MEAN(X)

Plot_traits_immunity_FREQ_ALLELES(X,trait=1)
Plot_traits_immunity_FREQ_ALLELES(X,trait=2)
Plot_P_recovery_innate_immu(X)
Plot_P_recovery_acquired_immu(X)

Plot_Covariation_fitn_evol_traits_infection_success(X,trait=1,Time=2:1000)
Plot_Covariation_fitn_evol_traits_infection_success(X,trait=2)
Plot_Covariation_fitn_evol_traits_immunity(X,trait=1,Time=2:1000)
Plot_Covariation_fitn_evol_traits_immunity(X,trait=2)
