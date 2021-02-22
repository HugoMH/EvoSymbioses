for sim in 1:Nsimulations
    global X = deepcopy(metaCommunity) ;
    EvolveOver_Or_UntilEquilibrium!(X, dt=dt, MaxDuration=MaxTimeToEquilibrium, RecordEvery=max(RecordEvery,MaxTimeToEquilibrium/5000.0), SpPop_WhereParasitesMustNotDisapear = :All)
    if OnlyRecordDynamicAfterEquilibrium
        X.Storage.TimeEvolved[] = metaCommunity.Storage.TimeEvolved[]
        X.Storage.NextTimeRecord[] = metaCommunity.Storage.NextTimeRecord[]
        for f in fieldnames(typeof(metaCommunity.Storage.Recorded))
            empty!(getfield(X.Storage.Recorded, f))
            append!(getfield(X.Storage.Recorded, f),  deepcopy(getfield(metaCommunity.Storage.Recorded, f)))
        end
    end
    println(X.Storage.EquilibriumAllPop[])
    if (!DiscardeSimulationIfEquilibriumNoReached) || (DiscardeSimulationIfEquilibriumNoReached & X.Storage.EquilibriumAllPop[])
        println(X.Storage.EquilibriumAllPop[])
        EvolveOver!(                X, dt=dt, Duration=TimeAfterEquilibrium   , RecordEvery=RecordEvery                                 , SpPop_WhereParasitesMustNotDisapear = :All )
        println(X.Storage.EquilibriumAllPop[])
        if !isnothing(RecordedStatistics)
            open(RecordedStatistics, "a") do f
                write(f, "\n"*join(RecordStatistics(X),"\t"))
            end
        end
    end
    println(X.Storage.EquilibriumAllPop[])
end
