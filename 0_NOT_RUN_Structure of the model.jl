"""
This file is not intended to be run
It represent a simplified version of the structure of the model
"""

:MetaCommunity
             |
             |_ :MetaPop
             |
             |
             |_ ...  "one 'MetaPop' per host species"
             |
             |
             |_ :MetaPop
                        |
                        |_ :HostSexS
                        |
                        |
                        |_ ... "one 'HostSexS' per population"
                        |
                        |
                        |_ :HostSexS
                                   |
                                   |_ :ParasitesPop
                                   |              |
                                   |              2×6 DataFrame
                                   |              │ Row │ IDp_HostShiftHistoryS         │ IDpS  │ IDpSphen │ N_    │ posiInParasitGenotypeS_List │ ParasitFreq │
                                   |              │     │ Symbol                        │ Int64 │ Int64    │ Int64 │ Int64                       │ Float16     │
                                   |              ├─────┼───────────────────────────────┼───────┼──────────┼───────┼─────────────────────────────┼─────────────┤
                                   |              │ 1   │ IDp_4__Sp°1°°Time°0.0°°IDh°26 │ 4     │ 4        │ 1017  │ 2                           │ 0.553       │
                                   |              │ 2   │ IDp_6__Sp°1°°Time°0.0°°IDh°26 │ 6     │ 6        │ 821   │ 4                           │ 0.4468      │
                                   |              ⋮
                                   |
                                   |
                                   |_ :HostsGenotypeS   "(males   : always first  position)"
                                   |_ :HostsGenotypeS   "(females : always second position)"
                                                    |
                                                    151×14 DataFrame
                                                    │ Row │ IDhS  │ IDhSphen │ IDpS  │ IDpSphen │ IDp_HostShiftHistoryS         │ IDiS             │ IDcategorieS                                                                                 │ virulenceS │ Precoveryinnateimmu │ Precoveryacquiredimmu │ N_    │ dN_   │ posiInHostsGenotypeS_List │ posiInParasitGenotypeS_List │
                                                    │     │ Int64 │ Int64    │ Int64 │ Int64    │ Symbol                        │ Symbol           │ Symbol                                                                                       │ Float64    │ Float64             │ Float64               │ Int64 │ Int64 │ Int64                     │ Int64                       │
                                                    ├─────┼───────┼──────────┼───────┼──────────┼───────────────────────────────┼──────────────────┼──────────────────────────────────────────────────────────────────────────────────────────────┼────────────┼─────────────────────┼───────────────────────┼───────┼───────┼───────────────────────────┼─────────────────────────────┤
                                                    │ 1   │ 26    │ 26       │ 0     │ 0        │ IDparasitCat_0_0              │ _0               │ IDcat_26_26_IDparasitCat_0_0                                                                 │ 0.0        │ 0.0                 │ 0.0                   │ 29    │ 0     │ 1                         │ 0                           │
                                                    │ 2   │ 26    │ 26       │ 4     │ 4        │ IDp_4__Sp°1°°Time°0.0°°IDh°26 │ _0               │ IDcat_26_26_IDparasitCat_4_4_HostShiftHistory_IDp_4__Sp°1°°Time°0.0°°IDh°26__0               │ 5.0        │ 0.8                 │ 0.0                   │ 96    │ 0     │ 1                         │ 2                           │
                                                    │ 3   │ 18    │ 18       │ 0     │ 0        │ IDparasitCat_0_0              │ _0               │ IDcat_18_18_IDparasitCat_0_0                                                                 │ 0.0        │ 0.0                 │ 0.0                   │ 28    │ 0     │ 5                         │ 0                           │
                                                    ⋮
                                                    │ 149 │ 26    │ 26       │ 0     │ 0        │ IDparasitCat_0_0              │ _0°2°2°2°2°2°4   │ IDcat_26_26_IDparasitCat_0_0__0°2°2°2°2°2°4                                                  │ 0.0        │ 0.0                 │ 0.0                   │ 1     │ 0     │ 1                         │ 0                           │
                                                    │ 150 │ 18    │ 18       │ 0     │ 0        │ IDparasitCat_0_0              │ _0°2°4           │ IDcat_18_18_IDparasitCat_0_0__0°2°4                                                          │ 0.0        │ 0.0                 │ 0.0                   │ 1     │ 0     │ 5                         │ 0                           │
                                                    │ 151 │ 26    │ 26       │ 0     │ 0        │ IDparasitCat_0_0              │ _0°1°2°2°2°4°4   │ IDcat_26_26_IDparasitCat_0_0__0°1°2°2°2°4°4                                                  │ 0.0        │ 0.0                 │ 0.0                   │ 1     │ 0     │ 1                         │ 0                           │
