export Pollu

struct PolluCoupler
    sys::Any
end

"""
    Pollu()

This atmospheric chemical system model is built based on this paper:

"GAUSS-SEIDEL ITERATION FOR STIFF ODES FROM CHEMICAL KINETICS" (1994), J. G. VERWER.

The model implements a 20-species, 25-reaction atmospheric chemistry ODE system from
the Appendix of the paper. Species concentrations are in ppb and time is in seconds.
The original rate constants are given in min⁻¹ (first order) and ppm⁻¹min⁻¹ (second order)
and are converted to s⁻¹ and (ppb⁻¹s⁻¹) respectively.

If the keyword argument `rxn_sys` is set to `true`, the function will return a reaction system instead of an ODE system.

# Example

```
using GasChem, EarthSciMLBase, DifferentialEquations, Plots
rs = Pollu()
sol = solve(ODEProblem(mtkcompile(rs), [], (0,360), []), AutoTsit5(Rosenbrock23()), saveat=10.0)
plot(sol)
```
"""
function Pollu(; name = :Pollu, rxn_sys = false)
    rx_sys = @reaction_network Pollu begin
        @ivs t [unit = u"s"]
        @parameters begin
            jNO2_O3P = 0.35 / 60.0, [unit = u"s^-1"]                     # Rxn 1: .350E+00 min⁻¹
            k2 = 0.266e2 / (60.0 * 1e3), [unit = u"1/(ppb*s)"]           # Rxn 2: .266E+02 ppm⁻¹min⁻¹
            k3 = 0.120e5 / (60.0 * 1e3), [unit = u"1/(ppb*s)"]           # Rxn 3: .120E+05 ppm⁻¹min⁻¹
            jH2COa = 0.86e-3 / 60.0, [unit = u"s^-1"]                    # Rxn 4: .860E-03 min⁻¹
            jH2COb = 0.82e-3 / 60.0, [unit = u"s^-1"]                    # Rxn 5: .820E-03 min⁻¹
            k6 = 0.15e5 / (60.0 * 1e3), [unit = u"1/(ppb*s)"]            # Rxn 6: .150E+05 ppm⁻¹min⁻¹
            jALD = 0.13e-3 / 60.0, [unit = u"s^-1"]                      # Rxn 7: .130E-03 min⁻¹
            k8 = 0.24e5 / (60.0 * 1e3), [unit = u"1/(ppb*s)"]            # Rxn 8: .240E+05 ppm⁻¹min⁻¹
            k9 = 0.165e5 / (60.0 * 1e3), [unit = u"1/(ppb*s)"]           # Rxn 9: .165E+05 ppm⁻¹min⁻¹
            k10 = 0.9e4 / (60.0 * 1e3), [unit = u"1/(ppb*s)"]            # Rxn 10: .900E+04 ppm⁻¹min⁻¹
            jPAN = 0.22e-1 / 60.0, [unit = u"s^-1"]                      # Rxn 11: .220E-01 min⁻¹
            k12 = 0.12e5 / (60.0 * 1e3), [unit = u"1/(ppb*s)"]           # Rxn 12: .120E+05 ppm⁻¹min⁻¹
            k13 = 1.88 / 60.0, [unit = u"1/s"]                           # Rxn 13: .188E+01 min⁻¹
            k14 = 0.163e5 / (60.0 * 1e3), [unit = u"1/(ppb*s)"]          # Rxn 14: .163E+05 ppm⁻¹min⁻¹
            k15 = 0.48e7 / 60.0, [unit = u"s^-1"]                        # Rxn 15: .480E+07 min⁻¹
            jO3_O1D = 0.35e-3 / 60.0, [unit = u"s^-1"]                   # Rxn 16: .350E-03 min⁻¹
            jO3_O3P = 0.175e-1 / 60.0, [unit = u"s^-1"]                  # Rxn 17: .175E-01 min⁻¹
            k18 = 0.1e9 / 60.0, [unit = u"s^-1"]                         # Rxn 18: .100E+09 min⁻¹
            k19 = 0.444e12 / 60.0, [unit = u"s^-1"]                      # Rxn 19: .444E+12 min⁻¹
            k20 = 0.124e4 / (60.0 * 1e3), [unit = u"1/(ppb*s)"]          # Rxn 20: .124E+04 ppm⁻¹min⁻¹
            jNO3_NO = 0.21e1 / 60.0, [unit = u"s^-1"]                    # Rxn 21: .210E+01 min⁻¹
            jNO3_NO2 = 0.578e1 / 60.0, [unit = u"s^-1"]                  # Rxn 22: .578E+01 min⁻¹
            k23 = 0.474e-1 / (60.0 * 1e3), [unit = u"1/(ppb*s)"]         # Rxn 23: .474E-01 ppm⁻¹min⁻¹
            k24 = 0.178e4 / (60.0 * 1e3), [unit = u"1/(ppb*s)"]          # Rxn 24: .178E+04 ppm⁻¹min⁻¹
            jN2O5 = 0.312e1 / 60.0, [unit = u"s^-1"]                     # Rxn 25: .312E+01 min⁻¹
        end

        @species begin
            NO2(t) = 4e-4, [unit = u"ppb"]
            NO(t) = 4e-4, [unit = u"ppb"]
            O3P(t) = 0, [unit = u"ppb", description = "ground-state atomic oxygen"]
            O3(t) = 40.0, [unit = u"ppb"]
            HO2(t) = 4e-6, [unit = u"ppb"]
            OH(t) = 4e-6, [unit = u"ppb"]
            CH2O(t) = 4e-6, [unit = u"ppb"]
            CO(t) = 100, [unit = u"ppb"]
            ALD(t) = 1e-11, [unit = u"ppb", description = "lumped non-formaldehyde aldehydes"]
            MEO2(t) = 4e-6, [unit = u"ppb", description = "CH3O2 (methyl peroxy radical)"]
            C2O3(t) = 0, [unit = u"ppb", description = "CH3C(O)O2 (acetyl peroxy radical)"]
            CO2(t) = 3.55e5, [unit = u"ppb"]
            PAN(t) = 1e-11, [unit = u"ppb", description = "CH3C(O)O2NO2 (peroxy nitrate)"]
            CH3O(t) = 0, [unit = u"ppb", description = "CH3O (methoxy radical)"]
            HNO3(t) = 4e-6, [unit = u"ppb"]
            O1D(t) = 0, [unit = u"ppb", description = "excited-state atomic oxygen"]
            SO2(t) = 1e-11, [unit = u"ppb"]
            SO4(t) = 1e-11, [unit = u"ppb", description = "SO4 (sulfate)"]
            NO3(t) = 4e-6, [unit = u"ppb"]
            N2O5(t) = 4e-6, [unit = u"ppb"]
        end

        # Gas-phase reactions (Appendix, Verwer 1994)
        k2, NO + O3 --> NO2                              # Rxn 2
        k3, HO2 + NO --> NO2 + OH                        # Rxn 3
        k6, CH2O + OH --> HO2 + CO                       # Rxn 6
        k8, ALD + OH --> C2O3                             # Rxn 8
        k9, C2O3 + NO --> NO2 + MEO2 + CO2               # Rxn 9
        k10, C2O3 + NO2 --> PAN                           # Rxn 10
        k12, MEO2 + NO --> CH3O + NO2                     # Rxn 12
        k13, CH3O --> CH2O + HO2                          # Rxn 13
        k14, NO2 + OH --> HNO3                            # Rxn 14
        k15, O3P --> O3                                   # Rxn 15
        k18, O1D --> 2OH                                  # Rxn 18
        k19, O1D --> O3P                                  # Rxn 19
        k20, SO2 + OH --> SO4 + HO2                       # Rxn 20
        k23, NO2 + O3 --> NO3                             # Rxn 23
        k24, NO3 + NO2 --> N2O5                           # Rxn 24

        # Photolysis reactions
        jNO2_O3P, NO2 --> NO + O3P                       # Rxn 1
        jH2COa, CH2O --> 2HO2 + CO                       # Rxn 4
        jH2COb, CH2O --> CO                               # Rxn 5
        jALD, ALD --> MEO2 + HO2 + CO                    # Rxn 7
        jPAN, PAN --> C2O3 + NO2                          # Rxn 11
        jO3_O1D, O3 --> O1D                               # Rxn 16
        jO3_O3P, O3 --> O3P                               # Rxn 17
        jNO3_NO, NO3 --> NO                               # Rxn 21
        jNO3_NO2, NO3 --> NO2 + O3P                       # Rxn 22
        jN2O5, N2O5 --> NO3 + NO2                         # Rxn 25

    end
    rxns = rx_sys
    if rxn_sys
        return rxns
    end
    # We set `combinatoric_ratelaws=false` because we are modeling macroscopic rather than microscopic behavior.
    # See [here](https://docs.juliahub.com/ModelingToolkit/Qmdqu/3.14.0/systems/ReactionSystem/#ModelingToolkit.oderatelaw)
    # and [here](https://github.com/SciML/Catalyst.jl/issues/311).
    convert(
        Catalyst.ReactionRateSystem,
        complete(rxns);
        combinatoric_ratelaws = false,
        name = name,
        metadata = Dict(CoupleType => PolluCoupler)
    )
end
