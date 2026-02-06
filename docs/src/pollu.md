# POLLU Atmospheric Chemistry System

## Overview

The POLLU system is a 20-species, 25-reaction atmospheric chemistry ODE system originally introduced as a benchmark problem for stiff ODE solvers.

**Reference**: Verwer, J. G. (1994). "Gauss-Seidel Iteration for Stiff ODEs from Chemical Kinetics." *SIAM Journal on Scientific Computing*, 15(5), 1243-1250. [CWI Report](https://ir.cwi.nl/pub/2262/2262D.pdf)

```@docs
Pollu
```

## Implementation

The model includes 20 chemical species and 25 reactions covering NO-NO2-O3 photochemistry, formaldehyde photolysis, aldehyde oxidation, PAN chemistry, SO2 oxidation, and nighttime NO3/N2O5 chemistry. Rate constants from the original paper (in min⁻¹ and ppm⁻¹min⁻¹) are converted to SI-compatible units (s⁻¹ and ppb⁻¹s⁻¹).

### State Variables

```@example pollu
using GasChem, ModelingToolkit, DataFrames, Symbolics, DynamicQuantities

model = Pollu()
vars = unknowns(model)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [ModelingToolkit.get_unit(v) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars],
    :Default => [ModelingToolkit.getdefault(v) for v in vars])
```

### Parameters

```@example pollu
pars = parameters(model)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in pars],
    :Units => [ModelingToolkit.get_unit(v) for v in pars],
    :Description => [ModelingToolkit.getdescription(v) for v in pars],
    :Default => [ModelingToolkit.getdefault(v) for v in pars])
```

### Equations

```@example pollu
eqs = equations(model)
```

## Analysis

### Reproducing Verwer (1994) Reference Solutions

The paper provides reference solutions at t = 60 s and t = 3600 s for a specific set of initial conditions. Here we reproduce those results.

```@example pollu
using OrdinaryDiffEqRosenbrock, Plots

sys = mtkcompile(model)
tspan = (0.0, 3600.0)

# Initial conditions from the paper (Verwer 1994, Appendix), converted from ppm to ppb
initial_conditions = [
    sys.NO2 => 0,
    sys.NO => 0.2 * 1e3,
    sys.O3P => 0,
    sys.O3 => 0.04 * 1e3,
    sys.HO2 => 0,
    sys.OH => 0,
    sys.CH2O => 0.1 * 1e3,
    sys.CO => 0.3 * 1e3,
    sys.ALD => 0.01 * 1e3,
    sys.MEO2 => 0,
    sys.C2O3 => 0,
    sys.CO2 => 0,
    sys.PAN => 0,
    sys.CH3O => 0,
    sys.HNO3 => 0,
    sys.O1D => 0,
    sys.SO2 => 0.007 * 1e3,
    sys.SO4 => 0,
    sys.NO3 => 0,
    sys.N2O5 => 0,
]

sol = solve(
    ODEProblem(sys, initial_conditions, tspan),
    Rosenbrock23(),
    saveat = 1.0,
    abstol = 1e-12,
    reltol = 1e-12,
)

p = plot(sol, idxs=[sys.NO2, sys.NO, sys.O3, sys.CH2O, sys.CO],
    xlabel = "Time (s)", ylabel = "Concentration (ppb)",
    title = "POLLU major species evolution",
    legend = :outerright)
p
```

### Radical and Minor Species

```@example pollu
p = plot(sol, idxs=[sys.OH, sys.HO2, sys.MEO2, sys.C2O3, sys.CH3O, sys.O1D, sys.O3P],
    xlabel = "Time (s)", ylabel = "Concentration (ppb)",
    title = "POLLU radical species evolution",
    legend = :outerright)
p
```

### Nitrogen and Sulfur Species

```@example pollu
p = plot(sol, idxs=[sys.HNO3, sys.PAN, sys.NO3, sys.N2O5, sys.SO2, sys.SO4],
    xlabel = "Time (s)", ylabel = "Concentration (ppb)",
    title = "POLLU nitrogen and sulfur species",
    legend = :outerright)
p
```
