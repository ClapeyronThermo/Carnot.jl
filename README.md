# Carnot.jl

[![Build Status](https://github.com/ClapeyronThermo/Carnot.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ClapeyronThermo/Carnot.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://clapeyronthermo.github.io/Carnot.jl/dev/)



<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/src/assets/logo_dark.svg">
    <source media="(prefers-color-scheme: light)" srcset="docs/src/assets/logo.svg">
    <img src="docs/src/assets/logo.svg" alt="Carnot Logo" width="300">
  </picture>
</p>


This package aims to solve Heat Pump and ORC systems for given known temperature glides. For now it is robust for sub-critical cycles. 

The thermodynamic computations use Clapeyron.jl. 

# Installation
For the latest release, first type `]` and then:

```julia
pkg> add Carnot
```

For the developer version type:

```julia
pkg> add https://github.com/ClapeyronThermo/Carnot.jl
```

# Usage
Usage Heat Pump Example :

```julia
julia> using Clapeyron, Carnot

julia> fluid = cPR(["cyclopentane"],idealmodel = ReidIdeal);

julia> η_comp = 0.75; pp_cond = 2; pp_evap = 2;

julia> T_evap_in = 273.15 + 10; T_evap_out = 273.15 + 0; T_cond_in = 273.15 + 50;  T_cond_out = 273.15+60;

julia> ΔT_sc = 3; ΔT_sh = 10;

julia> hp = HeatPump(fluid=fluid,z=[1.0],T_evap_in=T_evap_in,T_evap_out = T_evap_out,T_cond_in = T_cond_in,T_cond_out=T_cond_out,η_comp=η_comp,pp_evap=pp_evap,pp_cond=pp_cond,ΔT_sc = ΔT_sc,ΔT_sh = ΔT_sh);

julia> sol_hp = solve(hp,ThermoCycleParameters(autodiff=false))
SolutionState{Float64, Int64}([0.12829257763094135, 1.4551588056942617], 20, 4, [1.1044676284654997e-10, 1.475086719437968e-10], [0.07660159441435545, 0.07660159441435545], [1.6566058479359296, 1.6566058479359296], false, 2, 2.859273217366616e-7, 1.8427505452964792e-10, :subcritical)
 
julia> COP(hp,sol_hp)
-3.735868783511875
```

To plot do the following;

```julia
julia> using Plots

julia> plot(hp,sol_hp,N = 100)
```

![HP_cyclopentane](docs/src/Images/hp_cyclopentane.png)


# Limitation

1. Fluid models are limited to the ones provided by default in Clapeyron.jl. Now restricted to `CubicModel` and `SingleFluid` models. 
2. For now the solver is stable for sub-critical parameters. So if incase the solver does not converge please check if the parameters provided allow the solution to be subcritical. 
3. For mixtures, it is recommended to use parameters sufficently below the critical point as sometimes near crictical zone the computation of dew and bubble points can fail.   
4. If for solving with `autodiff = true`, the first run will have significant compilation time. The subsequent runs will be faster.