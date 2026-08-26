````@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: Carnot.jl
  text: A non-linear solver for Heat Pump and Organic Rankine Cycle systems
  image:
    src: logo.png
    alt: Carnot
  tagline: Solve steady-state heat pump and organic Rankine cycle systems, with an optional cycle optimization framework
  actions:
    - theme: brand
      text: Getting started
      link: /examples
    - theme: alt
      text: View on GitHub
      link: https://github.com/ClapeyronThermo/Carnot.jl

features:
  - icon: 🔄
    title: Cycle Systems
    details: Model Organic Rankine Cycle and Heat Pump systems, with or without internal heat exchangers
    link: /examples

  - icon: 🧮
    title: Non-linear Solver
    details: Newton-Raphson with box bounds solves for evaporator and condenser pressures given pinch-point temperatures
    link: /examples

  - icon: 🎯
    title: Cycle Optimization
    details: Find optimal superheating and subcooling temperatures using Metaheuristics.jl
    link: /optimization
---
````

Documentation of Carnot.jl.

The goal of this package is to provide a non-linear solver for Heat Pump and Organic Rankine Cycle systems. It solves for pressures at the evaporator and the condensor for given pinch-point temperatures and provides a framework for plotting the solution.
As of now the package is robust for subcritical cycle parameters. For thermodynamic properties, [Clapeyron.jl](https://github.com/ClapeyronThermo/Clapeyron.jl) is used as backend.

For details of modeling see:
[Carnot batteries for heat and power coupling: Energy, Exergy, Economic and Environmental (4E) analysis - Laterre, Antoine](https://hdl.handle.net/2268/333259), which describes modeling for pure fluids. This package extends the method for mixtures.

The nonlinear solver chosen is Newton-Raphson with box bounds which is inspired by the implementation in [NLboxsolve.jl](https://github.com/RJDennis/NLboxsolve.jl).

This package only supports steady-state applications.

There are 4 systems provided:

1. Organic Rankine Cycle : `ORC`
2. Organic Rankine Cycle with internal heat exchanger : `ORCEconomizer`
3. Heat Pump : `HeatPump`
4. Heat Pump with internal heat exchanger : `HeatPumpRecuperator`

The implemented version of these systems consist of the following components:

1. Compressor : Modelling with isentropic efficiency
2. Expander : Modelling with isentropic efficiency
3. Valve : Modeled as isenthalpic process
4. Evaporator : Volumes of equal change in enthalpy.
5. Condenser : Volumes of equal change in enthalpy.
6. Heat Exchangers (no phase change) : using $\epsilon$ as effectiveness of heat exchanger.

## Related packages

````@raw html
<div class="related-pkg-grid">
  <a class="related-pkg-card" href="https://clapeyronthermo.github.io/Clapeyron.jl/" target="_blank" rel="noreferrer">
    <div class="related-pkg-logo-wrap">
      <img class="related-pkg-logo" src="/assets/related_clapeyron_logo.svg" alt="Clapeyron.jl" />
    </div>
    <h3 class="related-pkg-title">Clapeyron.jl</h3>
    <p class="related-pkg-details">Provides every bulk equation of state used within Carnot, and is required alongside Carnot for essentially all use.</p>
  </a>
  <a class="related-pkg-card" href="https://clapeyronthermo.github.io/GCIdentifier.jl/" target="_blank" rel="noreferrer">
    <div class="related-pkg-logo-wrap">
      <img class="related-pkg-logo" src="/assets/related_gcidentifier_logo.png" alt="GCIdentifier.jl" />
    </div>
    <h3 class="related-pkg-title">GCIdentifier.jl</h3>
    <p class="related-pkg-details">Group contribution identification from SMILES, used for building heterosegmented and group-contribution models.</p>
  </a>
  <a class="related-pkg-card" href="https://clapeyronthermo.github.io/Langmuir.jl/" target="_blank" rel="noreferrer">
    <div class="related-pkg-logo-wrap">
      <img class="related-pkg-logo" src="/assets/related_langmuir_logo.png" alt="Langmuir.jl" />
    </div>
    <h3 class="related-pkg-title">Langmuir.jl</h3>
    <p class="related-pkg-details">Single- and multi-component adsorption equilibrium models.</p>
  </a>
</div>

<style>
.related-pkg-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
  gap: 16px;
  margin: 16px 0 32px;
}

.related-pkg-card {
  display: flex;
  flex-direction: column;
  align-items: center;
  text-align: center;
  padding: 24px;
  border: 1px solid var(--vp-c-bg-soft);
  border-radius: 12px;
  background-color: var(--vp-c-bg-soft);
  text-decoration: none !important;
  transition: border-color 0.25s, background-color 0.25s;
}

.related-pkg-card:hover {
  border-color: var(--vp-c-brand-1);
}

.related-pkg-logo-wrap {
  display: flex;
  align-items: center;
  justify-content: center;
  width: 100%;
  height: 96px;
  margin-bottom: 8px;
}

.related-pkg-logo {
  max-height: 100%;
  max-width: 100%;
  width: auto;
  height: auto;
  object-fit: contain;
}

.related-pkg-title {
  margin: 0;
  line-height: 24px;
  font-size: 16px;
  font-weight: 600;
  color: var(--vp-c-text-1);
  border-top: none;
  padding-top: 0;
}

.related-pkg-details {
  flex-grow: 1;
  margin: 8px 0 0;
  line-height: 22px;
  font-size: 14px;
  font-weight: 500;
  color: var(--vp-c-text-2);
}
</style>
````
