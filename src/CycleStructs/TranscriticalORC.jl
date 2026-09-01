"""
    TranscriticalORC{E<:EoSModel,T<:Real} <: ThermoCycleProblem

Thermodynamic problem definition for a transcritical Organic Rankine Cycle (ORC).

# Fields

- `fluid::E`: Equation-of-state model describing the ORC working fluid.
- `T_evap_in::T`: Inlet temperature of the secondary fluid to the evaporator.
- `T_evap_out::T`: Outlet temperature of the secondary fluid from the evaporator.
- `ΔT_sh_min::T`: Minimum degree of superheating required at the evaporator outlet.
- `T_cond_in::T`: Inlet temperature of the secondary fluid to the condenser.
- `T_cond_out::T`: Outlet temperature of the secondary fluid from the condenser.
- `ΔT_sc_min::T`: Minimum degree of subcooling required at the condenser outlet.
- `η_pump::T`: Isentropic efficiency of the pump.
- `η_expander::T`: Isentropic efficiency of the expander.
- `pp_evap::T`: Minimum pinch-point temperature difference in the evaporator.
- `pp_cond::T`: Minimum pinch-point temperature difference in the condenser.
"""
mutable struct TranscriticalORC{E<:EoSModel,T<:Real} <: ThermoCycleProblem
    fluid::E
    T_evap_in::T
    T_evap_out::T
    ΔT_sh_min::T
    T_cond_in::T
    T_cond_out::T
    ΔT_sc_min::T
    η_pump::T
    η_expander::T
    pp_evap::T
    pp_cond::T
end

function TranscriticalORC(;
    fluid::EoSModel,
    T_evap_in,
    T_evap_out,
    T_cond_in,
    T_cond_out,
    η_pump,
    η_expander,
    pp_evap,
    pp_cond,
    ΔT_sh_min,
    ΔT_sc_min
    )
    @assert fluid isa CubicModel || fluid isa SingleFluid || fluid isa MultiFluid "The type of EOS provided is not supported as of now."
    #default assertions
    @assert length(fluid.components) == 1 "Trancritical cycle is modelled only for pure fluid models"

    Tcrit,pcrit,_ = crit_pure(fluid)
    @assert Tcrit < T_evap_in "This is trancritical ORC hence, the inlet secondry fluid temperature has to be more than critical point."
    # Heat source (secondary fluid) temperature drop
    @assert T_evap_in > T_evap_out "Evaporator secondary fluid must cool down (T_evap_in > T_evap_out)"

    # Heat sink (secondary fluid) temperature rise
    @assert T_cond_out > T_cond_in "Condenser secondary fluid must heat up (T_cond_out > T_cond_in)"

    # ORC thermodynamic requirement
    if T_evap_out <= T_cond_in
        @warn """
        Evaporator outlet temperature is below or equal to condenser inlet
        temperature. This may lead to a non-functional cycle. Hence the
        condenser inlet temperature is set equal to the evaporator outlet
        temperature.
        """
        T_cond_in = T_evap_out
    end

    @assert T_evap_out >= T_cond_in "Working fluid evaporation temperature must exceed condensation temperature"

    # Efficiency assertions
    @assert 0 < η_pump <= 1 "Pump efficiency must be in (0, 1]"

    @assert 0 < η_expander <= 1 "Expander efficiency must be in (0, 1]"

    # Pinch-point assertions
    @assert pp_evap > 0 "Evaporator pinch point must be positive"

    @assert pp_cond > 0 "Condenser pinch point must be positive"

    # Minimum temperature-difference assertions
    @assert ΔT_sh_min ≥ 0 "Minimum evaporator superheating temperature must be non-negative"

    @assert ΔT_sc_min ≥ 0 "Minimum condenser subcooling temperature must be non-negative"

    # Promote scalars
    type_promoted = promote_type(
        typeof(T_evap_in),
        typeof(T_evap_out),
        typeof(T_cond_in),
        typeof(T_cond_out),
        typeof(η_pump),
        typeof(η_expander),
        typeof(pp_evap),
        typeof(pp_cond),
        typeof(ΔT_sh_min),
        typeof(ΔT_sc_min)
    )

    T_evap_in_T  = convert(type_promoted, T_evap_in)
    T_evap_out_T = convert(type_promoted, T_evap_out)
    T_cond_in_T  = convert(type_promoted, T_cond_in)
    T_cond_out_T = convert(type_promoted, T_cond_out)

    η_pump_T     = convert(type_promoted, η_pump)
    η_expander_T = convert(type_promoted, η_expander)

    pp_evap_T = convert(type_promoted, pp_evap)
    pp_cond_T = convert(type_promoted, pp_cond)

    ΔT_sh_min_T = convert(type_promoted, ΔT_sh_min)
    ΔT_sc_min_T = convert(type_promoted, ΔT_sc_min)

    return TranscriticalORC(
        fluid,
        T_evap_in_T,
        T_evap_out_T,
        ΔT_sh_min_T,
        T_cond_in_T,
        T_cond_out_T,
        ΔT_sc_min_T,
        η_pump_T,
        η_expander_T,
        pp_evap_T,
        pp_cond_T
    )
end



struct TranscriticalParamters
    N::Int 
    p_crit_max_ratio::Float64
end

"""
    F(prob::TranscriticalORC, x::AbstractVector; N::Int = 20)

Evaluate the evaporator and condenser pinch-point constraints and the thermal
efficiency of a transcritical Organic Rankine Cycle (ORC).

The cycle pressures and the minimum degrees of superheating and subcooling are
specified through the decision vector `x`. The heat exchangers are discretised
using the working-fluid enthalpy as the spatial coordinate.

# Arguments

- `prob::TranscriticalORC`: Transcritical ORC problem definition containing the
  working-fluid model, secondary-fluid temperatures, pinch-point constraints,
  and component efficiencies.
- `x::AbstractVector`: Decision vector containing:
    - `x[1]`: Evaporation pressure as a fraction of the critical pressure,
      `p_evap / p_crit`.
    - `x[2]`: Condensation pressure as a fraction of atmospheric pressure,
      `p_cond / 101325`.
    - `x[3]`: Degree of superheating at the evaporator outlet, `ΔT_sh`.
    - `x[4]`: Degree of subcooling at the condenser outlet, `ΔT_sc`.
- `N::Int`: Number of discretisation points used for the evaporator
  enthalpy-based discretisation.

# Returns

A tuple `(Δevap, Δcond), η_orc`, where:

- `Δevap`: Minimum evaporator pinch-point temperature difference relative to
  `prob.pp_evap`.
- `Δcond`: Minimum condenser pinch-point temperature difference relative to
  `prob.pp_cond`.
- `η_orc`: ORC thermal efficiency, calculated as the net specific work
  produced by the expander and pump divided by the specific enthalpy supplied
  in the evaporator.

Positive values of `Δevap` and `Δcond` indicate that the corresponding
minimum pinch-point temperature difference is satisfied.

# Notes

The evaporator is discretised between the pump outlet and evaporator outlet
enthalpies. The condenser is evaluated at the expander outlet, saturated
vapour, saturated liquid, and condenser outlet states.

The evaporator outlet temperature is defined relative to the critical
temperature of the working fluid, while the condenser outlet temperature is
defined relative to the saturation temperature at the condensation pressure.
"""
function F(prob::TranscriticalORC,x::AbstractVector;N::Int = 20)
    T_crit,p_crit,_ = crit_pure(prob.fluid)
    p_evap = x[1]*p_crit
    p_cond = x[2]*101325;
    ΔT_sh = x[3]
    ΔT_sc = x[4]

    z = [1.0]

    T_cond_out = saturation_temperature(prob.fluid,p_cond)[1] - ΔT_sc
    h_cond_out = enthalpy(prob.fluid,p_cond,T_cond_out,z)
    
    h_pump_out = isentropic_pump(p_cond,p_evap,prob.η_pump,h_cond_out,z,prob.fluid)

    T_evap_out = ΔT_sh + T_crit
    h_evap_out = enthalpy(prob.fluid,p_evap,T_evap_out,z)

    h_evap_array = collect(range(h_pump_out,h_evap_out,N)) 
    T_evap_array = similar(h_evap_array)
    T_evap_sf_array = similar(h_evap_array)
    Δevap_array = similar(h_evap_array)
    for i in eachindex(T_evap_array)
        T_evap_sf_array[i] = prob.T_evap_out + (prob.T_evap_in - prob.T_evap_out)*i/N
        T_evap_array[i] = Clapeyron.PH.temperature(prob.fluid,p_evap,h_evap_array[i],z)
        Δevap_array[i] = T_evap_sf_array[i] -  T_evap_array[i] - prob.pp_evap
    end
    Δevap = minimum(Δevap_array)

    h_expander_out =  isentropic_expander(p_evap,p_cond,prob.η_expander,h_evap_out,z,prob.fluid)
    
    T_cond_sat = Clapeyron.saturation_temperature(prob.fluid, p_cond)[1]
    h_cond_sat_liquid = Clapeyron.enthalpy(prob.fluid,p_cond,T_cond_sat,z,phase = :liquid)
    h_cond_sat_vapour = Clapeyron.enthalpy(prob.fluid,p_cond,T_cond_sat,z,phase = :vapour)
    h_cond_array = [h_expander_out,h_cond_sat_vapour,h_cond_sat_liquid,h_cond_out]
    T_cond_array = Clapeyron.PH.temperature.(prob.fluid,p_cond,h_cond_array,z)
    T_cond_sf_f(h) = prob.T_cond_out - (h_expander_out - h)*(prob.T_cond_out - prob.T_cond_in)/(h_expander_out - h_cond_out)
    T_cond_sf_array = T_cond_sf_f.(h_cond_array)

    Δcond = minimum(T_cond_array .- T_cond_sf_array .- prob.pp_cond)

    Δh_exp = h_expander_out - h_evap_out
    Δh_pump = h_pump_out - h_cond_out
    Δh_evap = h_evap_out - h_pump_out
    η_orc = (Δh_exp - Δh_pump)/Δh_evap

    return [Δevap,Δcond],η_orc
end

function η(prob::TranscriticalORC,x::AbstractVector,param::TranscriticalParamters)
    @assert length(x) == 4 "Not pinch point solver only"
    # if residues not met return 0 
    residue,η_orc = F(prob,x,N = param.N)
    if residue[1] < 0
        return 0.0
    end
    if residue[2] < 0
        return 0.0
    end
    return η_orc
end


function generate_box(prob::TranscriticalORC,param::TranscriticalParamters)
    # return box of p_evap,p_cond,Tsh, Tsc
    T_crit,p_crit,_ = crit_pure(prob.fluid)
    
    
    ΔT_sh_max = prob.T_evap_in - T_crit - prob.pp_evap
    ΔT_sc_max = prob.T_cond_out - prob.T_cond_in - prob.pp_cond

    psat_min_cond = saturation_pressure(prob.fluid,prob.T_cond_in + prob.pp_cond + prob.ΔT_sc_min)[1]./101325
    psat_max_cond = saturation_pressure(prob.fluid,prob.T_cond_out + prob.pp_cond + ΔT_sc_max)[1]./101325

    lb = [1.0, psat_min_cond, prob.ΔT_sh_min, prob.ΔT_sc_min]
    ub = [param.p_crit_max_ratio,psat_max_cond,ΔT_sh_max,ΔT_sc_max]
    return lb,ub
end


export TranscriticalORC, TranscriticalParamters, generate_box, η