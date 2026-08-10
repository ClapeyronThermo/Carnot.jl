
"""
   ORC{T<:Real} <: ThermoCycleProblem

Defines an Organic Rankine Cycle (ORC) problem with thermodynamic and design 
parameters specified in Kelvin and dimensionless efficiencies.

# Fields
- `fluid::CubicModel`: Equation of State (EoS) model representing the working fluid. For now it has to be Cubic EoS.
- `z::AbstractVector{T}`: Mole fraction composition vector of the working fluid.
- `T_evap_in::T`: Inlet temperature of the evaporator [K].
- `T_evap_out::T`: Outlet temperature of the evaporator [K].
- `ΔT_sh::T`: Degree of superheating at the expander inlet [K].
- `T_cond_in::T`: Inlet temperature of the condenser [K].
- `T_cond_out::T`: Outlet temperature of the condenser [K].
- `ΔT_sc::T`: Degree of subcooling at the pump inlet [K].
- `η_pump::T`: Isentropic efficiency of the pump [-].
- `η_expander::T`: Isentropic efficiency of the expander [-].
- `pp_evap::T`: Minimum temperature difference (pinch point) at the evaporator [K].
- `pp_cond::T`: Minimum temperature difference (pinch point) at the condenser [K].

"""
mutable struct ORC{E<:EoSModel,T<:Real,Z<:AbstractVector{T}} <: ThermoCycleProblem
    fluid::E
    z::Z
    T_evap_in::T
    T_evap_out::T
    ΔT_sh::T
    T_cond_in::T
    T_cond_out::T
    ΔT_sc::T
    η_pump::T
    η_expander::T
    pp_evap::T
    pp_cond::T
end


function ORC(; fluid::EoSModel, z, T_evap_in, T_evap_out, T_cond_in, T_cond_out,
             η_pump, η_expander, pp_evap, pp_cond, ΔT_sh, ΔT_sc)

    #@assert fluid isa CubicModel || fluid isa SingleFluid || fluid isa MultiFluid "The type of EOS provided is not supported as of now."

    @assert length(z) > 0 "Composition vector z must not be empty"
    @assert length(fluid.components) == length(z) "Composition vector z must match number of fluid components"

    # Heat source (secondary fluid) temperature drop
    @assert T_evap_in > T_evap_out "Evaporator secondary fluid must cool down (T_evap_in > T_evap_out)"

    # Heat sink (secondary fluid) temperature rise
    @assert T_cond_out > T_cond_in "Condenser secondary fluid must heat up (T_cond_out > T_cond_in)"

    # ORC thermodynamic requirement
    if T_evap_out <= T_cond_in
        @warn "Evaporator outlet temperature is below or equal to condenser inlet temperature. 
        This may lead to a non-functional cycle. Hence the condesner inlet temperature is set equal to the evaporator outlet temperature."
        T_cond_in = T_evap_out
    end
    @assert T_evap_out >= T_cond_in "Working fluid evaporation temperature must exceed condensation temperature"

    # Efficiency assertions
    @assert 0 < η_pump <= 1 "Pump efficiency must be in (0, 1]"
    @assert 0 < η_expander <= 1 "Expander efficiency must be in (0, 1]"

    @assert pp_evap > 0 "Evaporator pinch point must be positive"
    @assert pp_cond > 0 "Condenser pinch point must be positive"

    @assert ΔT_sh ≥ 0 "Superheating must be non-negative"
    @assert ΔT_sc ≥ 0 "Subcooling must be non-negative"

    # Ensure subcritical evaporation
    # Tcrit, _, _ = crit_mix(fluid, z)
    # @assert T_evap_out < Tcrit - pp_evap "Evaporator outlet must be below critical temperature minus pinch"

    # Promote scalars
    type_promoted = promote_type(
        eltype(z),
        typeof(T_evap_in), typeof(T_evap_out),
        typeof(T_cond_in), typeof(T_cond_out),
        typeof(η_pump), typeof(η_expander),
        typeof(pp_evap), typeof(pp_cond),
        typeof(ΔT_sh), typeof(ΔT_sc)
    )

    z_T = map(zi -> convert(type_promoted, zi), z)
    T_evap_in_T  = convert(type_promoted, T_evap_in)
    T_evap_out_T = convert(type_promoted, T_evap_out)
    T_cond_in_T  = convert(type_promoted, T_cond_in)
    T_cond_out_T = convert(type_promoted, T_cond_out)
    η_pump_T     = convert(type_promoted, η_pump)
    η_expander_T = convert(type_promoted, η_expander)
    pp_evap_T    = convert(type_promoted, pp_evap)
    pp_cond_T    = convert(type_promoted, pp_cond)
    ΔT_sh_T      = convert(type_promoted, ΔT_sh)
    ΔT_sc_T      = convert(type_promoted, ΔT_sc)

    return ORC( 
        fluid,
        z_T,
        T_evap_in_T,
        T_evap_out_T,
        ΔT_sh_T,
        T_cond_in_T,
        T_cond_out_T,
        ΔT_sc_T,
        η_pump_T,
        η_expander_T,
        pp_evap_T,
        pp_cond_T
    )
end



"""
p[1] -> condensor pressure
p[2] -> evaporator pressure
"""
function η(prob::ORC,sol::AbstractVector{T}) where {T<:Real}
    @assert length(sol) == 2 "Pressure vector p must be of length 2"
    p_evap,p_cond = sol .* 101325 # convert to Pa
    if p_evap < p_cond
        return 0.0
    end
    Tsat_cond = Clapeyron.bubble_temperature(prob.fluid, p_cond, prob.z)[1]
    h1 = Clapeyron.enthalpy(prob.fluid, p_cond, Tsat_cond - prob.ΔT_sc, prob.z)
    h2 = isentropic_pump(p_cond, p_evap, prob.η_pump, h1, prob.z, prob.fluid)

    Tsat_evap = Clapeyron.dew_temperature(prob.fluid, p_evap, prob.z)[1]
    h3 = Clapeyron.enthalpy(prob.fluid, p_evap, Tsat_evap + prob.ΔT_sh, prob.z)
    h4 = isentropic_expander(p_evap, p_cond, prob.η_expander, h3, prob.z, prob.fluid)
    return ((h4-h3) + (h2-h1))/(h3-h2)
end

"""
`show_parameters(prob::ORC)` prints parameters in REPL
"""
function show_parameters(prob::ORC)
    println("ORC Parameters:")
    println("Fluid: ", prob.fluid)
    println("Composition: ", prob.z)
    println("Evaporator Inlet Temperature: ", prob.T_evap_in, " K")
    println("Evaporator Outlet Temperature: ", prob.T_evap_out, " K")
    println("Condenser Inlet Temperature: ", prob.T_cond_in, " K")
    println("Condenser Outlet Temperature: ", prob.T_cond_out, " K")
    println("Pump Efficiency: ", prob.η_pump)
    println("Expander Efficiency: ", prob.η_expander)
    println("Evaporator Pinch Point: ", prob.pp_evap)
    println("Condenser Pinch Point: ", prob.pp_cond)
    println("Superheating: ", prob.ΔT_sh)
    println("Subcooling: ", prob.ΔT_sc)
end
export ORC, η


function F(prob::ORC, x::AbstractVector{T}; N::Int64) where {T<:Real}
    @assert length(x) == 2 "x must be a vector of length 2"

    if length(prob.fluid.components) == 1
        return F_pure(prob, x)
    end
    
    p_evap, p_cond = x .* 101325 # convert to Pa
    
    flash_res0_cond = Clapeyron.qp_flash_impl(prob.fluid, 0.0, p_cond, prob.z, RRQXFlash(equilibrium = :vle))
    flash_res1_evap = Clapeyron.qp_flash_impl(prob.fluid, 1.0, p_evap, prob.z, RRQXFlash(equilibrium = :vle))

    T_evap_out = Clapeyron.temperature(prob.fluid, flash_res1_evap) + prob.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.fluid, p_evap, T_evap_out, prob.z)

    h_exp_in  = h_evap_out
    h_exp_out = Carnot.isentropic_expander(p_evap, p_cond, prob.η_expander, h_exp_in, prob.z, prob.fluid)

    h_cond_in  = h_exp_out
    T_cond_out = Clapeyron.temperature(prob.fluid, flash_res0_cond) - prob.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.fluid, p_cond, T_cond_out, prob.z)

    # -------------------------
    # Condenser pinch (NaN-safe)
    # -------------------------
    ΔTpp_cond = begin
        Δmin = typemax(T)
        for i in 0:N-1
            α = i / (N-1)

            h = (1-α) * h_cond_out + α * h_cond_in
            T_hx = Clapeyron.PH.temperature(prob.fluid, p_cond, h, prob.z)
            T_sf = (1-α) * prob.T_cond_in + α * prob.T_cond_out

            if !isnan(T_hx)
                Δ = T_hx - T_sf
                if Δ < Δmin
                    Δmin = Δ
                end
            end
        end
        Δmin - prob.pp_cond
    end

    # Pump
    h_pump_in  = h_cond_out
    h_pump_out = Carnot.isentropic_pump(p_cond, p_evap, prob.η_pump, h_pump_in, prob.z, prob.fluid)
    # -------------------------
    # Evaporator pinch (NaN-safe)
    # -------------------------
    ΔTpp_evap = begin
        Δmin = typemax(T)
        for i in 0:N-1
            α = i / (N-1)

            h = (1-α) * h_pump_out + α * h_evap_out
            T_hx = Clapeyron.PH.temperature(prob.fluid, p_evap, h, prob.z)
            T_sf = (1-α) * prob.T_evap_out + α * prob.T_evap_in

            if !isnan(T_hx)
                Δ = T_sf - T_hx
                if Δ < Δmin
                    Δmin = Δ
                end
            end
        end
        Δmin - prob.pp_evap
    end

    return [ΔTpp_evap, ΔTpp_cond]
end


function F_pure(prob::ORC,x::AbstractVector{T}) where T<:Real
    @assert length(prob.fluid.components) == 1 "Pure fluid dispatch only"
    p_evap,p_cond = x .* 101325 # convert to Pa
    T_evap_out = Clapeyron.saturation_temperature(prob.fluid, p_evap)[1] + prob.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.fluid, p_evap, T_evap_out, prob.z)
    h_exp_in = h_evap_out;
    h_exp_out = Carnot.isentropic_expander(p_evap, p_cond, prob.η_expander, h_exp_in, prob.z, prob.fluid)
    h_cond_in = h_exp_out
    T_cond_out = Clapeyron.saturation_temperature(prob.fluid, p_cond)[1] - prob.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.fluid, p_cond, T_cond_out, prob.z)
    T_cond_sat = Clapeyron.saturation_temperature(prob.fluid, p_cond)[1]
    h_cond_sat_liquid = Clapeyron.enthalpy(prob.fluid,p_cond,T_cond_sat,prob.z,phase = :liquid)
    h_cond_sat_vapour = Clapeyron.enthalpy(prob.fluid,p_cond,T_cond_sat,prob.z,phase = :vapour)
    h_cond_array = [h_cond_in,h_cond_sat_vapour,h_cond_sat_liquid,h_cond_out]
    T_cond_array = Clapeyron.PH.temperature.(prob.fluid,p_cond,h_cond_array,prob.z)
    T_cond_sf_f(h) = prob.T_cond_out - (h_cond_in - h)*(prob.T_cond_out - prob.T_cond_in)/(h_cond_in - h_cond_out)
    
    h_pump_in = h_cond_out
    h_pump_out = Carnot.isentropic_pump(p_cond, p_evap, prob.η_pump, h_pump_in, prob.z, prob.fluid)
    h_evap_in = h_pump_out
    T_evap_sat = Clapeyron.saturation_temperature(prob.fluid, p_evap)[1]
    h_evap_sat_liquid = Clapeyron.enthalpy(prob.fluid,p_evap,T_evap_sat,prob.z,phase = :liquid)
    h_evap_sat_vapour = Clapeyron.enthalpy(prob.fluid,p_evap,T_evap_sat,prob.z,phase = :vapour)
    h_evap_array = [h_evap_out,h_evap_sat_vapour,h_evap_sat_liquid,h_pump_out]
    T_evap_array = Clapeyron.PH.temperature.(prob.fluid,p_evap,h_evap_array,prob.z)
    T_evap_sf_f(h) = prob.T_evap_in - (h_evap_out - h)*(prob.T_evap_in - prob.T_evap_out)/(h_evap_out - h_evap_in)

    ΔT_evap = minimum(T_evap_sf_f.(h_evap_array) .- T_evap_array) - prob.pp_evap
    ΔT_cond = minimum(T_cond_array .- T_cond_sf_f.(h_cond_array)) - prob.pp_cond

    return [ΔT_evap,ΔT_cond]
end

"""
Function that gives specific power ratings for ORC by fixing outlet power of expander to equal 1.
"""
function power_ratings(prob::ORC,sol::AbstractVector{T}) where T
    p_evap,p_cond = sol.*101325
    
    T_out_evap = dew_temperature(prob.fluid,p_evap,prob.z)[1] + prob.ΔT_sh
    h_out_evap = enthalpy(prob.fluid,p_evap,T_out_evap,prob.z)
    h_in_exp = h_out_evap;
    h_out_exp = Carnot.isentropic_expander(p_evap, p_cond, prob.η_expander, h_in_exp, prob.z, prob.fluid)
    Δh_exp = h_out_exp - h_in_exp
    if Δh_exp > 0
        @warn "something wrong in the system. Change in enthalpy of the fluid after expansion should be negative"
    end

    h_in_cond = h_out_exp
    T_out_cond = bubble_temperature(prob.fluid,p_cond,prob.z)[1] - prob.ΔT_sc
    h_out_cond = enthalpy(prob.fluid,p_cond,T_out_cond,prob.z)
    ΔQ_cond = h_out_cond - h_in_cond
    if ΔQ_cond > 0
        @warn "something wrong in the system. Change in enthalpy of the fluid after condensation should be negative"
    end
    h_in_pump = h_out_cond
    h_out_pump = Carnot.isentropic_pump(p_cond, p_evap, prob.η_pump, h_in_pump, prob.z, prob.fluid)
    Δh_pump = h_out_pump - h_in_pump
    if Δh_pump < 0 
         @warn "something wrong in the system. Change in enthalpy of the fluid after pump should be positive"
    end
    h_in_evap = h_out_pump
    ΔQ_evap = h_out_evap - h_in_evap
    if ΔQ_evap < 0
        @warn "something wrong in the system. Change in enthalpy of the fluid after evaporator should be positive"
    end

    return [Δh_exp/Δh_exp , Δh_pump/Δh_exp , ΔQ_evap/Δh_exp, ΔQ_cond/Δh_exp]
end

function power_ratings(prob::ORC,sol::SolutionState)
    return power_ratings(prob,sol.x)
end

"""
        ORCEconomizer{T<:Real} <: ThermoCycleProblem

Defines an Organic Rankine Cycle (ORC) configuration with an economiser (regenerative heat exchanger),
extending the base `ORC` problem with a specified effectiveness.

# Fields
- `orc::ORC{T}`: Base ORC system definition containing the thermodynamic parameters.
- `ϵ::T`: Effectiveness of the economiser (regenerator) [-].
"""
mutable struct ORCEconomizer{orc_struct<:ORC,T<:Real} <: ThermoCycleProblem
    orc::orc_struct
    ϵ::T
end

function ORCEconomizer(; orc::ORC, ϵ::Real)
    @assert 0.0 ≤ ϵ < 1.0 "Economizer effectiveness must be in [0, 1)"
    type_promoted = promote_type(eltype(orc.z), typeof(ϵ))
    ϵ_T = convert(type_promoted, ϵ)
    return Carnot.ORCEconomizer(orc, ϵ_T)
end

function F_pure(prob::ORCEconomizer,x::AbstractVector{T}) where T<:Real
    @assert length(prob.orc.fluid.components) == 1 "Pure fluid dispatch only"
    p_evap,p_cond = x .* 101325 # convert to Pa
    T_sat_evap = Clapeyron.saturation_temperature(prob.orc.fluid, p_evap)[1]
    T_sat_cond = Clapeyron.saturation_temperature(prob.orc.fluid, p_cond)[1]
    T_evap_out = Clapeyron.saturation_temperature(prob.orc.fluid, p_evap)[1] + prob.orc.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.orc.fluid, p_evap, T_evap_out, prob.orc.z)
    h_exp_in = h_evap_out;
    h_exp_out = Carnot.isentropic_expander(p_evap, p_cond, prob.orc.η_expander, h_exp_in, prob.orc.z, prob.orc.fluid)
    T_exp_out = Clapeyron.PH.temperature(prob.orc.fluid, p_cond, h_exp_out, prob.orc.z)

    T_cond_out = Clapeyron.saturation_temperature(prob.orc.fluid, p_cond)[1] - prob.orc.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.orc.fluid, p_cond, T_cond_out, prob.orc.z)
    h_pump_in = h_cond_out
    h_pump_out = Carnot.isentropic_pump(p_cond, p_evap, prob.orc.η_pump, h_pump_in, prob.orc.z, prob.orc.fluid)
    T_pump_out = Clapeyron.PH.temperature(prob.orc.fluid, p_evap, h_pump_out, prob.orc.z)

    q_ihex = IHEX_Q(prob.orc.fluid, prob.ϵ,T_exp_out,p_cond,T_pump_out,p_evap,prob.orc.z)

    h_evap_in = h_pump_out + q_ihex
    h_cond_in = h_exp_out - q_ihex

    h_evap_sat_liquid = Clapeyron.enthalpy(prob.orc.fluid,p_evap,T_sat_evap,prob.orc.z,phase = :liquid)
    h_evap_sat_gas = Clapeyron.enthalpy(prob.orc.fluid,p_evap,T_sat_evap,prob.orc.z,phase = :vapour)
    h_evap_array = [h_evap_out,h_evap_sat_gas,h_evap_sat_liquid,h_evap_in]
    T_evap_array = Clapeyron.PH.temperature.(prob.orc.fluid,p_evap,h_evap_array,prob.orc.z)
    T_evap_sf_f(h) = prob.orc.T_evap_in - (h_evap_out - h)*(prob.orc.T_evap_in - prob.orc.T_evap_out)/(h_evap_out - h_evap_in)
    ΔT_evap = minimum(T_evap_sf_f.(h_evap_array) .- T_evap_array) - prob.orc.pp_evap

    h_cond_sat_liquid = Clapeyron.enthalpy(prob.orc.fluid,p_cond,T_sat_cond,prob.orc.z,phase = :liquid)
    h_cond_sat_gas = Clapeyron.enthalpy(prob.orc.fluid,p_cond,T_sat_cond,prob.orc.z,phase = :vapour)
    h_cond_array = [h_cond_in,h_cond_sat_gas,h_cond_sat_liquid,h_cond_out]
    T_cond_array = Clapeyron.PH.temperature.(prob.orc.fluid,p_cond  ,h_cond_array,prob.orc.z)
    T_cond_sf_f(h) = prob.orc.T_cond_out - (h_cond_in - h)*(prob.orc.T_cond_out - prob.orc.T_cond_in)/(h_cond_in - h_cond_out)
    ΔT_cond = minimum(T_cond_array .- T_cond_sf_f.(h_cond_array)) - prob.orc.pp_cond
    return [ΔT_evap,ΔT_cond]
end



function F(prob::ORCEconomizer,x::AbstractVector{T};N::Int64) where {T<:Real}
    @assert length(x) == 2 "x must be a vector of length 2"
    if length(prob.orc.fluid.components) == 1
        return F_pure(prob,x)
    end
    p_evap,p_cond = x .* 101325 # convert to Pa
    flash_res0_cond = Clapeyron.qp_flash_impl(prob.orc.fluid,0.0, p_cond, prob.orc.z, RRQXFlash(equilibrium=:vle))
    flash_res1_evap = Clapeyron.qp_flash_impl(prob.orc.fluid,1.0, p_evap, prob.orc.z, RRQXFlash(equilibrium=:vle))
    T_sat_evap = Clapeyron.temperature(prob.orc.fluid, flash_res1_evap)
    T_sat_cond = Clapeyron.temperature(prob.orc.fluid, flash_res0_cond)
    T_evap_out = T_sat_evap + prob.orc.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.orc.fluid, p_evap, T_evap_out, prob.orc.z)
    h_exp_in = h_evap_out;
    h_exp_out = Carnot.isentropic_expander(p_evap, p_cond, prob.orc.η_expander, h_exp_in, prob.orc.z, prob.orc.fluid)
    T_exp_out = Clapeyron.PH.temperature(prob.orc.fluid, p_cond, h_exp_out, prob.orc.z)
    T_cond_out = T_sat_cond - prob.orc.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.orc.fluid, p_cond, T_cond_out, prob.orc.z)
    h_pump_in = h_cond_out
    h_pump_out = Carnot.isentropic_pump(p_cond, p_evap, prob.orc.η_pump, h_pump_in, prob.orc.z, prob.orc.fluid)
    T_pump_out = Clapeyron.PH.temperature(prob.orc.fluid, p_evap, h_pump_out, prob.orc.z)
    q_ihex = IHEX_Q(prob.orc.fluid, prob.ϵ,T_exp_out,p_cond,T_pump_out,p_evap,prob.orc.z)
    h_evap_in = h_pump_out + q_ihex
    h_cond_in = h_exp_out - q_ihex
    h_evap_array = collect(range(h_evap_out, h_evap_in, length=N))
    T_evap(h) = Clapeyron.PH.temperature(prob.orc.fluid, p_evap, h, prob.orc.z)
    T_evap_array = T_evap.(h_evap_array)
    T_evap_sf(h) = prob.orc.T_evap_in - (h_evap_out - h)*(prob.orc.T_evap_in - prob.orc.T_evap_out)/(h_evap_out - h_evap_in)
    T_evap_sf_array = T_evap_sf.(h_evap_array)
    ΔTpp_evap = minimum(T_evap_sf_array .- T_evap_array) - prob.orc.pp_evap
    h_cond_array = collect(range(h_cond_out, h_cond_in, length=N))
    T_cond(h) = Clapeyron.PH.temperature(prob.orc.fluid, p_cond, h, prob.orc.z)
    T_cond_array = T_cond.(h_cond_array)
    T_cond_sf(h) = prob.orc.T_cond_in + (h - h_cond_out)*(prob.orc.T_cond_out - prob.orc.T_cond_in)/(h_cond_in - h_cond_out)
    T_cond_sf_array = T_cond_sf.(h_cond_array)
    ΔTpp_cond = minimum(T_cond_array .- T_cond_sf_array) - prob.orc.pp_cond
    return [ΔTpp_evap, ΔTpp_cond]#, T_cond_array,
end


function η(prob::Carnot.ORCEconomizer,sol::AbstractVector{T}) where {T<:Real}
    @assert length(sol) == 2 "Pressure vector p must be of length 2"
    p_evap,p_cond = sol .* 101325 # convert to Pa
    if p_evap < p_cond
        return 0.0
    end
    Tsat_cond = Clapeyron.bubble_temperature(prob.orc.fluid, p_cond, prob.orc.z)[1]
    h1 = Clapeyron.enthalpy(prob.orc.fluid, p_cond, Tsat_cond - prob.orc.ΔT_sc, prob.orc.z)
    h2 = Carnot.isentropic_pump(p_cond, p_evap, prob.orc.η_pump, h1, prob.orc.z, prob.orc.fluid)
    T2 = Clapeyron.PH.temperature(prob.orc.fluid, p_evap, h2, prob.orc.z)
    Tsat_evap = Clapeyron.dew_temperature(prob.orc.fluid, p_evap, prob.orc.z)[1]
    h3 = Clapeyron.enthalpy(prob.orc.fluid, p_evap, Tsat_evap + prob.orc.ΔT_sh, prob.orc.z)
    h4 = Carnot.isentropic_expander(p_evap, p_cond, prob.orc.η_expander, h3, prob.orc.z, prob.orc.fluid)
    T4 = Clapeyron.PH.temperature(prob.orc.fluid, p_cond, h4, prob.orc.z)
    q_ihex = Carnot.IHEX_Q(prob.orc.fluid, prob.ϵ,T4,p_cond,T2,p_evap,prob.orc.z)
    h2_new = h2 + q_ihex
    return ((h4-h3) + (h2-h1))/(h3-h2_new)
end

export ORCEconomizer


"""
    η(prob::Carnot.ThermoCycleProblem, sol::SolutionState) -> Float64

Computes the thermal efficiency of a thermodynamic cycle given a problem definition 
and its corresponding solution state.

# Arguments
- `prob::Carnot.ThermoCycleProblem`: The thermodynamic cycle problem 
  containing fluid properties, boundary conditions, and component parameters.
- `sol::SolutionState`: The solution state object containing the converged 
  state variables (`x`), residuals, and convergence information.

# Returns
- `Float64`: The computed cycle efficiency, defined as the ratio of net work 
  output to heat input.

# Notes
This method acts as a wrapper that extracts the solution vector `x` from 
`sol` and calls the lower-level `η(prob, x)` implementation.
"""
function η(prob::Carnot.ThermoCycleProblem,sol::SolutionState)
    return η(prob,sol.x)
end



function _F(prob::ORC, x::AbstractVector{T}) where {T<:Real}
    @assert length(x) == 2 "x must be a vector of length 2"
    @assert length(prob.z) != 1 "This implementation is for mixtures"

    p_evap,p_cond = x .* 101325 # convert to Pa
    T_evap_out = Clapeyron.dew_temperature(prob.fluid, p_evap,prob.z)[1] + prob.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.fluid, p_evap, T_evap_out, prob.z)
    h_exp_in = h_evap_out;
    h_exp_out = Carnot.isentropic_expander(p_evap, p_cond, prob.η_expander, h_exp_in, prob.z, prob.fluid)
    h_cond_in = h_exp_out
    T_cond_in = Clapeyron.PH.temperature(prob.fluid,p_cond,h_cond_in,prob.z)
    T_cond_out = Clapeyron.bubble_temperature(prob.fluid, p_cond,prob.z)[1] - prob.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.fluid, p_cond, T_cond_out, prob.z)
    T_cond_dew = Clapeyron.dew_temperature(prob.fluid, p_cond,prob.z)[1]
    T_cond_bub = Clapeyron.bubble_temperature(prob.fluid, p_cond,prob.z)[1]
    h_cond_sat_liquid = Clapeyron.enthalpy(prob.fluid,p_cond,T_cond_bub,prob.z,phase = :liquid)
    h_cond_sat_vapour = Clapeyron.enthalpy(prob.fluid,p_cond,T_cond_dew,prob.z,phase = :vapour)
    h_cond_array = [h_cond_in,h_cond_sat_vapour,h_cond_sat_liquid,h_cond_out]
    T_cond_array = [T_cond_in,T_cond_dew,T_cond_bub,T_cond_out]
    T_cond_sf_f(h) = prob.T_cond_out - (h_cond_in - h)*(prob.T_cond_out - prob.T_cond_in)/(h_cond_in - h_cond_out)
    
    h_pump_in = h_cond_out
    h_pump_out = Carnot.isentropic_pump(p_cond, p_evap, prob.η_pump, h_pump_in, prob.z, prob.fluid)
    h_evap_in = h_pump_out
    T_evap_in = Clapeyron.PH.temperature(prob.fluid,p_evap,h_evap_in,prob.z)[1]
    T_evap_bub = Clapeyron.bubble_temperature(prob.fluid, p_evap,prob.z)[1]
    T_evap_dew = Clapeyron.dew_temperature(prob.fluid, p_evap,prob.z)[1]
    T_evap_out = T_evap_dew + prob.ΔT_sh
    h_evap_bub = Clapeyron.enthalpy(prob.fluid,p_evap,T_evap_bub,prob.z)
    h_evap_dew = Clapeyron.enthalpy(prob.fluid,p_evap,T_evap_dew,prob.z)
    h_evap_out = Clapeyron.enthalpy(prob.fluid,p_evap,T_evap_out,prob.z)
    h_evap_array = [h_evap_out,h_evap_bub,h_evap_dew,h_evap_in]
    T_evap_array = zeros(length(h_evap_array))
    for i in eachindex(h_evap_array)
        T_evap_array[i] = Clapeyron.PH.temperature(prob.fluid,p_evap,h_evap_array[i],prob.z)
    end
    
    T_evap_sf_f(h) = prob.T_evap_in - (h_evap_out - h)*(prob.T_evap_in - prob.T_evap_out)/(h_evap_out - h_evap_in)

    ΔT_evap = minimum(T_evap_sf_f.(h_evap_array) .- T_evap_array) - prob.pp_evap
    ΔT_cond = minimum(T_cond_array .- T_cond_sf_f.(h_cond_array)) - prob.pp_cond
    return [ΔT_evap,ΔT_cond]
end

function get_states(prob::ORC,sol::SolutionState)
    p_evap,p_cond = sol.x .* 101325
    T_cond_out = Clapeyron.bubble_temperature(prob.fluid, p_cond, prob.z)[1] - prob.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.fluid, p_cond, T_cond_out, prob.z)
    h_cond_out_spec = h_cond_out./Clapeyron.molecular_weight(prob.fluid,prob.z)
    s_cond_out_spec = entropy(prob.fluid,p_cond,T_cond_out,prob.z)./Clapeyron.molecular_weight(prob.fluid,prob.z)

    h_pump_out = Carnot.isentropic_pump(p_cond, p_evap, prob.η_pump, h_cond_out, prob.z, prob.fluid)
    T_pump_out = Clapeyron.PH.temperature(prob.fluid, p_evap, h_pump_out, prob.z)
    h_pump_out_spec = h_pump_out./Clapeyron.molecular_weight(prob.fluid,prob.z)
    s_pump_out_spec = entropy(prob.fluid,p_evap,T_pump_out,prob.z)./Clapeyron.molecular_weight(prob.fluid,prob.z)

    T_evap_out = Clapeyron.dew_temperature(prob.fluid, p_evap, prob.z)[1] + prob.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.fluid, p_evap, T_evap_out, prob.z)
    h_evap_out_spec = h_evap_out./Clapeyron.molecular_weight(prob.fluid,prob.z)
    s_evap_out_spec = entropy(prob.fluid,p_evap,T_evap_out,prob.z)./Clapeyron.molecular_weight(prob.fluid,prob.z)

    h_exp_out = Carnot.isentropic_expander(p_evap, p_cond, prob.η_expander, h_evap_out, prob.z, prob.fluid)
    T_exp_out = Clapeyron.PH.temperature(prob.fluid, p_cond, h_exp_out, prob.z)
    h_exp_out_spec = h_exp_out./Clapeyron.molecular_weight(prob.fluid,prob.z)
    s_exp_out_spec = entropy(prob.fluid,p_cond,T_exp_out,prob.z)./Clapeyron.molecular_weight(prob.fluid,prob.z)


    dict = Dict(
        :p_evap => p_evap,
        :p_cond => p_cond,
        :T_cond_out => T_cond_out,
        :h_cond_out => h_cond_out_spec,
        :s_cond_out => s_cond_out_spec,
        :T_pump_out => T_pump_out,
        :h_pump_out => h_pump_out_spec,
        :s_pump_out => s_pump_out_spec,
        :T_evap_out => T_evap_out,
        :h_evap_out => h_evap_out_spec,
        :s_evap_out => s_evap_out_spec,
        :T_exp_out => T_exp_out,
        :h_exp_out => h_exp_out_spec,
        :s_exp_out => s_exp_out_spec
    )
end

function get_states(prob::ORCEconomizer, sol::SolutionState)
    p_evap,p_cond = sol.x .* 101325
    mw_ = molecular_weight(prob.orc.fluid,prob.orc.z)
    flash_res0_cond = Clapeyron.qp_flash_impl(prob.orc.fluid,0.0, p_cond, prob.orc.z, RRQXFlash(equilibrium=:vle))
    flash_res1_evap = Clapeyron.qp_flash_impl(prob.orc.fluid,1.0, p_evap, prob.orc.z, RRQXFlash(equilibrium=:vle))

    T_evap_sat = Clapeyron.temperature(flash_res1_evap)
    T_evap_out = T_evap_sat + prob.orc.ΔT_sh
    h_evap_out = enthalpy(prob.orc.fluid,p_evap,T_evap_out,prob.orc.z)
    h_evap_out_spec = h_evap_out./mw_
    s_evap_out_spec = entropy(prob.orc.fluid,p_evap,T_evap_out,prob.orc.z)./mw_

    h_exp_out = Carnot.isentropic_expander(p_evap, p_cond, prob.orc.η_expander, h_evap_out, prob.orc.z, prob.orc.fluid)
    h_exp_out_spec = h_exp_out/mw_
    T_exp_out = Clapeyron.PH.temperature(prob.orc.fluid, p_cond, h_exp_out, prob.orc.z)
    s_exp_out_spec = entropy(prob.orc.fluid,p_cond,T_exp_out,prob.orc.z)./mw_

    T_sat_cond = Clapeyron.temperature(prob.orc.fluid, flash_res0_cond)
    T_cond_out = T_sat_cond - prob.orc.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.orc.fluid, p_cond, T_cond_out, prob.orc.z)
    h_cond_out_spec = h_cond_out/mw_
    s_cond_out_spec = Clapeyron.entropy(prob.orc.fluid,p_cond,T_cond_out,prob.orc.z)/mw_

    h_pump_in = h_cond_out
    h_pump_out = Carnot.isentropic_pump(p_cond, p_evap, prob.orc.η_pump, h_pump_in, prob.orc.z, prob.orc.fluid)
    T_pump_out = Clapeyron.PH.temperature(prob.orc.fluid, p_evap, h_pump_out, prob.orc.z)
    h_pump_out_spec = h_pump_out./mw_
    s_pump_out_spec = entropy(prob.fluid,p_evap,T_pump_out,prob.z)./mw_

    q_ihex = IHEX_Q(prob.orc.fluid, prob.ϵ,T_exp_out,p_cond,T_pump_out,p_evap,prob.orc.z)
    
    h_evap_in = h_pump_out + q_ihex
    h_eco_out_pump_end_spec = h_evap_in/mw_
    T_eco_out_pump_end = Clapeyron.PH.temperature(prob.orc.fluid,p_evap,h_evap_in,prob.orc.z)
    s_eco_out_pump_end_spec = entropy(prob.orc.fluid,p_evap,T_eco_out_pump_end,prob.orc.z)./mw_

    h_cond_in = h_exp_out - q_ihex
    h_cond_in_spec = h_cond_in./mw_
    T_cond_in = Clapeyron.PH.temperature(prob.orc.fluid,p_cond,h_cond_in,prob.orc.z)
    s_cond_in_spec = entropy(prob.orc.fluid,p_cond,T_cond_in,prob.orc.z)./mw_

    dict = Dict(
        :p_evap => p_evap,
        :p_cond => p_cond,
        :T_evap_out => T_evap_out,
        :h_evap_out => h_evap_out_spec,
        :s_evap_out => s_evap_out_spec,
        :T_exp_out => T_exp_out,
        :h_exp_out => h_exp_out_spec,
        :s_exp_out => s_exp_out_spec,
        :T_cond_out => T_cond_out,
        :h_cond_out => h_cond_out_spec,
        :s_cond_out => s_cond_out_spec,
        :T_pump_out => T_pump_out,
        :h_pump_out => h_pump_out_spec,
        :s_pump_out => s_pump_out_spec,
        :T_eco_out_pump_end => T_eco_out_pump_end,
        :h_eco_out_pump_end => h_eco_out_pump_end_spec,
        :s_eco_out_pump_end => s_eco_out_pump_end_spec,
        :T_eco_out_exp_end => T_cond_in,
        :h_eco_out_exp_end => h_cond_in_spec,
        :s_eco_out_exp_end => s_cond_in_spec

    )
    return dict
end


mutable struct ORCVarEff{F<:Function,T<:Real} <: ThermoCycleProblem
    fluid::EoSModel
    z::AbstractVector{T}
    T_evap_in::T
    T_evap_out::T
    ΔT_sh::T
    T_cond_in::T
    T_cond_out::T
    ΔT_sc::T
    η_pump::T
    η_expander::F
    pp_evap::T
    pp_cond::T
end


function ORCVarEff(;
    fluid::EoSModel,
    z,
    T_evap_in,
    T_evap_out,
    T_cond_in,
    T_cond_out,
    η_pump,
    η_expander::F,
    pp_evap,
    pp_cond,
    ΔT_sh,
    ΔT_sc
    ) where {F<:Function}

    @assert fluid isa CubicModel || fluid isa SingleFluid || fluid isa MultiFluid "The type of EOS provided is not supported as of now."

    # Composition checks
    @assert length(z) > 0 "Composition vector z must not be empty"
    @assert length(fluid.components) == length(z) "Composition vector z must match number of fluid components"

    # Secondary fluid temperature logic
    @assert T_evap_in > T_evap_out "Evaporator secondary fluid must cool down (T_evap_in > T_evap_out)"
    @assert T_cond_out > T_cond_in "Condenser secondary fluid must heat up (T_cond_out > T_cond_in)"

    # ORC thermodynamic requirement
    @assert T_evap_out > T_cond_in "Evaporation temperature must exceed condensation temperature"

    # Efficiencies
    @assert 0 < η_pump ≤ 1 "Pump efficiency must be in (0, 1]"

    # Very light sanity check on expander efficiency function
    @assert applicable(η_expander, T_evap_out, T_cond_in) ||
            applicable(η_expander, T_evap_out) ||
            applicable(η_expander) "η_expander must be a callable function"

    # Pinch points
    @assert pp_evap > 0 "Evaporator pinch point must be positive"
    @assert pp_cond > 0 "Condenser pinch point must be positive"

    # Superheating / subcooling
    @assert ΔT_sh ≥ 0 "Superheating must be non-negative"
    @assert ΔT_sc ≥ 0 "Subcooling must be non-negative"

    # Optional subcritical evaporation check (kept soft, like your HP logic)
    # Tcrit, _, _ = crit_mix(fluid, z)
    # @assert T_evap_out < Tcrit - pp_evap "Evaporator outlet must be below critical temperature minus pinch"

    # Type promotion
    type_promoted = promote_type(
        eltype(z),
        typeof(T_evap_in),
        typeof(T_evap_out),
        typeof(T_cond_in),
        typeof(T_cond_out),
        typeof(η_pump),
        typeof(pp_evap),
        typeof(pp_cond),
        typeof(ΔT_sh),
        typeof(ΔT_sc)
    )

    z_T          = convert(Vector{type_promoted}, z)
    T_evap_in_T  = convert(type_promoted, T_evap_in)
    T_evap_out_T = convert(type_promoted, T_evap_out)
    T_cond_in_T  = convert(type_promoted, T_cond_in)
    T_cond_out_T = convert(type_promoted, T_cond_out)
    η_pump_T     = convert(type_promoted, η_pump)
    pp_evap_T    = convert(type_promoted, pp_evap)
    pp_cond_T    = convert(type_promoted, pp_cond)
    ΔT_sh_T      = convert(type_promoted, ΔT_sh)
    ΔT_sc_T      = convert(type_promoted, ΔT_sc)

    return ORCVarEff{F,type_promoted}(
        fluid,
        z_T,
        T_evap_in_T,
        T_evap_out_T,
        ΔT_sh_T,
        T_cond_in_T,
        T_cond_out_T,
        ΔT_sc_T,
        η_pump_T,
        η_expander,
        pp_evap_T,
        pp_cond_T
    )
end



function F(prob::ORCVarEff,x::AbstractVector;N::Int)
    η_exp_orc = prob.η_expander(x[1]/x[2])
    orc = ORC(
        fluid = prob.fluid,
        z = prob.z,
        T_evap_in = prob.T_evap_in,
        T_evap_out = prob.T_evap_out,
        ΔT_sh = prob.ΔT_sh,
        T_cond_in = prob.T_cond_in,
        T_cond_out = prob.T_cond_out,
        ΔT_sc = prob.ΔT_sc,
        η_pump = prob.η_pump,
        η_expander = η_exp_orc,
        pp_evap = prob.pp_evap,
        pp_cond = prob.pp_cond
    )
    return F(orc,x;N=N)
end