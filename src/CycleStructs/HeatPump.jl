"""
    HeatPump{T<:Real} <: ThermoCycleProblem

A mutable structure representing a vapour-compression heat pump thermodynamic problem.

# Fields
- `fluid::EoSModel`: The equation of state (EoS) model defining the working fluid thermodynamic properties. For now it has to be Cubic EoS.
- `z::AbstractVector{T}`: The composition vector of the working fluid (for mixtures; typically `[1.0]` for pure fluids).
- `T_evap_in::T`: Inlet temperature to the evaporator [K].
- `T_evap_out::T`: Outlet temperature from the evaporator [K].
- `ΔT_sh::T`: Degree of superheating at the evaporator outlet [K].
- `T_cond_in::T`: Inlet temperature to the condenser [K].
- `T_cond_out::T`: Outlet temperature from the condenser [K].
- `ΔT_sc::T`: Degree of subcooling at the condenser outlet [K].
- `η_comp::T`: Isentropic efficiency of the compressor [-].
- `pp_evap::T`: Pinch point temperature difference for evaporator [K].
- `pp_cond::T`: Pinch point temperature difference for condensor [K].

# Optional arguments
- `check_subcritical` = true: checks if the input conditions (inlet temperature minus pinch point) are actually below the mixture critical point.
"""
mutable struct HeatPump{E<:EoSModel,T<:Real,Z<:AbstractVector{T}} <: ThermoCycleProblem
    fluid::E
    z::Z
    T_evap_in::T
    T_evap_out::T
    ΔT_sh::T
    T_cond_in::T
    T_cond_out::T
    ΔT_sc::T
    η_comp::T
    pp_evap::T
    pp_cond::T
    crit::NTuple{3,T}
end

function crit_mix!(model::HeatPump{E,T}) where {E,T}
    #TODO: this is probably not thread safe. safeguard this function behind a lock
    crit0 = model.crit
    if crit0 == (T(-1),T(-1),T(-1))
        Tc,Pc,Vc = crit_mix(model.fluid,model.z)
        model.crit = (T(Tc),T(Pc),T(Vc))
        return model.crit
    else
        return crit0
    end
end

crit_mix!(model::EoSModel,z) = crit_mix(model,z)

function HeatPump(;fluid::EoSModel,z,T_evap_in,T_evap_out,T_cond_in,T_cond_out,η_comp,pp_evap,pp_cond,ΔT_sh,ΔT_sc,check_subcritical = true)
    @assert fluid isa CubicModel || fluid isa SingleFluid || fluid isa MultiFluid "The type of EOS provided is not supported as of now."
    #default assertions
    @assert length(z) > 0 "Composition vector z must not be empty"
    @assert length(fluid.components) == length(z) "Composition vector z must match the number of components in the fluid model"
    @assert T_evap_in > T_evap_out "Evaporator inlet temperature must be more than outlet temperature for the secondary fluid"
    @assert T_cond_in < T_cond_out "Condenser inlet temperature must be less than outlet temperature for the secondary fluid"
    @assert η_comp > 0 && η_comp <= 1 "Compressor efficiency must be between 0 and 1"
    @assert pp_evap > 0 "Evaporator pinch point must be positive"
    @assert pp_cond > 0 "Condenser pinch point must be positive"
    @assert ΔT_sh >= 0 "Superheating temperature must be non-negative"
    @assert ΔT_sc >= 0 "Subcooling temperature must be non-negative"

    # Thermodynamic assertions
    # For heat-pump the inlet temperature of the condensor should be higher than outlet temperature of the evaporator
    if  T_cond_in < T_evap_out
        @warn  "Condenser inlet temperature must be higher than evaporator outlet temperature for the heat pump to function properly. Fixing the evap outlet to condensor inlet"
        T_evap_out = T_cond_in
    end
    # inlet temperature of the condensor should be subcritical - pinch point
    if check_subcritical
        T_crit,P_crit,V_crit = crit_mix(fluid,z)
        @assert T_cond_in < T_crit - pp_cond "Condenser inlet temperature must be less than critical temperature ($Tcrit) minus pinch point ($pp_cond) for the heat pump to function properly"
    else
        T_crit = convert(Base.promote_eltype(fluid,z),-1)
        P_crit,V_crit = T_crit,T_crit
    end

    type_promoted = promote_type(eltype(z), typeof(T_evap_in), typeof(T_evap_out), typeof(T_cond_in), typeof(T_cond_out), typeof(η_comp), typeof(pp_evap), typeof(pp_cond), typeof(ΔT_sh), typeof(ΔT_sc))
    z_T = map(zi -> convert(type_promoted, zi), z)
    ΔT_sh_T = convert(type_promoted, ΔT_sh)
    ΔT_sc_T = convert(type_promoted, ΔT_sc)
    T_evap_in_T = convert(type_promoted, T_evap_in)
    T_evap_out_T = convert(type_promoted, T_evap_out)
    T_cond_in_T = convert(type_promoted, T_cond_in)
    T_cond_out_T = convert(type_promoted, T_cond_out)
    η_comp_T = convert(type_promoted, η_comp)
    pp_evap_T = convert(type_promoted, pp_evap)
    pp_cond_T = convert(type_promoted, pp_cond)
    crit_T = (convert(type_promoted, T_crit),convert(type_promoted, P_crit),convert(type_promoted, V_crit))
    return HeatPump(
    fluid,         # EoSModel
    z_T,             # Z<:AbstractVector{T}
    T_evap_in_T,   # T
    T_evap_out_T,  # T
    ΔT_sh_T,       # T
    T_cond_in_T,   # T
    T_cond_out_T,  # T
    ΔT_sc_T,       # T
    η_comp_T,      # T
    pp_evap_T,     # T
    pp_cond_T,     # T
    crit_T,        # T
)
end


function COP(prob::HeatPump,sol::AbstractVector{T}) where {T<:Real}
    @assert length(sol) == 2 "Pressure vector p must be of length 2"
    p_evap,p_cond = sol .* 101325 # convert to Pa
    T_evap_out = dew_temperature(prob.fluid, p_evap, prob.z)[1] + prob.ΔT_sh 
    h_evap_out = Clapeyron.enthalpy(prob.fluid, p_evap, T_evap_out, prob.z)

    h_comp_in = h_evap_out;
    crit = crit_mix!(prob)
    h_comp_out = isentropic_compressor(p_evap, p_cond, prob.η_comp, h_comp_in, prob.z, prob.fluid, crit)

    T_cond_out = Clapeyron.bubble_temperature(prob.fluid, p_cond, prob.z)[1] - prob.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.fluid, p_cond, T_cond_out, prob.z)

    return  (h_cond_out - h_comp_out)/(h_comp_out - h_comp_in) 
end

"""
    HeatPumpRecuperator{HP<:HeatPump,T<:Real} <: ThermoCycleProblem

A mutable structure representing a heat pump cycle with an internal recuperator (economiser or heat exchanger) between the discharge and suction sides.

# Fields
- `hp::HP`: The base heat pump configuration, containing fluid properties and cycle parameters.
- `ϵ::T`: Effectiveness of the recuperator (dimensionless, typically between 0 and 1).
"""
mutable struct HeatPumpRecuperator{HP<:HeatPump,T<:Real} <:ThermoCycleProblem
    hp::HP
    ϵ::T
end

function HeatPumpRecuperator(;hp::HeatPump,ϵ)
    @assert ϵ <= 1 "Effictiveness has to be less than 1"
    @assert ϵ >=0 "Effictiveness has to be greater than 0"
    return HeatPumpRecuperator(hp,ϵ)
end

crit_mix!(model::HeatPumpRecuperator) = crit_mix!(model.hp)

function show_parameters(prob::HeatPump)
    println("Heat Pump Parameters:")
    println("Fluid: ", prob.fluid)
    println("Composition: ", prob.z)
    println("Evaporator Inlet Temperature: ", prob.T_evap_in, " K")
    println("Evaporator Outlet Temperature: ", prob.T_evap_out, " K")
    println("Superheating Temperature: ", prob.ΔT_sh, " K")
    println("Condenser Inlet Temperature: ", prob.T_cond_in, " K")
    println("Condenser Outlet Temperature: ", prob.T_cond_out, " K")
    println("Subcooling Temperature: ", prob.ΔT_sc, " K")
    println("Compressor Efficiency: ", prob.η_comp)
    println("Evaporator Pinch Point: ", prob.pp_evap)
    println("Condenser Pinch Point: ", prob.pp_cond)
end

export HeatPump, HeatPumpRecuperator, COP, show_parameters

"""
`F(prob::HeatPump, x::AbstractVector{T}; N::Int)` function call to `HeatPump`
"""
function F(prob::HeatPump, x::AbstractVector{T}; N::Int) where {T<:Real}
    @assert length(x) == 2 "x must be a vector of length 2"

    if length(prob.fluid.components) == 1
        return F_pure(prob, x)
    end

    p_evap = x[1] * 101_325
    p_cond = x[2] * 101_325

    flash_res0_cond = Clapeyron.qp_flash_impl(prob.fluid,0.0, p_cond, prob.z, RRQXFlash(equilibrium=:vle)) 
    flash_res1_evap = Clapeyron.qp_flash_impl(prob.fluid,1.0, p_evap, prob.z, RRQXFlash(equilibrium=:vle))
    
    # evaporator outlet
    T_evap_out = Clapeyron.temperature(prob.fluid, flash_res1_evap) + prob.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.fluid, p_evap, T_evap_out, prob.z)

    # compressor
    crit = crit_mix!(prob)
    h_comp_out = isentropic_compressor(p_evap, p_cond, prob.η_comp,
                                       h_evap_out, prob.z, prob.fluid, crit)

    # condenser outlet
    T_cond_out = Clapeyron.temperature(prob.fluid, flash_res0_cond) - prob.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.fluid, p_cond, T_cond_out, prob.z)
    # ----------------------------------
    # Condenser pinch point
    # ----------------------------------
    ΔTpp_cond = begin
        Δmin = typemax(T)
        for i in 0:N-1
            α = i / (N-1)
            h = (1-α) * h_cond_out + α * h_comp_out
            T_hx  = Clapeyron.PH.temperature(prob.fluid, p_cond, h, prob.z)
            T_sf  = (1-α) * prob.T_cond_in + α * prob.T_cond_out
            Δ     = T_hx - T_sf
            if Δ < Δmin
                Δmin = Δ
            end
        end
        Δmin - prob.pp_cond
    end

    # ----------------------------------
    # Evaporator pinch point
    # ----------------------------------
    ΔTpp_evap = begin
        Δmin = typemax(T)
        for i in 0:N-1
            α = i / (N-1)
            h = (1-α) * h_cond_out + α * h_evap_out   # linear enthalpy spacing
            T_hx  = Clapeyron.PH.temperature(prob.fluid, p_evap, h, prob.z)
            T_sf  = (1-α) * prob.T_evap_out + α * prob.T_evap_in
            Δ     = T_sf - T_hx
            if Δ < Δmin
                Δmin = Δ
            end
        end
        Δmin - prob.pp_evap
    end

    return SVector(ΔTpp_evap, ΔTpp_cond)  # avoids heap allocations
end

function F_pure(prob::HeatPump{E,T,Z},x::AbstractVector{T2}) where {E,T,Z,T2<:Real}
    @assert length(x) == 2 "x must be a vector of length 2"
    TT = promote_type(T, T2)
    p_evap,p_cond = x .* 101325 # convert to Pa
    T_sat_evap = saturation_temperature(prob.fluid,p_evap)[1]
    T_sat_cond = saturation_temperature(prob.fluid,p_cond)[1]
    T_evap_out = T_sat_evap + prob.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.fluid, p_evap, T_evap_out, prob.z)
    h_comp_in = h_evap_out;
    crit = crit_mix!(prob)
    h_comp_out = Carnot.isentropic_compressor(p_evap, p_cond, prob.η_comp, h_comp_in, prob.z, prob.fluid, crit, T_sat_cond)
    T_cond_out = T_sat_cond - prob.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.fluid, p_cond, T_cond_out, prob.z)
    h_cond_in = h_comp_out
    h_cond_vapour = Clapeyron.enthalpy(prob.fluid, p_cond, T_sat_cond, prob.z,phase =:vapour)
    h_cond_liquid = Clapeyron.enthalpy(prob.fluid, p_cond, T_sat_cond, prob.z,phase =:liquid)
    h_cond_array = TT[h_cond_in,h_cond_vapour,h_cond_liquid,h_cond_out]
    T_cond_sf_f(h) = prob.T_cond_out - (h_cond_in - h)*(prob.T_cond_out - prob.T_cond_in)/(h_cond_in - h_cond_out)
    Δmin_cond = typemax(TT)
    for h in h_cond_array
        T_hx = Clapeyron.PH.temperature(prob.fluid, p_cond, h, prob.z)::TT
        Δ = (T_hx - T_cond_sf_f(h))::TT
        if Δ < Δmin_cond
            Δmin_cond = Δ
        end
    end
    ΔT_cond = (Δmin_cond - prob.pp_cond)::TT
    h_valve_in = h_cond_out;
    h_valve_out = h_valve_in # isenthalpic expansion

    h_evap_in = h_valve_out
    h_evap_sat_vapour = Clapeyron.enthalpy(prob.fluid, p_evap, T_sat_evap, prob.z,phase =:vapour)
    h_evap_array = reverse(TT[h_evap_in,h_evap_sat_vapour,h_evap_out])
    T_evap_sf_f(h) = prob.T_evap_in - (h_evap_out - h)*(prob.T_evap_in - prob.T_evap_out)/(h_evap_out - h_evap_in)
    Δmin_evap = typemax(TT)
    for h in h_evap_array
        T_hx = Clapeyron.PH.temperature(prob.fluid, p_evap, h, prob.z)::TT
        Δ = (T_evap_sf_f(h) - T_hx)::TT
        if Δ < Δmin_evap
            Δmin_evap = Δ
        end
    end
    ΔT_evap = (Δmin_evap - prob.pp_evap)::TT
    [ΔT_cond,ΔT_evap]
end

"""
Function that gives specific power ratings for HP by fixing outlet power of compressor to equal 1.
"""
function power_ratings(prob::HeatPump,sol::AbstractVector{T}) where T
    p_evap,p_cond = sol.*101325
    
    T_out_evap = dew_temperature(prob.fluid,p_evap,prob.z)[1] + prob.ΔT_sh
    h_out_evap = enthalpy(prob.fluid,p_evap,T_out_evap,prob.z)
    h_in_comp = h_out_evap;
    crit = crit_mix!(prob)
    h_out_comp = Carnot.isentropic_compressor(p_evap, p_cond, prob.η_comp, h_in_comp, prob.z, prob.fluid, crit)
    Δh_comp = h_out_comp - h_in_comp
    if Δh_comp < 0
        @warn "something wrong in the system. Change in enthalpy of the fluid after compressor should be positive"
    end

    h_in_cond = h_out_comp
    T_out_cond = bubble_temperature(prob.fluid,p_cond,prob.z)[1] - prob.ΔT_sc
    h_out_cond = enthalpy(prob.fluid,p_cond,T_out_cond,prob.z)
    ΔQ_cond = h_out_cond - h_in_cond
    if ΔQ_cond > 0
        @warn "something wrong in the system. Change in enthalpy of the fluid after condensation should be negative"
    end
    h_in_valve = h_out_cond
    h_out_valve = h_in_valve
    Δh_valve = h_in_valve - h_out_valve

    h_in_evap = h_out_valve
    ΔQ_evap = h_out_evap - h_in_evap
    if ΔQ_evap < 0
        @warn "something wrong in the system. Change in enthalpy of the fluid after evaporator should be positive"
    end

    return [Δh_comp/Δh_comp , Δh_valve/Δh_comp , ΔQ_evap/Δh_comp, ΔQ_cond/Δh_comp]
end

export power_ratings



function F_pure(prob::HeatPumpRecuperator,x::AbstractVector{T}) where {T<:Real}
    @assert length(x) == 2 "x must be a vector of length 2"
    p_evap,p_cond = x .* 101325 # convert to Pa
    T_sat_evap = saturation_temperature(prob.hp.fluid,p_evap)[1]
    T_sat_cond = saturation_temperature(prob.hp.fluid,p_cond)[1]
    T_evap_out = T_sat_evap + prob.hp.ΔT_sh

    T_cond_out = T_sat_cond - prob.hp.ΔT_sc
    # @show T_evap_out, T_cond_out
    q_ihex = Carnot.IHEX_Q(prob.hp.fluid,prob.ϵ,T_cond_out, p_cond, T_evap_out, p_evap, prob.hp.z)
    h_evap_out = Clapeyron.enthalpy(prob.hp.fluid, p_evap, T_evap_out, prob.hp.z) 
    h_recup_out_comp_end = h_evap_out + q_ihex
    h_comp_in = h_recup_out_comp_end;
    crit = crit_mix!(prob)
    h_comp_out = Carnot.isentropic_compressor(p_evap, p_cond, prob.hp.η_comp, h_comp_in, prob.hp.z, prob.hp.fluid, crit, T_sat_cond)
    h_cond_in = h_comp_out
    h_cond_vapour = Clapeyron.enthalpy(prob.hp.fluid, p_cond, T_sat_cond, prob.hp.z,phase =:vapour)
    h_cond_liquid = Clapeyron.enthalpy(prob.hp.fluid, p_cond, T_sat_cond, prob.hp.z,phase =:liquid)
    h_cond_out = Clapeyron.enthalpy(prob.hp.fluid, p_cond, T_cond_out, prob.hp.z)
    h_cond_array = [h_cond_in,h_cond_vapour,h_cond_liquid,h_cond_out]
    T_cond_array = Clapeyron.PH.temperature.(prob.hp.fluid,p_cond,h_cond_array,prob.hp.z)
    T_cond_sf_f(h) = prob.hp.T_cond_out - (h_cond_in - h)*(prob.hp.T_cond_out - prob.hp.T_cond_in)/(h_cond_in - h_cond_out)
    ΔT_cond = minimum(T_cond_array .- T_cond_sf_f.(h_cond_array)) - prob.hp.pp_cond

    h_recup_out_valve_end = h_cond_out - q_ihex
    h_valve_in = h_recup_out_valve_end;
    h_valve_out = h_valve_in # isenthalpic expansion
    h_evap_in = h_valve_out
    h_evap_sat_vapour = Clapeyron.enthalpy(prob.hp.fluid, p_evap, T_sat_evap, prob.hp.z,phase =:vapour)
    h_evap_array = reverse([h_evap_in,h_evap_sat_vapour,h_evap_out])
    T_evap_array = Clapeyron.PH.temperature.(prob.hp.fluid,p_evap,h_evap_array,prob.hp.z)
    T_evap_sf_f(h) = prob.hp.T_evap_in - (h_evap_out - h)*(prob.hp.T_evap_in - prob.hp.T_evap_out)/(h_evap_out - h_evap_in)
    ΔT_evap = minimum(T_evap_sf_f.(h_evap_array) .- T_evap_array) - prob.hp.pp_evap
    return [ΔT_cond,ΔT_evap]
end

"""
`F(prob::HeatPumpRecuperator, x::AbstractVector{T}; N::Int)` function call to `HeatPump`
"""
function F(prob::HeatPumpRecuperator,x::AbstractVector{T};N::Int64) where {T<:Real}
    @assert length(x) == 2 "x must be a vector of length 2"
    if length(prob.hp.fluid.components) == 1
        return F_pure(prob,x)
    end    
    p_evap,p_cond = x .* 101325 # convert to Pa
    flash_res0_cond = Clapeyron.qp_flash_impl(prob.fluid,0.0, p_cond, prob.z, RRQXFlash(equilibrium=:vle)) 
    flash_res1_evap = Clapeyron.qp_flash_impl(prob.fluid,1.0, p_evap, prob.z, RRQXFlash(equilibrium=:vle))
    
    T_sat_evap = Clapeyron.temperature(prob.fluid, flash_res1_evap)
    T_sat_cond = Clapeyron.temperature(prob.fluid, flash_res0_cond)
    T_evap_out = T_sat_evap + prob.hp.ΔT_sh

    T_cond_out = T_sat_cond - prob.hp.ΔT_sc
    # @show T_evap_out, T_cond_out
    q_ihex = Carnot.IHEX_Q(prob.hp.fluid,prob.ϵ,T_cond_out, p_cond, T_evap_out, p_evap, prob.hp.z)
    h_evap_out = Clapeyron.enthalpy(prob.hp.fluid, p_evap, T_evap_out, prob.hp.z) 
    h_recup_out_comp_end = h_evap_out + q_ihex
    h_comp_in = h_recup_out_comp_end;
    crit = crit_mix!(prob)
    h_comp_out = Carnot.isentropic_compressor(p_evap, p_cond, prob.hp.η_comp, h_comp_in, prob.hp.z, prob.hp.fluid, crit)
    h_cond_in = h_comp_out
    T_cond(h) = Clapeyron.PH.temperature(prob.hp.fluid, p_cond, h, prob.hp.z)
    T_evap(h) = Clapeyron.PH.temperature(prob.hp.fluid, p_evap, h, prob.hp.z)
    h_cond_out = Clapeyron.enthalpy(prob.hp.fluid, p_cond, T_cond_out, prob.hp.z)
    h_cond_array = collect(range(h_cond_out, h_cond_in, length=N))
    T_cond_array = T_cond.(h_cond_array)
    # fix_nan!(T_cond_array)
    T_cond_sf_array = collect(range(prob.hp.T_cond_in, prob.hp.T_cond_out, length=N))
    ΔTpp_cond = minimum(T_cond_array .- T_cond_sf_array) - prob.hp.pp_cond


    T_cond_sf_f(h) = prob.hp.T_cond_out - (h_cond_in - h)*(prob.hp.T_cond_out - prob.hp.T_cond_in)/(h_cond_in - h_cond_out)
    ΔT_cond = minimum(T_cond_array .- T_cond_sf_f.(h_cond_array)) - prob.hp.pp_cond

    h_recup_out_valve_end = h_cond_out - q_ihex
    h_valve_in = h_recup_out_valve_end;
    h_valve_out = h_valve_in # isenthalpic expansion
    h_evap_in = h_valve_out
    h_evap_array = collect(range(h_evap_in, h_evap_out, length=N))
    T_evap_array = T_evap.(h_evap_array)
    # fix_nan!(T_evap_array)
    T_evap_sf_array = collect(range(prob.hp.T_evap_out, prob.hp.T_evap_in, length=N))
    ΔTpp_evap = minimum(T_evap_sf_array .- T_evap_array) - prob.hp.pp_evap
    return [ΔTpp_evap, ΔTpp_cond]
end


function show_parameters(prob::HeatPumpRecuperator)
    println("Heat Pump with Recuperator Parameters:")
    show_parameters(prob.hp)
    println("Recuperator Effectiveness: ", prob.ϵ)
end


function COP(prob::HeatPumpRecuperator,sol::AbstractVector{T}) where {T<:Real}
    @assert length(sol) == 2 "Pressure vector p must be of length 2"
    p_evap,p_cond = sol .* 101325 # convert to Pa
    T_evap_out = dew_temperature(prob.hp.fluid, p_evap, prob.hp.z)[1] + prob.hp.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.hp.fluid, p_evap, T_evap_out, prob.hp.z)
    T_cond_out = bubble_temperature(prob.hp.fluid, p_cond, prob.hp.z)[1] - prob.hp.ΔT_sc
    q_ihex = Carnot.IHEX_Q(prob.hp.fluid,prob.ϵ,T_cond_out, p_cond, T_evap_out, p_evap, prob.hp.z)
    h_comp_in = h_evap_out + q_ihex;
    crit = crit_mix!(prob)
    h_comp_out = isentropic_compressor(p_evap, p_cond, prob.hp.η_comp, h_comp_in, prob.hp.z, prob.hp.fluid, crit)

    h_cond_out = Clapeyron.enthalpy(prob.hp.fluid, p_cond, T_cond_out, prob.hp.z)

    return  (h_cond_out - h_comp_out)/(h_comp_out - h_comp_in) 
end


"""
Function that gives specific power ratings for HP by fixing outlet power of compressor to equal 1.
"""
function power_ratings(prob::HeatPumpRecuperator,sol::AbstractVector{T}) where T
    p_evap,p_cond = sol.*101325

    T_cond_out = bubble_temperature(prob.hp.fluid,p_cond,prob.hp.z)[1] - prob.hp.ΔT_sc
    T_evap_out = dew_temperature(prob.hp.fluid,p_evap,prob.hp.z)[1] + prob.hp.ΔT_sh
    h_out_evap = enthalpy(prob.hp.fluid,p_evap,T_evap_out,prob.hp.z)
    q_ihex = Carnot.IHEX_Q(prob.hp.fluid,prob.ϵ,T_cond_out, p_cond, T_evap_out, p_evap, prob.hp.z)

    if q_ihex < 0
        @warn "Recuperator is cooling the hot stream. T_cond_out < T_evap_out: ($T_cond_out < $T_evap_out)."
    end

    h_recup_out_comp_end = Clapeyron.enthalpy(prob.hp.fluid, p_evap, T_evap_out, prob.hp.z) + q_ihex
    h_in_comp = h_recup_out_comp_end;
    crit = crit_mix!(prob)
    h_out_comp = Carnot.isentropic_compressor(p_evap, p_cond, prob.hp.η_comp, h_in_comp, prob.hp.z, prob.hp.fluid, crit)
    Δh_comp = h_out_comp - h_in_comp
    if Δh_comp < 0
        @warn "something wrong in the system. Change in enthalpy of the fluid after compressor should be positive"
    end
    h_in_cond = h_out_comp
    T_out_cond = bubble_temperature(prob.hp.fluid,p_cond,prob.hp.z)[1] - prob.hp.ΔT_sc
    h_out_cond = enthalpy(prob.hp.fluid,p_cond,T_out_cond,prob.hp.z)
    ΔQ_cond = h_out_cond - h_in_cond
    if ΔQ_cond > 0
        @warn "something wrong in the system. Change in enthalpy of the fluid after condensation should be negative"
    end
    h_in_valve = h_out_cond - q_ihex
    h_out_valve = h_in_valve
    Δh_valve = h_in_valve - h_out_valve
    h_in_evap = h_out_valve
    ΔQ_evap = h_out_evap - h_in_evap
    if ΔQ_evap < 0
        @warn "something wrong in the system. Change in enthalpy of the fluid after evaporator should be positive"
    end
    return [Δh_comp/Δh_comp , Δh_valve/Δh_comp , ΔQ_evap/Δh_comp, ΔQ_cond/Δh_comp, q_ihex/Δh_comp]
end

"""
    COP(prob::Carnot.ThermoCycleProblem, sol::SolutionState) -> Float64

Computes the coefficient of performance (COP) of a thermodynamic cycle given the 
problem definition and its corresponding solution state.

# Arguments
- `prob::Carnot.ThermoCycleProblem`: The thermodynamic cycle problem 
  containing fluid models, boundary conditions, and component parameters.
- `sol::SolutionState`: The solution state object containing the converged 
  state variables (`x`), residuals, and convergence information.
"""
function COP(prob::Carnot.ThermoCycleProblem,sol::SolutionState)
    return COP(prob,sol.x)
end

function F_super(prob::HeatPump,x::AbstractVector,pcrit::Real,Tcrit::Real;N::Int = 30)
    p_evap = x[1] * 101_325
    p_cond = pcrit*x[2]

    T_evap_out = bubble_temperature(prob.fluid,p_evap,prob.z)[1] + prob.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.fluid,p_evap,T_evap_out,prob.z,phase = :vapour)
    
    crit = crit_mix!(prob)
    h_comp_out = Carnot.isentropic_compressor(p_evap,p_cond,prob.η_comp,h_evap_out,prob.z,prob.fluid,crit)
    T_comp_out = Clapeyron.PH.temperature(prob.fluid,p_cond,h_comp_out,prob.z)
    
    T_cond_in = T_comp_out
    h_cond_in = h_comp_out

    T_cond_out = Tcrit - prob.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.fluid,p_cond,T_cond_out,prob.z,phase = :liquid)

    h_cond_array = collect(range(h_cond_in,h_cond_out,length = N))
    T_cond_array = Clapeyron.PH.temperature.(Ref(prob.fluid),Ref(p_cond),h_cond_array,Ref(prob.z))

    T_cond_sf_array = collect(range(prob.T_cond_out,prob.T_cond_in,length = N))

    h_valve_in = h_cond_out
    h_valve_out = h_valve_in
    h_evap_in = h_valve_out
    T_sat_evap = dew_temperature(prob.fluid,p_evap,prob.z)[1]
    h_evap_sat = Clapeyron.enthalpy(prob.fluid,p_evap,T_sat_evap,prob.z,phase = :vapour)
    h_evap_array = [h_evap_out,h_evap_sat,h_evap_in]

    T_evap_array = Clapeyron.PH.temperature.(Ref(prob.fluid),Ref(p_evap),h_evap_array,Ref(prob.z))
    T_evap_sf_array_f(h) = prob.T_evap_in - (h_evap_out - h)*(prob.T_evap_in - prob.T_evap_out)/(h_evap_out - h_evap_in)
    T_evap_sf_array = T_evap_sf_array_f.(h_evap_array)

    return [minimum(T_evap_sf_array .- T_evap_array ) - prob.pp_evap ,minimum(T_cond_array .- T_cond_sf_array) - prob.pp_cond]
end


function power_ratings(prob::ThermoCycleProblem,sol::SolutionState)
    return power_ratings(prob,sol.x)
end



function _F(prob::HeatPump, x::AbstractVector{T}) where {T<:Real}
    @assert length(x) == 2 "x must be a vector of length 2"
    @assert length(prob.z) != 1 "This implementation is for mixtures"

    p_evap = x[1] * 101_325
    p_cond = x[2] * 101_325

    #dew temp, used in isentropic compressor
    T_cond_dew = dew_temperature(prob.fluid,p_cond,prob.z)[1]

    # evaporator outlet
    T_evap_out = dew_temperature(prob.fluid, p_evap, prob.z)[1] + prob.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.fluid, p_evap, T_evap_out, prob.z)

    # compressor
    crit = crit_mix!(prob)
    h_comp_out = isentropic_compressor(p_evap, p_cond, prob.η_comp,
                                       h_evap_out, prob.z, prob.fluid, crit, T_cond_dew)
    h_cond_in = h_comp_out
    T_cond_in = Clapeyron.PH.temperature(prob.fluid,p_cond,h_cond_in,prob.z)
    
    T_cond_bub = bubble_temperature(prob.fluid,p_cond,prob.z)[1]
    h_cond_vapour = Clapeyron.enthalpy(prob.fluid,p_cond,T_cond_dew,prob.z)
    h_cond_liquid = Clapeyron.enthalpy(prob.fluid,p_cond,T_cond_bub,prob.z)
    # condenser outlet
    T_cond_out = Clapeyron.bubble_temperature(prob.fluid, p_cond, prob.z)[1] - prob.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.fluid, p_cond, T_cond_out, prob.z)

    # ----------------------------------
    # Condenser pinch point
    # ----------------------------------
    h_cond_array = [h_cond_in,h_cond_vapour,h_cond_liquid,h_cond_out]
    T_cond_array = [T_cond_in,T_cond_dew,T_cond_bub,T_cond_out] 
    T_cond_sf_f(h) = prob.T_cond_out - (h_cond_in - h)*(prob.T_cond_out - prob.T_cond_in)/(h_cond_in - h_cond_out)
    ΔTpp_cond = minimum(T_cond_array .- T_cond_sf_f.(h_cond_array)) - prob.pp_cond
    # ----------------------------------
    # Evaporator pinch point
    # ----------------------------------
    h_valve_in = h_cond_out;
    h_valve_out = h_valve_in # isenthalpic expansion

    h_evap_in = h_valve_out
    T_evap_dew = dew_temperature(prob.fluid,p_evap,prob.z)[1]
    h_evap_sat_vapour = Clapeyron.enthalpy(prob.fluid, p_evap, T_evap_dew, prob.z,phase =:vapour)
    h_evap_array = reverse([h_evap_in,h_evap_sat_vapour,h_evap_out])
    T_evap_array = similar(h_evap_array)
    for i in eachindex(h_evap_array)
        T_evap_array[i] = Clapeyron.PH.temperature(prob.fluid,p_evap,h_evap_array[i],prob.z)
    end
    
    T_evap_sf_f(h) = prob.T_evap_in - (h_evap_out - h)*(prob.T_evap_in - prob.T_evap_out)/(h_evap_out - h_evap_in)
    ΔTpp_evap = minimum(T_evap_sf_f.(h_evap_array) .- T_evap_array) - prob.pp_evap
    return [ΔTpp_evap, ΔTpp_cond]  # avoids heap allocations
end


function get_states(prob::HeatPump,sol::SolutionState)
    p_evap,p_cond = 101325.0 .* sol.x
    T_evap_out = dew_temperature(prob.fluid,p_evap,prob.z)[1] + prob.ΔT_sh
    h_evap_out = enthalpy(prob.fluid,p_evap,T_evap_out,prob.z)

    h_evap_out_spec = enthalpy(prob.fluid,p_evap,T_evap_out,prob.z)./Clapeyron.molecular_weight(prob.fluid,prob.z)
    s_evap_out_spec = entropy(prob.fluid,p_evap,T_evap_out,prob.z)./Clapeyron.molecular_weight(prob.fluid,prob.z)

    crit = crit_mix!(prob)
    h_comp_out = Carnot.isentropic_compressor(p_evap,p_cond,prob.η_comp,h_evap_out,prob.z,prob.fluid,crit)
    T_comp_out = Clapeyron.PH.temperature(prob.fluid,p_cond,h_comp_out,prob.z)
    s_comp_out_spec = entropy(prob.fluid,p_cond,T_comp_out,prob.z)./Clapeyron.molecular_weight(prob.fluid,prob.z)
    h_comp_out_spec = enthalpy(prob.fluid,p_cond,T_comp_out,prob.z)./Clapeyron.molecular_weight(prob.fluid,prob.z)

    T_cond_out = bubble_temperature(prob.fluid,p_cond,prob.z)[1] - prob.ΔT_sc
    h_cond_out = enthalpy(prob.fluid,p_cond,T_cond_out,prob.z)
    h_cond_out_spec = enthalpy(prob.fluid,p_cond,T_cond_out,prob.z)./Clapeyron.molecular_weight(prob.fluid,prob.z)
    s_cond_out_spec = entropy(prob.fluid,p_cond,T_cond_out,prob.z)./Clapeyron.molecular_weight(prob.fluid,prob.z)

    h_valve_out = h_cond_out
    T_valve_out = Clapeyron.PH.temperature(prob.fluid,p_cond,h_valve_out,prob.z)
    h_valve_out_spec = h_cond_out./Clapeyron.molecular_weight(prob.fluid,prob.z)
    s_valve_out_spec = Clapeyron.PH.entropy(prob.fluid,p_evap,h_valve_out,prob.z)./Clapeyron.molecular_weight(prob.fluid,prob.z)

    return Dict(
        :p_evap => p_evap,
        :T_evap_out => T_evap_out,
        :h_evap_out => h_evap_out_spec,
        :s_evap_out => s_evap_out_spec,
        :T_comp_out => T_comp_out,
        :h_comp_out => h_comp_out_spec,
        :s_comp_out => s_comp_out_spec,
        :p_cond => p_cond,
        :T_cond_out => T_cond_out,
        :h_cond_out => h_cond_out_spec,
        :s_cond_out => s_cond_out_spec,
        :T_valve_out => T_valve_out,
        :h_valve_out => h_valve_out_spec,
        :s_valve_out => s_valve_out_spec
    )
end


function get_states(prob::HeatPumpRecuperator,sol::SolutionState)
    p_evap,p_cond = 101325.0 .* sol.x
    T_evap_out = dew_temperature(prob.hp.fluid,p_evap,prob.hp.z)[1] + prob.hp.ΔT_sh
    h_evap_out = enthalpy(prob.hp.fluid,p_evap,T_evap_out,prob.hp.z)
    h_evap_out_spec = enthalpy(prob.hp.fluid,p_evap,T_evap_out,prob.hp.z)./Clapeyron.molecular_weight(prob.hp.fluid,prob.hp.z)
    s_evap_out_spec = entropy(prob.hp.fluid,p_evap,T_evap_out,prob.hp.z)./Clapeyron.molecular_weight(prob.hp.fluid,prob.hp.z)

    q_ihex = Carnot.IHEX_Q(prob.hp.fluid,prob.ϵ,T_cond_out, p_cond, T_evap_out, p_evap, prob.hp.z)
    h_recup_out_comp_end = h_evap_out + q_ihex
    T_recup_out_comp_end = Clapeyron.PH.temperature(prob.hp.fluid,p_evap,h_recup_out_comp_end,prob.hp.z)
    s_recup_out_comp_end_spec = entropy(prob.hp.fluid,p_evap,T_recup_out_comp_end,prob.hp.z)./Clapeyron.molecular_weight(prob.hp.fluid,prob.hp.z)
    h_recup_out_comp_end_spec = h_recup_out_comp_end./Clapeyron.molecular_weight(prob.hp.fluid,prob.hp.z)
    crit = crit_mix!(prob)
    h_comp_out = Carnot.isentropic_compressor(p_evap,p_cond,prob.hp.η_comp,h_recup_out_comp_end,prob.hp.z,prob.hp.fluid,crit)
    T_comp_out = Clapeyron.PH.temperature(prob.hp.fluid,p_cond,h_comp_out,prob.hp.z)
    s_comp_out_spec = entropy(prob.hp.fluid,p_cond,T_comp_out,prob.hp.z)./Clapeyron.molecular_weight(prob.hp.fluid,prob.hp.z)
    h_comp_out_spec = h_comp_out./Clapeyron.molecular_weight(prob.hp.fluid,prob.hp.z)

    T_cond_out = bubble_temperature(prob.hp.fluid,p_cond,prob.hp.z)[1] - prob.hp.ΔT_sc
    h_cond_out = enthalpy(prob.hp.fluid,p_cond,T_cond_out,prob.hp.z)
    h_cond_out_spec = enthalpy(prob.hp.fluid,p_cond,T_cond_out,prob.hp.z)./Clapeyron.molecular_weight(prob.hp.fluid,prob.hp.z)
    s_cond_out_spec = entropy(prob.hp.fluid,p_cond,T_cond_out,prob.hp.z)./Clapeyron.molecular_weight(prob.hp.fluid,prob.hp.z)

    h_recup_out_valve_end = h_cond_out - q_ihex
    T_recup_out_valve_end = Clapeyron.PH.temperature(prob.hp.fluid,p_cond,h_recup_out_valve_end,prob.hp.z)
    h_recup_out_valve_end_spec = h_recup_out_valve_end./Clapeyron.molecular_weight(prob.hp.fluid,prob.hp.z)
    s_recup_out_valve_end_spec = entropy(prob.hp.fluid,p_cond,T_recup_out_valve_end,prob.hp.z)./Clapeyron.molecular_weight(prob.hp.fluid,prob.hp.z)

    h_valve_out = h_recup_out_valve_end
    T_valve_out = Clapeyron.PH.temperature(prob.hp.fluid,p_evap,h_valve_out,prob.hp.z)
    h_valve_out_spec = h_valve_out./Clapeyron.molecular_weight(prob.hp.fluid,prob.hp.z)
    s_valve_out_spec = Clapeyron.PH.entropy(prob.hp.fluid,p_evap,h_valve_out,prob.hp.z)./Clapeyron.molecular_weight(prob.hp.fluid,prob.hp.z)

    return Dict(
        :p_evap => p_evap,
        :T_evap_out => T_evap_out,
        :h_evap_out => h_evap_out,
        :T_comp_out => T_comp_out,
        :h_comp_out => h_comp_out,
        :p_cond => p_cond,
        :T_cond_out => T_cond_out,
        :h_cond_out => h_cond_out,
        :T_valve_out => T_valve_out,
        :h_valve_out => h_valve_out,
        :h_recup_out_comp_end => h_recup_out_comp_end,
        :h_recup_out_valve_end => h_recup_out_valve_end,
        :T_recup_out_valve_end => T_recup_out_valve_end,
        :h_recup_out_valve_end_spec => h_recup_out_valve_end_spec,
        :s_recup_out_valve_end_spec => s_recup_out_valve_end_spec
    )
    
end


mutable struct HeatPumpVarEff{E,F<:Function,T<:Real} <: ThermoCycleProblem
    η_comp::F
    fluid::E
    z::AbstractVector{T}
    T_evap_in::T
    T_evap_out::T
    ΔT_sh::T
    T_cond_in::T
    T_cond_out::T
    ΔT_sc::T
    pp_evap::T
    pp_cond::T
    crit::NTuple{3,T}
end


function HeatPumpVarEff(;
    η_comp::F,
    fluid::EoSModel,
    z,
    T_evap_in,
    T_evap_out,
    T_cond_in,
    T_cond_out,
    pp_evap,
    pp_cond,
    ΔT_sh,
    ΔT_sc,
    check_subcritical = true,
    ) where {F<:Function}

    @assert fluid isa CubicModel || fluid isa SingleFluid || fluid isa MultiFluid "The type of EOS provided is not supported as of now."

    # default assertions
    @assert length(z) > 0 "Composition vector z must not be empty"
    @assert length(fluid.components) == length(z) "Composition vector z must match the number of components in the fluid model"
    @assert T_evap_in > T_evap_out "Evaporator inlet temperature must be more than outlet temperature for the secondary fluid"
    @assert T_cond_in < T_cond_out "Condenser inlet temperature must be less than outlet temperature for the secondary fluid"
    @assert pp_evap > 0 "Evaporator pinch point must be positive"
    @assert pp_cond > 0 "Condenser pinch point must be positive"
    @assert ΔT_sh ≥ 0 "Superheating temperature must be non-negative"
    @assert ΔT_sc ≥ 0 "Subcooling temperature must be non-negative"

    # very light sanity check on compressor efficiency function
    @assert applicable(η_comp, T_evap_out, T_cond_in) ||
            applicable(η_comp, T_evap_out) ||
            applicable(η_comp) "η_comp must be a callable function (e.g. η(T_evap, T_cond))"

    # Thermodynamic assertions
    if T_cond_in < T_evap_out
        @warn "Condenser inlet temperature must be higher than evaporator outlet temperature for the heat pump to function properly. Fixing the evap outlet to condenser inlet"
        T_evap_out = T_cond_in
    end

    # inlet temperature of the condenser should be subcritical - pinch point
    if check_subcritical
        T_crit, P_crit, V_crit = crit_mix(fluid, z)
        @assert T_cond_in < T_crit - pp_cond "Condenser inlet temperature must be less than critical temperature ($Tcrit) minus pinch point ($pp_cond)"
    else
        T_crit = convert(Base.promote_eltype(fluid,z),-1)
        P_crit,V_crit = T_crit,T_crit
    end
    # type promotion
    type_promoted = promote_type(
        eltype(z),
        typeof(T_evap_in),
        typeof(T_evap_out),
        typeof(T_cond_in),
        typeof(T_cond_out),
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
    pp_evap_T    = convert(type_promoted, pp_evap)
    pp_cond_T    = convert(type_promoted, pp_cond)
    ΔT_sh_T      = convert(type_promoted, ΔT_sh)
    ΔT_sc_T      = convert(type_promoted, ΔT_sc)
    crit_T = (convert(type_promoted, T_crit),convert(type_promoted, P_crit),convert(type_promoted, V_crit))
    return HeatPumpVarEff{F,type_promoted}(
        η_comp,        # Function
        fluid,         # EoSModel
        z_T,           # AbstractVector{T}
        T_evap_in_T,   # T
        T_evap_out_T,  # T
        ΔT_sh_T,       # T
        T_cond_in_T,   # T
        T_cond_out_T,  # T
        ΔT_sc_T,       # T
        pp_evap_T,     # T
        pp_cond_T,     # T
        crit_T,        # NTuple{3,T}
    )
end

function crit_mix!(model::HeatPumpVarEff{F,E,T}) where {F,E,T}
    #TODO: this is probably not thread safe. safeguard this function behind a lock
    crit0 = model.crit
    if crit0 == (T(-1),T(-1),T(-1))
        Tc,Pc,Vc = crit_mix(model.fluid,model.z)
        model.crit = (T(Tc),T(Pc),T(Vc))
        return model.crit
    else
        return crit0
    end
end

function F_pure(prob::HeatPumpVarEff,x::AbstractVector{T}) where {T<:Real}
    @assert length(x) == 2 "x must be a vector of length 2"
    p_evap,p_cond = x .* 101325 # convert to Pa
    η_val = prob.η_comp(p_cond/p_evap)
    T_sat_evap = saturation_temperature(prob.fluid,p_evap)[1]
    T_sat_cond = saturation_temperature(prob.fluid,p_cond)[1]
    T_evap_out = T_sat_evap + prob.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.fluid, p_evap, T_evap_out, prob.z)
    h_comp_in = h_evap_out;

    crit = crit_mix!(prob)
    h_comp_out = Carnot.isentropic_compressor(p_evap, p_cond, η_val, h_comp_in, prob.z, prob.fluid, crit, T_sat_cond)
    T_cond_out = T_sat_cond - prob.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.fluid, p_cond, T_cond_out, prob.z)
    h_cond_in = h_comp_out
    h_cond_vapour = Clapeyron.enthalpy(prob.fluid, p_cond, T_sat_cond, prob.z,phase =:vapour)
    h_cond_liquid = Clapeyron.enthalpy(prob.fluid, p_cond, T_sat_cond, prob.z,phase =:liquid)
    h_cond_array = [h_cond_in,h_cond_vapour,h_cond_liquid,h_cond_out]
    T_cond_array = Clapeyron.PH.temperature.(prob.fluid,p_cond,h_cond_array,prob.z)
    T_cond_sf_f(h) = prob.T_cond_out - (h_cond_in - h)*(prob.T_cond_out - prob.T_cond_in)/(h_cond_in - h_cond_out)
    ΔT_cond = minimum(T_cond_array .- T_cond_sf_f.(h_cond_array)) - prob.pp_cond
    h_valve_in = h_cond_out;
    h_valve_out = h_valve_in # isenthalpic expansion

    h_evap_in = h_valve_out
    h_evap_sat_vapour = Clapeyron.enthalpy(prob.fluid, p_evap, T_sat_evap, prob.z,phase =:vapour)
    h_evap_array = reverse([h_evap_in,h_evap_sat_vapour,h_evap_out])
    T_evap_array = Clapeyron.PH.temperature.(prob.fluid,p_evap,h_evap_array,prob.z)
    T_evap_sf_f(h) = prob.T_evap_in - (h_evap_out - h)*(prob.T_evap_in - prob.T_evap_out)/(h_evap_out - h_evap_in)
    ΔT_evap = minimum(T_evap_sf_f.(h_evap_array) .- T_evap_array) - prob.pp_evap
    [ΔT_cond,ΔT_evap]
end

function F(prob::HeatPumpVarEff,x::AbstractVector{T};N::Int) where {T<:Real}
   if length(prob.fluid.components) == 1
        return F_pure(prob,x)
   end
    p_evap = x[1] * 101_325
    p_cond = x[2] * 101_325
    η_val = prob.η_comp(p_cond/p_evap)
    flash_res0_cond = Clapeyron.qp_flash_impl(prob.fluid,0.0, p_cond, prob.z, RRQXFlash(equilibrium=:vle)) 
    flash_res1_evap = Clapeyron.qp_flash_impl(prob.fluid,1.0, p_evap, prob.z, RRQXFlash(equilibrium=:vle))
    
    # evaporator outlet
    T_evap_out = Clapeyron.temperature(prob.fluid, flash_res1_evap) + prob.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.fluid, p_evap, T_evap_out, prob.z)
    T_dew_cond = Clapeyron.temperature(flash_res0_cond)
    # compressor
    crit = crit_mix!(prob)
    h_comp_out = isentropic_compressor(p_evap, p_cond, η_val,
                                       h_evap_out, prob.z, prob.fluid, crit, T_dew_cond)

    # condenser outlet
    T_cond_out = Clapeyron.temperature(prob.fluid, flash_res0_cond) - prob.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.fluid, p_cond, T_cond_out, prob.z)

    # ----------------------------------
    # Condenser pinch point
    # ----------------------------------
    ΔTpp_cond = begin
        Δmin = typemax(T)
        for i in 0:N-1
            α = i / (N-1)
            h = (1-α) * h_cond_out + α * h_comp_out
            T_hx  = Clapeyron.PH.temperature(prob.fluid, p_cond, h, prob.z)
            T_sf  = (1-α) * prob.T_cond_in + α * prob.T_cond_out
            Δ     = T_hx - T_sf
            if Δ < Δmin
                Δmin = Δ
            end
        end
        Δmin - prob.pp_cond
    end

    # ----------------------------------
    # Evaporator pinch point
    # ----------------------------------
    ΔTpp_evap = begin
        Δmin = typemax(T)
        for i in 0:N-1
            α = i / (N-1)
            h = (1-α) * h_cond_out + α * h_evap_out   # linear enthalpy spacing
            T_hx  = Clapeyron.PH.temperature(prob.fluid, p_evap, h, prob.z)
            T_sf  = (1-α) * prob.T_evap_out + α * prob.T_evap_in
            Δ     = T_sf - T_hx
            if Δ < Δmin
                Δmin = Δ
            end
        end
        Δmin - prob.pp_evap
    end

    return SVector(ΔTpp_evap, ΔTpp_cond)  # avoids heap allocations
end

export HeatPumpVarEff


function HeatPump_from_VarHeatPump(prob::HeatPumpVarEff,sol::SolutionState)
    η_val = prob.η_comp(sol.x[2]/sol.x[1])
    hp = HeatPump(
        η_comp = η_val,
        fluid = prob.fluid,
        z = prob.z,
        T_evap_in = prob.T_evap_in,
        T_evap_out = prob.T_evap_out,
        ΔT_sh = prob.ΔT_sh,
        T_cond_in = prob.T_cond_in,
        T_cond_out = prob.T_cond_out,
        ΔT_sc = prob.ΔT_sc,
        pp_evap = prob.pp_evap,
        pp_cond = prob.pp_cond
    )
    return hp
end

function get_states(prob::HeatPumpVarEff,sol::SolutionState)
    hp = HeatPump_from_VarHeatPump(prob,sol)
    return get_states(hp,sol)    
end

function COP(prob::HeatPumpVarEff,sol::SolutionState)
    hp = HeatPump_from_VarHeatPump(prob,sol)
    return COP(hp,sol)
end

export HeatPump_from_VarHeatPump, get_states

"""
`VHC(prob::HeatPump,sol::SolutionState)`

Returns the volumetric heating capacity [kJ/m³]. 
"""
function VHC(prob::HeatPump,sol::SolutionState)
    states = get_states(prob,sol)
    h_spec_comp = states[:h_comp_out]
    h_cond_out = states[:h_cond_out]
    ρ_suction = mass_density(prob.fluid,states[:p_evap],states[:T_evap_out],prob.z)
    return (h_spec_comp - h_cond_out)*ρ_suction/1e3
end

export VHC



mutable struct HeatPumpTranscritical <:ThermoCycleProblem
    fluid::EoSModel
    z::AbstractVector{<:Real}
    T_evap_in::Real
    T_evap_out::Real
    T_cond_in::Real
    T_cond_out::Real
    η_comp::Real
    pp_evap::Real
    pp_cond::Real
    ΔT_sh::Real
end


function HeatPumpTranscritical(;fluid::EoSModel,z,T_evap_in,T_evap_out,T_cond_in,T_cond_out,η_comp,pp_evap,pp_cond,ΔT_sh)
    @assert fluid isa CubicModel || fluid isa SingleFluid || fluid isa MultiFluid "The type of EOS provided is not supported as of now."
    #default assertions
    @assert length(z) == 1 "Composition vector z must have exactly one element. implementation for pure fluids only"
    @assert length(fluid.components) == length(z) "Composition vector z must match the number of components in the fluid model"
    @assert T_evap_in > T_evap_out "Evaporator inlet temperature must be more than outlet temperature for the secondary fluid"
    @assert T_cond_in < T_cond_out "Condenser inlet temperature must be less than outlet temperature for the secondary fluid"
    @assert η_comp > 0 && η_comp <= 1 "Compressor efficiency must be between 0 and 1"
    @assert pp_evap > 0 "Evaporator pinch point must be positive"
    @assert pp_cond > 0 "Condenser pinch point must be positive"


    # Thermodynamic assertions
    # For heat-pump the inlet temperature of the condensor should be higher than outlet temperature of the evaporator
    if  T_cond_in < T_evap_out
        @warn  "Condenser inlet temperature must be higher than evaporator outlet temperature for the heat pump to function properly. Fixing the evap outlet to condensor inlet"
        T_evap_out = T_cond_in
    end
    # inlet temperature of the condensor should be subcritical - pinch point
    Tcrit,_,_ = crit_mix(fluid,z)
    @assert T_cond_in < Tcrit - pp_cond "Condenser inlet temperature must be less than critical temperature ($Tcrit) minus pinch point ($pp_cond) for the heat pump to function properly"



    type_promoted = promote_type(eltype(z), typeof(T_evap_in), typeof(T_evap_out), typeof(T_cond_in), typeof(T_cond_out), typeof(η_comp), typeof(pp_evap), typeof(pp_cond))
    z_T = convert(Vector{type_promoted}, z)
    T_evap_in_T = convert(type_promoted, T_evap_in)
    T_evap_out_T = convert(type_promoted, T_evap_out)
    T_cond_in_T = convert(type_promoted, T_cond_in)
    T_cond_out_T = convert(type_promoted, T_cond_out)
    η_comp_T = convert(type_promoted, η_comp)
    pp_evap_T = convert(type_promoted, pp_evap)
    pp_cond_T = convert(type_promoted, pp_cond)
    ΔT_sh_T = convert(type_promoted,ΔT_sh)
    return HeatPumpTranscritical(
    fluid,         # EoSModel
    z_T,             # AbstractVector{T}
    T_evap_in_T,   # T
    T_evap_out_T,  # T
    T_cond_in_T,   # T
    T_cond_out_T,  # T
    η_comp_T,      # T
    pp_evap_T,     # T
    pp_cond_T,      # T
    ΔT_sh_T
)
end

export HeatPumpTranscritical

function F(prob::HeatPumpTranscritical,x::AbstractVector{TT};N::Int) where {TT<:Real}
    # x is the vector of unknowns (evaporator and condenser pressures)
    # p_cond is supercritical and p_evap is subcritical
    # T is the vector of temperatures (evaporator and condenser outlet temperatures) also unkown

    T_crit,p_crit,_ = crit_mix(prob.fluid,prob.z)
    p_evap = x[1] * 101325
    p_cond = x[2] * p_crit

    ΔT_sh = prob.ΔT_sh
    ΔT_sc = T_crit - (prob.T_cond_in + prob.pp_cond) #pinch will be met by default. So check for other pinch part from this? 

    T_evap_out = dew_temperature(prob.fluid,p_evap,prob.z)[1] + ΔT_sh
    h_evap_out = enthalpy(prob.fluid,p_evap,T_evap_out,prob.z)
    h_comp_out = Carnot.isentropic_compressor(p_evap, p_cond, prob.η_comp, h_evap_out, prob.z, prob.fluid)
    T_cond_out = T_crit - ΔT_sc
    h_cond_out = enthalpy(prob.fluid,p_cond,T_cond_out,prob.z)

    # ----------------------------------
    # Condenser pinch point
    # ----------------------------------
    ΔTpp_cond = begin
        Δmin = typemax(TT)
        for i in 0:N-1
            α = i / (N-1)
            h = (1-α) * h_cond_out + α * h_comp_out
            T_hx  = Clapeyron.PH.temperature(prob.fluid, p_cond, h, prob.z)
            T_sf  = (1-α) * prob.T_cond_in + α * prob.T_cond_out
            Δ     = T_hx - T_sf
            if Δ < Δmin
                Δmin = Δ
            end
        end
        Δmin - prob.pp_cond
    end

    h_valve_in = h_cond_out;
    h_valve_out = h_valve_in # isenthalpic expansion
    # ----------------------------------
    # Evaporator pinch point
    # ----------------------------------
    ΔTpp_evap = begin
        Δmin = typemax(TT)
        for i in 0:N-1
            α = i / (N-1)
            h = (1-α) * h_cond_out + α * h_evap_out   # linear enthalpy spacing
            T_hx  = Clapeyron.PH.temperature(prob.fluid, p_evap, h, prob.z)
            T_sf  = (1-α) * prob.T_evap_out + α * prob.T_evap_in
            Δ     = T_sf - T_hx
            if Δ < Δmin
                Δmin = Δ
            end
        end
        Δmin - prob.pp_evap
    end

    return SVector(ΔTpp_evap, ΔTpp_cond)
end


function COP(prob::HeatPumpTranscritical,x::AbstractVector,T::AbstractVector)
    T_crit,p_crit,_ = crit_mix(prob.fluid,prob.z)
    p_evap = x[1] * 101325
    p_cond = x[2] * p_crit
    ΔT_sh = T[1]
    ΔT_sc = T[2]
    T_evap_out = dew_temperature(prob.fluid,p_evap,prob.z)[1] + ΔT_sh
    h_evap_out = enthalpy(prob.fluid,p_evap,T_evap_out,prob.z)
    h_comp_out = Carnot.isentropic_compressor(p_evap, p_cond, prob.η_comp, h_evap_out, prob.z, prob.fluid)
    T_cond_out = T_crit - ΔT_sc
    h_cond_out = enthalpy(prob.fluid,p_cond,T_cond_out,prob.z)
    Q_evap = h_cond_out - h_comp_out
    Q_comp = h_comp_out - h_evap_out
    return Q_evap/Q_comp
end

function show(io::IO,prob::HeatPumpTranscritical)
    println(io, "Transcritical Heat Pump:")
    println(io, "  Fluid: $(prob.fluid)")
    println(io, "  Composition: $(prob.z)")
    println(io, "  Evaporator Inlet Temperature: $(prob.T_evap_in)")
    println(io, "  Evaporator Outlet Temperature: $(prob.T_evap_out)")
    println(io, "  Condenser Inlet Temperature: $(prob.T_cond_in)")
    println(io, "  Condenser Outlet Temperature: $(prob.T_cond_out)")
    println(io, "  Compressor Efficiency: $(prob.η_comp)")
    println(io, "  Evaporator Pressure Drop: $(prob.pp_evap)")
    println(io, "  Condenser Pressure Drop: $(prob.pp_cond)")
end
