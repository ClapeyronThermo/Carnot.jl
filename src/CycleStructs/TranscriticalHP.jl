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
    ΔT_sh_min::Real
    ΔT_sc_min::Real
end


function HeatPumpTranscritical(;fluid::EoSModel,z,T_evap_in,T_evap_out,T_cond_in,T_cond_out,η_comp,pp_evap,pp_cond,ΔT_sh_min,ΔT_sc_min)
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
    @assert T_cond_out > Tcrit + pp_cond "Condenser outlet temperature must be more than critical temperature ($Tcrit) plus pinch point ($pp_cond) for the heat pump to be transcritical"


    type_promoted = promote_type(eltype(z), typeof(T_evap_in), typeof(T_evap_out), typeof(T_cond_in), typeof(T_cond_out), typeof(η_comp), typeof(pp_evap), typeof(pp_cond), typeof(ΔT_sh_min), typeof(ΔT_sc_min))
    z_T = convert(Vector{type_promoted}, z)
    T_evap_in_T = convert(type_promoted, T_evap_in)
    T_evap_out_T = convert(type_promoted, T_evap_out)
    T_cond_in_T = convert(type_promoted, T_cond_in)
    T_cond_out_T = convert(type_promoted, T_cond_out)
    η_comp_T = convert(type_promoted, η_comp)
    pp_evap_T = convert(type_promoted, pp_evap)
    pp_cond_T = convert(type_promoted, pp_cond)
    ΔT_sh_min_T = convert(type_promoted, ΔT_sh_min)
    ΔT_sc_min_T = convert(type_promoted, ΔT_sc_min)
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
    ΔT_sh_min_T,    # T
    ΔT_sc_min_T     # T
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
    ΔT_sh = x[3]
    ΔT_sc = x[4]
 

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
    cop = (h_cond_out - h_comp_out) / (h_comp_out - h_evap_out)

    return SVector(ΔTpp_evap, ΔTpp_cond), cop
end


function COP(prob::HeatPumpTranscritical,x::AbstractVector,param::TranscriticalParamters)
    residue, cop = F(prob,x,N = param.N)
    if residue[1] < 0
        return 0.0
    end
    if residue[2] < 0
        return 0.0
    end
    return cop
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


function generate_box(prob::HeatPumpTranscritical,param::TranscriticalParamters)
    Tcrit, pcrit, _ = crit_pure(prob.fluid)
    ΔT_sh_max = prob.T_evap_in - prob.T_evap_out
    ΔT_sc_max = Tcrit - prob.T_cond_in

    p_evap_min = saturation_pressure(prob.fluid, prob.T_evap_out - prob.pp_evap - ΔT_sh_max)[1]./101325
    p_evap_max = saturation_pressure(prob.fluid, prob.T_evap_in - prob.pp_evap - prob.ΔT_sh_min)[1]./101325

    lb = [p_evap_min, 1.0, prob.ΔT_sh_min, prob.ΔT_sc_min]
    ub = [p_evap_max, param.p_crit_max_ratio, ΔT_sh_max, ΔT_sc_max]

    return lb, ub
end

export F, COP, HeatPumpTranscritical, generate_box