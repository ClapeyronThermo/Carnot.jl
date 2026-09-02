mutable struct OptHeatPump{E<:EoSModel,T<:Real,Z<:AbstractVector{T}} <: ThermoCycleProblem
    fluid::E
    z::Z
    T_evap_in::T
    T_evap_out::T
    T_cond_in::T
    T_cond_out::T
    η_comp::T
    pp_evap::T
    pp_cond::T
    crit::NTuple{3,T}
end

struct DirectOptParameters{T<:Real}
    N::Int
    ΔT_sh_min::T
    ΔT_sc_min::T
end

function DirectOptParameters(;N::Int, ΔT_sh_min, ΔT_sc_min)
    @assert N > 1 "N must be greater than 1"
    @assert ΔT_sh_min >= 0 "Minimum superheat temperature must be non-negative"
    @assert ΔT_sc_min >= 0 "Minimum subcool temperature must be non-negative"
    return DirectOptParameters(N, ΔT_sh_min, ΔT_sc_min)
end

function crit_mix!(model::OptHeatPump{E,T}) where {E,T}
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

function OptHeatPump(;fluid::EoSModel,z,T_evap_in,T_evap_out,T_cond_in,T_cond_out,η_comp,pp_evap,pp_cond,check_subcritical = true)
    @assert fluid isa CubicModel || fluid isa SingleFluid || fluid isa MultiFluid "The type of EOS provided is not supported as of now."
    #default assertions
    @assert length(z) > 0 "Composition vector z must not be empty"
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
    if check_subcritical
        T_crit,P_crit,V_crit = crit_mix(fluid,z)
        @assert T_cond_in < T_crit - pp_cond "Condenser inlet temperature must be less than critical temperature ($Tcrit) minus pinch point ($pp_cond) for the heat pump to function properly"
    else
        T_crit = convert(Base.promote_eltype(fluid,z),-1)
        P_crit,V_crit = T_crit,T_crit
    end

    type_promoted = promote_type(eltype(z), typeof(T_evap_in), typeof(T_evap_out), typeof(T_cond_in), typeof(T_cond_out), typeof(η_comp), typeof(pp_evap), typeof(pp_cond))
    z_T = map(zi -> convert(type_promoted, zi), z)
    T_evap_in_T = convert(type_promoted, T_evap_in)
    T_evap_out_T = convert(type_promoted, T_evap_out)
    T_cond_in_T = convert(type_promoted, T_cond_in)
    T_cond_out_T = convert(type_promoted, T_cond_out)
    η_comp_T = convert(type_promoted, η_comp)
    pp_evap_T = convert(type_promoted, pp_evap)
    pp_cond_T = convert(type_promoted, pp_cond)
    crit_T = (convert(type_promoted, T_crit),convert(type_promoted, P_crit),convert(type_promoted, V_crit))
    return OptHeatPump(
    fluid,         # EoSModel
    z_T,             # Z<:AbstractVector{T}
    T_evap_in_T,   # T
    T_evap_out_T,  # T
    T_cond_in_T,   # T
    T_cond_out_T,  # T
    η_comp_T,      # T
    pp_evap_T,     # T
    pp_cond_T,     # T
    crit_T,        # T
)
end

function F_pure(prob::OptHeatPump{E,T,Z},x::AbstractVector{T2}) where {E,T,Z,T2<:Real}
    @assert length(x) == 4 "x must be a vector of length 4"
    TT = promote_type(T, T2)
    p_evap = x[1] .* 101325 # convert to Pa
    p_cond = x[2] .* 101325 # convert to Pa
    ΔT_sh = x[3]
    ΔT_sc = x[4]
    T_sat_evap = saturation_temperature(prob.fluid,p_evap)[1]
    T_sat_cond = saturation_temperature(prob.fluid,p_cond)[1]
    T_evap_out = T_sat_evap + ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.fluid, p_evap, T_evap_out, prob.z)
    h_comp_in = h_evap_out;
    crit = crit_mix!(prob)
    h_comp_out = Carnot.isentropic_compressor(p_evap, p_cond, prob.η_comp, h_comp_in, prob.z, prob.fluid, crit, T_sat_cond)
    T_cond_out = T_sat_cond - ΔT_sc
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

    cop = (h_cond_out - h_cond_in)/(h_comp_out - h_comp_in)
    return [ΔT_cond,ΔT_evap], cop
end

function F(prob::OptHeatPump, x::AbstractVector{T}; N::Int) where {T<:Real}
    @assert length(x) == 4 "x must be a vector of length 4"

    if length(prob.fluid.components) == 1
        return F_pure(prob, x)
    end

    p_evap = x[1] * 101_325
    p_cond = x[2] * 101_325
    ΔT_sh = x[3]
    ΔT_sc = x[4]

    flash_res0_cond = Clapeyron.qp_flash_impl(prob.fluid,0.0, p_cond, prob.z, RRQXFlash(equilibrium=:vle)) 
    flash_res1_evap = Clapeyron.qp_flash_impl(prob.fluid,1.0, p_evap, prob.z, RRQXFlash(equilibrium=:vle))
    
    # evaporator outlet
    T_evap_out = Clapeyron.temperature(prob.fluid, flash_res1_evap) + ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.fluid, p_evap, T_evap_out, prob.z)

    # compressor
    crit = crit_mix!(prob)
    h_comp_out = Carnot.isentropic_compressor(p_evap, p_cond, prob.η_comp,
                                       h_evap_out, prob.z, prob.fluid, crit)

    # condenser outlet
    T_cond_out = Clapeyron.temperature(prob.fluid, flash_res0_cond) - ΔT_sc
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

    cop = (h_cond_out - h_comp_out) / (h_comp_out - h_evap_out)
    return [ΔTpp_evap, ΔTpp_cond], cop  # avoids heap allocations
end

function check_feasibility(prob::OptHeatPump,x::AbstractVector;N::Int)
    ΔTpp, cop = try F(prob,x;N=N)
    catch
        [-10,-10], 0.0
    end
    feasible = all(ΔTpp .>= 0)
    return feasible, ΔTpp, cop

end



function generate_box(prob::OptHeatPump,param::DirectOptParameters)
    @assert param.ΔT_sh_min >= 0 "Minimum superheat temperature must be non-negative"
    @assert param.ΔT_sc_min >= 0 "Minimum subcool temperature must be non-negative"
    @assert param.ΔT_sh_min ≤ prob.T_evap_in - prob.T_evap_out "Minimum superheat temperature must be less than the difference between evaporator inlet and outlet temperatures"
    @assert param.ΔT_sc_min ≤ prob.T_cond_out - prob.T_cond_in "Minimum subcool temperature must be less than the difference between condenser outlet and inlet temperatures"
    Δcond = prob.T_cond_out - prob.T_cond_in
    Δevap = prob.T_evap_in - prob.T_evap_out
    p_min = dew_pressure(prob.fluid,prob.T_evap_out - Δevap - prob.pp_evap, prob.z)[1]./101_325
    p_max = crit_mix!(prob)[2]/101_325
    lb = [p_min, p_min, param.ΔT_sh_min, param.ΔT_sc_min]
    ub = [p_max, p_max, prob.T_evap_in - prob.T_evap_out, prob.T_cond_out - prob.T_cond_in]
    return lb, ub
end

function COP(prob::OptHeatPump,x::AbstractVector)
    feasible, _, cop = check_feasibility(prob,x;N=10)
    if feasible
        return cop
    else
        return 0.0
    end
end

function convert_solution(prob::OptHeatPump,sol::SolutionState)
    @assert length(sol.x) == 4 "Solution vector must be of length 4"
    hp = HeatPump(
        fluid = prob.fluid,
        z = prob.z,
        T_evap_in = prob.T_evap_in,
        T_evap_out = prob.T_evap_out,
        T_cond_in = prob.T_cond_in,
        T_cond_out = prob.T_cond_out,
        η_comp = prob.η_comp,
        pp_evap = prob.pp_evap,
        pp_cond = prob.pp_cond,
        ΔT_sh = sol.x[3],
        ΔT_sc = sol.x[4]
    )
    sol = solve(hp, ThermoCycleParameters(N = 20, autodiff = false, max_iters = 10))

    return hp,sol
end

export generate_box, check_feasibility, COP, convert_solution, OptHeatPump, DirectOptParameters