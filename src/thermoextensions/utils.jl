
function isentropic_ideal_gas_out_temp(p_in,T_in,p_out,Cv = 1.5*Clapeyron.Rgas())
    R = Clapeyron.Rgas()
    #=
    Cv*(log(Tout/Tin)) = R*(Vout/Vin) #V = RT/p
    Cv*(log(Tout/Tin)) = R*(Tout/Tin) - R*log(pout/pin)
    (Cv - R)*log(Tout/Tin) = -R*log(pout/pin)
    =#
    dlogT = (-R/(Cv - R))*log(p_out/p_in)
    log_T_out = dlogT + log(T_in)
    return exp(log_T_out)
end


"""
Returns the outlet enthalpy in vapour phase for constant isentropic efficiency. Use in compression of gas. 
"""
function isentropic_compressor(p_in::T1, p_out::T2, η_isen::T3, h_in::T4, z::AbstractArray{TZ}, fluid::EoSModel, crit = crit_mix(fluid,z),Tdew = nothing) where {T1<:Real, T2<:Real, T3<:Real, T4<:Real, TZ<:Real}
    TT = promote_type(T1, T2, T3, T4, TZ)
    T_in = Clapeyron.Tproperty(fluid, p_in, h_in, z,enthalpy, phase = :vapour)::TT
    s_isen = Clapeyron.entropy(fluid, p_in, T_in, z, phase = :vapour)::TT
    # T_in = Tproperty(fluid,p_in,h_in,z,enthalpy,phase = :vapour)
    
    if isnothing(Tdew) #|| (Tdew isa Number && isnan(Tdew))
        T_dew = dew_temperature(fluid,p_out,z)[1]::TT
    else
        T_dew = TT(Tdew)
    end

    T0_isen_out = isnan(T_dew) ? isentropic_ideal_gas_out_temp(p_in,T_in,p_out) : T_dew
    
    #T_isen_out = Clapeyron.PS.temperature(fluid,p_out,s_isen,z, phase = :vapour,T0 = T0_isen_out)::TT
    T_isen_out = Clapeyron.Tproperty(fluid,p_out,s_isen,z,entropy, phase = :vapour,T0 = T0_isen_out)::TT
    h_isen = enthalpy(fluid,p_out,T_isen_out,z, phase = :vapour)::TT

    # h_isen = Clapeyron.PS.enthalpy(fluid, p_out, s_isen, z)

    ha =  (h_in + ((h_isen - h_in)/η_isen))::TT
    #  T_out = Clapeyron.PH.temperature(fluid,p_out,ha,z)
    Tcrit,pcrit,_ = crit
    if p_out < pcrit && pcrit > 0 && !isnan(T_dew)
        h_dew = enthalpy(fluid,p_out,T_dew,z,phase=:vapour)::TT
        if ha < h_dew
        # @warn "Fixing outlet of compressor at saturation temperature"
            return h_dew
        end
    end

    return ha
end

"""
Returns the outlet enthalpy in liquid phase for constant isentropic efficiency. Use in pump of ORC. 
"""
function isentropic_pump(p_in::T1, p_out::T2, η_isen::T3, h_in::T4, z::AbstractArray{TZ}, fluid::EoSModel) where {T1<:Real, T2<:Real, T3<:Real, T4<:Real, TZ<:Real}
    TT = promote_type(T1, T2, T3, T4, TZ)
    s_isen = Clapeyron.PH.entropy(fluid, p_in, h_in, z,phase = :liquid) ::TT
    T_isen_out = Tproperty(fluid,p_out,s_isen,z,entropy,phase = :liquid)::TT
    h_isen = enthalpy(fluid,p_out,T_isen_out,z,phase = :liquid)::TT
    return (h_in + (h_isen - h_in) / η_isen) ::TT
end


"""
Returns the outlet enthalpy in vapour phase for constant isentropic efficiency for turbine application. 
"""
function isentropic_expander(p_in::T1,p_out::T2,η_isen::T3,h_in::T4,z::AbstractVector{TZ},fluid::EoSModel) where {T1<:Real, T2<:Real, T3<:Real, T4<:Real, TZ<:Real}
    TT = promote_type(T1, T2, T3, T4, TZ)
    s_isen = Clapeyron.PH.entropy(fluid, p_in, h_in,z,phase = :vapour) ::TT
    dt = Clapeyron.dew_temperature(fluid, p_out, z)[1]::TT
    T_isen_out = Tproperty(fluid,p_out,s_isen,z,entropy,phase = :vapour,T0 = dt)::TT
    h_isen = enthalpy(fluid,p_out,T_isen_out,z,phase = :vapour)::TT
    h_out = (h_in - (h_in - h_isen) * η_isen) ::TT
    # force outlet to be gaseous
    h_gas = Clapeyron.enthalpy(fluid, p_out, dt, z,phase = :gas)::TT
    # @show   values(dt), h_gas, h_out, p_in, p_out
    if h_out < h_gas
        # @warn "The outlet enthalpy is below the gas enthalpy. Adjusting to gas phase."
        h_out = h_gas
    end

    return h_out::TT
end


"""
Computes the heat transfer in the recuperator. Assumes fluid entering from the left to be higher in temperature. Based over ϵ - effictiveness of heat exchanger. 
"""
function IHEX_Q(fluid::EoSModel,ϵ::T1,T_in_left::T2,p_in_left::T3,T_in_right::T4,p_in_right::T5,z::AbstractVector{TZ}) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real,T5<:Real, TZ<:Real}
    TT = promote_type(T1, T2, T3, T4, T5, TZ)
    T_in_left <= T_in_right && @warn  "Function assume left fluid is hotter. This didn't happen. Now heat transfer will be negative."
    c_in_left = Clapeyron.isobaric_heat_capacity(fluid,p_in_left,T_in_left,z) ::TT
    c_in_right = Clapeyron.isobaric_heat_capacity(fluid,p_in_right,T_in_right,z) ::TT
    C = min(c_in_left,c_in_right)
    Qmax = C*(T_in_left - T_in_right)
    Q = Qmax*ϵ
    return Q
end

"""
Computes the glide match coefficient for a given thermodynamic problem and solution state.
"""
function glide_match_coeff(prob::ThermoCycleProblem,sol::SolutionState;N::Int = 20)
 
end

function glide_match_val(wf::AbstractVector,sf::AbstractVector)

end



compute_isentropic_exponent(fluid::EoSModel,p,T,z) = Clapeyron.adiabatic_index(fluid,p_ref,T_ref,z)

function ideal_compression_relation(κ,V_ratio,η,x)
    πi = V_ratio^κ
    γ = (κ - 1)/κ
    return η*((x)^(γ) - 1)/(πi^γ - (γ*πi^(-1/κ))*(πi - x) - 1)
end

"""
https://data.dtu.dk/articles/dataset/Numerical_models_for_the_design_and_analysis_of_heat_pumps_with_zeotropic_mixtures/6825443?file=13709117
"""
function off_design_compressor_relation(fluid::EoSModel,z,η_isen_design,built_in_volume_ratio;p_ref = 101325.0, T_ref = 300,check = true)
    if check
        @assert 0.0 <= η_isen_design <= 1 "design point isentropic efficiency should be between (0,1)"
        @assert dew_temperature(fluid,p_ref,z)[1] ≤ T_ref "Reference state should be gas. Change p_ref or T_ref"
    end
    κ = Clapeyron.adiabatic_index(fluid,p_ref,T_ref,z,phase = :v)
    compression_relation(x) = ideal_compression_relation(κ,built_in_volume_ratio,η_isen_design,x)
    return compression_relation
end


"""
Similar to `off_design_compressor_relation` for modelling of expander
"""
function off_design_expander_relation(fluid::EoSModel,z,η_isen_design,built_in_volume_ratio;p_ref = 101325.0, T_ref = 300.0, check = true)
    if check
        @assert 0.0 <= η_isen_design <= 1 "design point isentropic efficiency should be between (0,1)"
        @assert dew_temperature(fluid,p_ref,z)[1] ≤ T_ref "Reference state should be gas. Change p_ref or T_ref"
    end
    κ = Clapeyron.adiabatic_index(fluid,p_ref,T_ref,z,phase = :v)
    expander_relation(x) = ideal_compression_relation(κ,built_in_volume_ratio,η_isen_design,1/x)
    return expander_relation
end

export off_design_compressor_relation, compute_isentropic_exponent
export off_design_expander_relation
