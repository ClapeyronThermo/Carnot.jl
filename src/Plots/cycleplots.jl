function plotting_data(fluid::EoSModel,z::AbstractVector;N = 30,p_min = nothing,nanfix = true)
    @assert length(fluid.components) == length(z) "Components and composition vector length mismatch"
    if isnothing(p_min)
        p_min = 101325*0.4
    end
    _,p_crit,_ = crit_mix(fluid,z)
    p_array = collect(range(p_min, p_crit, length = N))
    Tdew(x) = dew_temperature(fluid,x,z)[1]
    Tbub(x) = bubble_temperature(fluid,x,z)[1]
    Td = Tdew.(p_array);
    Tb = Tbub.(p_array);
    Td[end] = crit_mix(fluid,z)[1];
    Tb[end] = crit_mix(fluid,z)[1];
    if nanfix
        Carnot.fix_nan!(Td)
        Carnot.fix_nan!(Tb)
    end
    s_bubble = similar(Tb);
    for i in eachindex(Tb)
        s_bubble[i] = entropy(fluid,p_array[i],Tb[i],z,phase = :liquid)./Clapeyron.molecular_weight(fluid,z)
    end
    s_dew = similar(Td);
    for i in eachindex(Td)
        s_dew[i] = entropy(fluid, p_array[i], Td[i],z,phase = :vapor)./Clapeyron.molecular_weight(fluid,z)
    end
    return Dict(
        :Tb => Tb,
        :Td => Td,
        :s_bubble => s_bubble,
        :s_dew => s_dew
    )
end

function plotting_data(prob::HeatPump,sol::SolutionState;N = 30,p_min = nothing)
    if isnothing(p_min)
        p_min = 0.5*sol.x[1]*101325
    end
    p_evap , p_cond = sol.x .* 101325
    T_evap_out = dew_temperature(prob.fluid, p_evap, prob.z)[1] + prob.ΔT_sh
    h_comp_in = enthalpy(prob.fluid, p_evap, T_evap_out,prob.z, phase = :vapor)
    crit = crit_mix!(prob)
    h_comp_out = Carnot.isentropic_compressor(p_evap, p_cond, prob.η_comp, h_comp_in, prob.z, prob.fluid, crit)
    
    p_comp_array = collect(range(p_evap, p_cond, length = N))
    f_h(p_out) = isentropic_compressor(p_evap,p_out,prob.η_comp,h_comp_in,prob.z,prob.fluid,crit)
    h_comp_array = f_h.(p_comp_array)
    T_ph(p,h) = Clapeyron.PH.temperature(prob.fluid, p, h, prob.z)
    T_comp_array = T_ph.(p_comp_array, h_comp_array)
    s_ph_vapour(p,h) = Clapeyron.PH.entropy(prob.fluid, p, h, prob.z)
    s_comp_array = s_ph_vapour.(p_comp_array, h_comp_array)./molecular_weight(prob.fluid,prob.z)
    T_cond_out = Clapeyron.bubble_temperature(prob.fluid, p_cond, prob.z)[1] - prob.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.fluid, p_cond, T_cond_out, prob.z)
    h_cond_array = collect(range(h_cond_out, h_comp_out, length = N))
    T_cond_array = T_ph.(p_cond, h_cond_array)
    s_ph(p,h) = Clapeyron.PH.entropy(prob.fluid, p, h, prob.z)
    s_cond_array = s_ph.(p_cond, h_cond_array)./molecular_weight(prob.fluid,prob.z)
    h_valve_in = h_cond_out
    h_valve_out = h_valve_in
    h_valve_array = collect(range(h_valve_in, h_valve_out, length = N))
    p_valve_array = collect(range(p_cond, p_evap, length = N))
    T_valve_array = T_ph.(p_valve_array, h_valve_array)
    s_valve_array = s_ph.(p_valve_array, h_valve_array)./molecular_weight(prob.fluid,prob.z)

    h_evap_in = h_valve_out
    T_evap_out = dew_temperature(prob.fluid, p_evap, prob.z)[1] + prob.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.fluid, p_evap, T_evap_out, prob.z)
    h_evap_array = collect(range(h_evap_in, h_evap_out, length = N))
    T_evap_array = T_ph.(p_evap, h_evap_array)
    s_evap_array = s_ph.(p_evap, h_evap_array)./molecular_weight(prob.fluid,prob.z)

    # secondary fluids
    T_cond_sf_array = collect(range(prob.T_cond_in, prob.T_cond_out, length = N))
    T_evap_sf_array = collect(range(prob.T_evap_out, prob.T_evap_in, length = N))

    return Dict(
        :s_comp_array => s_comp_array,
        :T_comp_array => T_comp_array,
        :s_cond_array => s_cond_array,
        :T_cond_array => T_cond_array,
        :s_valve_array => s_valve_array,
        :T_valve_array => T_valve_array,
        :s_evap_array => s_evap_array,
        :T_evap_array => T_evap_array,
        :T_cond_sf_array => T_cond_sf_array,
        :T_evap_sf_array => T_evap_sf_array
    ), p_min
end

function plotting_data(prob::ORC,sol::SolutionState;N = 30,p_min = nothing)
    if isnothing(p_min)
      p_min = 0.5*sol.x[2]*101325
    end
    p_evap, p_cond = sol.x .* 101325 # convert to Pa

    T_pump_in = bubble_temperature(prob.fluid, p_cond, prob.z)[1] - prob.ΔT_sc
    h_pump_in = Clapeyron.enthalpy(prob.fluid, p_cond, T_pump_in, prob.z)
    h_pump_out = Carnot.isentropic_pump(p_cond, p_evap, prob.η_pump, h_pump_in, prob.z, prob.fluid)
    h_pump_array = collect(range(h_pump_in, h_pump_out, length = N))
    T_ph(p,h) = Clapeyron.PH.temperature(prob.fluid, p, h, prob.z)
    T_pump_array = T_ph.(p_cond, h_pump_array)
    s_pt(p,t) = entropy(prob.fluid, p, t, prob.z)
    s_ph(p,h) = Clapeyron.PH.entropy(prob.fluid, p, h, prob.z)
    s_pump_array = s_pt.(p_cond, T_pump_array)./molecular_weight(prob.fluid,prob.z)
   

    T_evap_out = dew_temperature(prob.fluid, p_evap, prob.z)[1] + prob.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.fluid, p_evap, T_evap_out, prob.z)
    h_evap_in = h_pump_out
    h_evap_array = collect(range(h_evap_in, h_evap_out, length = N))
    T_evap_array = T_ph.(p_evap, h_evap_array)
    s_evap_array = s_ph.(p_evap, h_evap_array)./molecular_weight(prob.fluid,prob.z)
 

    h_exp_in = h_evap_out
    h_exp_out = Carnot.isentropic_expander(p_evap, p_cond, prob.η_expander, h_exp_in, prob.z, prob.fluid)

    p_exp_array = collect(range(p_evap, p_cond, length = N))
    f_h(p_out) = isentropic_expander(p_evap,p_out,prob.η_expander,h_exp_in,prob.z,prob.fluid)
    h_exp_array = f_h.(p_exp_array)
    T_exp_array = T_ph.(p_exp_array, h_exp_array) 
    s_exp_array = s_ph.(p_exp_array, h_exp_array)./molecular_weight(prob.fluid,prob.z)


    h_cond_in = h_exp_out
    T_cond_out = bubble_temperature(prob.fluid, p_cond, prob.z)[1] - prob.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.fluid, p_cond, T_cond_out, prob.z)
    h_cond_array = collect(range(h_cond_in, h_cond_out, length = N))
    T_cond_array = T_ph.(p_cond, h_cond_array)
    s_cond_array = s_ph.(p_cond, h_cond_array)./molecular_weight(prob.fluid,prob.z)
 
    # secondary fluids
    T_evap_sf_array = collect(range(prob.T_evap_out, prob.T_evap_in, length = N))

    T_cond_sf_array = collect(range(prob.T_cond_out, prob.T_cond_in, length = N))
    return Dict(
        :s_pump_array => s_pump_array,
        :T_pump_array => T_pump_array,
        :s_evap_array => s_evap_array,
        :T_evap_array => T_evap_array,
        :s_exp_array => s_exp_array,
        :T_exp_array => T_exp_array,
        :s_cond_array => s_cond_array,
        :T_cond_array => T_cond_array,
        :T_evap_sf_array => T_evap_sf_array,
        :T_cond_sf_array => T_cond_sf_array
    ) ,  p_min
end

function plotting_data(prob::HeatPumpRecuperator,sol::SolutionState;N = 30,p_min = nothing)
    if isnothing(p_min)
      p_min = 0.5*sol.x[1]*101325
  end
    p_evap, p_cond = sol.x .* 101325 # convert to Pa


    T_evap_out = Clapeyron.dew_temperature(prob.hp.fluid, p_evap, prob.hp.z)[1] + prob.hp.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.hp.fluid, p_evap, T_evap_out, prob.hp.z)
    T_cond_out = Clapeyron.bubble_temperature(prob.hp.fluid, p_cond, prob.hp.z)[1] - prob.hp.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.hp.fluid, p_cond, T_cond_out, prob.hp.z)

    q_ihex = Carnot.IHEX_Q(prob.hp.fluid,prob.ϵ,T_cond_out, p_cond, T_evap_out, p_evap, prob.hp.z)
    
    crit = crit_mix!(prob)
    h_comp_in =  q_ihex + h_evap_out  #enthalpy(prob.fluid, p_evap, T_evap_out,prob.z, phase = :vapor)
    h_comp_out = Carnot.isentropic_compressor(p_evap, p_cond, prob.hp.η_comp, h_comp_in, prob.hp.z, prob.hp.fluid, crit)

    p_comp_array = collect(range(p_evap, p_cond, length = N))
    f_h(p_out) = isentropic_compressor(p_evap,p_out,prob.hp.η_comp,h_comp_in,prob.hp.z,prob.hp.fluid,crit)
    h_comp_array = f_h.(p_comp_array)
    T_ph(p,h) = Clapeyron.PH.temperature(prob.hp.fluid, p, h, prob.hp.z)
    T_comp_array = T_ph.(p_comp_array, h_comp_array)
    s_ph_vapour(p,h) = Clapeyron.PH.entropy(prob.hp.fluid, p, h, prob.hp.z)
    s_comp_array = s_ph_vapour.(p_comp_array, h_comp_array)./molecular_weight(prob.hp.fluid,prob.hp.z)


    # Recuperator
    h_recoup_evap_in = h_evap_out
    h_recout_evap_out = h_recoup_evap_in + q_ihex
    h_recoup_evap_array = collect(range(h_recoup_evap_in, h_recout_evap_out, length = N))
    T_recoup_comp = T_ph.(p_evap, h_recoup_evap_array)
    s_recoup_comp = s_ph_vapour.(p_evap, h_recoup_evap_array)./molecular_weight(prob.hp.fluid,prob.hp.z)


    h_recoup_cond_in = h_cond_out
    h_recout_cond_out = h_recoup_cond_in - q_ihex
    h_recoup_cond_array = collect(range(h_recoup_cond_in, h_recout_cond_out, length = N))
    T_recoup_cond = T_ph.(p_cond, h_recoup_cond_array)
    s_recoup_cond = s_ph_vapour.(p_cond, h_recoup_cond_array)./molecular_weight(prob.hp.fluid,prob.hp.z)

    T_cond_out = Clapeyron.bubble_temperature(prob.hp.fluid, p_cond, prob.hp.z)[1] - prob.hp.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.hp.fluid, p_cond, T_cond_out, prob.hp.z)
    h_cond_array = collect(range(h_cond_out, h_comp_out, length = N))
    T_cond_array = T_ph.(p_cond, h_cond_array)
    s_ph(p,h) = Clapeyron.PH.entropy(prob.hp.fluid, p, h, prob.hp.z)
    s_cond_array = s_ph.(p_cond, h_cond_array)./molecular_weight(prob.hp.fluid,prob.hp.z)

    h_valve_in = h_cond_out - q_ihex
    h_valve_out = h_valve_in
    h_valve_array = collect(range(h_valve_in, h_valve_out, length = N))
    p_valve_array = collect(range(p_cond, p_evap, length = N))
    T_valve_array = T_ph.(p_valve_array, h_valve_array)
    s_valve_array = s_ph.(p_valve_array, h_valve_array)./molecular_weight(prob.hp.fluid,prob.hp.z)

    h_evap_in = h_valve_out
    T_evap_out = Clapeyron.dew_temperature(prob.hp.fluid, p_evap, prob.hp.z)[1] + prob.hp.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.hp.fluid, p_evap, T_evap_out, prob.hp.z)
    h_evap_array = collect(range(h_evap_in, h_evap_out, length = N))
    T_evap_array = T_ph.(p_evap, h_evap_array)
    s_evap_array = s_ph.(p_evap, h_evap_array)./molecular_weight(prob.hp.fluid,prob.hp.z)


    # secondary fluids
    T_cond_sf_array = collect(range(prob.hp.T_cond_in, prob.hp.T_cond_out, length = N))


    T_evap_sf_array = collect(range(prob.hp.T_evap_out, prob.hp.T_evap_in, length = N))
    dict = Dict(
        :s_comp_array => s_comp_array,
        :T_comp_array => T_comp_array,
        :s_recoup_comp => s_recoup_comp,
        :T_recoup_comp => T_recoup_comp,
        :s_recoup_cond => s_recoup_cond,
        :T_recoup_cond => T_recoup_cond,
        :s_cond_array => s_cond_array,
        :T_cond_array => T_cond_array,
        :s_valve_array => s_valve_array,
        :T_valve_array => T_valve_array,
        :s_evap_array => s_evap_array,
        :T_evap_array => T_evap_array,
        :T_cond_sf_array => T_cond_sf_array,
        :T_evap_sf_array => T_evap_sf_array
    ) , p_min
    return dict
end

function plotting_data(prob::ORCEconomizer,sol::SolutionState;N = 30, p_min = nothing)
    if isnothing(p_min)
      p_min = 0.5*sol.x[2]*101325
    end
    
    p_evap, p_cond = sol.x .* 101325 # convert to Pa
    T_pump_in = bubble_temperature(prob.orc.fluid, p_cond, prob.orc.z)[1] - prob.orc.ΔT_sc
    h_pump_in = Clapeyron.enthalpy(prob.orc.fluid, p_cond, T_pump_in, prob.orc.z)
    h_pump_out = Carnot.isentropic_pump(p_cond, p_evap, prob.orc.η_pump, h_pump_in, prob.orc.z, prob.orc.fluid)
    T_pump_out = Clapeyron.PH.temperature(prob.orc.fluid, p_evap, h_pump_out, prob.orc.z)

    h_pump_array = collect(range(h_pump_in, h_pump_out, length = N))
    T_ph(p,h) = Clapeyron.PH.temperature(prob.orc.fluid, p, h, prob.orc.z)
    T_pump_array = T_ph.(p_cond, h_pump_array)
    s_pt(p,t) = entropy(prob.orc.fluid, p, t, prob.orc.z)
    s_ph(p,h) = Clapeyron.PH.entropy(prob.orc.fluid, p, h, prob.orc.z)
    s_pump_array = s_pt.(p_cond, T_pump_array)./molecular_weight(prob.orc.fluid,prob.orc.z)


    T_evap_out = dew_temperature(prob.orc.fluid, p_evap, prob.orc.z)[1] + prob.orc.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.orc.fluid, p_evap, T_evap_out, prob.orc.z)
    h_exp_in = h_evap_out
    h_exp_out = Carnot.isentropic_expander(p_evap, p_cond, prob.orc.η_expander, h_exp_in, prob.orc.z, prob.orc.fluid)

    p_exp_array = collect(range(p_evap, p_cond, length = N))
    f_h(p_out) = isentropic_expander(p_evap,p_out,prob.orc.η_expander,h_exp_in,prob.orc.z,prob.orc.fluid)
    h_exp_array = f_h.(p_exp_array)
    T_exp_array = T_ph.(p_exp_array, h_exp_array)
    s_exp_array = s_ph.(p_exp_array, h_exp_array)./molecular_weight(prob.orc.fluid,prob.orc.z)

    
    # Economizer
    h_econ_in = h_exp_out
    T_exp_out = T_ph(p_cond,h_econ_in)
    q_ihex = Carnot.IHEX_Q(prob.orc.fluid,prob.ϵ,T_exp_out, p_cond, T_pump_out, p_evap, prob.orc.z)
    h_cond_in = h_econ_in - q_ihex
    h_evap_in = h_pump_out + q_ihex
    h_econ_cond_array = collect(range(h_exp_out, h_cond_in, length = N))
    T_econ_cond = T_ph.(p_cond, h_econ_cond_array)
    s_econ_cond = s_ph.(p_cond, h_econ_cond_array)./molecular_weight(prob.orc.fluid,prob.orc.z)

    h_econ_evap_array = collect(range(h_pump_out, h_evap_in, length = N))
    T_econ_evap = T_ph.(p_evap, h_econ_evap_array)
    s_econ_evap = s_ph.(p_evap, h_econ_evap_array)./molecular_weight(prob.orc.fluid,prob.orc.z)



    h_evap_array = collect(range(h_evap_in, h_evap_out, length = N))
    T_evap_array = T_ph.(p_evap, h_evap_array)
    s_evap_array = s_ph.(p_evap, h_evap_array)./molecular_weight(prob.orc.fluid,prob.orc.z)



    T_cond_out = bubble_temperature(prob.orc.fluid, p_cond, prob.orc.z)[1] - prob.orc.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.orc.fluid, p_cond, T_cond_out, prob.orc.z)
    h_cond_array = collect(range(h_cond_in, h_cond_out, length = N))
    T_cond_array = T_ph.(p_cond, h_cond_array)
    s_cond_array = s_ph.(p_cond, h_cond_array)./molecular_weight(prob.orc.fluid,prob.orc.z)

    # secondary fluids
    T_evap_sf_array = collect(range(prob.orc.T_evap_out, prob.orc.T_evap_in, length = N))

    T_cond_sf_array = collect(range(prob.orc.T_cond_out, prob.orc.T_cond_in, length = N))

    # @show minimum(T_evap_sf_array .- T_evap_array)
    dict = Dict(
        :s_pump_array => s_pump_array,
        :T_pump_array => T_pump_array,
        :s_econ_evap => s_econ_evap,
        :T_econ_evap => T_econ_evap,
        :s_econ_cond => s_econ_cond,
        :T_econ_cond => T_econ_cond,
        :s_cond_array => s_cond_array,
        :T_cond_array => T_cond_array,
        :s_exp_array => s_exp_array,
        :T_exp_array => T_exp_array,
        :s_evap_array => s_evap_array,
        :T_evap_array => T_evap_array,
        :T_cond_sf_array => T_cond_sf_array,
        :T_evap_sf_array => T_evap_sf_array
    ) , p_min
    return dict
end

@recipe function f_phase(fluid::EoSModel,z::AbstractVector;N = 30,p_min = nothing,nanfix = true) 
  @assert length(fluid.components) == length(z) "Components and composition vector length mismatch"
  if isnothing(p_min)
      p_min = 101325*0.4
  end
  _,p_crit,_ = crit_mix(fluid,z)
  p_array = collect(range(p_min, p_crit, length = N))
  Tdew(x) = dew_temperature(fluid,x,z)[1]
  Tbub(x) = bubble_temperature(fluid,x,z)[1]
  Td = Tdew.(p_array);
  Tb = Tbub.(p_array);
  Td[end] = crit_mix(fluid,z)[1];
  Tb[end] = crit_mix(fluid,z)[1];
    if nanfix
        Carnot.fix_nan!(Td)
        Carnot.fix_nan!(Tb)
    end
    linewidth := 2
  @series begin
    linestyle := :solid
    label := false
    ylabel := "Temperature (K)"
    xlabel := "Specific Entropy (J/K/kg)"
    title := "$(fluid.components)"
      sb = similar(Tb);
      for i in eachindex(Tb)
        sb[i] = entropy(fluid,p_array[i],Tb[i],z,phase = :liquid)./Clapeyron.molecular_weight(fluid,z)
      end
        (sb,Tb)
  end

  
    @series begin
    #   seriestype := :auto
    linestyle := :solid
    label := false
    ylabel := "Temperature (K)"
    xlabel := "Specific Entropy (J/K/kg)"
    title := "$(fluid.components)"
    
    sd = similar(Td);
      for i in eachindex(Td)
        sd[i] = entropy(fluid, p_array[i], Td[i],z,phase = :vapor)./Clapeyron.molecular_weight(fluid,z)
      end
        (sd,Td)
  end
end

@recipe function f_plot(prob::HeatPump,sol::SolutionState;N = 30,p_min = nothing,nanfix = true)
    hpdata,_p_min = plotting_data(prob,sol,N=N,p_min=p_min)
    phasedata = plotting_data(prob.fluid,prob.z;N=N,p_min=_p_min,nanfix=nanfix)
    

    # phase envelope - dew
    @series begin
        
        linestyle := :solid
        linewidth := 2
        markercolor := :red
        label := false
        (phasedata[:s_dew], phasedata[:Td])
    end
    @series begin
        # phase envelope - bubble
        linewidth := 2
        linestyle := :solid
        markercolor := :blue
        label := false
        (phasedata[:s_bubble], phasedata[:Tb])
    end

    @series begin
        # compressor
        linewidth := 2
        linestyle := :solid
        linecolor := :black
        label := false
        (hpdata[:s_comp_array], hpdata[:T_comp_array])
    end
    @series begin
        # condensor
        linewidth := 2
        linestyle := :solid
        linecolor := :black
        label := false
        (hpdata[:s_cond_array], hpdata[:T_cond_array])
    end
    @series begin
        # valve
        linewidth := 2
        linestyle := :solid
        linecolor := :black
        label := false
        (hpdata[:s_valve_array], hpdata[:T_valve_array])
    end
    @series begin
        # evaporator
        linewidth := 2
        linestyle := :solid
        linecolor := :black
        label := false
        (hpdata[:s_evap_array], hpdata[:T_evap_array])
    end
    @series begin
        # sf evaporator
        linewidth := 2
        linestyle := :dash
        linecolor := :blue
        label := "Secondary Fluid Evaporator"
        (hpdata[:s_evap_array], hpdata[:T_evap_sf_array])
    end
    @series begin
        # sf condenser
        linewidth := 2
        linestyle := :dash
        linecolor := :red
        label := "Secondary Fluid Condenser"
        ylabel := "Temperature (K)"
        xlabel := "Specific Entropy (J/K/kg)"
        title := "$(prob.fluid.components)"
        (hpdata[:s_cond_array], hpdata[:T_cond_sf_array])
    end
end
@recipe function f_plot(prob::ORC,sol::SolutionState;N = 30,p_min = nothing,nanfix = true)
    orcdata,_p_min = plotting_data(prob,sol,N=N,p_min=p_min)
    phasedata = plotting_data(prob.fluid,prob.z;N=N,p_min=_p_min,nanfix=nanfix)
    

    # phase envelope - dew
    @series begin
        
        linestyle := :solid
        linewidth := 2
        markercolor := :red
        label := false
        (phasedata[:s_dew], phasedata[:Td])
    end
    @series begin
        # phase envelope - bubble
        linewidth := 2
        linestyle := :solid
        markercolor := :blue
        label := false
        (phasedata[:s_bubble], phasedata[:Tb])
    end

    @series begin
        # pump
        linewidth := 2
        linestyle := :solid
        linecolor := :black
        label := false
        (orcdata[:s_pump_array], orcdata[:T_pump_array])
    end
    @series begin
        # evaporator
        linewidth := 2
        linestyle := :solid
        linecolor := :black
        label := false
        (orcdata[:s_evap_array], orcdata[:T_evap_array])
    end
    @series begin
        # expander
        linewidth := 2
        linestyle := :solid
        linecolor := :black
        label := false
        (orcdata[:s_exp_array], orcdata[:T_exp_array])
    end
    @series begin
        # condensor
        linewidth := 2
        linestyle := :solid
        linecolor := :black
        label := false
        (orcdata[:s_cond_array], orcdata[:T_cond_array])
    end
    @series begin
        # sf evaporator
        linewidth := 2
        linestyle := :dash
        linecolor := :blue
        label := "Secondary Fluid Evaporator"
        (orcdata[:s_evap_array], orcdata[:T_evap_sf_array])
    end
    @series begin
        # sf condenser
        linewidth := 2
        linestyle := :dash
        linecolor := :red
        label := "Secondary Fluid Condenser"
        ylabel := "Temperature (K)"
        xlabel := "Specific Entropy (J/K/kg)"
        title := "$(prob.fluid.components)"
        (orcdata[:s_cond_array], orcdata[:T_cond_sf_array])
    end
    
end

@recipe function f_plot(prob::HeatPumpRecuperator,sol::SolutionState;N = 30,p_min = nothing,nanfix = true)
    hpdata,_p_min = plotting_data(prob,sol,N=N,p_min=p_min)
    phasedata = plotting_data(prob.hp.fluid,prob.hp.z;N=N,p_min=_p_min,nanfix=nanfix) 

    # phase envelope - dew
    @series begin
        
        linestyle := :solid
        linewidth := 2
        markercolor := :red
        label := false
        (phasedata[:s_dew], phasedata[:Td])
    end
    @series begin
        # phase envelope - bubble
        linewidth := 2
        linestyle := :solid
        markercolor := :blue
        label := false
        (phasedata[:s_bubble], phasedata[:Tb])
    end

    @series begin
        # compressor
        linewidth := 2
        linestyle := :solid
        linecolor := :black
        label := false
        (hpdata[:s_comp_array], hpdata[:T_comp_array])
    end
    @series begin
        # recuperator - comp side
        linewidth := 2
        linestyle := :dot
        linecolor := :green
        label := false
        (hpdata[:s_recoup_comp], hpdata[:T_recoup_comp])
    end
    @series begin
        # recuperator - cond side
        linewidth := 2
        linestyle := :dot
        linecolor := :green
        label := false
        (hpdata[:s_recoup_cond], hpdata[:T_recoup_cond])
    end
    @series begin
        # condensor
        linewidth := 2
        linestyle := :solid
        linecolor := :black
        label := false
        (hpdata[:s_cond_array], hpdata[:T_cond_array])
    end
    @series begin
        # valve
        linewidth := 2
        linestyle := :solid
        linecolor := :black
        label := false
        (hpdata[:s_valve_array], hpdata[:T_valve_array])
    end
    @series begin 
        # sf evaporator
        linewidth := 2
        linestyle := :dash
        linecolor := :blue
        label := "Secondary Fluid Evaporator"
        (hpdata[:s_evap_array], hpdata[:T_evap_sf_array])
    end

    @series begin 
        # sf condenser
        linewidth := 2
        linestyle := :dash
        linecolor := :red
        label := "Secondary Fluid Condenser"
        (hpdata[:s_cond_array], hpdata[:T_cond_sf_array])
    end
    @series begin
        # evaporator
        linewidth := 2
        linestyle := :solid
        linecolor := :black
        label := false
        ylabel := "Temperature (K)"
        xlabel := "Specific Entropy (J/K/kg)"
        title := "$(prob.hp.fluid.components)"
        (hpdata[:s_evap_array], hpdata[:T_evap_array])
    end
end


@recipe function f_plot(prob::ORCEconomizer,sol::SolutionState;N = 30,p_min = nothing,nanfix = true)
    orcdata,_p_min = plotting_data(prob,sol,N=N,p_min=p_min)
    phasedata = plotting_data(prob.orc.fluid,prob.orc.z;N=N,p_min=_p_min,nanfix=nanfix)


    # phase envelope - dew
    @series begin
        
        linestyle := :solid
        linewidth := 2
        markercolor := :red
        label := false
        (phasedata[:s_dew], phasedata[:Td])
    end
    @series begin
        # phase envelope - bubble
        linewidth := 2
        linestyle := :solid
        markercolor := :blue
        label := false
        (phasedata[:s_bubble], phasedata[:Tb])
    end

    @series begin
        # pump
        linewidth := 2
        linestyle := :solid
        linecolor := :black
        label := false
        (orcdata[:s_pump_array], orcdata[:T_pump_array])
    end
    @series begin
        # evaporator
        linewidth := 2
        linestyle := :solid
        linecolor := :black
        label := false
        (orcdata[:s_evap_array], orcdata[:T_evap_array])
    end
    @series begin
        # expander
        linewidth := 2
        linestyle := :solid
        linecolor := :black
        label := false
        (orcdata[:s_exp_array], orcdata[:T_exp_array])
    end
    @series begin
        # condensor
        linewidth := 2
        linestyle := :solid
        linecolor := :black
        label := false
        (orcdata[:s_cond_array], orcdata[:T_cond_array])
    end
    @series begin
        # sf evaporator
        linewidth := 2
        linestyle := :dash
        linecolor := :blue
        label := "Secondary Fluid Evaporator"
        (orcdata[:s_evap_array], orcdata[:T_evap_sf_array])
    end
    @series begin
        # sf condenser
        linewidth := 2
        linestyle := :dash
        linecolor := :red
        label := "Secondary Fluid Condenser"
        (orcdata[:s_cond_array], orcdata[:T_cond_sf_array])
    end

    @series begin
        #eco evap
        linewidth := 2
        linestyle := :dot
        linecolor := :green
        label := false
        (orcdata[:s_econ_evap], orcdata[:T_econ_evap])
    end

    @series begin
        #eco cond
        linewidth := 2
        linestyle := :dot
        linecolor := :green
        label := false
        ylabel := "Temperature (K)"
        xlabel := "Specific Entropy (J/K/kg)"
        title := "$(prob.orc.fluid.components)"
        (orcdata[:s_econ_cond], orcdata[:T_econ_cond])
    end
end

export f_plot

# function plotting_data(prob::HeatPumpTranscritical,x::AbstractVector,T::AbstractVector;N = 30,p_min = nothing, nanfix = true)
#     if isnothing(p_min)
#       p_min = 0.5*x[1]*101325
#   end
#     fluid = prob.fluid
#     z = prob.z
#     T_crit,p_crit,_ = crit_mix(fluid,z)

#     T_evap_out = dew_temperature(fluid, x[1]*101325, z)[1] + T[1]
#     T_cond_out = T_crit - T[2]

#     p_evap = x[1]*101325
#     p_cond = x[2]*p_crit
#     p_array = 
# end