module CarnotMetaheuristicsExt

using Carnot
using Metaheuristics
using Carnot.CommonSolve
import Carnot: optimize

function _build_residual(prob::HeatPump, N::Int)
    f = let prob = prob, N = N
        x -> begin 
            nc = length(prob.fluid.components)
            return nc == 1 ? F_pure(prob, x) : F(prob, x, N = N)
        end
    end
end

function objective(prob::HeatPump,param::ThermoCycleParameters,x::AbstractVector)
    @assert length(x) == 2 "Only super and sub cool temperatures"
    prob.ΔT_sh = x[1]
    prob.ΔT_sc = x[2] 
    sol =  solve(prob,param)
    if Carnot.norm(sol.residuals) > 1e-3
        return 0.0
    else
        return COP(prob,sol)
    end
end

function objective(prob::ORC,param::ThermoCycleParameters,x::AbstractVector)
    @assert length(x) == 2 "Only super and sub cool temperatures"
    prob.ΔT_sh = x[1]
    prob.ΔT_sc = x[2] 
    sol =  solve(prob,param)
    if Carnot.norm(sol.residuals) > 1e-3
        return 0.0
    else
        return η(prob,sol)
    end
end

function _build_objective(prob::HeatPump,param::ThermoCycleParameters)
    f = let prob = prob
        x -> begin
            objective(prob,param,x)
        end 
    end
end
function _build_objective(prob::ORC,param::ThermoCycleParameters)
    f = let prob = prob
        x -> begin
            objective(prob,param,x)
        end 
    end
end
export _build_objective
function generate_optimization_bounds(prob::HeatPump)
    ΔT_sh_min = 0.0
    ΔT_sh_max = prob.T_evap_in - prob.T_evap_out
    ΔT_sc_min = 0.0
    ΔT_sc_max = prob.T_cond_out - prob.T_cond_in
    lb = [ΔT_sh_min,ΔT_sc_min]
    ub = [ΔT_sh_max,ΔT_sc_max]
    return lb,ub
end

function generate_optimization_bounds(prob::ORC)
    ΔT_sh_min = 3.0
    ΔT_sh_max = prob.T_evap_in - prob.T_evap_out
    ΔT_sc_min = 3.0
    ΔT_sc_max = prob.T_cond_out - prob.T_cond_in
    lb = [ΔT_sh_min,ΔT_sc_min]
    ub = [ΔT_sh_max,ΔT_sc_max]
    return lb,ub
end


function optimize(prob::HeatPump,
    alg::Metaheuristics.AbstractAlgorithm,param::ThermoCycleParameters)

    @time "Building Objective function..." begin
    ℓ = _build_objective(prob,param)
    end
    @time "Generating bounds ..." begin
    lb,ub = generate_optimization_bounds(prob)
    end
    bounds = Metaheuristics.boxconstraints(lb = lb, ub = ub)
    opt_result = Metaheuristics.optimize(ℓ,bounds,alg)
    
    x_best = Metaheuristics.minimizer(opt_result)

    loss_opt_M = Metaheuristics.minimum(opt_result)

    hp_opt = prob
    hp_opt.ΔT_sh = x_best[1]
    hp_opt.ΔT_sc = x_best[2]
    sol_best = Carnot.solve(hp_opt,param)

    return x_best,sol_best
end

function optimize(prob::ORC,
    alg::Metaheuristics.AbstractAlgorithm,param::ThermoCycleParameters)

    @time "Building Objective function..." begin
    ℓ = _build_objective(prob,param)
    end
    @time "Generating bounds ..." begin
    lb,ub = generate_optimization_bounds(prob)
    end
    bounds = Metaheuristics.boxconstraints(lb = lb, ub = ub)
    opt_result = Metaheuristics.optimize(ℓ,bounds,alg)
    
    x_best = Metaheuristics.minimizer(opt_result)

    loss_opt_M = Metaheuristics.minimum(opt_result)

    orc_opt = prob
    orc_opt.ΔT_sh = x_best[1]
    orc_opt.ΔT_sc = x_best[2]
    sol_best = Carnot.solve(orc_opt,param)

    return x_best,sol_best
end


export optimize


end #module