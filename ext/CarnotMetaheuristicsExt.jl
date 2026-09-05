module CarnotMetaheuristicsExt

using Carnot
using Metaheuristics
using Carnot.CommonSolve
import Carnot: optimize

function _build_residual(prob::HeatPump, N::Int)
    f = let prob = prob, N = N
        x -> begin 
            nc = length(prob.fluid.components)s
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

function _build_objective(
    prob::HeatPump,
    param::ThermoCycleParameters,
    algo::Metaheuristics.AbstractAlgorithm,
)

    if algo.options.parallel_evaluation
        return let prob = prob, param = param
            X -> begin
                fitness = zeros(size(X, 1))

                Threads.@threads for i in axes(X, 1)
                    fitness[i] = objective(prob, param, X[i, :])
                end

                fitness
            end
        end
    else
        return let prob = prob, param = param
            x -> objective(prob, param, x)
        end
    end
end
function _build_objective(
    prob::ORC,
    param::ThermoCycleParameters,
    algo::Metaheuristics.AbstractAlgorithm,
    )   

    if algo.options.parallel_evaluation
        return let prob = prob, param = param
            X -> begin
                fitness = zeros(size(X, 1))

                Threads.@threads for i in axes(X, 1)
                    fitness[i] = objective(prob, param, X[i, :])
                end

                fitness
            end
        end
    else
        return let prob = prob, param = param
            x -> objective(prob, param, x)
        end
    end
end

function _build_objective(
    prob::TranscriticalORC,
    param::TranscriticalParamters,
    algo::Metaheuristics.AbstractAlgorithm,
    )

    if algo.options.parallel_evaluation
        return let prob = prob, param = param
            X -> begin
                fitness = zeros(size(X, 1))

                Threads.@threads for i in axes(X, 1)
                    fitness[i] = η(prob, X[i, :], param)
                end

                fitness
            end
        end
    else
        return let prob = prob, param = param
            x -> η(prob, x, param)
        end
    end
end


function _build_objective(
    prob::HeatPumpTranscritical,
    param::TranscriticalParamters,
    algo::Metaheuristics.AbstractAlgorithm,
    )

    if algo.options.parallel_evaluation
        return let prob = prob, param = param
            X -> begin
                fitness = zeros(size(X, 1))

                Threads.@threads for i in axes(X, 1)
                    fitness[i] = COP(prob, X[i, :], param)
                end

                fitness
            end
        end
    else
        return let prob = prob, param = param
            x -> COP(prob, x, param)
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
    ℓ = _build_objective(prob,param,alg)
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
    ℓ = _build_objective(prob,param,alg)
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


function optimize(prob::TranscriticalORC,alg::Metaheuristics.AbstractAlgorithm,param::TranscriticalParamters)
    ℓ = _build_objective(prob,param,alg)
    # generate box 
    lb,ub = generate_box(prob,param)
    bounds = Metaheuristics.boxconstraints(lb = lb, ub = ub)
    opt_result = Metaheuristics.optimize(ℓ,bounds,alg)

    x_best = Metaheuristics.minimizer(opt_result)
    loss_opt_M = Metaheuristics.minimum(opt_result)
    Δ,_ = Carnot.F(prob,x_best, N = param.N)

    sol = SolutionState(x_best,opt_result.f_calls,opt_result.iteration,Δ,lb,ub,false,2,NaN,NaN,:transcritical_optimal)
    return sol,loss_opt_M
end


function optimize(prob::HeatPumpTranscritical,alg::Metaheuristics.AbstractAlgorithm,param::TranscriticalParamters)
    ℓ = _build_objective(prob,param,alg)
    # generate box 
    lb,ub = generate_box(prob,param)
    bounds = Metaheuristics.boxconstraints(lb = lb, ub = ub)
    opt_result = Metaheuristics.optimize(ℓ,bounds,alg)

    x_best = Metaheuristics.minimizer(opt_result)
    loss_opt_M = Metaheuristics.minimum(opt_result)
    Δ,_ = Carnot.F(prob,x_best, N = param.N)

    sol = SolutionState(x_best,opt_result.f_calls,opt_result.iteration,Δ,lb,ub,false,2,NaN,NaN,:transcritical_optimal)
    return sol,loss_opt_M
end

function _build_objective(
    prob::OptHeatPump,
    param::DirectOptParameters,
    algo::Metaheuristics.AbstractAlgorithm,
    )

    f = let prob = prob, param = param
        x -> COP(prob, x, param)
    end

    if algo.options.parallel_evaluation

        f_parallel = let f = f
            X -> begin
                fitness = zeros(size(X, 1))

                Threads.@threads for i in axes(X, 1)
                    fitness[i] = f(X[i, :])
                end

                fitness
            end
        end

        return f_parallel
    end

    return f
end

function optimize(prob::OptHeatPump,alg::Metaheuristics.AbstractAlgorithm,param::DirectOptParameters)
    @info "Building objective function..."
    @time ℓ = _build_objective(prob,param,alg)
    @info "Objective function built."
    # generate box 
    @info "Generating box constraints..."
    lb,ub = generate_box(prob,param)
    bounds = Metaheuristics.boxconstraints(lb = lb, ub = ub)
    opt_result = Metaheuristics.optimize(ℓ,bounds,alg)
    x_best = Metaheuristics.minimizer(opt_result)
    loss_opt_M = Metaheuristics.minimum(opt_result)
    Δ,_ = Carnot.F(prob,x_best, N = 20)

    sol = SolutionState(x_best,opt_result.f_calls,opt_result.iteration,Δ,lb,ub,false,2,NaN,NaN,:optimal_no_solve_solution)
    return sol,loss_opt_M
end


export optimize


end #module