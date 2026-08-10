# function F(prob::ThermoCycleProblem,x::AbstractVector{T},N::Int,internal_pinch::Bool) where T<:Real
#     return F(prob,x,N=N,internal_pinch = internal_pinch)
# end

"""
`generate_initial_point(prob::HeatPump,lb::AbstractVector{T},ub::AbstractVector{T})` 
Generates initial point for solving the system.
"""
function generate_initial_point(prob::HeatPump,lb::AbstractVector{T},ub::AbstractVector{T},x0_init::Symbol) where T<: Real
    if x0_init == :default
        return [
            minimum(lb), maximum(ub)
        ]
    end
    if x0_init == :average 
        return (lb + ub)./2
    end    
end

"""
`generate_initial_point(prob::ORC,lb::AbstractVector{T},ub::AbstractVector{T})` 
Generates initial point for solving the system.
"""
function generate_initial_point(prob::ORC,lb::AbstractVector{T},ub::AbstractVector{T},x0_init::Symbol) where T<: Real
    if x0_init == :default
        return [
            minimum(lb), maximum(ub)
        ]
    end
    if x0_init == :average 
        return (lb + ub)./2
    end  
end

"""
`generate_initial_point(prob::ORCEconomizer,lb::AbstractVector{T},ub::AbstractVector{T})` 
Generates initial point for solving the system.
"""
function generate_initial_point(prob::ORCEconomizer,lb::AbstractVector{T},ub::AbstractVector{T},x0_init::Symbol) where T<: Real
    if x0_init == :default
        return [
            minimum(lb), maximum(ub)
        ]
    end
    if x0_init == :average 
        return (lb + ub)./2
    end     
end

"""
`generate_initial_point(prob::HeatPumpRecuperator,lb::AbstractVector{T},ub::AbstractVector{T})` 
Generates initial point for solving the system.
"""
function generate_initial_point(prob::HeatPumpRecuperator,lb::AbstractVector{T},ub::AbstractVector{T},x0_init::Symbol) where T<: Real
    if x0_init == :default
        return [
            minimum(lb), maximum(ub)
        ]
    end
    if x0_init == :average 
        return (lb + ub)./2
    end       
end

"""
`generate_initial_point(prob::HeatPumpVarEff,lb::AbstractVector{T},ub::AbstractVector{T})` 
Generates initial point for solving the system.
"""
function generate_initial_point(prob::HeatPumpVarEff,lb::AbstractVector{T},ub::AbstractVector{T},x0_init::Symbol) where T<: Real
    if x0_init == :default
        return [
            minimum(lb), maximum(ub)
        ]
    end
    if x0_init == :average 
        return (lb + ub)./2
    end       
end

function generate_initial_point(prob::HeatPumpTranscritical,lb::AbstractVector{T},ub::AbstractVector{T},x0_init::Symbol) where T<: Real
    if x0_init == :default
        return [
            minimum(lb), maximum(ub)
        ]
    end
    if x0_init == :average 
        return @error "not implemented"
    end       
end


"""
`generate_box_solve_bounds(prob::HeatPump) -> lb, ub`
Generates lower and upper bounds for the heat pump problem based on its parameters.
"""
function generate_box_solve_bounds(prob::HeatPump)
    Tcrit,pcrit,_ = crit_mix(prob.fluid, prob.z)
    lb = zeros(eltype(prob.z), 2)
    ub = zeros(eltype(prob.z), 2)
    if prob.T_cond_out > Tcrit
        psat_max = 0.99*pcrit 
    else 
        psat_max = bubble_pressure(prob.fluid,prob.T_cond_out + prob.pp_cond + prob.ΔT_sc,prob.z)[1] 
        if !isfinite(psat_max)
            psat_max = 0.99*pcrit 
        end
    end

    psat_min = dew_pressure(prob.fluid,prob.T_evap_out - prob.pp_evap - prob.ΔT_sh,prob.z)[1]
    
    ub[1] = psat_max
    lb[1] = psat_min

    ub[2] = psat_max
    lb[2] = psat_min
    if !isfinite(psat_max)
        throw(error("The upper bound on pressure is not finite. Check cycle parameters. Possible error on bubble_pressure calculation."))
    end
    return lb./101325, ub./101325 # normalize to 101325 Pa
end


function generate_box_solve_bounds(prob::HeatPumpTranscritical)
    #Tcrit,pcrit,_ = crit_mix(prob.fluid, prob.z)
    psat_min_evap = dew_pressure(prob.fluid,prob.T_evap_out - prob.pp_evap - prob.ΔT_sh,prob.z)[1]
    psat_max_evap = bubble_pressure(prob.fluid,prob.T_evap_in ,prob.z)[1]

    p_cond_min = 1.0 
    p_cond_max = 2.0

    lb = zeros(eltype(prob.z), 2)
    ub = zeros(eltype(prob.z), 2)
    lb[1] = psat_min_evap./101325; lb[2]= p_cond_min
    ub[1] = psat_max_evap ./101325; ub[2] = p_cond_max
    return lb, ub
end

"""
`generate_box_solve_bounds(prob::HeatPumpRecuperator) -> lb, ub`
Generates lower and upper bounds for the heat pump problem based on its parameters.
"""
function generate_box_solve_bounds(prob::HeatPumpRecuperator)
    return generate_box_solve_bounds(prob.hp)
end

"""
`generate_box_solve_bounds(prob::ORC) -> lb, ub`
Generates lower and upper bounds for the heat pump problem based on its parameters.
"""
function generate_box_solve_bounds(prob::ORC)
    Tcrit,pcrit,_ = crit_mix(prob.fluid, prob.z)
    lb = zeros(eltype(prob.z), 2)
    ub = zeros(eltype(prob.z), 2)
    if prob.T_evap_in > Tcrit
         psat_max = 0.95*pcrit
    end
    if prob.T_evap_out > Tcrit 
        throw(error("For now we handel subcritical ORC. The outlet temperature of the evap : 
        $(prob.T_evap) is higher than critical temperature of the fluid $Tcrit, this will not allow to meet the pinch points.
        "))
    end
    psat_min_evap = dew_pressure(prob.fluid,prob.T_evap_out - prob.pp_evap- prob.ΔT_sh,prob.z)[1]
    psat_max_evap = bubble_pressure(prob.fluid,prob.T_evap_in - prob.pp_evap- prob.ΔT_sh,prob.z)[1]
    if !isfinite(psat_max_evap)
        psat_max_evap = 0.95*pcrit
    end
    psat_min_cond = dew_pressure(prob.fluid,prob.T_cond_in + prob.pp_cond + prob.ΔT_sc,prob.z)[1]
    psat_max_cond = bubble_pressure(prob.fluid,prob.T_cond_out + prob.pp_cond + prob.ΔT_sc,prob.z)[1]
    @assert psat_min_cond <= psat_max_cond "Box generation min has to be less than max"
    @assert psat_min_evap <= psat_max_evap "Box generation min has to be less than max"
    # psat_min = dew_pressure(prob.fluid,prob.T_cond_in,prob.z)[1]
    # psat_max = bubble_pressure(prob.fluid,prob.T_evap_in - prob.ΔT_sh,prob.z)[1] #pcrit*0.9#dew_pressure(prob.fluid,prob.T_evap_in,prob.z)[1]
    ub[1] = psat_max_evap#dew_pressure(prob.fluid,prob.T_evap_in - prob.pp_evap - prob.ΔT_sh,prob.z)[1] # evaporator pressure
    lb[1] = psat_min_evap#bubble_pressure(prob.fluid,prob.T_evap_out - prob.pp_evap - prob.ΔT_sh,prob.z)[1] # evaporator pressure
    ub[2] = psat_max_cond
    lb[2] = psat_min_cond 

    return lb./101325, ub./101325 # normalize to 101325 Pa
end

"""
`generate_box_solve_bounds(prob::ORCEconomizer) -> lb, ub`
Generates lower and upper bounds for the heat pump problem based on its parameters.
"""
function generate_box_solve_bounds(prob::ORCEconomizer)
    return generate_box_solve_bounds(prob.orc)
end

"""
`generate_box_solve_bounds(prob::HeatPumpVarEff) -> lb, ub`
Generates lower and upper bounds for the heat pump problem based on its parameters.
"""
function generate_box_solve_bounds(prob::HeatPumpVarEff)
    Tcrit,pcrit,_ = crit_mix(prob.fluid, prob.z)
    lb = zeros(eltype(prob.z), 2)
    ub = zeros(eltype(prob.z), 2)
    if prob.T_cond_out > Tcrit
        psat_max = 0.99*pcrit 
    else 
        psat_max = bubble_pressure(prob.fluid,prob.T_cond_out + prob.pp_cond + prob.ΔT_sc,prob.z)[1] 
        if !isfinite(psat_max)
            psat_max = 0.99*pcrit 
        end
    end

    psat_min = dew_pressure(prob.fluid,prob.T_evap_out - prob.pp_evap - prob.ΔT_sh,prob.z)[1]
    
    ub[1] = psat_max
    lb[1] = psat_min

    ub[2] = psat_max
    lb[2] = psat_min
    if !isfinite(psat_max)
        throw(error("The upper bound on pressure is not finite. Check cycle parameters. Possible error on bubble_pressure calculation."))
    end
    return lb./101325, ub./101325 # normalize to 101325 Pa
end

"""
`ThermoCycleParameters` -  A struct to hold the solver parameters of the nonlinear solver.

- `N::Int`: Heat Exchanger discretization.
- `max_iters::Int`: Maximum number of iterations
- `autodiff::Bool`: A flag indicating whether automatic differentiation was used.
- `fd_order::Int`: The order of finite difference used if autodiff is false.
- `xtol::Float64`: convergece criteria on `x`.
- `ftol::Float64`: convergece criteria on `f`.
- `restart_TOL::Float64`: Restrat strategy with tolerace.
- `internal_pinch::Bool` : Check for interal pinch for mixtures
"""
mutable struct ThermoCycleParameters{AD}
    N::Int
    autodiff::Val{AD}
    internal_pinch::Bool
    fd_order::Int
    xtol::Float64
    ftol::Float64
    restart_TOL::Float64
    max_iters::Int
    x0_init::Symbol
    verbose::Bool
end

function ThermoCycleParameters(; 
    N::Int = 20,
    autodiff::Union{Bool,Val{true},Val{false}} = Val{true}(),
    internal_pinch::Bool = true,
    fd_order::Int = 2,
    xtol::Real = 1e-6,
    ftol::Real = 1e-6,
    restart_TOL::Real = 1e-3,
    max_iters::Int = 100,
    x0_init::Symbol = :default,
    verbose::Bool = false
)
    N > 0 || error("N must be positive")
    is_ad = autodiff == false || autodiff == Val{false}()
    if !is_ad && fd_order < 2
        error("If autodiff is false, fd_order must be ≥ 2 (higher-order finite differences required).")
    end

    ad_typed = autodiff isa Bool ? Val(autodiff) : autodiff

    @assert x0_init == :default || x0_init == :average "x0 initizer has two versions :default and :average" 
    return ThermoCycleParameters(N, ad_typed, internal_pinch, fd_order, xtol, ftol, restart_TOL, max_iters,x0_init, verbose)
end

function switch_x0(x0_init::Symbol)
    if x0_init == :default
        return :average
    end
    if x0_init == :average
        return :default
    end
    return :invalid_switch_x0
end

function _build_residual(prob::HeatPump, N::Int)
    f = let prob = prob, N = N
        x -> begin 
            nc = length(prob.fluid.components)
            return nc == 1 ? F_pure(prob, x) : F(prob, x, N = N)
        end
    end
end

function _build_residual(prob::ORC, N::Int)
    f = let prob = prob, N = N
        x -> begin 
            nc = length(prob.fluid.components)
            return nc == 1 ? F_pure(prob, x) : F(prob, x, N = N)
        end
    end
end

function _build_residual(prob::ThermoCycleProblem, N::Int)
    return x -> F(prob, x, N = N)
end

function _build_least_squares_objective(prob, N)
    f = let prob = prob, N = N
        x -> begin 
            nc = length(prob.fluid.components)
            θ = nc == 1 ? F_pure(prob, x) : F(prob, x, N = N)
            return sum(abs2,θ)/2
        end
    end
end

function _build_least_squares_objective(prob::ThermoCycleProblem, N)
    f = let prob = prob, N = N
        x -> begin 
            nc = length(prob.fluid.components)
            θ = F(prob, x, N = N)
            return sum(abs2,θ)/2
        end
    end
end


function solve_ad(prob::ThermoCycleProblem,lb::AbstractVector,ub::AbstractVector;N::Int64 = 20,restart_TOL = 1e-3,xtol = 1e-8,ftol = 1e-8,max_iters = 1000,x0_init::Symbol = :default,verbose::Bool = false)
    f = _build_residual(prob, N)
    x0 = generate_initial_point(prob,lb,ub,x0_init)
    sol = constrained_newton_ad(f, x0, lb, ub; xtol = xtol, ftol = ftol, max_iters = max_iters,verbose = verbose)
    sol.soltype = :subcritical
    if norm(sol.residuals)  > restart_TOL
        # f(x::AbstractVector{T}) where {T<:Real} = F_transcritical(prob, x,N = N)
        # x0 = generate_bounds(prob,lb,ub)
        # sol = constrained_newton_ad(f, x0, lb, ub; xtol = xtol, ftol = ftol, max_iters = max_iters)
        # sol.soltype = :transcritical
        x0_init = switch_x0(x0_init)
        x0 = generate_initial_point(prob,lb,ub,x0_init)
        sol = constrained_newton_ad(f, x0, lb, ub; xtol = xtol, ftol = ftol, max_iters = max_iters) 
        sol.soltype = :subcritical_restart_mode
    end
    return sol
end

function solve_fd(prob::ThermoCycleProblem,lb::AbstractVector,ub::AbstractVector;N::Int64 = 20,restart_TOL = 1e-3,fd_order = 2,xtol = 1e-8,ftol = 1e-8,max_iters = 1000,x0_init::Symbol = :default,verbose::Bool = false)
    f = _build_residual(prob, N)
    x0 = generate_initial_point(prob,lb,ub,x0_init)
    sol = constrained_newton_fd(f, x0, lb, ub; xtol = xtol, ftol = ftol, max_iters = max_iters,fd_order = fd_order,verbose = verbose)
    sol.soltype = :subcritical
    if norm(sol.residuals) > restart_TOL
        x0_init = switch_x0(x0_init)
        x0 = generate_initial_point(prob,lb,ub,x0_init)
        sol = constrained_newton_fd(f, x0, lb, ub; xtol = xtol, ftol = ftol, max_iters = max_iters,fd_order = fd_order)
        sol.soltype = :subcritical_restart_mode
    end
    return sol
end


"""
Solves for pressure values in HP and ORC cycles for the given glide and problem parameters. 
Define those problems in the respective structs. 
For now the default box-nonlinear solver is newton-raphson, but this can be changed to other solvers in the future.
"""
function CommonSolve.solve(prob::ThermoCycleProblem; autodiff::Union{Bool,Val{true},Val{false}} = Val(true), fd_order =2 , N::Int64 = 20,restart_TOL = 1e-3,xtol = 1e-6,ftol = 1e-6,max_iters = 1000,x0_init::Symbol=:default,verbose::Bool = false)
    alg = ThermoCycleParameters(;autodiff,fd_order,N,restart_TOL,xtol,ftol,max_iters,x0_init,verbose)
    return CommonSolve.solve(prob,alg)
end

"""
Solves for pressure values in HP and ORC cycles for the given glide and problem parameters. 
Define those problems in the respective structs. 
For now the default box-nonlinear solver is newton-raphson, but this can be changed to other solvers in the future.
"""
function CommonSolve.solve(prob::ThermoCycleProblem,param::ThermoCycleParameters)
    Base.@nospecialize prob
    lb,ub = generate_box_solve_bounds(prob)
    return _solve_with_params(prob, param, lb, ub)
end

@noinline function _solve_with_params(prob, param::ThermoCycleParameters{false}, lb, ub)
    return solve_fd(
        prob,
        lb,
        ub;
        N = param.N,
        fd_order = param.fd_order,
        restart_TOL = param.restart_TOL,
        xtol = param.xtol,
        ftol = param.ftol,
        max_iters = param.max_iters,
        x0_init = param.x0_init,
        verbose = param.verbose,
    )
end

@noinline function _solve_with_params(prob, param::ThermoCycleParameters{true}, lb, ub)
    return solve_ad(
        prob,
        lb,
        ub;
        N = param.N,
        restart_TOL = param.restart_TOL,
        xtol = param.xtol,
        ftol = param.ftol,
        max_iters = param.max_iters,
        x0_init = param.x0_init,
        verbose = param.verbose,
    )
end

export solve, ThermoCycleParameters

function show(io::IO,params::ThermoCycleParameters)
    println(io, "ThermoCycleParameters:")
    println(io, "  N               = ", params.N)
    println(io, "  autodiff        = ", params.autodiff)
    println(io, "  internal_pinch  = ", params.internal_pinch)
    println(io, "  fd_order        = ", params.fd_order)
    println(io, "  xtol            = ", params.xtol)
    println(io, "  ftol            = ", params.ftol)
    println(io, "  restart_TOL     = ", params.restart_TOL)
    println(io, "  max_iters       = ", params.max_iters)
    println(io, "  x0_init         = ", params.x0_init)
    println(io, "  verbose         = ", params.verbose)
end

function optimize end