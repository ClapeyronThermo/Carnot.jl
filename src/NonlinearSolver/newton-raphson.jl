
function constrained_newton_verbose_step(xn,fxn,iter,lenx,lenf)
    @info "Iteration: $iter, x: $xn, f(x): $fxn, lenx: $lenx, lenf: $lenf"
end

"""
`SolutionState` -  A struct to hold the solution state of the nonlinear solver.

- `x::Vector{T}`: The solution vector.
- `f_calls::I`: The number of function calls made during the solving process.
- `iterations::I`: The number of iterations taken to converge.
- `residuals::Vector{T}`: The residuals at the solution.
- `lb::Vector{T}`: The lower bounds used in the solver.
- `ub::Vector{T}`: The upper bounds used in the solver.
- `autodiff::Bool`: A flag indicating whether automatic differentiation was used.
- `fd_order::I`: The order of finite difference used if autodiff is false.
- `lenx::T`: The final change in the solution vector.
- `lenf::T`: The final change in the residuals.
- `soltype::Symbol`: A symbol indicating the type of cycle (:unknown,:subcritical,:transcritical). This will be updated with the cycle type after solving.
"""
mutable struct SolutionState{T<:Real,I<:Integer}
    x::Vector{T}
    f_calls::I
    iterations::I
    residuals::Vector{T}
    lb::Vector{T}
    ub::Vector{T}
    autodiff::Bool
    fd_order::I
    lenx::T
    lenf::T
    soltype::Symbol
    function SolutionState(x, f_calls, iterations, residuals,
                           lb, ub, autodiff, fd_order, lenx, lenf, soltype)
        T = promote_type(eltype(x), eltype(residuals), eltype(lb), eltype(ub),
                         typeof(lenx), typeof(lenf))
        I = promote_type(typeof(f_calls), typeof(iterations), typeof(fd_order))
        new{T,I}(Vector{T}(x), convert(I, f_calls), convert(I, iterations),
                  Vector{T}(residuals), Vector{T}(lb), Vector{T}(ub),
                  autodiff, convert(I, fd_order),
                  convert(T, lenx), convert(T, lenf), soltype)
    end
end

function show_parameters(sol::SolutionState)
    println("Iterations: ", sol.iterations)
    println("Function calls: ", sol.f_calls)
    println("Final residuals: ", sol.residuals)
    println("Final residual norm: ", norm(sol.residuals))
    println("Final x: ", sol.x)
    println("Final lenx ||xn - xk||/||xk||: ", sol.lenx)
    println("Final lenf ||f(xn)||: ", sol.lenf)
    println("Lower bounds: ", sol.lb)
    println("Upper bounds: ", sol.ub)
    println("Autodiff: ", sol.autodiff)
    if !sol.autodiff
        println("Finite difference order: ", sol.fd_order)
    end
    println("Solution type: ", sol.soltype)
    return nothing
end

function NaNSolution(lb,ub,autodiff,fd_order)
    return SolutionState([NaN,NaN],NaN,NaN,[NaN,NaN],lb,ub,autodiff,fd_order,NaN,NaN,:NaNSolution)
end

function box_projection(x::Array{T,1},lb::Array{TT,1},ub::Array{TTT,1}) where {T <: Real, TT <: Real, TTT <: Real}
    y = copy(x)
    for i in eachindex(x)
        if y[i] < lb[i]
            y[i] = lb[i]
        elseif y[i] > ub[i]
            y[i] = ub[i]
        end
    end
    return y
end

function box_projection!(x::Array{T,1},lb::Array{TT,1},ub::Array{TTT,1}) where {T <: Real, TT <: Real, TTT <: Real}
    for i in eachindex(x)
        if x[i] < lb[i]
            x[i] = lb[i]
        elseif x[i] > ub[i]
            x[i] = ub[i]
        end
    end
end


function constrained_newton_fd(f::Function,x::Array{T,1},
    lb::Array{TT,1},ub::Array{TTT,1};xtol::TOL=1e-16,ftol::TOL=1e-16,max_iters::S=100,
    fd_order::M = 2,verbose::Bool=false) where  {T <: Real, S <: Integer, TT <: Real, TTT <: Real, M <: Int,TOL<:Real}

    fd_method  = central_fdm(fd_order,1)

    n = length(x)
    xk = copy(x)
    xn = similar(x)
    jk = Array{T,2}(undef,n,n)

    lenx = one(T)
    lenf = one(T)

    iter_ = 0

    f_calls = 0
    if verbose
        @info "Starting constrained Newton-Raphson with finite difference Jacobian (order $fd_order)"
    end
    for iter in 1:max_iters
        jk .= first(FiniteDifferences.jacobian(fd_method,f,xk))
        f_calls += fd_order*n # approx number of function calls for finite difference jacobian
        if !all(isfinite,jk)
            error("The jacobian has non-finite elements")
        end
           
        # if rank(jk) == n # Compute a Newton step
            xn .= xk - jk\f(xk)
        # else # Compute a Levenberg-Marquardt step
        #     λ = 0.001*norm(f(xk))^2
         
        #     xn .= xk - (jk'jk + λ*I)\(jk'f(xk))
        # end
        f_calls += 1 # for f(xk) call above
        box_projection!(xn,lb,ub)

        lenx = norm((xn.-xk))/norm(xk)
        lenf = norm(f(xn))

        xk .= xn

        iter_ += 1
        if iter_ >= max_iters || (lenx <= xtol || lenf <= ftol)
            break
        end
        if verbose
            constrained_newton_verbose_step(xn, f(xn), iter_, lenx, lenf)
        end
    end
  
    return SolutionState(xn,f_calls,iter_,f(xn),lb,ub,false,fd_order,lenx,lenf,:unknown)
end

function constrained_newton_ad(f::Function,x::Array{T,1},
    lb::Array{TT,1},ub::Array{TTT,1};xtol::TOL=1e-16,ftol::TOL=1e-16,max_iters::S=100,verbose::Bool=false) where 
    {T <: Real, S <: Integer, TT <: Real, TTT <: Real,TOL<:Real}

    type_promoted = promote_type(eltype(x), eltype(lb), eltype(ub))
    xk = copy(x)
    n = length(x)
    xn = similar(x)
    jk = similar(x, type_promoted, n, n)
    # @show x, jn, 
    iter_ = 0
    f_calls = 0
    lenx = one(T)
    lenf = one(T)
    if verbose
        @info "Starting constrained Newton-Raphson with automatic differentiation"
    end
    for iter in 1:max_iters
        jk = ForwardDiff.jacobian(f, xk)
        f_calls += 1
        fx = f(xk)
        xn = xk - jk \ fx
        f_calls += 1

        # simple box projection (not differentiable, but works if outside AD path)
        for i in eachindex(xn)
            xn[i] = clamp(xn[i], lb[i], ub[i])
        end
        lenx = norm((xn.-xk))/norm(xk)
        lenf = norm(f(xn))
        iter_ += 1
        if iter >= max_iters || (lenx <= xtol || lenf <= ftol)
            break
        end

        xk = xn
        if verbose
            constrained_newton_verbose_step(xn, fx, iter_, length(x), length(fx))
        end
    end

    return SolutionState(xn,f_calls,iter_,f(xn),lb,ub,true,0,lenx,lenf,:unknown)
end




export SolutionState