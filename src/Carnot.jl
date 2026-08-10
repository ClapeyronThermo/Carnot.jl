module Carnot

using Clapeyron, Polynomials, Interpolations, FiniteDifferences
using ForwardDiff#, LinearAlgebra
using StaticArrays
using RecipesBase
using CommonSolve
using CommonSolve: solve
# import LinearAlgebra: norm, rank
import Base: show, length, copy, similar, promote_type

import Clapeyron: molecular_weight

abstract type ThermoCycleProblem end

#NonlinearSolver - NR
norm(x) = sqrt(sum(abs2,x))
include("NonlinearSolver/newton-raphson.jl")

# Thermo-fixes
include("thermoextensions/utils.jl")
include("thermoextensions/fix_instabilites.jl")

#Data map file
include("DataMap/DataMap.jl")

# Cycle Structs
include("CycleStructs/ORC.jl")
include("CycleStructs/HeatPump.jl")

#solve
include("Solve/solve.jl")



#Plotting
include("Plots/cycleplots.jl")



    function show(prob::ThermoCycleProblem)
        show_parameters(prob)
    end
    function show(sol::SolutionState)
        show_parameters(sol)
    end

    export show

    export get_states, ThermoCycleProblem
end
