using Documenter, Carnot

makedocs(
sitename="Carnot.jl",
modules = [Carnot],
doctest = true,
pages = [
    "Home" => "index.md",
    "Examples" => "examples.md",
    "Optimization" => "optimization.md",
    # "Cycle Optimization" => "Optimization.md",
    "References" => "reference.md"
]
)

deploydocs(repo = "github.com/ClapeyronThermo/Carnot.jl",push_preview = true)