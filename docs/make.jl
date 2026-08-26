using Documenter, DocumenterVitepress, Carnot

makedocs(
    sitename = "Carnot.jl",
    modules = [Carnot],
    doctest = true,
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/ClapeyronThermo/Carnot.jl",
    ),
    warnonly = Documenter.except(),
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Optimization" => "optimization.md",
        "References" => "reference.md",
    ],
)

DocumenterVitepress.deploydocs(
    repo = "github.com/ClapeyronThermo/Carnot.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true,
)
