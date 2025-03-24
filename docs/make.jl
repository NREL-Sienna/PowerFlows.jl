using Documenter, PowerSystems, DocStringExtensions, PowerFlows, DataStructures

pages = OrderedDict(
    "Welcome Page" => "index.md",
    "How-to-Guides" => "how-tos/stub.md",
    "Tutorials" => "tutorials/stub.md",
    "Explanation" => "explanation/stub.md",
    "Reference" => Any[
        "Code Base Developer Guide" => "reference/developers/developer.md",
        "Public API Reference" => "reference/api/public.md",
        "Internal API Reference" => "reference/api/internal.md"],
)

makedocs(;
    modules = [PowerFlows],
    format = Documenter.HTML(; mathengine = Documenter.MathJax()),
    sitename = "PowerFlows.jl",
    pages = Any[p for p in pages],
)

deploydocs(;
    repo = "github.com/NREL-Sienna/PowerFlows.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "main",
    devurl = "dev",
    push_preview = true,
    versions = ["stable" => "v^", "v#.#"],
)
