using Documenter, PowerSystems, DocStringExtensions, PowerFlows, DataStructures
using DocumenterInterLinks

links = InterLinks(
    "DocumenterInterLinks" => "http://juliadocs.org/DocumenterInterLinks.jl/stable/",
    "PowerSystems" => "https://nrel-sienna.github.io/PowerSystems.jl/stable/",
    "PowerNetworkMatrices" => "https://nrel-sienna.github.io/PowerNetworkMatrices.jl/stable/",
)

pages = OrderedDict(
    "Welcome Page" => "index.md",
    "How-to-Guides" => "how-tos/stub.md",
    "Tutorials" => "tutorials/stub.md",
    "Explanation" => "explanation/stub.md",
    "Reference" => Any[
        "Code Base Developer Guide" => "reference/developers/developer.md",
        "LCC Model Implementation" => "reference/developers/lcc_model.md",
        "Running Tests" => "reference/developers/running_tests.md",
        "Public API Reference" => "reference/api/public.md",
        "Internal API Reference" => "reference/api/internal.md",
    ],
)

makedocs(;
    modules = [PowerFlows],
    format = Documenter.HTML(; mathengine = Documenter.MathJax()),
    plugins = [links],
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
