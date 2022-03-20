using Documenter, PowerSystems, DocStringExtensions, PowerFlows

makedocs(
    modules=[PowerFlows],
    format=Documenter.HTML(mathengine=Documenter.MathJax()),
    sitename="PowerFlows.jl",
    pages=Any[
        "Welcome Page" => "index.md",
        "Code Base Developer Guide" => Any["Developer Guide" => "code_base_developer_guide/developer.md",],
        "Public API Reference" => "api/public.md",
        "Internal API Reference" => "api/internal.md",
    ],
)

deploydocs(
    repo="github.com/NREL-SIIP/PowerFlows.jl.git",
    target="build",
    branch="gh-pages",
    devbranch="main",
    devurl="dev",
    push_preview=true,
    versions=["stable" => "v^", "v#.#"],
)
