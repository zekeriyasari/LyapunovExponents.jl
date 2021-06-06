using LyapunovExponents
using Documenter

DocMeta.setdocmeta!(LyapunovExponents, :DocTestSetup, :(using LyapunovExponents); recursive=true)

makedocs(;
    modules=[LyapunovExponents],
    authors="Zekeriya Sarı",
    repo="https://github.com/zekeriyasari/LyapunovExponents.jl/blob/{commit}{path}#{line}",
    sitename="LyapunovExponents.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://zekeriyasari.github.io/LyapunovExponents.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/zekeriyasari/LyapunovExponents.jl",
)
