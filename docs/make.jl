using Documenter, StructuralIdentifiability

makedocs(
    sitename="StructuralIdentifiability.jl",
    authors="SciML",
    modules=[StructuralIdentifiability],
    clean=true,doctest=false,
    format = Documenter.HTML(#analytics = "UA-90474609-3",
                             assets = ["assets/favicon.ico"],
                             canonical="https://si.sciml.ai/stable/"),
    pages=[
        "Home" => "index.md",
        "Tutorials" => Any[
        ],
        "Basics" => Any[
        ],
    ]
)

deploydocs(
   repo = "github.com/SciML/StructuralIdentifiability.jl.git";
   push_preview = true
)