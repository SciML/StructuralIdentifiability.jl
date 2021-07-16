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
            "input/input.md",
            "identifiability/identifiability.md",
        ],
        "Library" => Any[
            "Local Identifiability Tools" => "utils/local_identifiability.md",
            "Elimination"=>"utils/elimination.md",
            "ODE Tools" => "utils/ode.md",
            "Power Series Tools" => "utils/power_series_utils.md",
            "Primality Chekcs" => "utils/primality.md",
            "Wronskian Tools" => "utils/wronskian.md",
            "Input-Output Equation tools"=>"ioequations/ioequations.md",
            "Other Utilities"=>"utils/util.md"
        ],
        "Export" => Any[
            "export/export.md"
        ],
        
    ]
)

deploydocs(
   repo = "github.com/SciML/StructuralIdentifiability.jl.git";
   push_preview = true
)
