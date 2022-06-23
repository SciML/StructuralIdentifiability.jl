using Documenter, StructuralIdentifiability
using BenchmarkTools

include("pages.jl")

makedocs(; sitename = "StructuralIdentifiability.jl",
         authors = "SciML",
         modules = [StructuralIdentifiability],
         clean = true, doctest = false,
         format = Documenter.HTML(; analytics = "UA-90474609-3",
                                  assets = ["assets/favicon.ico"],
                                  canonical = "https://si.sciml.ai/stable/"),
         pages = pages)

deploydocs(; repo = "github.com/SciML/StructuralIdentifiability.jl.git",
           push_preview = true)
