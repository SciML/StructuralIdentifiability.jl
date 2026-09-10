@safetestset "Assessing identifiability" begin
    # Nemo fpMatrix Julia-owned row pointers raise InexactError on i686 (Nemo#2358).
    # Gate the include here: top-level `return` inside the body file does not stop
    # SafeTestsets evaluation of the rest of the included file.
    if Sys.WORD_SIZE != 64
        @info "Skipping identifiability.jl on $(Sys.WORD_SIZE)-bit (Nemo#2358)"
    else
        include(joinpath(@__DIR__, "bodies", "identifiability.jl"))
    end
end
