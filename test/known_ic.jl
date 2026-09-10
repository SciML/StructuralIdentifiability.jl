@safetestset "Identifiable functions with known generic initial conditions" begin
    # Nemo fpMatrix Julia-owned row pointers raise InexactError on i686 (Nemo#2358).
    if Sys.WORD_SIZE != 64
        @info "Skipping known_ic.jl on $(Sys.WORD_SIZE)-bit (Nemo#2358)"
    else
        include(joinpath(@__DIR__, "bodies", "known_ic.jl"))
    end
end
