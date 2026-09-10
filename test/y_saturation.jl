@safetestset "Output saturation" begin
    # Nemo fpMatrix Julia-owned row pointers raise InexactError on i686 (Nemo#2358).
    if Sys.WORD_SIZE != 64
        @info "Skipping y_saturation.jl on $(Sys.WORD_SIZE)-bit (Nemo#2358)"
    else
        include(joinpath(@__DIR__, "bodies", "y_saturation.jl"))
    end
end
