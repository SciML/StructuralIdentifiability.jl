using Logging

# Adapted from https://github.com/JuliaLang/julia/blob/5023ee21d70b734edf206aab3cac7c202ee0235a/stdlib/REPL/test/TerminalMenus/runtests.jl#L7
function simulate_input(keys...; kwargs...)
    keydict = Dict(:up => "\e[A", :down => "\e[B", :enter => "\r", :newline => "\n")

    new_stdin = Base.BufferStream()
    for key in keys
        if isa(key, Symbol)
            write(new_stdin, keydict[key])
        else
            write(new_stdin, "$key")
        end
    end

    return new_stdin
end

@testset "Interactive reparametrizations" begin
    ode = @ODEmodel(x'(t) = x*u(t), y(t) = x*a)
    new_stdin = simulate_input(:down, :enter, "d", "X", :enter, "\n")
    res = reparametrize_interactive(
        ode,
        input = new_stdin,
        loglevel = Logging.Error,
        output = devnull,
    )
    @test res[1] isa ODE
    @test string.(res[1].x_vars) == ["X"]
    @test isempty(res[1].parameters)
    old_var_to_new = Dict(v => k for (k, v) in res[2])
    @test old_var_to_new[x * a // 1] == res[1].x_vars[1]

    ode = @ODEmodel(x'(t) = x + b^2 + 1, y(t) = x)
    new_stdin = simulate_input(
        :enter,
        "d",
        "b",
        :enter,
        :newline,
        "bb",
        :enter,
        :newline,
        "b^2 + 1",
        :enter,
        :newline,
        "B",
        :enter,
        :newline,
        :down,
        :enter,
        "d",
        :enter,
        :newline,
    )
    res = reparametrize_interactive(
        ode,
        input = new_stdin,
        loglevel = Logging.Error,
        output = devnull,
    )
    @test string.(res[1].x_vars) == ["X1"]
    @test string.(res[1].parameters) == ["B"]
end
