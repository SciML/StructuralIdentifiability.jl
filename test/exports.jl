@testset "Exporting to other formats" begin
    ode = @ODEmodel(x'(t) = a * x(t) + b * u(t)^2, y(t) = 1 // x(t))
    @test print_for_maple(ode) ==
          "read \"../IdentifiabilityODE.mpl\";\n\nsys := [\ndiff(x(t), t) = x(t)*a + u(t)^2*b,\ny(t) = (1) / (x(t))\n];\nCodeTools[CPUTime](IdentifiabilityODE(sys, GetParameters(sys)));"
    @test print_for_maple(ode, :DifferentialAlgebra) ==
          "with(DifferentialAlgebra):\nring_diff := DifferentialRing(blocks = [[x], [y] , [u]], derivations = [t]):\nsys := [\ndiff(x(t), t) - (x(t)*a + u(t)^2*b),\ny(t) - ((1) / (x(t)))\n];\nres := CodeTools[CPUTime](RosenfeldGroebner(sys, ring_diff, singsol=none));"
    @test print_for_maple(ode, :DifferentialThomas) ==
          "with(DifferentialThomas):\nwith(Tools):\nRanking([t], [[x], [y] , [u]]):\nsys := [\ndiff(x(t), t) - (x(t)*a + u(t)^2*b),\ny(t) - ((1) / (x(t)))\n];\nres := CodeTools[CPUTime](ThomasDecomposition(sys));"
    @test print_for_DAISY(ode) ==
          "B_:={u, y, x}\$\nFOR EACH EL_ IN B_ DO DEPEND EL_,T\$\n\nB1_:={a, b}\$\n %NUMBER OF STATES\nNX_:=1\$\n %NUMBER OF INPUTS\nNU_:=1\$\n %NUMBER OF OUTPUTS\nNY_:=1\$\n\nC_:={df(x, t) = x*a + u^2*b,\ny = (1) / (x)}\$\nFLAG_:=1\$\nSHOWTIME\$\nDAISY()\$\nSHOWTIME\$\nEND\$\n"
    @test print_for_GenSSI(ode) ==
          "function model = SMTH()\nsyms x\nsyms a b\nsyms x0\nsyms u\nmodel.sym.p = [a; b; x0];\nmodel.sym.x = [x];\nmodel.sym.g = [u];\nmodel.sym.x0 = [x0];\nmodel.sym.xdot = [x*a + u^2*b];\nmodel.sym.y = [(1) / (x)];\nend"
    @test print_for_COMBOS(ode) == "dx1/dt = x1*a + u1^2*b;\ny1 = (1) / (x1)"
end
