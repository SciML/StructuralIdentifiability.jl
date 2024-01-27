using StructuralIdentifiability, Aqua
@testset "Aqua" begin
     Aqua.find_persistent_tasks_deps(StructuralIdentifiability)
     Aqua.test_ambiguities(StructuralIdentifiability, recursive = false)
     Aqua.test_deps_compat(StructuralIdentifiability)
     Aqua.test_piracies(StructuralIdentifiability,
         treat_as_own = [])
     Aqua.test_project_extras(StructuralIdentifiability)
     Aqua.test_stale_deps(StructuralIdentifiability)
     Aqua.test_unbound_args(StructuralIdentifiability)
     Aqua.test_undefined_exports(StructuralIdentifiability)
end
