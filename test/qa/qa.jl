using SciMLTesting, StructuralIdentifiability, Test

# ExplicitImports only walks an extension that is actually loaded (it resolves each
# `[extensions]` entry through `Base.get_extension`, which is `nothing` otherwise), so
# every weakdep is loaded here to bring `ModelingToolkitSIExt` into the check.
using ModelingToolkitBase, SymbolicUtils, Symbolics

# ExplicitImports silently skips an extension that fails to load, so assert the
# extension module actually exists rather than trusting a green run_qa.
@testset "Extensions loaded" begin
    @test Base.get_extension(StructuralIdentifiability, :ModelingToolkitSIExt) !== nothing
end

run_qa(
    StructuralIdentifiability;
    explicit_imports = true,
    # Heavy `using AbstractAlgebra`/`using Nemo`/... bring ~70 implicit imports
    # (core computer-algebra operators used pervasively); making them all explicit
    # is a large, risky refactor tracked in
    # https://github.com/SciML/StructuralIdentifiability.jl/issues/527.
    ei_broken = (:no_implicit_imports,),
    api_docs_kwargs = (; rendered = true),
    ei_kwargs = (;
        all_qualified_accesses_are_public = (;
            ignore = (
                :CoreLogging,             # Base internal
                :filter,                  # Base.Iterators (non-public)
                :isidentifier,            # Base internal
                :FracFieldElem,           # AbstractAlgebra.Generic (non-public)
                :GF,                      # Nemo.Native (non-public)
                :check_algebraicity_modp, # RationalFunctionFields (non-public)
                :dennums_to_fractions,    # RationalFunctionFields (non-public)
                :rational_function_cmp,   # RationalFunctionFields (non-public)
                :enable_progressbar,      # ParamPunPam (non-public)
                :postwalk,                # MacroTools (non-public)
                :seed!,                   # Random (non-public)
                :hasshift,                # ModelingToolkitBase (non-public)
                # SymbolicUtils (non-public). `Symbolics.scalarize`/`Symbolics.shape`
                # are public but SymbolicUtils owns both bindings, so the
                # owner-correct spelling is necessarily the non-public one.
                :scalarize,
                :shape,
                # StructuralIdentifiability's own internals, reached from
                # ModelingToolkitSIExt. An extension is its own module root, so
                # ExplicitImports' internal-access allowance does not cover it.
                :DDS,
                :_assess_identifiability,
                :_assess_identifiability_kic,
                :_assess_local_identifiability,
                :_assess_local_identifiability_discrete_aux,
                :_find_identifiable_functions,
                :_find_identifiable_functions_kic,
                :x_vars,
            ),
        ),
        all_explicit_imports_are_public = (;
            ignore = (
                # all RationalFunctionFields names (non-public, no public API yet)
                :_reduce_mod_p,
                :check_constructive_field_membership,
                :eval_at_dict,
                :fractions_to_dennums,
                :gen_tag_names,
                :is_rational_func_const,
                :parent_ring_change,
                :select_pivots,
                :str_to_var,
                :total_degree_frac,
                :unpack_fraction,
                :var_to_str,
                # StructuralIdentifiability's own internals, imported by
                # ModelingToolkitSIExt (see the note in the qualified-access list).
                :_si_logger,
                :_to,
                :nonrational_error,
                :reset_timings,
                :restart_logging,
            ),
        ),
    ),
)
