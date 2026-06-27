using SciMLTesting, StructuralIdentifiability, Test

run_qa(
    StructuralIdentifiability;
    explicit_imports = true,
    # Heavy `using AbstractAlgebra`/`using Nemo`/... bring ~70 implicit imports
    # (core computer-algebra operators used pervasively); making them all explicit
    # is a large, risky refactor tracked in
    # https://github.com/SciML/StructuralIdentifiability.jl/issues/527.
    ei_broken = (:no_implicit_imports,),
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
            ),
        ),
    ),
)
