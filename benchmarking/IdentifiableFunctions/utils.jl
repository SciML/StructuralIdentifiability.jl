# shared utilities for benchmarking

ID_TIME_CATEGORIES = [
    :id_io_time,
    :id_global_time,
    :id_inclusion_check,
    :id_inclusion_check_mod_p,
    :id_groebner_time,
    :id_beautifulization,
    :id_normalforms_time,
    :id_gbfan_time,
    :id_total,
]
ID_DATA_CATEGORIES = []
ALL_CATEGORIES = union(ID_TIME_CATEGORIES, ID_DATA_CATEGORIES)

HUMAN_READABLE_CATEGORIES = Dict(
    :id_io_time => "io",
    :id_primality_evaluate => "io/primality-evaluate",
    :id_uncertain_factorization => "io/uncertain-factor",
    :id_global_time => "global id.",
    :id_ideal_time => "gen. ideal",
    :id_inclusion_check => "inclusion",
    :id_inclusion_check_mod_p => "inclusion Zp",
    :id_groebner_time => "ParamPunPam.jl",
    :id_total => "total",
    :id_beautifulization => "beautifulization",
    :id_normalforms_time => "normal forms",
    :id_gbfan_time => "GB fan",
)

function parse_keywords(keywords)
    sets_of_keywords = map(strip, split(keywords, ";"))
    @assert !isempty(sets_of_keywords)
    kws_named_tuples = []
    for kwset in sets_of_keywords
        @info "" kwset
        nt = eval(Meta.parse(kwset))
        push!(kws_named_tuples, nt)
    end
    kws_named_tuples
end

function keywords_to_id(keywords)
    get(keywords, :strategy, :default)
end
