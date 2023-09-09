# shared utilities for benchmarking
using Printf

const BENCHMARK_RESULTS = "results"

const categories = [:implicit_relations, :dim_before, :dim_after]

const ID_TIME_CATEGORIES = [
    :id_total,
    #:id_io_time,
    #:id_global_time,
    #:id_inclusion_check,
    #:id_inclusion_check_mod_p,
    #:id_groebner_time,
    #:id_beautifulization,
    #:id_normalforms_time,
    #:id_gbfan_time,
    #:id_ranking,
]
const ID_DATA_CATEGORIES = []
const ALL_CATEGORIES = union(ID_TIME_CATEGORIES, ID_DATA_CATEGORIES)
const ALL_POSSIBLE_CATEGORIES = union(ALL_CATEGORIES, Symbol[])

const HUMAN_READABLE_CATEGORIES = Dict(
    :id_io_time => "io",
    :id_primality_evaluate => "io/primality-evaluate",
    :id_uncertain_factorization => "io/uncertain-factor",
    :id_global_time => "global id.",
    :id_ideal_time => "gen. ideal",
    :id_inclusion_check => "inclusion",
    :id_inclusion_check_mod_p => "inclusion Zp",
    :id_groebner_time => "ParamPunPam.jl",
    :id_total => "Runtime",
    :id_beautifulization => "beautifulization",
    :id_normalforms_time => "normal forms",
    :id_gbfan_time => "GB fan",
    :id_ranking => "Score",
    :implicit_relations => "Algebraic relations",
    :dim_before => "Dim. before",
    :dim_after => "Dim. after",
)

const CATEGORY_FORMAT = Dict()
for cat in ALL_POSSIBLE_CATEGORIES
    CATEGORY_FORMAT[cat] = (val) -> if val isa Real
        @sprintf("%.2f", val)
    else
        string(val)
    end
end

function parse_keywords(keywords)
    if isempty(keywords)
        return [[]]
    end
    if isempty(eval(Meta.parse(keywords)))
        return [[]]
    end
    sets_of_keywords = map(strip, split(keywords, ";"))
    @assert !isempty(sets_of_keywords)
    kws_named_tuples = []
    for kwset in sets_of_keywords
        @debug "" kwset
        nt = eval(Meta.parse(kwset))
        push!(kws_named_tuples, nt)
    end
    kws_named_tuples
end

function keywords_to_global_id(keywords)
    if isempty(keywords)
        return Symbol()
    end
    id = get(keywords, :strategy, :default)
    if haskey(keywords, :with_states)
        if keywords.with_states
            id = Symbol(id, :_with_states)
        end
    end
    if haskey(keywords, :adjoin_identifiable)
        if keywords.adjoin_identifiable
            id = Symbol(id, :_adjoin_identifiable)
        end
    end
    id
end

function timings_filename(kwid)
    generic_filename("timings", kwid)
end

function result_filename(kwid)
    generic_filename("result", kwid)
end

function data_filename(kwid)
    generic_filename("data", kwid)
end

function generic_filename(name, kwid)
    if kwid === Symbol("")
        "$name"
    else
        "$name_$kwid"
    end
end
