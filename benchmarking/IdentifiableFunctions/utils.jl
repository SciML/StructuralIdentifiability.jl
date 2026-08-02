# shared utilities for benchmarking

const BENCHMARK_RESULTS = "results"
const BENCHMARK_LOGS = "logs"

const ID_TIME_CATEGORIES = [
    :id_total,
    :id_io_time,
    :id_global_time,
    :id_inclusion_check,
    :id_inclusion_check_mod_p,
    :id_groebner_time,
    :id_beautifulization,
    :id_normalforms_time,
    :id_gbfan_time,
    :id_ranking,
]
const ID_DATA_CATEGORIES =
    [:id_npoints_degree, :id_npoints_interpolation, :are_id_funcs_polynomial]
const ALL_CATEGORIES = union(ID_TIME_CATEGORIES, ID_DATA_CATEGORIES)

function parse_keywords(input)
    isempty(strip(input)) && return [NamedTuple()]

    return map(split(input, ";")) do keyword_set
        expression = Meta.parse(strip(keyword_set))
        keywords = Core.eval(@__MODULE__, expression)
        keywords == () && return NamedTuple()
        keywords isa NamedTuple ||
            throw(ArgumentError("Expected a named tuple of keyword arguments, got $keywords"))
        return keywords
    end
end

function benchmark_comparator(label, ode)
    if label in (:default, :cmp_default)
        return StructuralIdentifiability.default_cmp(ode)
    elseif label === :cmp_lie
        return StructuralIdentifiability.cmp_lie(
            ode,
            StructuralIdentifiability.RationalFunctionFields.rational_function_cmp,
        )
    end
    throw(ArgumentError("Unknown comparator label: $label"))
end

function benchmark_keyword_label(value)
    value isa Symbol && return value
    label = value isa Function ? typeof(value) : value
    return Symbol(replace(string(label), r"[^A-Za-z0-9_]+" => "_"))
end

function keywords_to_global_id(keywords)
    isempty(keywords) && return Symbol()

    id = get(keywords, :strategy, Symbol())
    if get(keywords, :with_states, false)
        id = isempty(string(id)) ? :with_states : Symbol(id, :_with_states)
    end
    if haskey(keywords, :rational_interpolator)
        id = Symbol(id, :_, keywords.rational_interpolator)
    end
    if haskey(keywords, :cmp)
        comparator = benchmark_keyword_label(keywords.cmp)
        id = isempty(string(id)) ? comparator : Symbol(id, :_, comparator)
    end
    if get(keywords, :adjoin_identifiable, false)
        id = isempty(string(id)) ? :adjoin_identifiable : Symbol(id, :_adjoin_identifiable)
    end
    return id
end

timings_filename(id) = generic_filename("timings", id)
result_filename(id) = generic_filename("result", id)
data_filename(id) = generic_filename("data", id)

function generic_filename(name, id)
    filename = isempty(string(id)) ? name : "$(name)_$id"
    return replace(filename, ":" => "")
end
