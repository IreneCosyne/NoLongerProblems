module NoLongerProblems

using DataFrames
using LightGraphs
using Distributions
try using Statistics catch; end

export type_array
export exclude_rows, exclude_columns, include_rows, include_columns, nquantile, count_cells_per_sample
export z_to_pvalue_2sided
export find_out, common_noncommon_indexes, common_noncommon
export substitute_NA_for_1, substitute_NA_for_0, substitute_NA_for_this, substitute_that_for_this
export select_columns_that_contain
export obtain_keys_that_contain, obtain_keys_that_are_dataframes
export return_index_which_array_contains
export split_by, columns_containing
export equals, morethan, Dict, lessthan
export summary_factor

include("Functions.jl")

end # module