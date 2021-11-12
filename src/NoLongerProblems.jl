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
export reverse_boolean_vector, exclude_rows_with_NA
export get_values
export join_in_all_common_columns, common_keys, concat_values
export tss_tss_distance_bp
export mean, select_string_columns
export split_by_range,split_by_equally_sized
export subset_by_padj_when_NA

include("Functions.jl")
include("Functions2.jl")

end # module
