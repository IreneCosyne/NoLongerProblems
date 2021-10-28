
function nquantile(array, groups)
    return quantile(array, 1/groups:1/groups:1)
end

function type_array(TYPE, array)
    dimensions = size(array)

    if TYPE == "Float"
        new_array = Array{Float64}(dimensions)
        for i in eachindex(array)
            new_array[i] = convert(Float, array[i])
        end

        elseif TYPE == "Int"
        new_array = Array{Int}(dimensions)
        for i in eachindex(array)
            new_array[i] = convert(Int, array[i])
        end

    else
        println("Type not included yet in type_array function")
    end
    return new_array
end

function exclude_rows(list_to_exclude::Array, data::Array; return_excluded = false)

    # Obtain size of the new array to create
    number_of_rows_to_exclude = length(list_to_exclude)
    old_dimensions = size(data)
    new_dimensions = (old_dimensions[1]-number_of_rows_to_exclude, old_dimensions[2])

    # Create an empty new array with the new dimensions
    new_array = Array{Any}(new_dimensions)

    # Copy the non excluded rows to the new array
    count = 0
    for i in 1:old_dimensions[1]
        if in(i, list_to_exclude) != true
            count +=1
            new_array[count, :] = data[i, :]
        end
    end

    #This bit is in case we would like to have separately the excluded rows in other array
    if return_excluded == true
        excluded_array_dimensions = (number_of_rows_to_exclude + 1, old_dimensions[2])
        excluded_array = Array{Any}(excluded_array_dimensions)

        count = 1
        excluded_array[count, :] = data[1, :] # To copy the first line which contains which information is contained in tha column

        for i in 1:old_dimensions[1]
            if in(i, list_to_exclude) == true
                count +=1
                excluded_array[count, :] = data[i, :]
            end
        end
        return new_array, excluded_array
    else
        return new_array
    end
end

function exclude_columns(list_to_exclude::Array, data::Array; return_excluded = false)
    # This function probably rewuires an update because it assumes that the
    # first line of the array shoud be conserved when you are excluding rows
    # Obtain size of the new array to create
    number_of_columns_to_exclude = length(list_to_exclude)
    old_dimensions = size(data)
    new_dimensions = (old_dimensions[1], old_dimensions[2] - number_of_columns_to_exclude)

    # Create an empty new array with the new dimensions
    new_array = Array{Any}(new_dimensions)

    # Copy the non excluded columns to the new array
    count = 0
    for i in 1:old_dimensions[2]
        if in(i, list_to_exclude) != true
            count +=1
            new_array[!, count] = data[!, i]
        end
    end

    # This bit is in case we would like to have separately the excluded columns in other array
    # I do not think it is goint to be useful for an array[cells, features]  but one never knows
    if return_excluded == true
        excluded_array_dimensions = (old_dimensions[1], number_of_columns_to_exclude)
        excluded_array = Array{Any}(excluded_array_dimensions)

        count = 0 # here we are not copying the first column because it does not contain

        for i in 1:old_dimensions[2]
            if in(i, list_to_exclude) == true
               count +=1
                excluded_array[!, count] = data[!, i]
            end
        end
        return new_array, excluded_array
    else
        return new_array
    end
end

function include_rows(list_to_exclude::Array, data::Array)
    return exclude_rows(list_to_exclude, data; return_excluded = true)[2]
end

function include_columns(list_to_exclude::Array, data::Array)
    return exclude_columns(list_to_exclude, data; return_excluded = true)[2]
end


############ New NoLongerProblems functions (25/01/17)


function find_out(collectionA, collectionB)
    all_index = 1:length(collectionA)
    index_common = findin(collectionA, collectionB)

    index_out = []

    for index in all_index
        if in(index, index_common) == false
            push!(index_out, index)
        end
    end

    return index_out

end

function common_noncommon_indexes(collectionA, collectionB)
    index_common = findin(collectionA, collectionB)
    index_noncommon = find_out(collectionA, collectionB)
    return index_common, index_noncommon

end



function common_noncommon(collectionA, collectionB)
    index_com_noncom = common_noncommon_indexes(collectionA, collectionB)
    things_in_common = collectionA[index_com_noncom[1]]
    things_notin_common = collectionA[index_com_noncom[2]]
    return things_in_common, things_notin_common

end

export reverse_boolean_vector, exclude_rows_with_NA
function reverse_boolean_vector(boolean_vector)

    new_boolean_vector = Array{Bool, 1}()


    for x in boolean_vector
        if x == true
            push!(new_boolean_vector, false)
        else
            push!(new_boolean_vector, true)
        end
    end
    return new_boolean_vector
end


function exclude_rows_with_NA(dataframe, colname)
    whereareNA = Array(map(isna, eachcol(dataframe)))
    colnames = names(dataframe)
    n_col = 0
    for col in colnames
        n_col += 1
        if col == colname
            break
        end
    end
    booleanvector = reverse_boolean_vector(whereareNA[!, n_col])
    new_dataframe = dataframe[booleanvector, :]
end


function substitute_NA_for_this(df::DataFrame, column::Symbol, this)
  data = df[!,column]
    for i in eachindex(data)
        try
            if isna(data[i])
                data[i] = this
            end
        catch 
            if !ismissing(data[i])
                data[i] = this
            end
        end
    end
    df[!,column] = data
    return df
end

function substitute_NA_for_1(df::DataFrame, column::Symbol)
  substitute_NA_for_this(df::DataFrame, column::Symbol, 1)
end

function substitute_NA_for_0(df::DataFrame, column::Symbol)
  substitute_NA_for_this(df::DataFrame, column::Symbol, 0)
end

function substitute_that_for_this(df::DataFrame, column::Symbol, that, this)
  data = df[column]
    for i in eachindex(data)
        if data[i] == that
            data[i] = this
        end
    end
    df[column] = data
    return df
end



export get_values
function get_values(dict)
    values = []
    for key in keys(dict)
        push!(values, dict[key])
    end
    return values
end

export join_in_all_common_columns, common_keys, concat_values

function join_in_all_common_columns(dataframe1, dataframe2)
    cols1 = names(dataframe1)
    cols2 = names(dataframe2)
    common_names = findin(cols1, cols2)
    columns_to_join = cols1[common_names]
    new_dataframe = DataFrame()
    for common_name in columns_to_join
        new_dataframe[!,common_name] = vcat(dataframe1[!,common_name], dataframe2[!,common_name])
    end
    return new_dataframe
end

function join_in_all_common_columns(dfs::Array{DataFrames.DataFrame,1})
    n_df = length(dfs)
    df = dfs[1]
    for i in 2:n_df
        df = NoLongerProblems.join_in_all_common_columns(df, dfs[i])
    end
    return df
end

function join_in_all_common_columns(dfs...)
    df = dfs[1]
    
    for i in 2:length(dfs)
        df = join_in_all_common_columns(df, dfs[i])
    end
    
    return df
end

function findin(a, b)
    findall((in)(b), a)
end


function common_keys(dict1, dict2)
    keys1 = collect(keys(dict1))
    keys2 = collect(keys(dict2))
    indexes = findin(keys1, keys2)
    common_keys = keys1[indexes]
    return common_keys
end

function concat_values(dict::Dict)
    values = Array{Real, 1}()
    for key in keys(dict)
        values = vcat(values, dict[key])
    end
    return values
end

function Base.Dict(array1::Vector, array2::Vector)
    new_dict = Dict()
    if length(array1) == length(array2)
        for i in 1:length(array1)
            new_dict[array1[i]] = array2[i]
        end
    else
        print("ERROR, vector dimensions do not match")
    end
    return new_dict
end

function select_chr(dataframe, chromosome::String)
    f(x) = x == chromosome
    new_dataframe = dataframe[[f(i) for i in dataframe[!,:chr]], :]
    #println(size(new_dataframe)[1], " genes in chromosome $chromosome")
    return new_dataframe
end

function select_chr_and_order_by_tss(dataframe, chromosome::String)
    new_dataframe = select_chr(dataframe, chromosome)
    new_dataframe = sort(new_dataframe, cols = :tss, rev = true)
end
export select_chr, select_chr_and_order_by_tss, dict_chr

function dict_chr(data; order_by_TSS = false)
    new_dict = Dict()
    chromosomes = unique(data[!,:chr])
    for chrom in chromosomes
      if order_by_TSS == true
        new_dict[chrom] = select_chr_and_order_by_tss(data, chrom)
      else
        new_dict[chrom] = select_chr(data, chrom)
      end
    end
    return new_dict
end
export tss_tss_distance_bp

function tss_tss_distance_bp(genes_in_chromosome::DataFrame)
    count_genes = size(genes_in_chromosome)[1]
    distance_array = zeros(count_genes, count_genes)
    for gene1 in 1:count_genes
        gene1_dimension_values = genes_in_chromosome[gene1, :tss]
        for gene2 in 1:count_genes
            gene2_dimension_values = genes_in_chromosome[gene2, :tss]
            value = gene1_dimension_values - gene2_dimension_values            # Calculate euclidean distance between 2 nodes
            if value >= 0
              distance_array[gene1, gene2] = value    # Add the value into the array
            else
              distance_array[gene1, gene2] = value*(-1)
            end
        end
    end
    return distance_array
end



function select_columns_that_contain(dataframe::DataFrame, thing::String; extra_column = :NameofColumn)
    cols = names(dataframe)
    #cols_int = Array{Symbol, 1}()
    #for col in cols
    #    if contains(string(col), thing) == true
   
    cols_int = [] # 6th June 2019
      for col in cols
        if occursin(thing,string(col)) == true
            
            push!(cols_int, col)
            elseif col  == extra_column
            push!(cols_int, col)
        end
    end
    return dataframe[!, cols_int]
end



function obtain_keys_that_contain(dict::Dict, thing::String)
    cols = keys(dict)
    cols_int = Array{Symbol, 1}()
    for col in cols
        if contains(string(col), thing) == true
            push!(cols_int, col)
        end
    end
    return cols_int
end

function obtain_keys_that_are_dataframes(dict::Dict)
    cols = keys(dict)
    cols_int = Array{Symbol, 1}()
    for col in cols
        if typeof(dict[col]) == DataFrame
            push!(cols_int, col)
        end
    end
    return cols_int
end



function return_index_which_array_contains(listofarrays, containingarray)
    
    matches = []
    for i in eachindex(listofarrays)
        
        if issubset(containingarray, listofarrays[i])
            push!(matches, i)
        end
        
    end
    return matches
end



function split_by(df, thing)

    samples = unique(df[!,thing])
    dfs = Dict()
    for sample in samples
        f(x) = x == sample
         d = df[[f(x) for x in df[!,thing]],:]
        dfs[sample] = d
    end
    return dfs

end

function columns_containing(df, containing)
    list = names(df)
    bool = [occursin(uppercase(containing), uppercase(string(i))) for i in list]
    list[bool]
end


function rename_content_column(df::DataFrame, column::Symbol, replacement)
    n = nrow(df)
    df[column] = fill(replacement, n)
    df
end



function rename_column_and_merge_dfs(column::Symbol, dfs::Array{DataFrame, 1}, replacements::Array)
    n_df = length(dfs)
    n_rep = length(replacements)
    
    if n_df == n_rep
        for i in 1:n_df
            dfs[i] = rename_content_column(dfs[i], column, replacements[i])
        end
    else
        return error
    end
    
    if n_df == 1
        return dfs[1]
    else
        return dfs = NoLongerProblems.join_in_all_common_columns(dfs)
    end
    
    
end



function equals(dataarray, thingequal)
    f(x) = x == thingequal
    return [try f(i) catch; false end for i in collect(dataarray)]
end

function equals(df, column, value)
    f(x) = [i == value for i in x]
    df[f(df[column]), :]
end

function morethan(dataarray, thingequal)
    f(x) = x > thingequal
    return [f(i) for i in collect(dataarray)]
end

function morethan(df, column, value)
    f(x) = [i > value for i in x]
    df[f(df[!,column]), :]
end

function lessthan(dataarray, thingequal)
    f(x) = x > thingequal
    return [f(i) for i in collect(dataarray)]
end

function lessthan(df, column, value)
    f(x) = [i .< value for i in x]
    df[f(df[column]), :]
end

function Base.Dict(array1, array2)
    new_dict = Dict()
    if length(array1) == length(array2)
        for i in 1:length(array1)
            new_dict[array1[i]] = array2[i]
        end
    end
    return new_dict
end



function z_to_pvalue_2sided(x)
    dist = Normal(0, 1)
    pvalue = 2*cdf(dist, -abs(x))
end
                                            


function summary_factor(df; col = :ResultDiffDisp)
    df_grouped = split_by(df, col)
    
    new_df = DataFrame()
    new_df[!,:Factor] = collect(keys(df_grouped))
    new_df[!,:N] =  [nrow(df_grouped[i]) for i in new_df[!,:Factor] ]
    new_df
    
end

function remove_column(df::DataFrame, colsymbols::Symbol...)
        cols = names(df)
        futurecols =  filter!(e -> !in(e, colsymbols), cols)
        return df = df[!,[futurecols...]]
end

function drop_narows(df, col)
    df[.!isna.(df[col]), :]
end

function count_cells_per_sample(a_df_with_Sample_column)
    samples = unique(a_df_with_Sample_column[!,:Sample])
    df = DataFrames.DataFrame()
    df[!,:Samples] = samples
    df[!,:N_cells] = [count(x -> x !=0, a_df_with_Sample_column[!,:Sample] .== s) for s in samples]
    df

end

function transform_pvalue_in_stars(pval)
    if pval <= 0.0001; return "****"
        elseif pval <= 0.001; return "***"
        elseif pval <= 0.01; return "**"
        elseif pval <= 0.05; return "*"
        else; return "ns"
    end
        
end

function factortable(df, fcolumn1, fcolumn2) 
    fcolumn1 = Symbol(fcolumn1)
    fcolumn2 = Symbol(fcolumn2)
    
    f1 = ["$ii" for ii in unique(df[fcolumn1])]
    f2 = ["$ii" for ii in unique(df[fcolumn2])]
    
    new_df = DataFrames.DataFrame()
    new_df[!,:Feature] = f1
    
    for f_2 in unique(df[fcolumn2])
        array = []
        array_f = []
        bool2 = df[fcolumn2] .== f_2
        for f_1 in unique(df[fcolumn1])
            bool1 = df[fcolumn1] .== f_1
            n = count(x -> x != 0, bool1)
            bool = bool2 .* bool1
            push!(array, count(x -> x != 0, bool) )
            push!(array_f, count(x -> x != 0, bool)/n)
        end
        new_df[Symbol(f_2)] =  array
        new_df[Symbol("$f_2"*"_Fraction")] =  array_f
    end
    
    new_df
        
end

function calculate_positives(df, column, by; limit = 0)
    samples = []
    for ii in 1:nrow(df)
        s = ""
        for jj in df[ii, by]
            s = s*jj*"_"
        end
        s = s[1:(length(s)-1)]; push!(samples, s)
    end
    df[!,:Sample] = samples
    samplesunique = unique(samples)
    sp_df = split_by(df, :Sample)
    newdf = DataFrames.DataFrame(Sample=samplesunique)
    for ii in 1:length(by)
        newdf[!,by[ii]] = [split(jj, "_")[ii] for jj in newdf[!,:Sample]]
    end
    newdf[!, :N] = [sum(df[!,:Sample].==ii) for ii in samplesunique]
    newdf[!, :N_positive] = [sum(sp_df[ii][!,y].>limit) for ii in samplesunique]
    newdf[!, :F_positive] = newdf[!, :N_positive] ./ newdf[!, :N]
    newdf
end


