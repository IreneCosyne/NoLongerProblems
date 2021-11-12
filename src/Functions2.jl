

function Statistics.mean(a::Array, b::Array)
    if length(a) == length(b)
        means = [(a[i] + b[i])/2 for i in eachindex(a)]
        return means
    end
end

function select_string_columns(df::DataFrame)
    cols = names(df)
    new_cols = Array{Symbol, 1}()
    for col in cols
        if typeof(df[col]) == DataArrays.DataArray{String,1}
            push!(new_cols, col)
        end
    end
    df[new_cols]
end

export string_to_int, string_to_float

function string_to_int(strvec)
    a = map(x->(v = tryparse(Int,x); isnull(v) ? 0.0 : get(v)),strvec)
    return a
end

function string_to_float(strvec)
    a = map(x->(v = tryparse(Float64,x); isnull(v) ? 0.0 : get(v)),strvec)
    return a
end

export sentencecase

function sentencecase(word::String)
    capital = uppercase(word[1])
    lower = lowercase(word[2:end])
    return string(capital, lower)
end


function log2Mean(df, col1, col2)
    # log(x+1) to avoid errors
    [log2((df[u, col1]+df[u, col2])/2 + 1) for u in 1:nrow(df)]
end

function Meanlog2(df, col1, col2)
    # log(x+1) to avoid errors
    [(log2(df[u, col1] + 1)+log2(df[u, col2] + 1))/2 for u in 1:nrow(df)]
end



function split_by_range(da::DataFrame, by, range)
    new_data = []
    range = collect(range)
    count = 0
    for i in eachindex(range)
        if i !=1
            d = da[da[by] .< range[i], :]
            d = d[d[by] .> range[i-1], :]
            push!(new_data, d)
        end
    end
    return new_data
end

function split_by_equally_sized(da::DataFrame, by, groups)
new_data = []
    count = 0
    
    quantiles = nquantile(da[by], groups)
    
    for i in eachindex(quantiles)
        if i !=1
            d = da[da[by] .<= quantiles[i], :]
            d = d[d[by] .> quantiles[i-1], :]
            push!(new_data, d)
            elseif i == 1
            d = da[da[by] .<= quantiles[i], :]
            push!(new_data, d)
        end
    end
    return new_data
end

function add_quantile(df, columnsymbol; quantiles = [0, 0.5, 0.75, 0.9])
    
    quantiles = sort!(quantiles, rev = false)
    
    quants = quantile(df[columnsymbol], quantiles)
    
    ## Delete quantiles == to 0
    delpos = []
    
    for a in 2:length(quants)
        if quants[a] == 0
            push!(delpos, a)
        end
    end
    
    deleteat!(quants, delpos)
    deleteat!(quantiles, delpos)
    
    ## End of Delete quantiles == to 0
    
    
    dictq = Dict(quants, [try string(quantiles[i]*100)*"-"*string(quantiles[i+1]*100) catch; string(quantiles[i]*100)*"-100.0" end for i in 1:length(quantiles)])
    
    quant_list = []

    
    for a in df[columnsymbol]
        d = a
        for quant in quants
            if a >= quant
                d = dictq[quant]
            end
        end
         push!(quant_list, d)
    end
    
    df[Symbol(string(columnsymbol)*"__Quantile")] = quant_list
    
    df
end

##### Improvements on subsetting by padj, now there in no need to construct a new_array
# function everytime we join a new dataframe

function subset_by_padj_when_NA(dataframe, padj; dataset = 1)
    if dataset == 1
      new_dataframe = exclude_rows_with_NA(dataframe, :padj)
      new_dataframe = new_dataframe[new_dataframe[:padj] .<= padj, :]
    else
        dataset -= 1
        namecolumn = string("padj_", string(dataset))
        namecolumn = convert(Symbol, namecolumn)
        new_dataframe = exclude_rows_with_NA(dataframe, :padj)
        new_dataframe = new_dataframe[new_dataframe[namecolumn] .<= padj, :]
    end

    return new_dataframe
end