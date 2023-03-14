const DO = joinpath(homedir(), "Downloads")

const PK = dirname(@__DIR__)

const IN = joinpath(PK, "input")

const OU = joinpath(PK, "output")

# ----------------------------------------------------------------------------------------------- #

using DataFrames

using Test

using BioLab

# --------------------------------------------- #

function remove_all_equal_column(da)

    println(size(da))

    da = da[!, [!allequal(co) for co in eachcol(da)]]

    println(size(da))

    return da

end

# --------------------------------------------- #

function read(di, na)

    println("📖 $na")

    return remove_all_equal_column(BioLab.Table.read(joinpath(di, na)))

end
