const DO = joinpath(homedir(), "Downloads")

const PK = dirname(@__DIR__)

const IN = joinpath(PK, "input")

const OU = joinpath(PK, "output")

# ----------------------------------------------------------------------------------------------- #

using DataFrames

using Test

using BioLab

# --------------------------------------------- #

function read(di, na)

    return BioLab.Table.read(joinpath(di, na))

end
