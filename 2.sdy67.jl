include("_.jl")

# --------------------------------------------- #

da = "SDY67"

ind = joinpath(IN, da)

ou = mkpath(joinpath(OU, da))

validate = (true, true)

# --------------------------------------------- #

nap = "Participant ID"

nat = "Study Time Collected"

# --------------------------------------------- #

pe_x_ex_x_an = read(ind, "gene_expression_files.tsv")

sapet_ = [
    sa => "$pe.$(convert(Int,ti))" for (pe, sa, ti) in
    eachrow(sort!(select!(pe_x_ex_x_an, nap, "Biosample Accession", nat), [nap, nat]))
]

pet_ = [sapet[2] for sapet in sapet_]

@test length(unique(sapet[1] for sapet in sapet_)) == length(unique(pet_))

# --------------------------------------------- #

pe_x_de_x_an = read(ind, "demographics.tsv")

pe_x_ig_x_an = read(ind, "neut_ab_titer.tsv")

rename!(pe_x_ig_x_an, "Value Preferred" => "Antibody Titer")

pe_x_ig_x_an = subset(groupby(pe_x_ig_x_an, nap), nat => ti_ -> ti_ .== maximum(ti_))

pe_x_io_x_an = outerjoin(pe_x_de_x_an, pe_x_ig_x_an; on = nap, validate)

pe_id = Dict(pe => id for (id, pe) in enumerate(pe_x_io_x_an[!, nap]))

# --------------------------------------------- #

io_x_sa_x_an = DataFrame("Information" => names(pe_x_io_x_an))

for pet in pet_

    pe = join(split(pet, '.')[1:2], '.')

    io_x_sa_x_an[!, pet] = collect(pe_x_io_x_an[pe_id[pe], :])

end

# --------------------------------------------- #

BioLab.Table.write(joinpath(ou, "information_x_sample_x_anything.tsv"), io_x_sa_x_an)

# --------------------------------------------- #

fe_x_sa_x_nu = read(ind, "SDY67_PBMC_HealthyAdults.tsv")

select!(rename!(fe_x_sa_x_nu, sapet_), 1, pet_...)

# --------------------------------------------- #

fe_x_sa_x_nu[!, 1] =
    BioLab.Gene.rename(fe_x_sa_x_nu[!, 1], BioLab.Gene.map_ensembl(), BioLab.Gene.map_hgnc())[1]

nar = "Gene"

rename!(fe_x_sa_x_nu, 1 => nar)

fe_x_sa_x_nuc = BioLab.DataFrame.collapse(fe_x_sa_x_nu)

# --------------------------------------------- #

ts = joinpath(ou, "gene_x_sample_x_numbercollapsed.tsv")

BioLab.Table.write(ts, fe_x_sa_x_nuc)

# --------------------------------------------- #

ti_ = [parse(Int, BioLab.String.split_and_get(pet, '.', 3)) for pet in pet_]

# --------------------------------------------- #

BioLab.Plot.plot_heat_map(
    fe_x_sa_x_nuc;
    nar,
    nac = "Sample",
    grc_ = ti_,
    layout = Dict("title" => Dict("text" => da)),
    ht = BioLab.Path.replace_extension(ts, "html"),
)

# --------------------------------------------- #

for ti in sort(unique(ti_))

    id_ = findall(ti2 == ti for ti2 in ti_)

    idd_ = vcat(1, [id + 1 for id in id_])

    io_x_sat_x_an = io_x_sa_x_an[:, idd_]

    fe_x_sat_x_nuc = fe_x_sa_x_nuc[:, idd_]

    BioLab.Table.write(joinpath(ou, "information_x_sample$(ti)_x_anything.tsv"), io_x_sat_x_an)

    local ts = joinpath(ou, "gene_x_sample$(ti)_x_numbercollapsed.tsv")

    BioLab.Table.write(ts, fe_x_sat_x_nuc)

    BioLab.Plot.plot_heat_map(
        fe_x_sat_x_nuc;
        nar,
        nac = "Sample",
        layout = Dict("title" => Dict("text" => "$da 🗓️ $ti")),
        ht = BioLab.Path.replace_extension(ts, "html"),
    )

end
