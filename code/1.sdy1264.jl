include("_.jl")

# --------------------------------------------- #

da = "SDY1264"

ind = joinpath(IN, da)

ou = mkpath(joinpath(OU, da))

validate = (true, true)

# --------------------------------------------- #

pe_x_io_x_an = read(ind, "gene_expression_files_2023-03-12_12-02-58.tsv")

sort!(pe_x_io_x_an, [1, 3])

sapet_ = [sa => "$pe.$(Int(ti))" for (pe, sa, ti) in eachrow(pe_x_io_x_an)]

pet_ = (sapet[2] for sapet in sapet_)

@test length(unique(sapet[1] for sapet in sapet_)) == length(unique(pet_))

# --------------------------------------------- #

pe_x_io_x_an = outerjoin(
    read(ind, "demographics_2023-03-12_11-08-48.tsv"),
    read(ind, "cohort_membership_2023-03-12_11-06-40.tsv"),
    read(ind, "neut_ab_titer_2023-03-12_11-06-06.tsv");
    on = "Participant ID",
    validate,
)

rename!(pe_x_io_x_an, "Value Preferred" => "Antibody Titer")

# --------------------------------------------- #

io_x_sa_x_an = DataFrame("Information" => names(pe_x_io_x_an))

pe_id = Dict(pe => id for (id, pe) in enumerate(pe_x_io_x_an[!, 1]))

for pet in pet_

    pe = join(split(pet, '.')[1:2], '.')

    id = pe_id[pe]

    io_x_sa_x_an[!, pet] = collect(pe_x_io_x_an[id, :])

end

sa_ = names(io_x_sa_x_an)[2:end]

# --------------------------------------------- #

BioLab.Table.write(joinpath(ou, "information_x_sample_x_anything.tsv"), io_x_sa_x_an)

# --------------------------------------------- #

fe_x_sa_x_nu = outerjoin(
    read(ind, "SDY1264_PBMC_Trial1_Geo.tsv"),
    read(ind, "SDY1264_PBMC_Trial2_Geo.tsv");
    on = "feature_id",
    validate,
)

rename!(fe_x_sa_x_nu, sapet_...)

@test sa_ == names(fe_x_sa_x_nu)[2:end]

# --------------------------------------------- #

fe_ge = BioLab.DataFrame.map_to(
    read(ind, "FeatureAnnotation_2023-03-12_11-04-20.tsv"),
    BioLab.Dict.set_with_last!,
    ["Feature Id"],
    "Gene Symbol",
)

fe_ = fe_x_sa_x_nu[!, 1]

fe2_ = BioLab.Gene.rename(fe_, fe_ge)[1]

fe3_ = BioLab.Gene.rename(fe2_, BioLab.Gene.map_ensembl(), BioLab.Gene.map_hgnc())[1]

fe_x_sa_x_nu[!, 1] = fe3_

nar = "Gene"

rename!(fe_x_sa_x_nu, 1 => nar)

# --------------------------------------------- #

fe_x_sa_x_nuc = BioLab.DataFrame.collapse(fe_x_sa_x_nu)

# --------------------------------------------- #

ts = joinpath(ou, "gene_x_sample_x_numbercollapsed.tsv")

BioLab.Table.write(ts, fe_x_sa_x_nuc)

# --------------------------------------------- #

ti_ = [parse(Int, BioLab.String.split_and_get(sa, '.', 3)) for sa in sa_]

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

    ts = joinpath(ou, "gene_x_sample$(ti)_x_numbercollapsed.tsv")

    BioLab.Table.write(ts, fe_x_sat_x_nuc)

    BioLab.Plot.plot_heat_map(
        fe_x_sat_x_nuc;
        nar,
        nac = "Sample",
        layout = Dict("title" => Dict("text" => "$da 🗓️ $ti")),
        ht = BioLab.Path.replace_extension(ts, "html"),
    )

end
