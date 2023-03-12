include("_.jl")

# --------------------------------------------- #

da = "SDY1264"

ind = joinpath(IN, da)

ou = mkpath(joinpath(OU, da))

validate = (true, true)

# --------------------------------------------- #

pa_x_io_x_an = read(ind, "gene_expression_files_2023-03-12_12-02-58.tsv")

sort!(pa_x_io_x_an, [1, 3])

sa_pat__ = [sa => "$pa.$(Int(ti))" for (pa, sa, ti) in eachrow(pa_x_io_x_an)]

pat_ = (sa_pat[2] for sa_pat in sa_pat__)

@test length(unique(sa_pat[1] for sa_pat in sa_pat__)) == length(unique(pat_))


# --------------------------------------------- #

pa_x_io_x_an = outerjoin(
    read(ind, "demographics_2023-03-12_11-08-48.tsv"),
    read(ind, "cohort_membership_2023-03-12_11-06-40.tsv"),
    read(ind, "neut_ab_titer_2023-03-12_11-06-06.tsv");
    on = "Participant ID",
    validate,
)

rename!(pa_x_io_x_an, "Value Preferred" => "Antibody Titer")

# --------------------------------------------- #

io_x_sa_x_an = DataFrame("Information" => names(pa_x_io_x_an))

pa_id = Dict(pa => id for (id, pa) in enumerate(pa_x_io_x_an[!, 1]))

for pat in pat_

    pa = join(split(pat, '.')[1:2], '.')

    id = pa_id[pa]

    io_x_sa_x_an[!, pat] = collect(pa_x_io_x_an[id, :])

end

# --------------------------------------------- #

BioLab.Table.write(joinpath(ou, "information_x_sample_x_anything.tsv"), io_x_sa_x_an)

# --------------------------------------------- #

fe_x_sa_x_nu = outerjoin(
    read(ind, "SDY1264_PBMC_Trial1_Geo.tsv"),
    read(ind, "SDY1264_PBMC_Trial2_Geo.tsv");
    on = "feature_id",
    validate,
)

rename!(fe_x_sa_x_nu, sa_pat__...)

@test names(io_x_sa_x_an)[2:end] == names(fe_x_sa_x_nu)[2:end]

fe_ge = BioLab.DataFrame.map_to(
    read(ind, "FeatureAnnotation_2023-03-12_11-04-20.tsv"),
    BioLab.Dict.set_with_last!,
    ["Feature Id"],
    "Gene Symbol",
)

fe_ = fe_x_sa_x_nu[!, 1]

fe_ = BioLab.Gene.rename(fe_, fe_ge)[1]

fe_ = BioLab.Gene.rename(fe_, BioLab.Gene.map_ensembl(), BioLab.Gene.map_hgnc())[1]

fe_x_sa_x_nu[!, 1] = fe_

nar = "Gene"

rename!(fe_x_sa_x_nu, 1 => nar)

# --------------------------------------------- #

fe_x_sa_x_nuc = BioLab.DataFrame.collapse(fe_x_sa_x_nu)

# --------------------------------------------- #

ts = joinpath(ou, "gene_x_sample_x_numbercollapsed.tsv")

BioLab.Table.write(ts, fe_x_sa_x_nuc)

# --------------------------------------------- #

BioLab.Plot.plot_heat_map(
    fe_x_sa_x_nuc;
    nar,
    nac = "Sample",
    grc_ = collect(io_x_sa_x_an[4, 2:end]),
    layout = Dict("title" => Dict("text" => da)),
    ht = BioLab.Path.replace_extension(ts, "html"),
)
