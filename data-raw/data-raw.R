source('R/hotspot_utils.R')

data_raw_dir = file.path(
    here::here(),
    "data-raw"
)

gOncoX <- list()

gOncoX[['basic']] <- geneOncoX::get_basic(
    cache_dir = data_raw_dir
)

gOncoX[['basic']]$records <-
    gOncoX[['basic']]$records |>
    dplyr::mutate(entrezgene = as.integer(
        entrezgene
    ))

gOncoX[['alias']] <- geneOncoX::get_alias(
    cache_dir = data_raw_dir)$records |>
    dplyr::filter(n_primary_map == 1 &
                      alias != symbol) |>
    dplyr::select(
        alias, entrezgene
    ) |>
    dplyr::arrange(entrezgene) |>
    dplyr::mutate(property = "alias") |>
    dplyr::rename(value = alias) |>
    dplyr::mutate(entrezgene = as.integer(
        entrezgene
    ))

cancer_hotspots <- load_cancer_hotspots(
    data_raw_dir = data_raw_dir,
    gOncoX = gOncoX
)

usethis::use_data(cancer_hotspots, overwrite = T)

write.table(
    cancer_hotspots$wide, file =
        "~/project_data/data__misc/gvanno/data/grch38/cancer_hotspots/cancer_hotspots.tsv",
    col.names = T, row.names = F, quote = F, sep="\t")
