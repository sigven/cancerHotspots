load_cancer_hotspots <- function(
        data_raw_dir = NA,
        gOncoX = NULL){

    cancer_hotspots <- list()

    cancer_hotspots[['indel']] <-
        openxlsx::read.xlsx(
            xlsxFile = file.path(
                data_raw_dir,
                "hotspots_v2.xlsx"),
            sheet = 2,
            startRow = 1) |>
        janitor::clean_names() |>
        dplyr::select(
            hugo_symbol, qvalue,
            variant_amino_acid, samples, tm) |>
        dplyr::mutate(qvalue = stringr::str_trim(
            format(as.numeric(
                as.character(qvalue)),
                scientific = TRUE, digits = 2))) |>
        dplyr::mutate(var_aa = stringr::str_replace(
            variant_amino_acid,"\\*","X")) |>
        dplyr::mutate(hgvsp = paste0(
            "p.",stringr::str_replace(
                var_aa, ":[0-9]{1,}$",""))) |>
        dplyr::mutate(amino_acid_position = stringr::str_split_fixed(
            tm, " ",2)[,2]) |>
        dplyr::select(-c(var_aa,
                         variant_amino_acid, tm)) |>
        tidyr::separate_rows(samples, sep="\\|") |>
        dplyr::rename(tumor_type_freq = samples) |>
        rename_hotspot_tumor_types() |>
        tidyr::separate(tumor_type_freq,
                        into = c("ttype","freq"),
                        sep = ":",
                        remove = T) |>
        dplyr::mutate(freq = as.integer(freq)) |>
        dplyr::mutate(reference_amino_acid = as.character(NA),
                      codon = as.character(NA)) |>
        dplyr::select(
            hugo_symbol, qvalue, ttype, freq, hgvsp, amino_acid_position,
            reference_amino_acid, ttype) |>
        dplyr::distinct()

    cancer_hotspots[['indel']] <- resolve_gene_symbol(
        df = cancer_hotspots[['indel']], gOncoX = gOncoX) |>
        dplyr::mutate(MUTATION_HOTSPOT = paste0(
            symbol, "|", entrezgene, "|",
            amino_acid_position, "|",
            qvalue)) |>
        dplyr::mutate(MUTATION_HOTSPOT2 = stringr::str_replace_all(
            MUTATION_HOTSPOT, "\\*","X"
        )) |>
        dplyr::mutate(hgvsp2 = stringr::str_replace_all(
            hgvsp, "\\*","X"
        ))

    cancer_hotspots[['snv']] <-  openxlsx::read.xlsx(
        xlsxFile = file.path(
            data_raw_dir,
            "hotspots_v2.xlsx"),
        sheet = 1,
        startRow = 1) |>
        janitor::clean_names() |>
        dplyr::select(
            hugo_symbol, mutation_count, amino_acid_position,
            reference_amino_acid, qvalue, qvalue_pancan,
            qvaluect, detailed_cancer_types, variant_amino_acid,
            samples, total_samples) |>
        dplyr::mutate(reference_amino_acid = stringr::str_replace(
            reference_amino_acid,":[0-9]{1,}","")) |>
        dplyr::mutate(variant_amino_acid = stringr::str_replace(
            variant_amino_acid,":[0-9]{1,}","")) |>
        dplyr::mutate(qvalue = stringr::str_trim(
            format(as.numeric(
                as.character(qvalue)),
                scientific = TRUE, digits = 2))) |>
        dplyr::mutate(hgvsp = paste0(
            "p.", reference_amino_acid,
            amino_acid_position, variant_amino_acid)) |>
        dplyr::mutate(codon = paste0(
            "p.", reference_amino_acid,
            amino_acid_position)) |>
        tidyr::separate_rows(samples, sep="\\|") |>
        dplyr::rename(tumor_type_freq = samples) |>
        rename_hotspot_tumor_types() |>
        tidyr::separate(tumor_type_freq,
                        into = c("ttype","freq"),
                        sep = ":",
                        remove = T) |>
        dplyr::mutate(freq = as.integer(freq)) |>
        dplyr::select(
            hugo_symbol, qvalue, ttype, freq, hgvsp, codon,
            amino_acid_position,
            reference_amino_acid, variant_amino_acid, ttype) |>
        dplyr::distinct() |>
        dplyr::mutate(codon = dplyr::if_else(
            stringr::str_detect(hgvsp,"splice"),
            as.character(NA),
            as.character(codon)
        )) |>
        dplyr::mutate(reference_amino_acid = dplyr::if_else(
            stringr::str_detect(hgvsp,"splice"),
            as.character(NA),
            as.character(reference_amino_acid)
        )) |>
        dplyr::mutate(variant_amino_acid = dplyr::if_else(
            stringr::str_detect(hgvsp,"splice"),
            as.character(NA),
            as.character(variant_amino_acid)
        )) |>
        dplyr::mutate(hgvsp = dplyr::if_else(
            stringr::str_detect(hgvsp,"splice"),
            as.character(NA),
            as.character(hgvsp)
        )) |>
        dplyr::mutate(hgvsc = NA) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X187_splice" &
                hugo_symbol == "TP53",
            "c.559+1G,c.559+2T,c.560-1G,c.560-2A",
            as.character(hgvsc)
        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X187_splice" &
                hugo_symbol == "TP53",
            "c.559:c.560",
            as.character(amino_acid_position)
        )) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X307_splice" &
                hugo_symbol == "TP53",
            "c.919+1G,c.919+2T,c.920-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X307_splice" &
                hugo_symbol == "TP53",
            "c.919:c.920",
            as.character(amino_acid_position)
        )) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X126_splice" &
                hugo_symbol == "TP53",
            "c.376-1G,c.376-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X126_splice" &
                hugo_symbol == "TP53",
            "c.376",
            as.character(amino_acid_position)
        )) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X331_splice" &
                hugo_symbol == "TP53",
            "c.993+1G,c.993+2T",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X331_splice" &
                hugo_symbol == "TP53",
            "c.993",
            as.character(amino_acid_position)
        )) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X225_splice" &
                hugo_symbol == "TP53",
            "c.673-1G,c.673-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X225_splice" &
                hugo_symbol == "TP53",
            "c.673",
            as.character(amino_acid_position)
        )) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X261_splice" &
                hugo_symbol == "TP53",
            "c.782+1G,c.782+2T,c.783-1G,c.783-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X261_splice" &
                hugo_symbol == "TP53",
            "c.782:c.783",
            as.character(amino_acid_position)
        )) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X1010_splice" &
                hugo_symbol == "MET",
            "c.3028+1G,c.3028+2T",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X1010_splice" &
                hugo_symbol == "MET",
            "c.3028",
            as.character(amino_acid_position)
        )) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X70_splice" &
                hugo_symbol == "PTEN",
            "c.209+1G,c.209+2T,c.210-1G,c.210-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X70_splice" &
                hugo_symbol == "PTEN",
            "c.209:c.210",
            as.character(amino_acid_position)
        )) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X212_splice" &
                hugo_symbol == "PTEN",
            "c.634+1G,c.634+2T,c.635-1G,c.635-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X212_splice" &
                hugo_symbol == "PTEN",
            "c.634:c.635",
            as.character(amino_acid_position)
        )) |>
        dplyr::filter(
            !stringr::str_detect(amino_acid_position,"splice") |
                (stringr::str_detect(amino_acid_position,"splice") &
                     !is.na(hgvsc))
        ) |>
        tidyr::separate_rows(
            hgvsc, sep=","
        )


    cancer_hotspots[['snv']] <- resolve_gene_symbol(
        df = cancer_hotspots[['snv']], gOncoX = gOncoX) |>
        dplyr::mutate(MUTATION_HOTSPOT = dplyr::if_else(
            is.na(hgvsc),
            paste0(
                symbol, "|", entrezgene, "|",
                reference_amino_acid,
                amino_acid_position, "|",
                variant_amino_acid,"|",
                as.character(qvalue)),
            paste0(
                symbol, "|", entrezgene, "|",
                amino_acid_position, "|", hgvsc,"|",
                as.character(qvalue))
            )) |>
        dplyr::mutate(MUTATION_HOTSPOT2 = stringr::str_replace_all(
            MUTATION_HOTSPOT, "\\*","X"
        )) |>
        dplyr::mutate(hgvsp2 = stringr::str_replace_all(
            hgvsp, "\\*","X"
        ))


    site_freqs <- list()
    site_freqs[['snv']] <- as.data.frame(
        cancer_hotspots[['snv']] |>
            dplyr::select(
                symbol,
                entrezgene,
                amino_acid_position,
                reference_amino_acid,
                ttype,
                freq) |>
            dplyr::distinct() |>
            dplyr::group_by(
                symbol,
                entrezgene,
                amino_acid_position,
                reference_amino_acid,
                ttype) |>
            dplyr::summarise(
                ttype_site_freq = sum(as.integer(freq)),
                .groups = "drop")
    )

    site_freqs[['indel']] <- as.data.frame(
        cancer_hotspots[['indel']] |>
            dplyr::select(
                symbol,
                entrezgene,
                amino_acid_position,
                ttype,
                freq) |>
            dplyr::group_by(
                symbol,
                entrezgene,
                amino_acid_position,
                ttype) |>
            dplyr::summarise(
                ttype_site_freq = sum(as.integer(freq)),
                .groups = "drop")
    )

    for(t in c('snv','indel')){
        cancer_hotspots[[t]] <- cancer_hotspots[[t]] |>
            dplyr::left_join(site_freqs[[t]]) |>
            dplyr::mutate(vartype = t)
    }

    hotspots_long <- dplyr::bind_rows(
        cancer_hotspots$snv,
        cancer_hotspots$indel
    )

    hotspots_wide <-  as.data.frame(
        dplyr::bind_rows(
            cancer_hotspots$snv,
            cancer_hotspots$indel) |>
            dplyr::mutate(
                MUTATION_HOTSPOT_CANCERTYPE = dplyr::if_else(
                    is.na(hgvsc),
                     paste(
                         ttype, ttype_site_freq, freq, sep="|"),
                    paste(
                        ttype, ttype_site_freq, "", sep="|")
                    )) |>
            dplyr::group_by(
                symbol,
                entrezgene,
                amino_acid_position,
                reference_amino_acid,
                vartype,
                qvalue,
                codon,
                hgvsc,
                hgvsp,
                hgvsp2,
                MUTATION_HOTSPOT,
                MUTATION_HOTSPOT2) |>
            dplyr::summarise(
                MUTATION_HOTSPOT_CANCERTYPE = paste(
                    sort(MUTATION_HOTSPOT_CANCERTYPE), collapse=","
                ),
                .groups = "drop")
    )




    return(list('wide' = hotspots_wide, 'long' = hotspots_long))

}

resolve_gene_symbol <- function(df, gOncoX = NULL){

    df <- df |>
        dplyr::left_join(
            dplyr::select(
                gOncoX$basic$records, symbol, entrezgene),
            by = c("hugo_symbol" = "symbol")
        )

    df1 <- df |>
        dplyr::filter(!is.na(entrezgene)) |>
        dplyr::rename(symbol = hugo_symbol)

    df2 <- df |>
        dplyr::filter(is.na(entrezgene)) |>
        dplyr::select(-entrezgene) |>
        dplyr::left_join(
            dplyr::select(
                gOncoX$alias, value, entrezgene),
            by = c("hugo_symbol" = "value")
        ) |>
        dplyr::select(-hugo_symbol) |>
        dplyr::filter(!is.na(entrezgene)) |>
        dplyr::left_join(
            dplyr::select(
                gOncoX$basic$records, symbol, entrezgene),
            by = "entrezgene")

    df <- dplyr::bind_rows(
        df1, df2
    )

    return(df)
}


rename_hotspot_tumor_types <- function(df){

    df <- df |>
        dplyr::mutate(tumor_type_freq = stringr::str_to_title(tumor_type_freq)) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Adrenalgland","Adrenal_Gland"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Biliarytract","Biliary_Tract"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Bladder","Bladder@Urinary_Tract"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Blood","Myeloid"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Cnsbrain","CNS@Brain"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Bowel","Colon@Rectum"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Esophagusstomach","Esophagus@Stomach"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Ovaryfallopiantube","Ovary@Fallopian_Tube"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Headandneck","Head_and_Neck"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Ampullaofvater","Ampulla_of_Vater"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Vulvavagina","Vulva@Vagina"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Softtissue","Soft_Tissue"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Unk","Unknown"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Lymph","Lymphoid"
        ))

    return(df)
}


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

