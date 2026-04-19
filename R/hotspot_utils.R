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
            hugo_symbol, qvalue,
            ttype, freq, hgvsp, amino_acid_position,
            reference_amino_acid, ttype) |>
        dplyr::distinct()

    cancer_hotspots[['indel']] <- resolve_gene_symbol(
        df = cancer_hotspots[['indel']], gOncoX = gOncoX) |>
        dplyr::mutate(MUTATION_HOTSPOT = paste0(
            symbol, "|", entrezgene, "|",
            amino_acid_position, "||",
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
            hugo_symbol,  amino_acid_position,
            reference_amino_acid, qvalue,
            variant_amino_acid,
            samples, total_samples) |>
        #dplyr::mutate(amino_acid_position = as.numeric(
        #    amino_acid_position
        #)) |>
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
        dplyr::mutate(hgvsc = NA_character_) |>
        apply_splice_hotspot_curation(data_raw_dir = data_raw_dir) |>

        dplyr::filter(
            !stringr::str_detect(amino_acid_position,"splice") |
                (stringr::str_detect(amino_acid_position,"splice") &
                     !is.na(hgvsc))
        ) |>
        tidyr::separate_rows(
            hgvsc, sep=","
        ) |>
        dplyr::select(
            hugo_symbol, qvalue, ttype, freq,
            hgvsc, hgvsp, codon,
            amino_acid_position,
            reference_amino_acid,
            variant_amino_acid)


    hotspot_v3_qvalues <- readxl::read_excel(
        "data-raw/hotspots_v3_single_residue_and_indels.xlsx",
        sheet = "Sheet") |>
        janitor::clean_names() |>
        dplyr::select(hugo_symbol, codon,
                      codon_position, q_value) |>
        dplyr::mutate(qvalue = stringr::str_trim(
            format(as.numeric(
                as.character(q_value)),
                scientific = TRUE, digits = 2))) |>
        dplyr::mutate(amino_acid_position = as.character(codon_position)) |>
        dplyr::mutate(
            reference_amino_acid = stringr::str_replace(
                codon,"[0-9]{1,}","")) |>
        dplyr::select(hugo_symbol, reference_amino_acid,
                      amino_acid_position, qvalue, codon)

    cancer_hotspots[['snv_v3']] <-
        readxl::read_excel(
            "data-raw/hotspots_v3_single_residue_and_indels.xlsx",
            sheet = "SNV_Variants") |>
        janitor::clean_names() |>
        dplyr::filter(hotspot_version == "v3") |>
        dplyr::mutate(amino_acid_position = as.character(codon_position)) |>
        dplyr::select(-c("codon_position","hotspot_version")) |>
        dplyr::left_join(hotspot_v3_qvalues) |>
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
        dplyr::mutate(freq = as.integer(freq),
                      hgvsc = as.character(NA)) |>
        dplyr::select(
            hugo_symbol, qvalue, ttype, freq,
            hgvsc, hgvsp, codon,
            dplyr::everything())

    cancer_hotspots[['snv']] <- dplyr::bind_rows(
        cancer_hotspots[['snv']],
        cancer_hotspots[['snv_v3']]) |>
        dplyr::distinct()


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





    # cancer_hotspots[['snv_v3']] <-  openxlsx::read.xlsx(
    #     xlsxFile = file.path(
    #         data_raw_dir,
    #         "hotspots_v3.xlsx"),
    #     sheet = 1,
    #     startRow = 1) |>
    #     janitor::clean_names() |>
    #     dplyr::select(
    #         hugo_symbol, mutation_count, amino_acid_position,
    #         reference_amino_acid, qvalue, qvalue_pancan,
    #         qvaluect, detailed_cancer_types, variant_amino_acid,
    #         samples, total_samples) |>
    #     dplyr::mutate(reference_amino_acid = stringr::str_replace(
    #         reference_amino_acid,":[0-9]{1,}","")) |>
    #     dplyr::mutate(variant_amino_acid = stringr::str_replace(
    #         variant_amino_acid,":[0-9]{1,}","")) |>
    #     dplyr::mutate(qvalue = stringr::str_trim(
    #         format(as.numeric(
    #             as.character(qvalue)),
    #             scientific = TRUE, digits = 2))) |>
    #     dplyr::mutate(hgvsp = paste0(
    #         "p.", reference_amino_acid,
    #         amino_acid_position, variant_amino_acid)) |>
    #     dplyr::mutate(codon = paste0(
    #         "p.", reference_amino_acid,
    #         amino_acid_position)) |>
    #     tidyr::separate_rows(samples, sep="\\|") |>
    #     dplyr::rename(tumor_type_freq = samples) |>
    #     rename_hotspot_tumor_types() |>
    #     tidyr::separate(tumor_type_freq,
    #                     into = c("ttype","freq"),
    #                     sep = ":",
    #                     remove = T) |>
    #     dplyr::mutate(freq = as.integer(freq)) |>
    #     dplyr::select(
    #         hugo_symbol, qvalue, ttype, freq, hgvsp, codon,
    #         amino_acid_position,
    #         reference_amino_acid, variant_amino_acid, ttype)


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


    metadata_hotspots <-
        data.frame(
            'source' = 'cancerhotspots.org',
            'source_description' = paste0(
                'A resource for statistically significant mutations in cancer'),
            'source_url' = 'https://www.cancerhotspots.org/#/home',
            'source_citation' = 'Chang et al., Cancer Discov, 2018; 29247016 | Bandlamudi et al., Cancer Cell. 2026; 41895280',
            'source_version' = 'v2+v3',
            'source_abbreviation' = 'hotspots',
            'source_license' = 'ODbL v1.0',
            'source_license_url' = 'https://opendatacommons.org/licenses/odbl/1-0/'
        )


    return(list(
        'metadata' = metadata_hotspots,
        'wide' = hotspots_wide,
        'long' = hotspots_long))

}

read_splice_hotspot_curation <- function(data_raw_dir){

    utils::read.delim(
        file = file.path(
            data_raw_dir,
            "splice_hotspot_curation.tsv"
        ),
        sep = "\t",
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

}

apply_splice_hotspot_curation <- function(df, data_raw_dir){

    splice_hotspot_curation <- read_splice_hotspot_curation(
        data_raw_dir = data_raw_dir
    )

    df |>
        dplyr::left_join(
            splice_hotspot_curation,
            by = c("hugo_symbol", "amino_acid_position")
        ) |>
        dplyr::mutate(
            amino_acid_position = dplyr::coalesce(
                curated_amino_acid_position,
                amino_acid_position
            ),
            hgvsc = dplyr::coalesce(
                curated_hgvsc,
                hgvsc
            )
        ) |>
        dplyr::select(
            -curated_amino_acid_position,
            -curated_hgvsc
        )

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
            tumor_type_freq, "Adrenalgland|Adrenal_gland","Adrenal_Gland"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Biliarytract|Biliary_tract","Biliary_Tract"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Bladder","Bladder@Urinary_Tract"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Blood","Myeloid"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Cnsbrain|Brain","CNS@Brain"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Bowel","Colon@Rectum"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Esophagusstomach","Esophagus@Stomach"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Ovaryfallopiantube|Ovary","Ovary@Fallopian_Tube"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Headandneck|Head_neck","Head_and_Neck"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Ampullaofvater|Ampulla_of_vater","Ampulla_of_Vater"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Vulvavagina|Vulva","Vulva@Vagina"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Softtissue|Soft_tissue","Soft_Tissue"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Unk","Unknown"
        )) |>
        dplyr::mutate(tumor_type_freq = stringr::str_replace(
            tumor_type_freq, "Lymph","Lymphoid"
        ))

    return(df)
}
