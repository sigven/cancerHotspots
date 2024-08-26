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
            amino_acid_position == "X125_splice" &
                hugo_symbol == "TP53",
            "c.375+1G,c.375+2T",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X125_splice" &
                hugo_symbol == "TP53",
            "c.375",
            as.character(amino_acid_position)
        )) |>

        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X224_splice" &
                hugo_symbol == "TP53",
            "c.672+1G,c.672+2T",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X224_splice" &
                hugo_symbol == "TP53",
            "c.672",
            as.character(amino_acid_position)
        )) |>

        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X33_splice" &
                hugo_symbol == "TP53",
            "c.97-1G,c.97-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X33_splice" &
                hugo_symbol == "TP53",
            "c.97",
            as.character(amino_acid_position)
        )) |>

        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X332_splice" &
                hugo_symbol == "TP53",
            "c.994-1G,c.994-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X332_splice" &
                hugo_symbol == "TP53",
            "c.994",
            as.character(amino_acid_position)
        )) |>

        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X32_splice" &
                hugo_symbol == "TP53",
            "c.96+1G,c.96+2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X32_splice" &
                hugo_symbol == "TP53",
            "c.96",
            as.character(amino_acid_position)
        )) |>


        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X114_splice" &
                hugo_symbol == "VHL",
            "c.340+1G,c.340+2T,c.341-1G,c.341-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X114_splice" &
                hugo_symbol == "VHL",
            "c.340:c.341",
            as.character(amino_acid_position)
        )) |>

        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X155_splice" &
                hugo_symbol == "VHL",
            "c.463+1G,c.463+2T,c.464-1G,c.464-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X155_splice" &
                hugo_symbol == "VHL",
            "c.463:c.464",
            as.character(amino_acid_position)
        )) |>


        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X342_splice" &
                hugo_symbol == "PTEN",
            "c.1026+1G,c.1026+2T",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X342_splice" &
                hugo_symbol == "PTEN",
            "c.1026",
            as.character(amino_acid_position)
        )) |>

        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X267_splice" &
                hugo_symbol == "PTEN",
            "c.801+1G,c.801+2T",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X267_splice" &
                hugo_symbol == "PTEN",
            "c.801",
            as.character(amino_acid_position)
        )) |>


        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X1010_splice" &
                hugo_symbol == "MET",
            "c.3028+1G,c.3028+2T,c.3082+1G,c.3082+2T",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X1010_splice" &
                hugo_symbol == "MET",
            "c.3028,c.3082",
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
            amino_acid_position == "X343_splice" &
                hugo_symbol == "PTEN",
            "c.1027-1G,c.1027-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X343_splice" &
                hugo_symbol == "PTEN",
            "c.1027",
            as.character(amino_acid_position)
        )) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X342_splice" &
                hugo_symbol == "PTEN",
            "c.1026+1G,c.1026+2T",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X342_splice" &
                hugo_symbol == "PTEN",
            "c.1026",
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
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X268_splice" &
                hugo_symbol == "PTEN",
            "c.802-1G,c.802-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X268_splice" &
                hugo_symbol == "PTEN",
            "c.802",
            as.character(amino_acid_position)
        )) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X85_splice" &
                hugo_symbol == "PTEN",
            "c.253+1G,c.253+2T,c.254-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X85_splice" &
                hugo_symbol == "PTEN",
            "c.253:c.254",
            as.character(amino_acid_position)
        )) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X51_splice" &
                hugo_symbol == "CDKN2A",
            "c.151-1G,c.151-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X51_splice" &
                hugo_symbol == "CDKN2A",
            "c.151",
            as.character(amino_acid_position)
        )) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X153_splice" &
                hugo_symbol == "CDKN2A",
            "c.457+1G,c.457+2T,c.458-1G,c.458-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X153_splice" &
                hugo_symbol == "CDKN2A",
            "c.457:c.458",
            as.character(amino_acid_position)
        )) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X155_splice" &
                hugo_symbol == "STK11",
            "c.464+1G,c.465-1G,c.465-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X155_splice" &
                hugo_symbol == "STK11",
            "c.464:c.465",
            as.character(amino_acid_position)
        )) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X245_splice" &
                hugo_symbol == "STK11",
            "c.734+1G,c.734+2T,c.735-1G",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X245_splice" &
                hugo_symbol == "STK11",
            "c.734:c.735",
            as.character(amino_acid_position)
        )) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X307_splice" &
                hugo_symbol == "STK11",
            "c.920+1G,c.921-2A,c.921-1G",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X307_splice" &
                hugo_symbol == "STK11",
            "c.920:c.921",
            as.character(amino_acid_position)
        )) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X288_splice" &
                hugo_symbol == "STK11",
            "c.862+1G,c.863-1G,c.863-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X288_splice" &
                hugo_symbol == "STK11",
            "c.862:c.863",
            as.character(amino_acid_position)
        )) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X97_splice" &
                hugo_symbol == "STK11",
            "c.290+1G,c.291-1G,c.291-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X97_splice" &
                hugo_symbol == "STK11",
            "c.290:c.291",
            as.character(amino_acid_position)
        )) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X405_splice" &
                hugo_symbol == "RB1",
            "c.1215+1G",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X405_splice" &
                hugo_symbol == "RB1",
            "c.1215",
            as.character(amino_acid_position)
        )) |>


        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X406_splice" &
                hugo_symbol == "RB1",
            "c.1216-1G,c.1216-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X406_splice" &
                hugo_symbol == "RB1",
            "c.1216",
            as.character(amino_acid_position)
        )) |>

        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X463_splice" &
                hugo_symbol == "RB1",
            "c.1389+1G",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X463_splice" &
                hugo_symbol == "RB1",
            "c.1389",
            as.character(amino_acid_position)
        )) |>

        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X500_splice" &
                hugo_symbol == "RB1",
            "c.1498+1G,c.1498+2T,c.1499-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X500_splice" &
                hugo_symbol == "RB1",
            "c.1498:c.1499",
            as.character(amino_acid_position)
        )) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X474_splice" &
                hugo_symbol == "RB1",
            "c.1421+1G,c.1422-1G,c.1422-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X474_splice" &
                hugo_symbol == "RB1",
            "c.1421:c.1422",
            as.character(amino_acid_position)
        )) |>

        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X654_splice" &
                hugo_symbol == "RB1",
            "c.1960+1G,c.1960+2T,c.1961-1G,c.1961-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X654_splice" &
                hugo_symbol == "RB1",
            "c.1960:c.1961",
            as.character(amino_acid_position)
        )) |>


        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X830_splice" &
                hugo_symbol == "RB1",
            "c.2489+1G,c.2489+2T,c.2490-1G",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X830_splice" &
                hugo_symbol == "RB1",
            "c.2489:c.2490",
            as.character(amino_acid_position)
        )) |>

        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X445_splice" &
                hugo_symbol == "RB1",
            "c.1333-1G",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X445_splice" &
                hugo_symbol == "RB1",
            "c.1333",
            as.character(amino_acid_position)
        )) |>

        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X127_splice" &
                hugo_symbol == "RB1",
            "c.380+1G,c.381-1G,c.381-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X127_splice" &
                hugo_symbol == "RB1",
            "c.380:c.381",
            as.character(amino_acid_position)
        )) |>

        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X203_splice" &
                hugo_symbol == "RB1",
            "c.607+1G,c.608-1G",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X203_splice" &
                hugo_symbol == "RB1",
            "c.607:c.608",
            as.character(amino_acid_position)
        )) |>

        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X46_splice" &
                hugo_symbol == "RB1",
            "c.137+1G,c.137+2T,c.138-1G,c.138-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X46_splice" &
                hugo_symbol == "RB1",
            "c.137:c.138",
            as.character(amino_acid_position)
        )) |>

        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X702_splice" &
                hugo_symbol == "RB1",
            "c.2106+1G,c.2106+2T",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X702_splice" &
                hugo_symbol == "RB1",
            "c.2106",
            as.character(amino_acid_position)
        )) |>

        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X605_splice" &
                hugo_symbol == "RB1",
            "c.1814+1G,c.1814+2T,c.1815-1G",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X605_splice" &
                hugo_symbol == "RB1",
            "c.1814:c.1815",
            as.character(amino_acid_position)
        )) |>


        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X582_splice" &
                hugo_symbol == "PIK3R1",
            "c.1746-1G,c.1746-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X582_splice" &
                hugo_symbol == "PIK3R1",
            "c.1746",
            as.character(amino_acid_position)
        )) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X87_splice" &
                hugo_symbol == "SDHAF2",
            "c.261-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X87_splice" &
                hugo_symbol == "SDHAF2",
            "c.261",
            as.character(amino_acid_position)
        )) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X278_splice" &
                hugo_symbol == "CDH1",
            "c.832+1G,c.833-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X278_splice" &
                hugo_symbol == "CDH1",
            "c.832:c.833",
            as.character(amino_acid_position)
        )) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X177_splice" &
                hugo_symbol == "CDH1",
            "c.531+1G,c.531+2T",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X177_splice" &
                hugo_symbol == "CDH1",
            "c.531",
            as.character(amino_acid_position)
        )) |>
        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X440_splice" &
                hugo_symbol == "CDH1",
            "c.1320+1G",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X440_splice" &
                hugo_symbol == "CDH1",
            "c.1320",
            as.character(amino_acid_position)
        )) |>

        dplyr::mutate(hgvsc = dplyr::if_else(
            amino_acid_position == "X722_splice" &
                hugo_symbol == "CDH1",
            "c.2165-2A",
            as.character(hgvsc)

        )) |>
        dplyr::mutate(amino_acid_position = dplyr::if_else(
            amino_acid_position == "X722_splice" &
                hugo_symbol == "CDH1",
            "c.2165",
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


    metadata_hotspots <-
        data.frame(
            'source' = 'cancerhotspots.org',
            'source_description' = paste0(
                'A resource for statistically significant mutations in cancer'),
            'source_url' = 'https://www.cancerhotspots.org/#/home',
            'source_citation' = 'Chang et al., Cancer Discov, 2018; 29247016',
            'source_version' = 'v2',
            'source_abbreviation' = 'hotspots',
            'source_license' = 'ODbL v1.0',
            'source_license_url' = 'https://opendatacommons.org/licenses/odbl/1-0/'
        )


    return(list(
        'metadata' = metadata_hotspots,
        'wide' = hotspots_wide,
        'long' = hotspots_long))

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

write.table(
    cancer_hotspots$wide, file =
        "~/project_data/data__misc/gvanno/data/grch38/cancer_hotspots/cancer_hotspots.tsv",
    col.names = T, row.names = F, quote = F, sep="\t")
