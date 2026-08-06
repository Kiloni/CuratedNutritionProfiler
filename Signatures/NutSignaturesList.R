NutSignatures = list(
  "A_Retinol_GWAS" = c("FFAR4", "B4GALT6", "TTR", "ADAMTS3", "HNRNPA1P67", 
                       "RBP4", "GCKR", "PCCB", "RPL31P23", "SOCS2", "CRADD", 
                       "DHRS9", "UGT2B11", "DHRS4", "UGT2A1", "ADH1A", "ADH1B", 
                       "ADH1C", "ADH4", "ADH5", "ADH6", "ADH7", "RDH12", 
                       "CYP1A1", "CYP1A2", "CYP2A6", "CYP2A7", "CYP3A7", 
                       "CYP2A13", "CYP2B6", "CYP2C19", "CYP2C8", 
                       "CYP2C9", "CYP2C18", "RDH10", "CYP3A4", "CYP3A5", 
                       "CYP4A11", "AWAT2", "CYP26A1", "ALDH1A1", "CYP4A22", 
                       "DHRS4L2", "CYP26C1", "RDH8", "RDH11", "BCO1", 
                       "UGT2B28", "UGT1A10", "UGT1A8", "UGT1A7", "UGT1A6", 
                       "UGT1A5", "UGT1A9", "UGT1A4", "UGT1A1", "UGT1A3", 
                       "RETSAT", "CYP26B1", "RDH5", "RPE65", "CYP3A43", 
                       "UGT2B4", "UGT2B7", "UGT2B10", "UGT2B15", "UGT2B17", 
                       "UGT2A3", "PNPLA4", "DGAT2", "RDH16", "DGAT1", 
                       "ALDH1A2", "LRAT", "DHRS3"),
  "A_Retinol_MSigDB_Filtered" = c("PNPLA4", "LRAT", "RPE65", "ADH1A", 
                                  "DHRS4", "CYP2A6", "DHRS9", "RDH12", 
                                  "CYP1A1", "BCO1", "AWAT2", "ALDH1A1", 
                                  "CYP3A7", "CYP2W1", "ADH7"),
  "A_Retinol_GLM_43" = c("Npas2", "Coq10b", "Slc45a3", "1700016C15Rik",
                         "Fpgs", "Itgb6", "Hnrnpa3", "Slc50a1", "Ecm1", 
                         "Slc16a1", "Snx7", "Calb1", "Plppr1", "Phactr4", 
                         "Wdtc1", "Wscd2", "Atp6v1b1", "Cyp26b1", "Itpr2", 
                         "Nphs1", "Klk1b11", "Klk1", "Scnn1g", "Mtmr7", 
                         "Clpx", "Cgnl1", "Cish", "Dhx30", "Ipcef1", "Ftcd", 
                         "Avpr1a", "Hist3h2a", "Sez6", "Trib2", "Gmpr", "Nr1d2", 
                         "Slitrk1", "Angpt1", "Ndrg1", "Cldn8", "Tex16", "Rgag1", 
                         "Hccs"),
  "A_Retinol_Quiles_20" = c("Lrat", "Tnni3k", "Dancr", "Suv39h2", "Apol7c", "Tnni3", 
                            "Fndc5", "St6galnac5", "Neil3", "Cxcr6", "Slc47a1", 
                            "Myoz2", "Cenpi", "Slc17a7", "Lrch2", "Shcbp1l", 
                            "Cox6a2", "Hspb7", "Grik4", "Gm7511"),
  "B1_Thiamine" = c("CCDC82", "WRN", "TMEM16E", "SERPINI1", "RHOH",
                    "NETO1", "MUSK", "MAML2", "NRG1", "NELL1", 
                    "ANO5", "GOLIM4", "N4BP2", "LINC01899", 
                    "CBLN2", "LPAR1", "BCKDK", "SLC19A2", "TPK", 
                    "OGDH", "PDHA1", "BCKDHA", "SLC25A19", "TKT", 
                    "SLC19A3"),
  "B3_Niacin" = c("SLC22A3", "HERPUD1", "CETP"),
  "B5_Pantothenic Acid" = c("TMEM132D", "IGFBP7", "RTEL1", "RTEL1-TNFRSF6B", "LY96", 
                            "TSNAX", "SNRPD2P2", "MIR302F", "RNU6-408P", "RECK", 
                            "OR2AM1P", "RNA5SP517", "MIR891A", "SYT16", "LINC02198", 
                            "PIK3R1", "ERG", "LINC02731"),
  "B6_Pyridoxine" = c("ALPL", "NBPF3", "ADCYAP1R1", "ZNF608", "SMAD3-DT", "OSBPL3", "CYCS",
                        "PYROXD2", "CACNA1C"),
  "B9_Folate"      = c("MTHFR", "DHFR", "SLC19A1", "FOLR1", "SHMT1", "MTHFD1"),
  "B12_Cobalamin" = c("E2F1", "Rps3", "Rps6", "Rpl5", "Rpl7", "Rpl11", "Rpl13a", "Rps14",
                      "Rps27a", "Rpl22"),
  "Vitamin_C" = c("SLC23A1", "SLC23A3", "CHPT1", "BCAS3", "SNRPF", "RER1", "MAF", "GSTA5", 
                  "RGS14", "AKT1", "FADS1"),
  "Vitamin_D_Aigner_35" = c("CD24", "CDH1", "CDH11", "CDH3", "CDH4", "CLDN7", "CRB3", "CXADR", 
                            "DMKN", "DSC2", "DSP", "EPCAM", "EPPK1", "F11R", "GJB2", "GJB3", 
                            "MAL2", "MARVELD2", "MPZL2", "MUC1", "OCLN", "PATJ", "PCDH7", "PKP3", 
                            "PMEPA1", "PPL", "SCEL", "SFN", "SH3YL1", "SHROOM3", "SYTL1", 
                            "TACSTD2", "TMEM30B", "TSPAN1", "TSPAN15"),
  "Vitamin_E" = c("CD36", "SRD5A1", "GCLC", "PKC", "CCNB2", "CDC2", "CDC6", "TGM2", "MMP1", 
                  "MIR122", "MIR125B"),
  "Vitamin_K" = c("Gm35339", "Gm20559", "AC123061.1", "Gm7240", "Gm37494", "Btbd8",
                  "D430020J02Rik", "Zfp862-ps", "Gm9856", "AC117588.1", "Gm19605", "AC174678.1",
                  "Gm35315", "Gm10033", "BC024063", "Gm45250", "Tes3-ps", "Gm45871", "Nup62",
                  "Mpv17", "Gm43682", "Chmp1b", "Gm46430", "Samd4b", "AC131586.2", "Peg13",
                  "Gm44250", "Gm4130", "4933412E12Rik", "Gm18860", "Gm43597", "A830082N09Rik",
                  "Gm45133", "Gm43681", "Gm20075", "AC154636.2", "Gm32856", "Gm4366", "Snhg1",
                  "Gm45456", "AC163032.1", "Gtf3c2", "Hspa14", "Prmt1", "0610012G03Rik", 
                  "Gm9616", "AC121965.1", "Mrpl33", "Gm10095", "AC129931.1"),
  "Iron_Homeostasis" = c("TFRC", "FTH1", "FTL", "SLC11A2", "ACO1", "IREB2"),
  "Omega_3_FattyAcids" = c("FADS1", "FADS2", "ELOVL2", "ELOVL5", "SCD"),
  "Carnitine_Mito" = c("TMLHE", "SLC22A5", "CPT1A", "CPT2", "CRAT")
)

Human_NutSignatures <- list(
  "A_Retinol" = c("LRAT", "TNNI3K", "DANCR", "SUV39H2", "APOL7C", "TNNI3", 
                  "FNDC5", "ST6GALNAC5", "NEIL3", "CXCR6", "SLC47A1", 
                  "MYOZ2", "CENPI", "SLC17A7", "LRCH2", "SHCBP1L", 
                  "COX6A2", "HSPB7", "GRIK4"),
  "B1_Thiamine" = c("CCDC82", "WRN", "TMEM16E", "SERPINI1", "RHOH", "NETO1", 
                    "MUSK", "MAML2", "NRG1", "NELL1", "ANO5", "GOLIM4", 
                    "N4BP2", "LINC01899", "CBLN2", "LPAR1", "BCKDK", 
                    "SLC19A2", "TPK", "OGDH", "PDHA1", "BCKDHA", 
                    "SLC25A19", "TKT", "SLC19A3"),
  "B3_Niacin" = c("SLC22A3", "HERPUD1", "CETP"),
  "B5_Pantothenic_Acid" = c("TMEM132D", "IGFBP7", "RTEL1", "RTEL1-TNFRSF6B", "LY96", 
                            "TSNAX", "SNRPD2P2", "MIR302F", "RNU6-408P", "RECK", 
                            "OR2AM1P", "RNA5SP517", "MIR891A", "SYT16", "LINC02198", 
                            "PIK3R1", "ERG", "LINC02731"),
  "B6_Pyridoxine" = c("ALPL", "NBPF3", "ADCYAP1R1", "ZNF608", "SMAD3-DT", "OSBPL3", 
                      "CYCS", "PYROXD2", "CACNA1C"),
  "B9_Folate" = c("MTHFR", "DHFR", "SLC19A1", "FOLR1", "SHMT1", "MTHFD1"),
  "B12_Cobalamin" = c("E2F1", "RPS3", "RPS6", "RPL5", "RPL7", "RPL11", "RPL13A", 
                      "RPS14", "RPS27A", "RPL22"),
  "Vitamin_C" = c("SLC23A1", "SLC23A3", "CHPT1", "BCAS3", "SNRPF", "RER1", "MAF", 
                  "GSTA5", "RGS14", "AKT1", "FADS1"),
  "Vitamin_D" = c("CD24", "CDH1", "CDH11", "CDH3", "CDH4", "CLDN7", "CRB3", "CXADR", 
                  "DMKN", "DSC2", "DSP", "EPCAM", "EPPK1", "F11R", "GJB2", "GJB3", 
                  "MAL2", "MARVELD2", "MPZL2", "MUC1", "OCLN", "PATJ", "PCDH7", 
                  "PKP3", "PMEPA1", "PPL", "SCEL", "SFN", "SH3YL1", "SHROOM3", 
                  "SYTL1", "TACSTD2", "TMEM30B", "TSPAN1", "TSPAN15"),
  "Vitamin_E" = c("CD36", "SRD5A1", "GCLC", "PKC", "CCNB2", "CDC2", "CDC6", "TGM2", 
                  "MMP1", "MIR122", "MIR125B"),
  "Vitamin_K" = c("BTBD8", "NUP62", "MPV17", "CHMP1B", "SAMD4B", "PEG13", "SNHG1", 
                  "GTF3C2", "HSPA14", "PRMT1", "MRPL33"),
  "Iron_Homeostasis" = c("TFRC", "FTH1", "FTL", "SLC11A2", "ACO1", "IREB2"),
  "Omega_3_FattyAcids" = c("FADS1", "FADS2", "ELOVL2", "ELOVL5", "SCD"),
  "Carnitine_Mito" = c("TMLHE", "SLC22A5", "CPT1A", "CPT2", "CRAT")
)

###Run this code if you want to add the NutSignatures to TBsignatures

TBsignatures[["A_Retinol"]] = c("LRAT", "TNNI3K", "DANCR", "SUV39H2", "APOL7C", "TNNI3", 
                                "FNDC5", "ST6GALNAC5", "NEIL3", "CXCR6", "SLC47A1", 
                                "MYOZ2", "CENPI", "SLC17A7", "LRCH2", "SHCBP1L", 
                                "COX6A2", "HSPB7", "GRIK4")

TBsignatures[["B1_Thiamine"]] = c("CCDC82", "WRN", "TMEM16E", "SERPINI1", "RHOH",
                                  "NETO1", "MUSK", "MAML2", "NRG1", "NELL1", 
                                  "ANO5", "GOLIM4", "N4BP2", "LINC01899", 
                                  "CBLN2", "LPAR1", "BCKDK", "SLC19A2", "TPK", 
                                  "OGDH", "PDHA1", "BCKDHA", "SLC25A19", "TKT", 
                                  "SLC19A3")

TBsignatures[["B3_Niacin"]] = c("SLC22A3", "HERPUD1", "CETP")

TBsignatures[["B5_Pantothenic Acid"]] = c("TMEM132D", "IGFBP7", "RTEL1", "RTEL1-TNFRSF6B", "LY96", 
                                          "TSNAX", "SNRPD2P2", "MIR302F", "RNU6-408P", "RECK", 
                                          "OR2AM1P", "RNA5SP517", "MIR891A", "SYT16", "LINC02198", 
                                          "PIK3R1", "ERG", "LINC02731")

TBsignatures[["B6_Pyridoxine"]] = c("ALPL", "NBPF3", "ADCYAP1R1", "ZNF608", "SMAD3-DT", "OSBPL3", 
                                    "CYCS", "PYROXD2", "CACNA1C")
TBsignatures[["B9_Folate"]] = c("MTHFR", "DHFR", "SLC19A1", "FOLR1", "SHMT1", "MTHFD1")

TBsignatures[["B12_Cobalamin"]] = c("E2F1", "RPS3", "RPS6", "RPL5", "RPL7", "RPL11", "RPL13A", 
                                  "RPS14", "RPS27A", "RPL22")

TBsignatures[["Vitamin_C"]] = c("SLC23A1", "SLC23A3", "CHPT1", "BCAS3", "SNRPF", "RER1", "MAF", "GSTA5", 
                              "RGS14", "AKT1", "FADS1")

TBsignatures[["Vitamin_D"]] = c("CD24", "CDH1", "CDH11", "CDH3", "CDH4", "CLDN7", "CRB3", "CXADR", 
                                "DMKN", "DSC2", "DSP", "EPCAM", "EPPK1", "F11R", "GJB2", "GJB3", 
                                "MAL2", "MARVELD2", "MPZL2", "MUC1", "OCLN", "PATJ", "PCDH7", "PKP3", 
                                "PMEPA1", "PPL", "SCEL", "SFN", "SH3YL1", "SHROOM3", "SYTL1", 
                                "TACSTD2", "TMEM30B", "TSPAN1", "TSPAN15")

TBsignatures[["Vitamin_E"]] = c("CD36", "SRD5A1", "GCLC", "PKC", "CCNB2", "CDC2", "CDC6", "TGM2", "MMP1", 
                              "MIR122", "MIR125B")

TBsignatures[["Vitamin_K"]] = c("BTBD8", "NUP62", "MPV17", "CHMP1B", "SAMD4B", "PEG13", "SNHG1", 
                              "GTF3C2", "HSPA14", "PRMT1", "MRPL33")

TBsignatures[["Iron_Homeostasis"]] = c("TFRC", "FTH1", "FTL", "SLC11A2", "ACO1", "IREB2")

TBsignatures[["Omega_3_FattyAcids"]] = c("FADS1", "FADS2", "ELOVL2", "ELOVL5", "SCD")

TBsignatures[["Carnitine_Mito"]] = c("TMLHE", "SLC22A5", "CPT1A", "CPT2", "CRAT")
