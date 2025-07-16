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
  "Vitamin_D_Aigner_35" = c("CD24", "CDH1", "CDH11", "CDH3", "CDH4", "CLDN7", "CRB3", "CXADR", 
                         "DMKN", "DSC2", "DSP", "EPCAM", "EPPK1", "F11R", "GJB2", "GJB3", 
                         "MAL2", "MARVELD2", "MPZL2", "MUC1", "OCLN", "PATJ", "PCDH7", "PKP3", 
                         "PMEPA1", "PPL", "SCEL", "SFN", "SH3YL1", "SHROOM3", "SYTL1", 
                         "TACSTD2", "TMEM30B", "TSPAN1", "TSPAN15"), 
)




###Run this code if you want to add the NutSignatures to TBsignatures

TBsignatures[["A_Retinol"]] = c("Lrat", "Tnni3k", "Dancr", "Suv39h2", "Apol7c", "Tnni3", 
                                "Fndc5", "St6galnac5", "Neil3", "Cxcr6", "Slc47a1", 
                                "Myoz2", "Cenpi", "Slc17a7", "Lrch2", "Shcbp1l", 
                                "Cox6a2", "Hspb7", "Grik4", "Gm7511")

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

TBsignatures[["Vitamin_D"]] = c("CD24", "CDH1", "CDH11", "CDH3", "CDH4", "CLDN7", "CRB3", "CXADR", 
                                "DMKN", "DSC2", "DSP", "EPCAM", "EPPK1", "F11R", "GJB2", "GJB3", 
                                "MAL2", "MARVELD2", "MPZL2", "MUC1", "OCLN", "PATJ", "PCDH7", "PKP3", 
                                "PMEPA1", "PPL", "SCEL", "SFN", "SH3YL1", "SHROOM3", "SYTL1", 
                                "TACSTD2", "TMEM30B", "TSPAN1", "TSPAN15")

