# ╔═══════════════════════════════════════════════════════════════════╗
# ║                          rename_tree.R                            ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : diversity-cereal                                 ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2024-11-04                                       ║
# ║ Last Modified  : 2024-11-04                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/diversity-cereal  ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## I/O handling

input <- "/home/sergio/projects/trees_diversity/edited_fus_its2_msa_mafft_linsi.nex.con.tre"
output <- "/home/sergio/projects/trees_diversity/renamed_fus_its2_msa_mafft_linsi.nex.con.tre"

## Replacements

replacements <- c(
  "Fus.ASV.CON.1_1-165_F.sp_1e4911d" = "ASV_CON1_F.sp",
  "Fus.ASV.CON.2_1-147_F.sp_243070c" = "ASV_CON2_F.sp",
  "Fus.ASV.CON.3_1-149_F.venenatum_" = "ASV_CON3_F.venenatum",
  "Fus.ASV.CON.4_1-164_F.sp_ad8600b" = "ASV_CON4_F.sp",
  "Fus.ASV.CON.5_1-162_F.algeriense" = "ASV_CON5_F.algeriense",
  "Fus.ASV.CON.6_1-167_F.sp_60dcdd7" = "ASV_CON6_F.sp",
  "Fus.ASV.ECO.1_1-163_F.tricinctum" = "ASV_ECO1_F.tricinctum",
  "Fus.ASV.ECO.2_1-159_F.sp_29f950c" = "ASV_ECO2_F.sp",
  "AF414967.1_129-278_F.poae_FpII" = "AF414967.1_F.poae_FpII",
  "AF414966.1_129-280_F.poae_FpIV" = "AF414966.1_F.poae_FpIV",
  "MH865933.1_129-294_F.acuminatum_" = "MH865933.1_F.acuminatum_CBS131258",
  "HM068326.1_129-294_F.acuminatum_" = "HM068326.1_F.acuminatum_NRRL54218",
  "HM068323.1_129-294_F.acuminatum_" = "HM068323.1_F.acuminatum_NRRL54215",
  "MH858189.1_130-293_F.avenaceum_C" = "MH858189.1_F.avenaceum_CBS387.62",
  "MH860648.1_129-292_F.avenaceum_C" = "MH860648.1_F.avenaceum_CBS121.73",
  "OL832320.1_129-292_F.avenaceum_N" = "OL832320.1_F.avenaceum_NRRL54939",
  "JX162372.1_129-288_F.brachygibbo" = "JX162372.1_F.brachygibbosum_CBS131252",
  "AF006341.1_129-279_F.cerealis_NR" = "AF006341.1_F.cerealis_NRRL25805",
  "AF006340.1_129-279_F.cerealis_NR" = "AF006340.1_F.cerealis_NRRL25491",
  "AF006342.1_129-279_F.culmorum_NR" = "AF006342.1_F.culmorum_NRRL25475",
  "MH862856.1_129-279_F.culmorum_CB" = "MH862856.1_F.culmorum_CBS110262",
  "EU926221.1_129-292_F.domesticum_" = "EU926221.1_F.domesticum_CBS102407",
  "MH862468.1_129-279_F.equiseti_CB" = "MH862468.1_F.equiseti_CBS307.94",
  "GQ505733.1_129-279_F.equiseti_NR" = "GQ505733.1_F.equiseti_NRRL36136",
  "GQ505742.1_129-279_F.equiseti_NR" = "GQ505742.1_F.equiseti_NRRL36466",
  "MH855481.1_129-279_F.equiseti_CB" = "MH855481.1_F.equiseti_CBS185.34",
  "MZ890564.1_129-292_F.flocciferum" = "MZ890564.1_F.flocciferum_NL19-97008",
  "MZ890561.1_129-292_F.flocciferum" = "MZ890561.1_F.flocciferum_JW267001",
  "MW827608.1_129-292_F.fujikuroi_C" = "MW827608.1_F.fujikuroi_CBS221.76",
  "MH857023.1_129-292_F.fujikuroi_C" = "MH857023.1_F.fujikuroi_CBS257.52",
  "MH865804.1_129-279_F.graminearum" = "MH865804.1_F.graminearum_CBS130609",
  "MH865931.1_129-279_F.graminearum" = "MH865931.1_F.graminearum_CBS131262",
  "GQ505704.1_129-279_F.incarnatum_" = "GQ505704.1_F.incarnatum_NRRL32866",
  "GQ505705.1_129-279_F.incarnatum_" = "GQ505705.1_F.incarnatum_NRRL32867",
  "NR_121214.1_129-278_F.langsethia" = "NR_121214.1_F.langsethiae_CBS113234",
  "AB587022.1_129-278_F.langethiae_" = "AB587022.1_F.langethiae_FRC_T-0992",
  "AB587023.1_129-278_F.langethiae_" = "AB587023.1_F.langethiae_FRC_T-1000",
  "OK482905.1_129-288_F.lateritium_" = "OK482905.1_F.lateritium_NCP10",
  "KC453998.1_129-288_F.lateritium_" = "KC453998.1_F.lateritium_KCTC46029",
  "AB586998.1_129-289_F.merismoides" = "AB586998.1_F.merismoides_MAFF_236504",
  "MH866031.1_129-278_F.oxysporum_C" = "MH866031.1_F.oxysporum_CBS132473",
  "MH866025.1_129-278_F.oxysporum_C" = "MH866025.1_F.oxysporum_CBS132481",
  "KF255448.1_129-278_F.oxysporum_C" = "KF255448.1_F.oxysporum_CBS133023",
  "MH858407.1_129-278_F.poae_CBS_17" = "MH858407.1_F.poae_CBS_177.64",
  "MH859029.1_129-278_F.poae_CBS446" = "MH859029.1_F.poae_CBS446.67",
  "EF453150.1_129-292_F.proliferatu" = "EF453150.1_F.proliferatum_NRRL43667",
  "JX162388.1_129-291_F.proliferatu" = "JX162388.1_F.proliferatum_CBS131785",
  "MT435064.1_129-292_F.redolens_NR" = "MT435064.1_F.redolens_NRRL25600",
  "JX162380.1_129-298_F.solani_CBS1" = "JX162380.1_F.solani_CBS131775",
  "KF255440.1_129-298_F.solani_CBS1" = "KF255440.1_F.solani_CBS132898",
  "MH855269.1_129-278_F.sporotrichi" = "MH855269.1_F.sporotrichioides_CBS180.32",
  "AF006348.1_129-278_F.sporotrichi" = "AF006348.1_F.sporotrichioides_NRRL25479",
  "AF006347.1_129-278_F.sporotrichi" = "AF006347.1_F.sporotrichioides_NRRL25474",
  "AB587026.1_129-278_F.sporotrichi" = "AB587026.1_F.sporotrichioides_CBS119839",
  "HQ165919.1_129-280_F.subglutinan" = "HQ165919.1_F.subglutinans_PUF016",
  "AY898251.1_129-280_F.subglutinan" = "AY898251.1_F.subglutinans_ATCC38016",
  "MH856607.1_129-293_F.tricinctum_" = "MH856607.1_F.tricinctum_CBS253.50",
  "MH862424.1_129-294_F.tricinctum_" = "MH862424.1_F.tricinctum_CBS393.93",
  "M068317.1_129-294_F.tricinctum_N" = "M068317.1_F.tricinctum_NRRL25481",
  "F453174.1_129-280_F.verticillioi" = "F453174.1_F.verticillioides_NRRL43697",
  "MH861171.1_129-280_F.verticillio" = "MH861171.1_F.verticillioides_CBS576.78",
  "NR_121320.1_129-279_F.asiaticum_" = "NR_121320.1_F.asiaticum_NRRL26156",
  "DQ459836.1_129-279_F.asiaticum_N" = "DQ459836.1_F.asiaticum_NRRL28720",
  "DQ459834.1_129-279_F.asiaticum_N" = "DQ459834.1_F.asiaticum_NRRL13818",
  "DQ459841.1_129-279_F.meridionale" = "DQ459841.1_F.meridionale_NRRL29010",
  "DQ459840.1_129-279_F.meridionale" = "DQ459840.1_F.meridionale_NRRL28436",
  "NR_121203.1_129-279_F.boothii_NR" = "NR_121203.1_F.boothii_NRRL29011",
  "MH860690.1_129-279_F.boothii_CBS" = "MH860690.1_F.boothii_CBS316.73",
  "HQ322366.1_129-281_Outgroup_Ceph" = "HQ322366.1_Outgroup_Cephalosporium"
)

tree <- readLines(input, warn = FALSE)
tree <- paste(tree, collapse = "\n")

for (old_name in names(replacements)) {
  new_name <- replacements[old_name]
  tree <- gsub(old_name, new_name, tree, fixed = TRUE)
}

writeLines(tree, output)

