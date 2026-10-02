# IEEE Access 2021 mtDNA Panel (Continuity Dataset)

The 30 complete mitochondrial genomes used as the benchmark panel in the original
IEEE Access manuscript (`archive/manuscripts/ieee_access/genome_ieeeaccess_submitted.pdf`,
Section IV-B and Appendix A), re-downloaded from NCBI on 2026-08-18 for continuity
with the BMC rebuild.

30 species, 27 genera, 13 families, 6 orders (all Mammalia). The original paper
simulated shotgun reads directly from these sequences (sequencing depth 1X+); it
did not use real sequencing reads. This panel is kept
for continuity/reproducibility comparison against the original results, and as a
source of "divergent related reference" sequences within genera/families.

## Provenance

- Source: `archive/manuscripts/ieee_access/genome_ieeeaccess_submitted.pdf`, Table 1 / Appendix A.
- Downloaded via NCBI E-utilities (`efetch`, `db=nuccore`, `rettype=fasta`) on 2026-08-18.
- Every downloaded FASTA header was checked against the expected species name before
  being accepted; all 30 matched (one taxonomic synonym: NCBI still lists
  `Cephalopachus bancanus` under its older name `Tarsius bancanus`, same accession/sequence).

## Manifest

| # | Species | Genus | Family | Order | Accession | File |
|---|---|---|---|---|---|---|
| 1 | Lagothrix lagotricha | Lagothrix | Atelidae | Primates | NC_021951 | `1_Lagothrix_lagotricha_NC_021951.fasta` |
| 2 | Ateles belzebuth | Ateles | Atelidae | Primates | NC_019800 | `2_Ateles_belzebuth_NC_019800.fasta` |
| 3 | Alouatta seniculus | Alouatta | Atelidae | Primates | NC_027825 | `3_Alouatta_seniculus_NC_027825.fasta` |
| 4 | Aotus azarai | Aotus | Aotidae | Primates | NC_021939 | `4_Aotus_azarai_NC_021939.fasta` |
| 5 | Callithrix penicillata | Callithrix | Callitrichidae | Primates | NC_030788 | `5_Callithrix_penicillata_NC_030788.fasta` |
| 6 | Callimico goeldii | Callimico | Callitrichidae | Primates | NC_024628 | `6_Callimico_goeldii_NC_024628.fasta` |
| 7 | Leontopithecus rosalia | Leontopithecus | Callitrichidae | Primates | NC_021952 | `7_Leontopithecus_rosalia_NC_021952.fasta` |
| 8 | Saguinus oedipus | Saguinus | Callitrichidae | Primates | NC_021960 | `8_Saguinus_oedipus_NC_021960.fasta` |
| 9 | Saimiri boliviensis | Saimiri | Cebidae | Primates | NC_021966 | `9_Saimiri_boliviensis_NC_021966.fasta` |
| 10 | Pygathrix cinerea | Pygathrix | Cercopithecidae | Primates | NC_018063 | `10_Pygathrix_cinerea_NC_018063.fasta` |
| 11 | Procolobus verus | Procolobus | Cercopithecidae | Primates | NC_020666 | `11_Procolobus_verus_NC_020666.fasta` |
| 12 | Rhinopithecus strykeri | Rhinopithecus | Cercopithecidae | Primates | NC_018059 | `12_Rhinopithecus_strykeri_NC_018059.fasta` |
| 13 | Rhinopithecus brelichi | Rhinopithecus | Cercopithecidae | Primates | NC_018057 | `13_Rhinopithecus_brelichi_NC_018057.fasta` |
| 14 | Simias concolor | Simias | Cercopithecidae | Primates | NC_020667 | `14_Simias_concolor_NC_020667.fasta` |
| 15 | Cercocebus torquatus | Cercocebus | Cercopithecidae | Primates | NC_023964 | `15_Cercocebus_torquatus_NC_023964.fasta` |
| 16 | Allenopithecus nigroviridis | Allenopithecus | Cercopithecidae | Primates | NC_023965 | `16_Allenopithecus_nigroviridis_NC_023965.fasta` |
| 17 | Macaca silenus | Macaca | Cercopithecidae | Primates | NC_025221 | `17_Macaca_silenus_NC_025221.fasta` |
| 18 | Macaca arctoides | Macaca | Cercopithecidae | Primates | NC_025201 | `18_Macaca_arctoides_NC_025201.fasta` |
| 19 | Macaca tonkeana | Macaca | Cercopithecidae | Primates | NC_025222 | `19_Macaca_tonkeana_NC_025222.fasta` |
| 20 | Papio papio | Papio | Cercopithecidae | Primates | NC_020009 | `20_Papio_papio_NC_020009.fasta` |
| 21 | Homo sapiens | Homo | Hominidae | Primates | NC_012920 | `21_Homo_sapiens_NC_012920.fasta` |
| 22 | Gorilla gorilla | Gorilla | Hominidae | Primates | NC_001645 | `22_Gorilla_gorilla_NC_001645.fasta` |
| 23 | Cephalopachus bancanus | Cephalopachus | Tarsiidae | Primates | NC_002811 | `23_Cephalopachus_bancanus_NC_002811.fasta` |
| 24 | Varecia variegata | Varecia | Lemuridae | Primates | NC_012773 | `24_Varecia_variegata_NC_012773.fasta` |
| 25 | Lemur catta | Lemur | Lemuridae | Primates | NC_004025 | `25_Lemur_catta_NC_004025.fasta` |
| 26 | Galeopterus variegatus | Galeopterus | Cynocephalidae | Dermoptera | NC_004031 | `26_Galeopterus_variegatus_NC_004031.fasta` |
| 27 | Oryctolagus cuniculus | Oryctolagus | Leporidae | Lagomorpha | NC_001913 | `27_Oryctolagus_cuniculus_NC_001913.fasta` |
| 28 | Mus musculus | Mus | Muroidea | Rodentia | KF937876 | `28_Mus_musculus_KF937876.fasta` |
| 29 | Tupaia belangeri | Tupaia | Tupaiidae | Scandentia | AF217811 | `29_Tupaia_belangeri_AF217811.fasta` |
| 30 | Sphaerias blanfordi | Sphaerias | Pteropodidae | Chiroptera | NC_046933 | `30_Sphaerias_blanfordi_NC_046933.fasta` |
