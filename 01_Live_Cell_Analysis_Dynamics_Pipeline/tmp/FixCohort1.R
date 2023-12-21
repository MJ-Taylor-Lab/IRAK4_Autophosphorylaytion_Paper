
# Table <- fread("/Users/u_deliz/Desktop/PhasePortraitAnalysis/Big Essential.csv.gz")

library(dtplyr)

Table2 <-
  Table %>% 
  filter(
    COHORT != "MyD88 IRAK4 Inhibitor_PF06650833",
    COHORT != "MyD88 IRAK4 Kinase_Dead",
    COHORT != "MyD88 IRAK4 LV",
    COHORT != "MyD88 IRAK4 DMSO",
    IMAGE != "20210922  well7E_10nM 085_MyD88_IRAK4KI_DMSOcntrl_002",
    IMAGE != "20210922  well8F_5nM 085_MyD88_RAK4inhibPF06650833_001",
    IMAGE != "20210922  well8G_5nM 085_MyD88_RAK4inhibPF06650833_001"
  ) %>% 
  mutate(
    COHORT = ifelse(IMAGE == "20190619 EL4 MyD88gfp homozygous IRAK4mScar 001", "MyD88 IRAK4", COHORT),
    COHORT = ifelse(IMAGE == "20190620 CL854af 5nm MyD88gfp IRAK4mScar 001", "MyD88 IRAK4", COHORT),
    COMPLEMENTARY_PROTEIN_1 = ifelse(IMAGE == "20190620 CL854af 5nm MyD88gfp IRAK4mScar 001" & COMPLEMENTARY_PROTEIN_1 == "IRAK1", "IRAK4", COMPLEMENTARY_PROTEIN_1),
    COMPLEMENTARY_PROTEIN_1 = ifelse(IMAGE == "20190702 cl85 4AF il1646 20nm 001" & COMPLEMENTARY_PROTEIN_1 == "IRAK1", "IRAK4", COMPLEMENTARY_PROTEIN_1),
    COMPLEMENTARY_PROTEIN_1 = ifelse(IMAGE == "20190718 0,1nm CL82 3D IRAK1mScar myd88gfp well2 001" & COMPLEMENTARY_PROTEIN_1 == "IRAK4", "IRAK1", COMPLEMENTARY_PROTEIN_1),
    COMPLEMENTARY_PROTEIN_1 = ifelse(IMAGE == "20190718 0,25nm cl082 3d myd88gfp well2 " & COMPLEMENTARY_PROTEIN_1 == "IRAK4", "IRAK1", COMPLEMENTARY_PROTEIN_1),
    COMPLEMENTARY_PROTEIN_1 = ifelse(IMAGE == "20190718 0,25nm cl082 3d myd88gfp well3 " & COMPLEMENTARY_PROTEIN_1 == "IRAK4", "IRAK1", COMPLEMENTARY_PROTEIN_1),
    COMPLEMENTARY_PROTEIN_1 = ifelse(IMAGE == "20190619 EL4 MyD88gfp homozygous IRAK4mScar 001" & COMPLEMENTARY_PROTEIN_1 == "IRAK1", "IRAK4", COMPLEMENTARY_PROTEIN_1),
    
    PROTEIN = ifelse(IMAGE == "20190620 CL854af 5nm MyD88gfp IRAK4mScar 001" & PROTEIN == "IRAK1", "IRAK4", PROTEIN),
    PROTEIN = ifelse(IMAGE == "20190702 cl85 4AF il1646 20nm 001" & PROTEIN == "IRAK1", "IRAK4", PROTEIN),
    PROTEIN = ifelse(IMAGE == "20190718 0,1nm CL82 3D IRAK1mScar myd88gfp well2 001" & PROTEIN == "IRAK4", "IRAK1", PROTEIN),
    PROTEIN = ifelse(IMAGE == "20190718 0,25nm cl082 3d myd88gfp well2 " & PROTEIN == "IRAK4", "IRAK1", PROTEIN),
    PROTEIN = ifelse(IMAGE == "20190718 0,25nm cl082 3d myd88gfp well3 " & PROTEIN == "IRAK4", "IRAK1", PROTEIN),
    PROTEIN = ifelse(IMAGE == "20190619 EL4 MyD88gfp homozygous IRAK4mScar 001" & PROTEIN == "IRAK1", "IRAK4", PROTEIN)
  )

Table2 <- as_tibble(Table2)

Summary <-
  Table2 %>% 
  filter(
    PROTEIN == "MyD88",
    COMPLEMENTARY_PROTEIN_1 != ""
  ) %>% 
  select(
    COHORT,
    LIGAND_DENSITY_CAT,
    IMAGE,
    PROTEIN,
    # IMAGE
  ) %>% 
  distinct() %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    PROTEIN,
    # IMAGE
  ) %>% 
  summarize(
    N = n(),
    # TEST = grepl(PROTEIN, COHORT)
  ) 

