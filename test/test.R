load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
resultdir = "G:/lab/DCcomboNet/Rpackage/tryout_result/"
drugcandidate = data.frame(drug = "Sorafenib", drug = "DB00398")
manual_input = FALSE
drugtarget = NULL #read.table(paste0(load_dir,'data/drugtarget.csv'), sep =  ',',header = TRUE, stringsAsFactors = FALSE)
druggene = NULL

DComboNet(load_dir = load_dir,
          resultdir = resultdir,
          model = "L1", # To choose level one model
          manual_input = FALSE, # To shield manually input drug name
          drugcandidate = drugcandidate,
          drugnetWeight = TRUE, # Confirm if drug network should be weighted
          featuretype = 'integrated_score')

#2) runing level two model(L2)

DComboNet(load_dir = load_dir,
          resultdir = resultdir,
          model = "L2", # Choose level two model
          manual_input = FALSE, # You can use manual input function (see level one model example)
          drugcandidate = drugcandidate,
          drugnetWeight = TRUE,
          featuretype = 'integrated_score',
          dataset = "92742",
          cellline = "HEPG2",
          treatment_time = 6,
          # the absolute value of fold-change between drug-treated group and control group will be above 0.5
          foldchange_DEG = 0.5,
          # the p-value from the significent test (t.test) between drug-treated group and control group will be below 0.05
          pvalue_DEG = 0.05,
          foldchange_DEP = 0,
          pvalue_DEP = 0.05)

drugseed = 'Sorafenib'
drugcandidate = 'Vorinostat'
drugtarget = data.frame(Drug = c(rep(drugseed,10),rep(drugcandidate,5)),
                        target = c("BRAF", "FGFR1", "FLT1", "FLT3", "FLT4", "KDR", "KIT", "PDGFRB", "RAF1","RET", "HDAC1", "HDAC2", "HDAC3", "HDAC6", "HDAC8"))
cellline = 'HEPG2'
network_extract(drugseed = drugseed, drugcandidate = drugcandidate,cellline = cellline,
                drugtarget = drugtarget, generank_filter = 0.02, pathwayrank_filter = 0.1,
                model = "L2", load_dir = load_dir, resultdir = resultdir)
network_visualization(drugseed = drugseed, drugcandidate = drugcandidate,cellline = cellline,
                      drugtarget = drugtarget, model = "L2", load_dir = load_dir,
                      resultdir = resultdir)

outputdir = 'G:/lab/Projects/p1_DComboNet/Rpackage/input_data/test/OCILY3_data/'
load_dir = 'G:/lab/Projects/p1_DComboNet/Rpackage/OCL_LY3/'
resultdir = "G:/lab/Projects/p1_DComboNet/DComboNet_improve/OCI_LY3_tryout_results/"
dir.create(resultdir)
drugpair = read.csv(paste0(outputdir,'drugpair_oci_ly3.csv'),stringsAsFactors = F)
drugcandidate <-  read.csv(paste0(outputdir,"druglist_oci_ly3.csv"), sep = ',', header = T, stringsAsFactors = F)
CDK_FP  <- read.csv(paste0(outputdir,'structure_sim/fingerprint.csv'), sep = ',',header = T, stringsAsFactors = F)
pubchemFP <- read.csv(paste0(outputdir,'structure_sim/PubChem.csv'), sep = ',',header = T, stringsAsFactors = F)
MACCS_FP <- read.csv(paste0(outputdir,'structure_sim/MACCF.csv'), sep = ',',header = T, stringsAsFactors = F)

drugnetWeight = TRUE
featuretype = 'integrated_score'
drugtarget = read.csv(paste0(outputdir,'drug_target_oci_ly3.csv'), sep = ',',header = T, stringsAsFactors = F)
druggene = NULL
# dt_lib = read.csv(paste0(outputdir,'drug_gene/drugtarget.csv'))

cellline = 'OCILY3'
drugDEG = read.table('G:/lab/Projects/p1_DComboNet/Rpackage/input_data/test/OCILY3_data/drug_DEG_12hrs.txt',sep='\t',header=T,stringsAsFactors = F)
drugDEP = read.table('G:/lab/Projects/p1_DComboNet/Rpackage/input_data/test/OCILY3_data/drug_DEP_12hrs.txt',sep='\t',header=T,stringsAsFactors = F)

cancergene = read.table(paste0(outputdir,'OCILY3_genelist.txt'),sep='\t',header=T,stringsAsFactors = F)
dtweight = 2
dgweight = 1
dDEGweight = 1
model = 'L2'

drugnetWeight = 'T'
featuretype = 'integrated_score2'

DComboNet(load_dir = load_dir,
          resultdir = resultdir,
          model = 'L2',
          drugcandidate=drugcandidate,
          CDK_FP = CDK_FP,
          pubchemFP = pubchemFP,
          MACCS_FP = MACCS_FP,
          drugnetWeight = TRUE,
          featuretype = featuretype,
          drugtarget = drugtarget,
          dataset = NULL,
          cellline = 'OCILY3',
          treatment_time = NULL,
          drugDEG = drugDEG,
          drugDEP = drugDEP,
          cancergene = cancergene,
          dtweight = 2,
          dgweight = 1,
          dDEGweight = 1,
          dgpweight = 1,
          dDEpweight = 1)

network_extract(drugseed = 'Sorafenib', drugcandidate = 'Vorinostat',
                drugtarget = drugtarget, generank_filter = 0.01, pathwayrank_filter = 0.1,
                model = "L2", load_dir = load_dir, resultdir = resultdir)
network_visualization(drugseed = 'Sorafenib', drugcandidate = 'Vorinostat',
                      drugtarget = drugtarget, model = "L2", load_dir = load_dir,
                      resultdir = resultdir)
