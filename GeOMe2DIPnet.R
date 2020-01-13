#To install packages. Unfortunately seqinr Does not read fasta names with spaces. Have to use Biostrings
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biostrings")
# install.packages("XLConnect")
# install.packages("XLConnect")
# install.packages("dplyr")



geome2dipnet <- function(GeomePath=GeomePath, fastaPath=fastaPath) {
  ##############################
  #Adds sequence data and other required fields to a dataframe that is compatible
  # for analysis with the DIPnet_Stats_Functions. Reads from geOMe xls file in 
  # IraMoana template.
  ##############################
  require(Biostrings,quietly = T)
  require(XLConnect)
  # require(seqinr)
  require(dplyr)
 
  #Read xls and fasta
  geometest1 <- XLConnect::readWorksheetFromFile(file = GeomePath,sheet="Samples") #Can safely ignore 'Illegal reflective acces' warning. Java bug to be fixed with next release.
  
  #Make dataframe of fasta file
  fastaFile <- readDNAStringSet(fastaPath)
  seq_names = names(fastaFile)
  sequences = paste(fastaFile)
  fastadf <- data.frame(seq_names, sequences)
  fastadf$seq_name <- as.character(fastadf$seq_name) #Force into character so that dplyr inner_join() will work silently
  #Remove white spaces
  geometest1$materialSampleID <- gsub('[() ]','',geometest1$materialSampleID)
  fastadf$seq_name <- gsub('[() ]','',fastadf$seq_name)
  fastadf$seq_name
  
  #Check for any seqs not in geOMe metadata sheet
  missing_metadata <- setdiff(fastadf$seq_name,geometest1$materialSampleID)
  
  #Add sequence data to GEOME dataframe
  geometest1 <- geometest1 %>%
    inner_join(fastadf,
               by = c('materialSampleID' = 'seq_name'))
  
  #Missing data check
  #Print how many seqs had no metadata in geOMe sheet
  if (length(missing_metadata) > 0) {
    warning((paste("Error reading files:", length(missing_metadata), " sequences from fasta file are missing metadata. Check that names match in GeomePath and fastaPath")),noBreaks. = F)
    print("The following seqs had no metadata")
    print(missing_metadata)
    }# End of missing data check
  
  ###Add locus info
  #Add dummy locus
  # geometest1$locus <- rep("locus",nrow(geometest1))
  
  #Add locus from associatedSequences column in geOMe df
  # geometest1$locus <- sub(".*MTDNA_(.*)_.*", "\\1", geometest1$associatedSequences)
  # r <- regexpr("MTDNA_(.*)_", geometest1$associatedSequences)
  # regmatches(geometest1$associatedSequences, r)
  
  #Add locus from fasta file name.
  marker_list <- c("18S","S7","TMO","GnRH","ND2","16S","CR","CO1","CO2","CYB","RAG","A68")
  fasta_name <- basename(fastaPath)
  for (m in marker_list) {
    if(grepl(m,fasta_name))  locus = m}
  #print if no locus
  ########TO DO#########
  # print(locus)

  
  #Add dummy column $sample with '1' for testing. Not interested in breaking up to ecoregions ect at this time.
  geometest1$sample <- rep(1,nrow(geometest1))
  
  #Create $Genus_species_locus column
  geometest1$Genus_species_locus <- paste(geometest1$scientificName, locus)
  geometest1$Genus_species_locus <- gsub(' ','_',geometest1$Genus_species_locus) #Remove white space, not sure if it affects dipnet functions.
  return(geometest1)
} #end geome2dipnet()



# Convert from GeOMe/Fasta format to DIPnet analysis compatible dataframe
GeomePath <- '/home/iggy/BioinformatIg/IraMoana/Tool/ExampleGeome/IraMoana_GEOME_metadata_template_181219.xlsx'
fastaPath <- '/home/iggy/BioinformatIg/IraMoana/Tool/ExampleGeome/LR031Cumming2016_MTDNA_CO1_oncnig.fasta'
geometest1 <- geome2dipnet(GeomePath = GeomePath, fastaPath = fastaPath)



####TO TEST DIPNET functions##########

setwd('/home/iggy/BioinformatIg/IraMoana/Tool/DipnetOUT')
# source("~/git_clones/popgenDB/DIPnet_Stats_Functions.R") # This one is broken due to updates in dependencies
source("~/GitRepos/IraMoana/Updated_DIPnet_Stats_Functions.R")


# Various DIPnet popgen analysis examples.
diffstats<-pairwise.structure.mtDNA.db(ipdb=geometest1, gdist = "PhiST", minseqs = 6, minsamps = 3, mintotalseqs = 0, nrep = 0, num.cores = 1, ABGD = F, regionalization = "sample")
divstats<-genetic.diversity.mtDNA.db(ipdb=geometest1, basic_diversity = T, sequence_diversity = T, coverage_calc = F, coverage_correction = F, minseqs = 6, minsamps = 3, mintotalseqs = 0, ABGD=F,regionalization = "sample", keep_all_gsls=F, mincoverage = 0.4, hill.number = 0)
hierstats<-hierarchical.structure.mtDNA.db(ipdb = geometest1 ,level1 = "sample",level2="ECOREGION",model="raw",nperm=1)

for(g in c("WC Theta","PhiST", "Jost D")){
  diffstats<-pairwise.structure.mtDNA.db(ipdb=geometest1, gdist = g, regionalization = "sample", minseqs = 6, minsamps= 3, mintotalseqs= 0, num.cores = 2)
  write.stats(diffstats,filename=file.path(paste("DIPnet_structure_060315_", g,".csv",sep="")),structure=T) # for an excel-readable csv. Ignore warnings. Note this function will not overwrite, it will append to existing files
  }



#OLD CODE
# # # Shorten locality column?
# 
# 
# #Add this other stuff from ipdb paper to see if the columns match up with GEOME
# # spatial_path<-"/home/iggy/git_clones/IPDB/ipdb_spatial.txt"
# # spatial2_path<-"/home/iggy/git_clones/IPDB/ipdb_HypothesisRegionAssigned.txt"
# # abgd_path<-"/home/iggy/git_clones/IPDB/ABGD_group_assignment_JCdistance_defaults.txt"
# # spatial<-read.table(spatial_path, header=T, sep="\t",stringsAsFactors = F, na.strings=c("NA"," ",""), quote="")
# # spatial2<-read.table(spatial2_path, header=T,sep="\t", stringsAsFactors = F, na.strings=c("NA"," ",""), quote="")
# # abgd<-read.table(abgd_path, header=T, sep="\t", stringsAsFactors = F)
# # ipdb<-join(ipdb,spatial, by = "IPDB_ID",type = "left")
# # ipdb<-join(ipdb,spatial2[,c(2,18:24)], by = "IPDB_ID", type = "left")
# # ipdb<-join(ipdb,abgd[,c(1,3)], by = "IPDB_ID",type = "left")
# # ipdb<-ipdb[ipdb$IPDB_ID %in% drops == FALSE, ]
# # 
# # ipdb_full_col <- colnames(ipdb)
# ##Load ibdp dataset for testing comparisons
# ipdb <- read.table('/home/iggy/BioinformatIg/IraMoana/Tool/ExampleDIPnet/ipdb_sub1.txt',sep="\t",header=T,stringsAsFactors = F,quote="", na.strings=c("NA"," ","")) 
# #Add dummy sample column
# ipdb$sample <- rep(1,nrow(ipdb))
# 
# 
# ##To list columns that overlap between Geome and IBD
# Geome_col <- colnames(geometest1)
# sort(Geome_col)
# 
# ipdb_col <- colnames(ipdb)
# sort(ipdb_col)
# 
# identical_col <- intersect(Geome_col,ipdb_col)
# length(identical_col)
# sort(identical_col)
# 
# need_col <- setdiff(ipdb_col,Geome_col)
# length(need_col)
# sort(need_col)
# 
# avail_col <- setdiff(Geome_col,ipdb_col)
# length(avail_col)
# sort(avail_col)

# #Using seqinr to load fasta into df
# sequence <- c(sequence=c(unlist(getSequence(fastatest1,as.string=T))))
# seq_name <- c(seq_name=c(getName(fastatest1)))
# as.character(seq_name)
# fastadf <- data.frame(cbind(seq_name,sequence))

