# Analyses for Neutral Theory and the Somatic Evolution of Cancer
# By: Vincent L. Cannataro and Jeffrey P. Townsend


library(reshape2)
library(ggplot2)

# import table from Miyata T, Miyazawa S, Yasunaga T. 1979. Two types of amino acid substitutions in protein evolution. J. Mol. Evol. 12:219–236.

miyata.table <- matrix(data = NA,nrow = 20,ncol = 20)
colnames(miyata.table) <- c("CYS","PRO","ALA","GLY","SER","THR","GLN","GLU","ASN","ASP","HIS","LYS","ARG","VAL","LEU","ILE","MET","PHE","TYR","TRP")
rownames(miyata.table) <- c("CYS","PRO","ALA","GLY","SER","THR","GLN","GLU","ASN","ASP","HIS","LYS","ARG","VAL","LEU","ILE","MET","PHE","TYR","TRP")

miyata.table[,"PRO"] <- c(1.33,rep(NA,19))
miyata.table[,"ALA"] <- c(1.39,0.06,rep(NA,18))
miyata.table[,"GLY"] <- c(2.22,0.97,0.91,rep(NA,17))
miyata.table[,"SER"] <- c(1.84,0.56,0.51,0.85,rep(NA,16))
miyata.table[,"THR"] <- c(1.45,0.87,0.90,1.70,0.89,rep(NA,15))
miyata.table[,"GLN"] <- c(2.48,1.92,1.92,2.48,1.65,1.12,rep(NA,14))
miyata.table[,"GLU"] <- c(3.26,2.48,2.46,2.78,2.06,1.83,0.84,rep(NA,13))
miyata.table[,"ASN"] <- c(2.83,1.80,1.78,1.96,1.31,1.40,0.99,0.85,rep(NA,12))
miyata.table[,"ASP"] <- c(3.48,2.40,2.37,2.37,1.87,2.05,1.47,0.90,0.65,rep(NA,11))
miyata.table[,"HIS"] <- c(2.56,2.15,2.17,2.78,1.94,1.32,0.32,0.96,1.29,1.72,rep(NA,10))
miyata.table[,"LYS"] <- c(3.27,2.94,2.96,3.54,2.71,2.10,1.06,1.14,1.84,2.05,0.79,rep(NA,9))
miyata.table[,"ARG"] <- c(3.06,2.90,2.92,3.58,2.74,2.03,1.13,1.45,2.04,2.34,0.82,0.40,rep(NA,8))
miyata.table[,"VAL"] <- c(0.86,1.79,1.85,2.76,2.15,1.42,2.13,2.97,2.76,3.40,2.11,2.70,2.43,rep(NA,7))
miyata.table[,"LEU"] <- c(1.65,2.70,2.76,3.67,3.04,2.25,2.70,3.53,3.49,4.10,2.59,2.98,2.62,0.91,rep(NA,6))
miyata.table[,"ILE"] <- c(1.63,2.62,2.69,3.60,2.95,2.14,2.57,3.39,3.37,3.98,2.45,2.84,2.49,0.85,0.14,rep(NA,5))
miyata.table[,"MET"] <- c(1.46,2.36,2.42,3.34,2.67,1.86,2.30,3.13,3.08,3.69,2.19,2.63,2.29,0.62,0.41,0.29,rep(NA,4))
miyata.table[,"PHE"] <- c(2.24,3.17,3.23,4.14,3.45,2.60,2.81,3.59,3.70,4.27,2.63,2.85,2.47,1.43,0.63,0.61,0.82,rep(NA,3))
miyata.table[,"TYR"] <- c(2.38,3.12,3.18,4.08,3.33,2.45,2.48,3.22,3.42,3.95,2.27,2.42,2.02,1.52,0.94,0.86,0.93,0.48,NA,NA)
miyata.table[,"TRP"] <- c(3.34,4.17,4.23,5.13,4.38,3.50,3.42,4.08,4.39,4.88,3.16,3.11,2.72,2.51,1.73,1.72,1.89,1.11,1.06,NA)

# import table from McLachlan AD. 1972. Repeating sequences and gene duplication in proteins. J. Mol. Biol. 64:417–437.
McLachlan.1972.vec.order <- c("VAL","LEU","ILE","MET","PHE","TRP","TYR","GLY","ALA","PRO","SER","THR","CYS","HIS","ARG","LYS","GLN","GLU","ASN","ASP")

McLachlan.matrix <- matrix(NA,nrow=length(McLachlan.1972.vec.order),ncol=length(McLachlan.1972.vec.order))

rownames(McLachlan.matrix) <- McLachlan.1972.vec.order
colnames(McLachlan.matrix) <- McLachlan.1972.vec.order

McLachlan.matrix["VAL",] <- c(662,rep(NA,19))
McLachlan.matrix["LEU",] <- c(82,545,rep(NA,18))
McLachlan.matrix["ILE",] <- c(102,72,411,rep(NA,17))
McLachlan.matrix["MET",] <- c(28,44,19,227,rep(NA,16))
McLachlan.matrix["PHE",] <- c(27,43,18,12,249,rep(NA,15))
McLachlan.matrix["TRP",] <- c(5,7,5,1,10,85,rep(NA,14))
McLachlan.matrix["TYR",] <- c(19,18,14,5,49,16,247,rep(NA,13))
McLachlan.matrix["GLY",] <- c(35,18,13,5,2,4,6,535,rep(NA,12))
McLachlan.matrix["ALA",] <- c(79,43,31,20,17,4,13,78,878,rep(NA,11))
McLachlan.matrix["PRO",] <- c(25,12,11,3,5,0,3,26,52,331,rep(NA,10))
McLachlan.matrix["SER",] <- c(50,38,31,17,18,9,24,73,135,46,1039,rep(NA,9))
McLachlan.matrix["THR",] <- c(46,35,33,20,12,5,13,37,87,27,144,751,rep(NA,8))
McLachlan.matrix["CYS",] <- c(6,2,2,3,0,1,2,5,6,2,12,8,68,rep(NA,7))
McLachlan.matrix["HIS",] <- c(12,11,7,6,10,2,10,11,19,8,25,22,3,249,rep(NA,6))
McLachlan.matrix["ARG",] <- c(18,17,7,3,5,3,7,20,24,12,48,25,3,18,359,rep(NA,5))
McLachlan.matrix["LYS",] <- c(41,34,13,10,6,4,13,46,68,27,87,67,2,25,65,701,rep(NA,4))
McLachlan.matrix["GLN",] <- c(21,21,6,9,1,3,5,21,43,16,53,37,1,14,30,46,450,rep(NA,3))
McLachlan.matrix["GLU",] <- c(27,12,10,6,4,3,9,41,68,28,62,47,1,8,20,52,56,568,rep(NA,2))
McLachlan.matrix["ASN",] <- c(18,18,10,8,4,1,11,46,44,8,98,43,4,20,21,53,33,35,534,NA)
McLachlan.matrix["ASP",] <- c(21,18,7,8,6,2,10,48,47,20,69,43,5,18,13,42,34,79,59,549)



# import amino acid codon sequences and names 
translations <- read.csv(file = "translations.csv",header = T,stringsAsFactors = F)

# Import data from our selection intensity analysis (https://www.biorxiv.org/content/early/2018/02/27/229724) 
load("combined_selection_output.RData")


# Clean selection data such that it does not have STOP codons, data from non-coding regions, etc., to match
# analysis by McLachlan and then Kimura
combined_all_data$distance <- NA

unique(toupper(translations$AA_short)) %in% colnames(miyata.table)
translations$AA_short <- toupper(translations$AA_short)
if(length(which(combined_all_data$AA_Change=="*"))>0){
  combined_all_data <- combined_all_data[-which(combined_all_data$AA_Change=="*"),]
}

if(length(which(combined_all_data$AA_Ref=="*"))>0){
  combined_all_data <- combined_all_data[-which(combined_all_data$AA_Ref=="*"),] 
}

if(length(which(is.na(combined_all_data$AA_Ref)))>0){
  combined_all_data <- combined_all_data[-which(is.na(combined_all_data$AA_Ref)),]
}


combined_all_data$distance[which(combined_all_data$AA_Ref==combined_all_data$AA_Change)] <- 0

nucleotides <- c("A","T","C","G")



# Make a matrix of the potential AA mutations given all nucleotide mutations
transition.matrix <- matrix(data=0, nrow=nrow(translations), ncol=nrow(translations))
rownames(transition.matrix) <- paste(translations$Nucs,translations$AA_short,translations$AA_letter,sep=":")
colnames(transition.matrix) <- paste(translations$Nucs,translations$AA_short,translations$AA_letter,sep=":")



for(i in 1:nrow(transition.matrix)){
  this.codon <- translations[i,"Nucs"]
  this.codon.split <- unlist(strsplit(this.codon,split = ""))
  
  for(j in 1:3){
    these.changes <- nucleotides[which(!(nucleotides %in% this.codon.split[j]))]
    
    for(k in 1:3){
      new.codon <- this.codon.split
      new.codon[j] <- these.changes[k]
      
      transition.matrix[i,which(translations[,"Nucs"]== paste(new.codon,collapse = ""))] <- transition.matrix[i,which(translations[,"Nucs"]== paste(new.codon,collapse = ""))]+1
      
    }
    
  }
  
}



# Make a data structure that stores frequency and distance 
miyata.melt <- melt(miyata.table)
miyata.melt$value[which(miyata.melt$Var1 == miyata.melt$Var2)] <- 0
miyata.melt <- miyata.melt[-which(is.na(miyata.melt$value)),]
miyata.melt$potential_changes <- NA
for(i in 1:nrow(miyata.melt)){
  if(miyata.melt$Var1[i]==miyata.melt$Var2[i]){
    miyata.melt$potential_changes[i] <- sum(transition.matrix[which(translations[,"AA_short"]==miyata.melt$Var1[i]),which(translations[,"AA_short"]==miyata.melt$Var2[i])])/2
  }else{
    miyata.melt$potential_changes[i] <- sum(transition.matrix[which(translations[,"AA_short"]==miyata.melt$Var1[i]),which(translations[,"AA_short"]==miyata.melt$Var2[i])])
  }
}

# Only want substitutions 
miyata.melt <- miyata.melt[-which(miyata.melt$potential_changes==0),]

miyata.melt$proportion_changes_expected <- miyata.melt$potential_changes/sum(miyata.melt$potential_changes)

miyata.melt <- cbind(miyata.melt,matrix(data = NA,nrow = nrow(miyata.melt),ncol = (length(unique(combined_all_data$tumor_type))+1)))
colnames(miyata.melt) <- c("AA_Ref","AA_Change","Distance","Potential_changes","Potential_changes_proportion",unique(combined_all_data$tumor_type),"Between species")

miyata.melt <- melt(miyata.melt, id.vars = c("AA_Ref","AA_Change","Distance","Potential_changes","Potential_changes_proportion"))
tail(miyata.melt)
colnames(miyata.melt) <- c(colnames(miyata.melt)[1:5],"Tumor_type","Frequency")

# Fill in the data from the selection analysis

# takes a while, so saved output: 
load("miyata_melt_initialized.Rdata")
# for(i in 1:nrow(miyata.melt)){
#   
#   if(miyata.melt$Tumor_type[i] != "Between species"){
#     miyata.melt$Frequency[i] <- sum(combined_all_data$freq[which(combined_all_data$tumor_type==miyata.melt$Tumor_type[i] & 
#                                                                    combined_all_data$AA_Ref== translations$AA_letter[which(translations$AA_short== miyata.melt$AA_Ref[i])[1]] & 
#                                                                    combined_all_data$AA_Change== translations$AA_letter[which(translations$AA_short== miyata.melt$AA_Change[i])[1]])],combined_all_data$freq[which(combined_all_data$tumor_type==miyata.melt$Tumor_type[i] & 
#                                                                                                                                                                                                                      combined_all_data$AA_Ref== translations$AA_letter[which(translations$AA_short== miyata.melt$AA_Change[i])[1]] & 
#                                                                                                                                                                                                                      combined_all_data$AA_Change== translations$AA_letter[which(translations$AA_short== miyata.melt$AA_Ref[i])[1]])])
#   }else{
#     if(miyata.melt$AA_Ref[i]==miyata.melt$AA_Change[i]){
#       miyata.melt$Frequency[i] <- sum(McLachlan.matrix[which(rownames(McLachlan.matrix)==miyata.melt$AA_Ref[i]),which(colnames(McLachlan.matrix)==miyata.melt$AA_Change[i])],McLachlan.matrix[which(rownames(McLachlan.matrix)==miyata.melt$AA_Change[i]),which(colnames(McLachlan.matrix)==miyata.melt$AA_Ref[i])],na.rm = T)/2
#     }else{
#       miyata.melt$Frequency[i] <- sum(McLachlan.matrix[which(rownames(McLachlan.matrix)==miyata.melt$AA_Ref[i]),which(colnames(McLachlan.matrix)==miyata.melt$AA_Change[i])],McLachlan.matrix[which(rownames(McLachlan.matrix)==miyata.melt$AA_Change[i]),which(colnames(McLachlan.matrix)==miyata.melt$AA_Ref[i])],na.rm = T)  
#     }
#     
#   }
#   print(i/nrow(miyata.melt))
# }


# tail(miyata.melt)


miyata.melt$proportion_observed <- NA

for(i in 1:length(unique(miyata.melt$Tumor_type))){
  miyata.melt$proportion_observed[which(miyata.melt$Tumor_type== unique(miyata.melt$Tumor_type)[i])] <- miyata.melt$Frequency[which(miyata.melt$Tumor_type== unique(miyata.melt$Tumor_type)[i])]/sum(miyata.melt$Frequency[which(miyata.melt$Tumor_type== unique(miyata.melt$Tumor_type)[i])])
}






miyata.melt$relative_freq <- miyata.melt$proportion_observed/miyata.melt$Potential_changes_proportion



McLachlan.matrix.diag <- McLachlan.matrix[lower.tri(McLachlan.matrix,diag=F)]

sum(McLachlan.matrix.diag,na.rm = T)*2

positions_i <- c(394,321,250,115,154,53,159,353,493,203,502,370,89,122,173,382,218,255,282,304)
McLachlan.ni <- c(394,321,250,115,154,53,159,353,493,203,502,370,89,122,173,382,218,255,282,304)

multiplied.mat <- positions_i %o% positions_i

sum(multiplied.mat[lower.tri(multiplied.mat,diag = F)])*2


McLachlan.alpha <- (sum(McLachlan.matrix.diag,na.rm = T)*2)/(sum(multiplied.mat[lower.tri(multiplied.mat,diag = F)])*2)

# 9438/25290138


# Example for LEU and ILE
72/(321*250* McLachlan.alpha)







###
# the read method ---- 
###




# need to find the number of amino acids in each of 
length(unique(combined_all_data$Gene))

# load table with total number of each amino acid in the sequenced exome from our selection intensity analysis
load("combined_amino_acid_tally.RData")

total_positions <- colSums(amino.acid.tally_combined,na.rm = T)

multiplied.mat <- total_positions %o% total_positions

N2 <- sum(multiplied.mat) # (sum(multiplied.mat[lower.tri(multiplied.mat,diag = F)])*2)+ sum(diag(multiplied.mat))
N2.vec <- c(rep(N2,23),25290138)

tally_array <- array(data=NA,dim = c(length(McLachlan.1972.vec.order),length(McLachlan.1972.vec.order),length(unique(combined_all_data$tumor_type))+1,2),dimnames = list(McLachlan.1972.vec.order,McLachlan.1972.vec.order,c(unique(combined_all_data$tumor_type),"Between species"),c("Frequency","Relative_frequency")))

N1.vec <- rep(NA,length(unique(combined_all_data$tumor_type))+1)

# dimnames(tally_array)[[4]]

# fill in the tally array with frequency data from the tumor types
for(i in 1:length(unique(miyata.melt$Tumor_type))){
  this.miyata <- miyata.melt[which(miyata.melt$Tumor_type==unique(miyata.melt$Tumor_type)[i]),]
  for(j in 1:nrow(this.miyata)){
    tally_array[as.character(this.miyata$AA_Ref[j]),as.character(this.miyata$AA_Change[j]),i,1] <- this.miyata$Frequency[j] 
    tally_array[as.character(this.miyata$AA_Change[j]),as.character(this.miyata$AA_Ref[j]),i,1] <- this.miyata$Frequency[j] 
  }
}

# calculating alpha
for(i in 1:length(unique(miyata.melt$Tumor_type))){
  if(i < 24){
    this.diag <- tally_array[,,i,1][upper.tri(tally_array[,,i,1],diag=F)]
    
    N1.vec[i] <-  (sum(this.diag,na.rm = T)*2)+sum(diag(tally_array[,,i,1]),na.rm = T)
  }else{
    this.diag <- tally_array[,,i,1][upper.tri(tally_array[,,i,1],diag=F)]
    N1.vec[i] <- sum(this.diag,na.rm = T)*2
  }
}
# do not have full N1 because deleted substitutions not possible with more than 1 SNV. 
N1.vec[24] <- 9438

alpha.vec <- N1.vec/N2.vec
# filling in the relative frequency 
for(z in 1:23){
  for(i in 1:20){
    for(j in 1:20){
      tally_array[i,j,z,2] <- tally_array[i,j,z,1]/(total_positions[i]*total_positions[j]*alpha.vec[z])
    }
  }
}

for(i in 1:20){
  for(j in 1:20){
    tally_array[i,j,24,2] <- tally_array[i,j,24,1]/(McLachlan.ni[i]*McLachlan.ni[j]*alpha.vec[24])
  }
}

tally_array[,,1,2]


mean(tally_array[,,1,2][upper.tri(tally_array[,,1,2],diag = F)],na.rm=T)

# filling this relative frequency back into the miyata.melt dataframe 
miyata.melt$McLachlan_rel_freq <- NA

for(i in 1:nrow(miyata.melt)){
  miyata.melt$McLachlan_rel_freq[i] <- tally_array[as.character(miyata.melt$AA_Ref[i]),as.character(miyata.melt$AA_Change[i]),as.character(miyata.melt$Tumor_type[i]),2]  
}

miyata.melt$McLachlan_rel_freq[which(miyata.melt$Tumor_type=="Between species" & miyata.melt$AA_Ref==miyata.melt$AA_Change)] <- NA




###
# normalizing by trinucleotide rates per tumor ----- 
###





# flips the nucleotide to the other strand

flip.function <- function(nucleotide){
  if(nucleotide=="A"){
    return("T")
  }
  if(nucleotide=="T"){
    return("A")
  }
  if(nucleotide=="C"){
    return("G")
  }
  if(nucleotide=="G"){
    return("C")
  }
  if(nucleotide=="N"){
    return("N")
  }
}



nucs <- c("A","C","G","T")

load("combined_expanded_codon_tally.RData")

sum(colSums(expanded.codon.tally_combined,na.rm = T))

expanded.codon.sums <- colSums(expanded.codon.tally_combined,na.rm = T)

expanded.codon <- apply(expand.grid(c("A","C","G","T"),translations[,"Nucs"],c("A","C","G","T")),1,paste,collapse=".")

expanded.codon.df <- expand.grid(c("A","C","G","T"),translations[,"Nucs"],c("A","C","G","T"))

load("combined_selection_output.RData")

tumor.list <- unique(combined_all_data$tumor_type)

trinuc.mutation.df <- expand.grid(rep(translations[,"Nucs"],9),tumor.list)
colnames(trinuc.mutation.df) <- c("codon_1","tumor_type")

trinuc.mutation.df$AA_1 <- NA
trinuc.mutation.df$codon_2 <- NA
trinuc.mutation.df$AA_2 <- NA
trinuc.mutation.df$mutation_rate <- NA

all.muts.function <- function(original.codon){
  original.codon.split <- strsplit(original.codon,"")
  new.codons <- NULL
  for(k in 1:3){
    for(l in 1:3){
      this.change <- original.codon.split[[1]][k]
      this.change <- nucs[!(nucs %in% this.change)][l]
      this.new.codon <- original.codon.split[[1]]
      this.new.codon[k] <- this.change
      new.codons <- c(new.codons,paste(this.new.codon,collapse = ""))
    }
  }
  return(new.codons)
}


proportion.finder <- function(trinuc,mutant){
  this.trinuc.split <- strsplit(trinuc,"")
  if(this.trinuc.split[[1]][2] %in% c("C","T")){
    return(trinuc.mutation_data[which(trinuc.mutation_data$Upstream==this.trinuc.split[[1]][1] &
                                        trinuc.mutation_data$mutated_from==this.trinuc.split[[1]][2] &
                                        trinuc.mutation_data$Downstream==this.trinuc.split[[1]][3] &
                                        trinuc.mutation_data$mutated_to==mutant),"proportion"])
  }else{
    return(trinuc.mutation_data[which(trinuc.mutation_data$Upstream==flip.function(this.trinuc.split[[1]][3]) &
                                        trinuc.mutation_data$mutated_from==flip.function(this.trinuc.split[[1]][2]) &
                                        trinuc.mutation_data$Downstream==flip.function(this.trinuc.split[[1]][1]) &
                                        trinuc.mutation_data$mutated_to==flip.function(mutant)),"proportion"])
    
  }
  
}


rownames(translations) <- translations[,"Nucs"]

# load in trinucleotide contexts

for(tumor in 1:length(tumor.list)){
  load(paste("trinuc_signature_data/trinuc_data_",tumor.list[tumor],".RData",sep=""))
  # load(paste("~/Documents/Selection_analysis/",tumor.list[tumor],"/trinuc_output/trinuc_data_",tumor.list[tumor],".RData",sep=""))
  
  for(codon in 1:nrow(translations)){
    this.position.list <- which(trinuc.mutation.df$codon_1==translations[codon,"Nucs"] & trinuc.mutation.df$tumor_type==tumor.list[tumor])
    
    trinuc.mutation.df[this.position.list,"AA_1"] <- translations[codon,"AA_short"]
    trinuc.mutation.df[this.position.list,"codon_2"] <- all.muts.function(original.codon = translations[codon,"Nucs"])
    trinuc.mutation.df[this.position.list,"AA_2"] <- translations[trinuc.mutation.df[this.position.list,"codon_2"],"AA_short"]
    
    
    for(position in 1:length(this.position.list)){
      if(position <=3){
        
        
        this.codon2 <- trinuc.mutation.df[this.position.list[position],"codon_2"]
        these.sums <- NULL
        for(this.upstream in 1:4){
          these.sums <- c(these.sums,sum(expanded.codon.sums[which(expanded.codon.df$Var2== translations[codon,"Nucs"] & expanded.codon.df$Var1 == nucs[this.upstream])]))
        }
        # these.sums <- these.sums/sum(these.sums)
        
        these.rates <- NULL
        
        for(this.rate in 1:4){
          these.rates <- c(these.rates,proportion.finder(trinuc = paste(c(nucs[this.rate],
                                                                          strsplit(translations[codon,"Nucs"],"")[[1]][1:2]),collapse = ""),mutant = strsplit(trinuc.mutation.df[this.position.list[position],"codon_2"],"")[[1]][1]))
        }
        
        trinuc.mutation.df[this.position.list[position],"mutation_rate"] <- sum(these.sums*these.rates)
        
      }
      if(position>=4 & position <=6){
        
        
        this.codon2 <- trinuc.mutation.df[this.position.list[position],"codon_2"]
        
        trinuc.mutation.df[this.position.list[position],"mutation_rate"] <- proportion.finder(trinuc = as.character(trinuc.mutation.df[this.position.list[position],"codon_1"]),
                                                                                              mutant = strsplit(trinuc.mutation.df[this.position.list[position],"codon_2"],"")[[1]][2]) * sum(expanded.codon.sums[which(expanded.codon.df$Var2==as.character(trinuc.mutation.df[this.position.list[position],"codon_1"]))])
        
        
      }
      if(position>=7){
        
        
        
        this.codon2 <- trinuc.mutation.df[this.position.list[position],"codon_2"]
        these.sums <- NULL
        for(this.downstream in 1:4){
          these.sums <- c(these.sums,sum(expanded.codon.sums[which(expanded.codon.df$Var2== translations[codon,"Nucs"] & expanded.codon.df$Var3 == nucs[this.downstream])]))
        }
        # these.sums <- these.sums/sum(these.sums)
        
        these.rates <- NULL
        
        for(this.rate in 1:4){
          these.rates <- c(these.rates,proportion.finder(trinuc = paste(c(strsplit(translations[codon,"Nucs"],"")[[1]][2:3],nucs[this.rate]),collapse = ""),mutant = strsplit(trinuc.mutation.df[this.position.list[position],"codon_2"],"")[[1]][3]))
        }
        
        trinuc.mutation.df[this.position.list[position],"mutation_rate"] <- sum(these.sums*these.rates)
        
        
      }
      
    }
    
    
  }
  
}

miyata.melt$McLachlan_rel_freq_prop <- NA
for(i in 1:length(tumor.list)){
  these.tumors <- which(miyata.melt$Tumor_type==tumor.list[i]) 
  miyata.melt$McLachlan_rel_freq_prop[these.tumors] <- miyata.melt$McLachlan_rel_freq[these.tumors]/sum(miyata.melt$McLachlan_rel_freq[these.tumors],na.rm = T)
}

head(miyata.melt)


miyata.melt$trinuc_rate <- NA

for(i in 1:nrow(miyata.melt)){
  miyata.melt$trinuc_rate[i] <- sum(trinuc.mutation.df$mutation_rate[which(
    ((trinuc.mutation.df$AA_1==as.character(miyata.melt$AA_Ref[i]) & 
        trinuc.mutation.df$AA_2==as.character(miyata.melt$AA_Change[i])) |
       (trinuc.mutation.df$AA_1==as.character(miyata.melt$AA_Change[i]) & 
          trinuc.mutation.df$AA_2==as.character(miyata.melt$AA_Ref[i])))  & trinuc.mutation.df$tumor_type == as.character(miyata.melt$Tumor_type[i]))] )
}

head(miyata.melt)
miyata.melt$trinuc_rate_prop <- NA
for(i in 1:length(tumor.list)){
  these.tumors <- which(miyata.melt$Tumor_type==tumor.list[i]) 
  miyata.melt$trinuc_rate_prop[these.tumors] <- miyata.melt$trinuc_rate[these.tumors]/sum(miyata.melt$trinuc_rate[these.tumors],na.rm = T)
}



# miyata.melt$McLachlan_rel_freq_prop/miyata.melt$trinuc_rate_prop
tail(miyata.melt)


# miyata.melt

miyata.melt$trinuc_normalized <- miyata.melt$proportion_observed/miyata.melt$trinuc_rate_prop

miyata.melt$trinuc_normalized[which(miyata.melt$Tumor_type=="Between species")] <- miyata.melt$McLachlan_rel_freq[which(miyata.melt$Tumor_type=="Between species")]


library(ggplot2)

head(miyata.melt)
unique(miyata.melt$Tumor_type)
library(magrittr)
library(dplyr)

max.trinuc_norm <- as_tibble(miyata.melt[-which(miyata.melt$Tumor_type=="Between species"),]) %>%
  group_by(Tumor_type) %>%
  summarize(max = max(trinuc_normalized)) %>%
  arrange(desc(max))


miyata.melt$Tumor_type <- factor(miyata.melt$Tumor_type,levels = c(as.character(max.trinuc_norm$Tumor_type),"Between species"),labels = c(as.character(max.trinuc_norm$Tumor_type)[1:3],
                                                                                                                                          expression(HPV^{"+"}~HNSC),
                                                                                                                                          as.character(max.trinuc_norm$Tumor_type)[5:19],
                                                                                                                                          expression(HPV^{"−"}~HNSC),
                                                                                                                                          as.character(max.trinuc_norm$Tumor_type)[21:23],"Between~species"))

frequency.plot_trinuc <- ggplot(data = subset(miyata.melt, Frequency>0), aes(x=Distance, y= trinuc_normalized)) + geom_point(alpha=0.8,size=.5) + geom_smooth(method="glm",method.args=list(family=gaussian(link="log")),se = F) + facet_wrap(~Tumor_type, scales = "free", ncol=4, labeller = "label_parsed") + labs(y="Normalized frequency",x="Molecular distance") + theme_classic() 
frequency.plot_trinuc
ggsave(filename = "normalized_freq_vs_distace_trinuc.png",plot = frequency.plot_trinuc,height = 8,width = 6)


# density plots ------


combined_all_data$synonymous <-  combined_all_data$AA_Ref==combined_all_data$AA_Change

combined_all_data.noNA <- combined_all_data[which(!is.na(combined_all_data$AA_Ref)),]


library(scales)
# density.plot <- ggplot(data = subset(combined_all_data.noNA,freq>0), aes(x=gamma_epistasis)) + geom_density(aes(fill=synonymous),position = "stack") + facet_wrap(~tumor_type,scales = "free") + scale_x_continuous(trans = "log2") + theme_classic()
mysqrt_trans <- function() {
  domain <- c(0, Inf)
  transform <- function(x) x^(1/3)
  range <- transform(domain)
  trans_new("mysqrt", 
            transform = transform,
            inverse = function(x) squish(x, range=range)^3,
            domain = domain)
}



density.plot <- ggplot(data = subset(combined_all_data.noNA,freq>0), aes(x=gamma_epistasis)) + geom_density(aes(fill=synonymous),position = "stack") + facet_wrap(~tumor_type,scales = "free",ncol=4) + 
  scale_x_continuous(trans="mysqrt",expand = c(0,0)) + labs(x="Selection intensity", y="Density")+ 
  scale_y_continuous(expand=c(0,0)) + theme_classic()

density.plot

# ggsave(filename = "density_plot_cuberoot.png",plot = density.plot,height = 8,width = 6)



library(scales)

# adding names
combined_all_data.noNA$name <- NA

combined_all_data.noNA$name <- paste(combined_all_data.noNA$Gene,rep(" ",nrow(combined_all_data.noNA)),combined_all_data.noNA$AA_Ref,combined_all_data.noNA$AA_Pos,combined_all_data.noNA$AA_Change,sep="")


# BLCA density ---- 

this.tumor <- "BLCA"
these.data <- subset(combined_all_data.noNA,tumor_type==this.tumor & freq>1)
this.size <- 10
sub.1.y <- 0.15
sub.1.x <- max(these.data$gamma_epistasis[which(these.data$freq>1)])-2.2e5
text.1.df <- data.frame(sub.1.x,sub.1.y)

arrow.pointer.x <-  max(these.data$gamma_epistasis[which(these.data$freq>1)])-1e3
arrow.pointer.y <- 0.01
arrow.1 <- data.frame(sub.1.x,sub.1.y,arrow.pointer.x,arrow.pointer.y)


these.breaks <- c(1,1e4,1e5)
these.labels <- c(expression(10^0),expression(10^4),expression(10^5))

these.limits <- c(1,max(combined_all_data.noNA$gamma_epistasis[which(combined_all_data.noNA$tumor_type==this.tumor & combined_all_data.noNA$freq>1)]))
this.max.selection <- which.max(combined_all_data.noNA$gamma_epistasis[which(combined_all_data.noNA$tumor_type==this.tumor & 
                                                                               combined_all_data.noNA$freq>1)])
this.max.label <- combined_all_data.noNA$name[which(combined_all_data.noNA$tumor_type==this.tumor & combined_all_data.noNA$freq>1 & combined_all_data.noNA$gamma_epistasis==max(combined_all_data.noNA$gamma_epistasis[which(combined_all_data.noNA$tumor_type==this.tumor& combined_all_data.noNA$freq>1)]))]


max.freq <- these.data[order(these.data$freq,decreasing = T),][1,]
label.2.coords <- data.frame(x1=max.freq$gamma_epistasis,x2 = max.freq$gamma_epistasis,
                             y1=0.01,y2=.1,lab=max.freq$name)

density.blca <- ggplot(data = subset(combined_all_data.noNA,freq>0 & tumor_type==this.tumor), aes(x=gamma_epistasis)) +
  geom_density(aes(fill=synonymous),position = "stack") + 
  geom_point(data = subset(combined_all_data.noNA,freq>1 & tumor_type==this.tumor), 
             aes(x=max(gamma_epistasis),
                 y=c(0)),
             shape="|",size=10) +
  geom_text(data=text.1.df,aes(x=sub.1.x,
                               y=sub.1.y),label=this.max.label,size=this.size,vjust=0) +
  geom_segment(data=arrow.1,aes(x=sub.1.x,y=sub.1.y,
                                xend=arrow.pointer.x,yend=arrow.pointer.y),
               arrow = arrow(length = unit(0.03, "npc"))) + 
  geom_point(data = label.2.coords, 
             aes(x=x1,
                 y=c(0)),
             shape="|",size=10) + 
  geom_text(data=label.2.coords,aes(x=x1,
                                    y=y2),label=label.2.coords$lab,size=this.size,vjust=0) +
  geom_segment(data=label.2.coords,aes(x=x2,y=y2,
                                       xend=x1,yend=y1),
               arrow = arrow(length = unit(0.03, "npc"))) + 
  scale_x_continuous(breaks=these.breaks,labels=these.labels,trans="mysqrt",expand = c(0,0)) + 
  coord_cartesian(xlim = these.limits) + guides(fill=FALSE)+
  scale_y_continuous(expand=c(0,0)) + 
  theme(legend.position="none") +
  theme_classic() + labs(title="BLCA")
density.blca

ggsave(filename = "figures/BLCA_density.png",plot = density.blca)




library(cowplot)

# test.plot <- plot_grid(density.blca,density.blca,align = 'h')
# save_plot(filename = "figures/test.png",plot = test.plot,base_height = 8.24,base_width = 9.25*2)

# function to create different density plots ---- 
density_plotter_function <- function(main_df, 
                                     this_tumor,
                                     main_size=10,
                                     sub_1 = "PIK3CA E545K",
                                     sub_1_lab_y = 0.1,
                                     sub_1_lab_x = 1,
                                     sub_1_arrow_x = 1,
                                     sub_1_arrow_y = 0.01,
                                     sub_1_lab_lab = sub_1,
                                     sub_1_arrow_alpha=1,
                                     sub_2,
                                     sub_2_lab_y = 0.1,
                                     sub_2_lab_x = 1,
                                     sub_2_arrow_x = 1,
                                     sub_2_arrow_y = 0.01,
                                     sub_2_lab_lab = sub_2,
                                     sub_2_arrow_alpha = 1,
                                     point_shape = "|",
                                     point_size = 10,
                                     x_axis_breaks = c(1,1e4,1e5),
                                     x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5)),
                                     y_axis_breaks = c(0,.1,.2),
                                     y_axis_labels = c(0,.1,.2),
                                     y_axis_limits = c(0,.25),
                                     figure_title=this_tumor,
                                     extender_percent=0.01,
                                     minimum_maximum=1.1e6,
                                     margins=c(2,2,2,2)
){
  # recurrent mutations in this tumor
  data_for_this_tumor_recur <- subset(main_df, tumor_type==this_tumor & freq>1)
  
  # data for substitution 1 
  sub_1_data <- subset(data_for_this_tumor_recur, name==sub_1)
  message(paste("Selection intensity of substitution 1:",round(sub_1_data$gamma_epistasis,2)))
  sub_1_df <- data.frame(sub_1_lab_x,
                         sub_1_lab_y,
                         sub_1_arrow_x=sub_1_arrow_x,
                         sub_1_arrow_y,
                         sub_1_lab_lab)
  
  # data for substitution 2 
  sub_2_data <- subset(data_for_this_tumor_recur, name==sub_2)
  message(paste("Selection intensity of substitution 2:",round(sub_2_data$gamma_epistasis,2)))
  sub_2_df <- data.frame(sub_2_lab_x,
                         sub_2_lab_y,
                         sub_2_arrow_x=sub_2_arrow_x,
                         sub_2_arrow_y,
                         sub_2_lab_lab)
  
  
  # create the plot 
  this_density_plot <- ggplot(data = subset(main_df, tumor_type==this_tumor), 
                              aes(x=gamma_epistasis)) + 
    geom_density(aes(fill=synonymous),position = "stack") + 
    
    # first substitution 
    geom_point(data = sub_1_data, 
               aes(x=gamma_epistasis,
                   y=c(0)),
               shape=point_shape,size=point_size) + 
    geom_text(data=sub_1_df,
              aes(x=sub_1_lab_x,
                  y=sub_1_lab_y),
              label=sub_1_lab_lab,
              size=main_size*(5/14),
              vjust=0) +
    geom_segment(data=sub_1_df,
                 aes(x=sub_1_lab_x,
                     y=sub_1_lab_y,
                     xend=sub_1_arrow_x,
                     yend=sub_1_arrow_y),
                 arrow = arrow(length = unit(0.05, "npc")),
                 alpha=sub_1_arrow_alpha) + 
    
    
    # second substitution 
    geom_point(data = sub_2_data, 
               aes(x=gamma_epistasis,
                   y=c(0)),
               shape=point_shape,size=point_size) + 
    geom_text(data=sub_2_df,
              aes(x=sub_2_lab_x,
                  y=sub_2_lab_y),
              label=sub_2_lab_lab,
              size=main_size*(5/14),
              vjust=0) +
    geom_segment(data=sub_2_df,
                 aes(x=sub_2_lab_x,
                     y=sub_2_lab_y,
                     xend=sub_2_arrow_x,
                     yend=sub_2_arrow_y),
                 arrow = arrow(length = unit(0.05, "npc")),
                 alpha=sub_2_arrow_alpha) + 
    # title in the plot
    
    # annotate("text",label=this_tumor,x=(max(data_for_this_tumor_recur$gamma_epistasis)+(max(data_for_this_tumor_recur$gamma_epistasis)*extender_percent))/2,y=)
    
    # scaling, etc. 
    scale_x_continuous(breaks=x_axis_breaks,
                       labels=x_axis_labels,
                       trans="mysqrt",
                       expand = c(0,0)) + 
    coord_cartesian(xlim = c(1,ifelse(minimum_maximum>max(data_for_this_tumor_recur$gamma_epistasis),minimum_maximum+(minimum_maximum*extender_percent),
                                      (max(data_for_this_tumor_recur$gamma_epistasis)+(max(data_for_this_tumor_recur$gamma_epistasis)*extender_percent))))) + 
    guides(fill=FALSE)+
    scale_y_continuous(expand=c(0,0), 
                       limits=y_axis_limits,
                       breaks=y_axis_breaks,
                       labels = y_axis_labels) + 
    theme(  axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text = element_text(size=main_size),
            text = element_text(size = main_size),
            plot.title = element_text(size=main_size,face = "plain"),
            plot.margin=unit(c(0,0,0,0), "mm")) #+
  # theme_classic() + 
  # labs(title=figure_title) 
  
  
  return(this_density_plot)
}



# test_plot <- density_plotter_function(main_df = combined_all_data.noNA,this_tumor = "LUAD",sub_1 = "KRAS G12C")
# 
# 
# test_plot <- density_plotter_function(main_df = combined_all_data.noNA,this_tumor = "HNSC_HPVneg",sub_1 = "PIK3CA E545K",figure_title = "HPV negative HNSCC",sub_1_lab_x = 5e4,sub_2 = "TP53 R283P",sub_2_lab_x = 1e5,sub_2_lab_y = .2,main_size = 10,sub_2_arrow_y = .04,sub_1_arrow_y = 0.04,extender_percent = 0.02,sub_2_arrow_x = 358905.26-40000,sub_1_arrow_x = 13434.61-1000)
# 
# 
# library(cowplot)
# 
# test_plot_combine <- plot_grid(test_plot,test_plot,test_plot,test_plot,test_plot,test_plot,ncol=3,align = 'hv')
# 
# save_plot(filename = "figures/test_combined.png",base_width = 6.5,plot = test_plot_combine,dpi=600)






# building plots 

# common parameters 

common.size <-7

# BLCA ----
BLCA_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "BLCA",
                                         main_size = common.size,
                                         sub_1 = "FBXW7 R425G",
                                         sub_2 = "HRAS G13R",
                                         sub_1_lab_y = .2,
                                         sub_1_lab_x = 4e5,
                                         sub_1_arrow_x = 372257.26,
                                         sub_1_arrow_y = 0.07,
                                         sub_2_lab_y = .1,
                                         sub_2_lab_x = 7e4,
                                         sub_2_arrow_x = 106564.92,
                                         sub_2_arrow_y = .07,
                                         extender_percent = 0.06,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),
                                         y_axis_limits = c(0,.4),
                                         y_axis_breaks = c(0,0.1,0.2,0.3),
                                         y_axis_labels = c(0,0.1,0.2,0.3),
                                         point_size = 4)

BLCA_density <- ggdraw(BLCA_density) + draw_label("BLCA",x = .5,y=.85,size = common.size)

save_plot(filename = "figures/BLCA_density.png",plot = BLCA_density,base_height = 5/6,base_width = 6.5/4)

# BRCA ----
BRCA_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "BRCA",
                                         main_size = common.size,
                                         sub_1 = "PIK3CA H1047R",
                                         sub_2 = "TP53 Y220S",
                                         sub_1_lab_y = .1,
                                         sub_1_lab_x = 6e5,
                                         sub_1_arrow_x =226324.32,
                                         sub_1_arrow_y = 0.07,
                                         sub_2_lab_y = .2,
                                         sub_2_lab_x = 7e4,
                                         sub_2_arrow_x = 81082.88,
                                         sub_2_arrow_y = .08,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),
                                         x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),
                                         y_axis_limits = c(0,.4),
                                         y_axis_breaks = c(0,0.1,0.2,0.3),
                                         y_axis_labels = c(0,0.1,0.2,0.3),
                                         point_size = 4)

BRCA_density <- ggdraw(BRCA_density) + draw_label("BRCA",x = .5,y=.9,size = common.size)

save_plot(filename = "figures/BRCA_density.png",plot = BRCA_density,base_height = 5/6,base_width = 6.5/4)


# CESC ---- 

CESC_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "CESC",
                                         main_size = common.size,
                                         sub_1 = "FBXW7 R505G",
                                         sub_2 = "MED12 D23Y",
                                         sub_1_lab_y = .25,
                                         sub_1_lab_x = 4e5,
                                         sub_1_arrow_x =531100.12,
                                         sub_1_arrow_y = 0.08,
                                         sub_2_lab_y = .15,
                                         sub_2_lab_x = 6e4,
                                         sub_2_arrow_x = 68681.85,
                                         sub_2_arrow_y = .08,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),
                                         x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),extender_percent = 0.06,
                                         y_axis_limits = c(0,.4),
                                         y_axis_breaks = c(0,0.1,0.2,0.3),
                                         y_axis_labels = c(0,0.1,0.2,0.3),
                                         point_size = 4)

CESC_density <- ggdraw(CESC_density) + draw_label("CESC",x = .5,y=.9,size = common.size)

save_plot(filename = "figures/CESC_density.png",plot = CESC_density,base_height = 5/6,base_width = 6.5/4)


# COAD ----

COAD_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "COAD",
                                         main_size = common.size,
                                         sub_1 = "BRAF V600E",
                                         sub_2 = "TP53 A159P",
                                         sub_1_lab_y = .4,
                                         sub_1_lab_x = 4e5,
                                         sub_1_arrow_x =1200000.5,
                                         sub_1_arrow_y = 0.12,
                                         sub_2_lab_y = .25,
                                         sub_2_lab_x = 6e4,
                                         sub_2_arrow_x = 399823.27,
                                         sub_2_arrow_y = .12,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),
                                         x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),
                                         extender_percent = 0.02,
                                         y_axis_limits = c(0,.65),
                                         y_axis_breaks = c(0,0.2,0.4,0.6),
                                         y_axis_labels = c(0,0.2,0.4,0.6),
                                         point_size = 4)

COAD_density <- ggdraw(COAD_density) + draw_label("COAD",x = .5,y=.9,size = common.size)

save_plot(filename = "figures/COAD_density.png",plot = COAD_density,base_height = 5/6,base_width = 6.5/4)

# ESCA ---- 


ESCA_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "ESCA",
                                         main_size = common.size,
                                         sub_1 = "TP53 R175G",
                                         sub_2 = "PIK3CA H1047L",
                                         sub_1_lab_y = .13,
                                         sub_1_lab_x = 7e5,
                                         sub_1_arrow_x = 1707238.91,
                                         sub_1_arrow_y = 0.04,
                                         sub_2_lab_y = .08,
                                         sub_2_lab_x = 240937.09,
                                         sub_2_arrow_x = 240937.09,
                                         sub_2_arrow_y = .04,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),extender_percent = 0.02,
                                         y_axis_limits = c(0,.22),
                                         y_axis_breaks = c(0,0.1,0.2),
                                         y_axis_labels = c(0,0.1,0.2),
                                         point_size = 4)

ESCA_density <- ggdraw(ESCA_density) + draw_label("ESCA",x = .5,y=.9,size = common.size)

save_plot(filename = "figures/ESCA_density.png",plot = ESCA_density,base_height = 5/6,base_width = 6.5/4)


# GBM ----

GBM_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                        this_tumor = "GBM",
                                        main_size = common.size,
                                        sub_1 = "BRAF V600E",
                                        sub_2 = "TP53 Y220C",
                                        sub_1_lab_y = .05,
                                        sub_1_lab_x = 5e5,
                                        sub_1_arrow_x = 55377.811,
                                        sub_1_arrow_y = 0.025,
                                        sub_2_lab_y = .1,
                                        sub_2_lab_x = 5e4,
                                        sub_2_arrow_x =  32406.12,
                                        sub_2_arrow_y = .03,
                                        x_axis_breaks = c(1,1e4,1e5,1e6),x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),extender_percent = 0.06,sub_1_arrow_alpha = 1,sub_2_arrow_alpha = 1,
                                        y_axis_limits = c(0,.22),
                                        y_axis_breaks = c(0,0.1,0.2),
                                        y_axis_labels = c(0,0.1,0.2),
                                        point_size = 4)

GBM_density <- ggdraw(GBM_density) + draw_label("GBM",x = .5,y=.9,size = common.size)

save_plot(filename = "figures/GBM_density.png",plot = GBM_density,base_height = 5/6,base_width = 6.5/4)



# HNSC_HPVpos ----

HNSC_HPVpos_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                                this_tumor = "HNSC_HPVpos",
                                                main_size = common.size,
                                                sub_1 = "FBXW7 R505G",
                                                sub_2 = "FGFR3 S249C",
                                                sub_1_lab_y = .07,
                                                sub_1_lab_x = 4e5,
                                                sub_1_arrow_x = 1140442.02,
                                                sub_1_arrow_y = 0.025,
                                                sub_2_lab_y = .04,
                                                sub_2_lab_x = 2e5,
                                                sub_2_arrow_x =  145184.52,
                                                sub_2_arrow_y = .025,
                                                x_axis_breaks = c(1,1e4,1e5,1e6),
                                                x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),extender_percent = 0.02,
                                                sub_1_arrow_alpha = 1,
                                                sub_2_arrow_alpha = 1,
                                                y_axis_limits = c(0,.13),
                                                y_axis_breaks = c(0,0.1,0.2),
                                                y_axis_labels = c(0,0.1,0.2),
                                                point_size = 4)

HNSC_HPVpos_density <- ggdraw(HNSC_HPVpos_density) + draw_label(label = expression(HPV^{"+"}~HNSC),x = .55,y=.9,size = common.size)

save_plot(filename = "figures/HNSC_HPVpos_density.png",plot = HNSC_HPVpos_density,base_height = 5/6,base_width = 6.5/4)


# HNSC_HPVneg ----

HNSC_HPVneg_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                                this_tumor = "HNSC_HPVneg",
                                                main_size = common.size,
                                                sub_1 = "TP53 R283P",
                                                sub_2 = "TP53 E298*",
                                                sub_1_lab_y = .14,
                                                sub_1_lab_x = 5e5,
                                                sub_1_arrow_x = 358905.26,
                                                sub_1_arrow_y = 0.05,
                                                sub_2_lab_y = .1,
                                                sub_2_lab_x = 6e4,
                                                sub_2_arrow_x =  200186.43,
                                                sub_2_arrow_y = .04,
                                                x_axis_breaks = c(1,1e4,1e5,1e6),
                                                x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),
                                                extender_percent = 0.06,
                                                sub_1_arrow_alpha = 1,
                                                sub_2_arrow_alpha = 1,
                                                y_axis_limits = c(0,.3),
                                                y_axis_breaks = c(0,0.1,0.2),
                                                y_axis_labels = c(0,0.1,0.2),
                                                point_size = 4)

HNSC_HPVneg_density <- ggdraw(HNSC_HPVneg_density) + draw_label(label = expression(HPV^{"−"}~HNSC),x = .55,y=.9,size = common.size)

save_plot(filename = "figures/HNSC_HPVneg_density.png",plot = HNSC_HPVneg_density,base_height = 5/6,base_width = 6.5/4)


# KIRC ---- 


KIRC_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "KIRC",
                                         main_size = common.size,
                                         sub_1 = "VHL S68*",
                                         sub_2 = "VHL S65*",
                                         sub_1_lab_y = .1,
                                         sub_1_lab_x = 5e5,
                                         sub_1_arrow_x =  185360.85,
                                         sub_1_arrow_y = 0.025,
                                         sub_2_lab_y = .05,
                                         sub_2_lab_x = 8.5e4,
                                         sub_2_arrow_x =  118487.13,
                                         sub_2_arrow_y = .025,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),
                                         x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),extender_percent = 0.06,
                                         sub_1_arrow_alpha = 1,
                                         sub_2_arrow_alpha = 1,
                                         y_axis_limits = c(0,.22),
                                         y_axis_breaks = c(0,0.1,0.2),
                                         y_axis_labels = c(0,0.1,0.2),
                                         point_size = 4)

KIRC_density <- ggdraw(KIRC_density) + draw_label(label = "KIRC",x = .5,y=.9,size = common.size)

save_plot(filename = "figures/KIRC_density.png",plot = KIRC_density,base_height = 5/6,base_width = 6.5/4)


# LAML ----

LAML_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "LAML",
                                         main_size = common.size,
                                         sub_1 = "KIT D816V",
                                         sub_2 = "KRAS G12V",
                                         sub_1_lab_y = .075,
                                         sub_1_lab_x = 8e5,
                                         sub_1_arrow_x =  1969488.61,
                                         sub_1_arrow_y = 0.025,
                                         sub_2_lab_y = .05,
                                         sub_2_lab_x = 3e5,
                                         sub_2_arrow_x =   474810.3,
                                         sub_2_arrow_y = .025,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),
                                         x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),extender_percent = 0.02,
                                         sub_1_arrow_alpha = 1,
                                         sub_2_arrow_alpha = 1,
                                         y_axis_limits = c(0,.13),
                                         y_axis_breaks = c(0,0.1,0.2),
                                         y_axis_labels = c(0,0.1,0.2),
                                         point_size = 4)

LAML_density <- ggdraw(LAML_density) + draw_label(label = "LAML",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/LAML_density.png",plot = LAML_density,base_height = 5/6,base_width = 6.5/4)


# LGG ---- 
LGG_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                        this_tumor = "LGG",
                                        main_size = common.size,
                                        sub_1 = "IDH1 R132G",
                                        sub_2 = "IDH1 R132S",
                                        sub_1_lab_y = .13,
                                        sub_1_lab_x = 6e6,
                                        sub_1_arrow_x =  13233018.63,
                                        sub_1_arrow_y = 0.025,
                                        sub_2_lab_y = .08,
                                        sub_2_lab_x = 1e6,
                                        sub_2_arrow_x =   2848965.8,
                                        sub_2_arrow_y = .025,
                                        x_axis_breaks = c(1,1e5,1e6,1e7),
                                        x_axis_labels = c(expression(10^0),expression(10^5),expression(10^6),expression(10^7)),extender_percent = 0.02,sub_1_arrow_alpha = 1,sub_2_arrow_alpha = 1,
                                        y_axis_limits = c(0,.22),
                                        y_axis_breaks = c(0,0.1,0.2),
                                        y_axis_labels = c(0,0.1,0.2),
                                        point_size = 4)

LGG_density <- ggdraw(LGG_density) + draw_label(label = "LGG",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/LGG_density.png",plot = LGG_density,base_height = 5/6,base_width = 6.5/4)



# LIHC ---- 
LIHC_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "LIHC",
                                         main_size = common.size,
                                         sub_1 = "BAZ2A R1642L",
                                         sub_2 = "CTNNB1 I35S",
                                         sub_1_lab_y = .15,
                                         sub_1_lab_x = 3e5,
                                         sub_1_arrow_x =  949298.36,
                                         sub_1_arrow_y = 0.025,
                                         sub_2_lab_y = .065,
                                         sub_2_lab_x = 1.5e5,
                                         sub_2_arrow_x =   228601.05,
                                         sub_2_arrow_y = .025,
                                         x_axis_breaks = c(1,1e5,1e6,1e7),
                                         x_axis_labels = c(expression(10^0),expression(10^5),expression(10^6),expression(10^7)),extender_percent = 0.06,
                                         sub_1_arrow_alpha = 1
                                         ,sub_2_arrow_alpha = 1,
                                         y_axis_limits = c(0,.22),
                                         y_axis_breaks = c(0,0.1,0.2),
                                         y_axis_labels = c(0,0.1,0.2),
                                         point_size = 4)

LIHC_density <- ggdraw(LIHC_density) + draw_label(label = "LIHC",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/LIHC_density.png",plot = LIHC_density,base_height = 5/6,base_width = 6.5/4)

# LUAD ---- 

LUAD_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "LUAD",
                                         main_size = common.size,
                                         sub_1 = "EGFR L858R",
                                         sub_2 = "CTNNB1 S37F",
                                         sub_1_lab_y = .15,
                                         sub_1_lab_x = 5e5,
                                         sub_1_arrow_x =  280236.22,
                                         sub_1_arrow_y = 0.05,
                                         sub_2_lab_y = .08,
                                         sub_2_lab_x = 6e4,
                                         sub_2_arrow_x =    120122.92,
                                         sub_2_arrow_y = .05,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),
                                         x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),extender_percent = 0.06,
                                         sub_1_arrow_alpha = 1,
                                         sub_2_arrow_alpha = 1,
                                         y_axis_limits = c(0,.3),
                                         y_axis_breaks = c(0,0.1,0.2),
                                         y_axis_labels = c(0,0.1,0.2),
                                         point_size = 4)

LUAD_density <- ggdraw(LUAD_density) + draw_label(label = "LUAD",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/LUAD_density.png",plot = LUAD_density,base_height = 5/6,base_width = 6.5/4)




# LUSC ---- 

LUSC_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "LUSC",
                                         main_size = common.size,
                                         sub_1 = "TP53 Y234S",
                                         sub_2 = "NFE2L2 R34G",
                                         sub_1_lab_y = .1,
                                         sub_1_lab_x = 5e5,
                                         sub_1_arrow_x =  90448.39,
                                         sub_1_arrow_y = 0.03,
                                         sub_2_lab_y = .16,
                                         sub_2_lab_x =5e4,
                                         sub_2_arrow_x =    69242.49,
                                         sub_2_arrow_y = .05,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),
                                         x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),
                                         extender_percent = 0.02,
                                         sub_1_arrow_alpha = 1,
                                         sub_2_arrow_alpha = 1,
                                         y_axis_limits = c(0,.3),
                                         y_axis_breaks = c(0,0.1,0.2),
                                         y_axis_labels = c(0,0.1,0.2),
                                         point_size = 4)

LUSC_density <- ggdraw(LUSC_density) + draw_label(label = "LUSC",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/LUSC_density.png",plot = LUSC_density,base_height = 5/6,base_width = 6.5/4)




# OV ---- 

OV_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                       this_tumor = "OV",
                                       main_size = common.size,
                                       sub_1 = "MGRN1 G64A",
                                       sub_2 = "TP53 V157F",
                                       sub_1_lab_y = .13,
                                       sub_1_lab_x = 5e5,
                                       sub_1_arrow_x =   409080.25,
                                       sub_1_arrow_y = 0.025,
                                       sub_2_lab_y = .065,
                                       sub_2_lab_x = 8e4,
                                       sub_2_arrow_x =    184993.66,
                                       sub_2_arrow_y = .025,
                                       x_axis_breaks = c(1,1e4,1e5,1e6),
                                       x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),extender_percent = 0.06,
                                       sub_1_arrow_alpha = 1,
                                       sub_2_arrow_alpha = 1,
                                       y_axis_limits = c(0,.22),
                                       y_axis_breaks = c(0,0.1,0.2),
                                       y_axis_labels = c(0,0.1,0.2),
                                       point_size = 4)

OV_density <- ggdraw(OV_density) + draw_label(label = "OV",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/OV_density.png",plot = OV_density,base_height = 5/6,base_width = 6.5/4)



# PAAD ---- 

PAAD_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "PAAD",
                                         main_size = common.size,
                                         sub_1 = "KRAS G12R",
                                         sub_2 = "KRAS G12V",
                                         sub_1_lab_y = .12,
                                         sub_1_lab_x = 2e6,
                                         sub_1_arrow_x =   9513627.87,
                                         sub_1_arrow_y = 0.025,
                                         sub_2_lab_y = .065,
                                         sub_2_lab_x = 8e5,
                                         sub_2_arrow_x =    2801115.02,
                                         sub_2_arrow_y = .025,
                                         x_axis_breaks = c(1,1e5,1e6),
                                         x_axis_labels = c(expression(10^0),expression(10^5),expression(10^6)),extender_percent = 0.02,sub_1_arrow_alpha = 1,
                                         sub_2_arrow_alpha = 1,
                                         y_axis_limits = c(0,.22),
                                         y_axis_breaks = c(0,0.1,0.2),
                                         y_axis_labels = c(0,0.1,0.2),
                                         point_size = 4)

PAAD_density <- ggdraw(PAAD_density) + draw_label(label = "PAAD",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/PAAD_density.png",plot = PAAD_density,base_height = 5/6,base_width = 6.5/4)


# PRAD ---- 

PRAD_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "PRAD",
                                         main_size = common.size,
                                         sub_1 = "SPOP Y87S",
                                         sub_2 = "SPOP W131G",
                                         sub_1_lab_y = .13,
                                         sub_1_lab_x = 2e5,
                                         sub_1_arrow_x =   683195.22,
                                         sub_1_arrow_y = 0.025,
                                         sub_2_lab_y = .065,
                                         sub_2_lab_x = 1e5,
                                         sub_2_arrow_x =    552894.51,
                                         sub_2_arrow_y = .025,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),
                                         x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),extender_percent = 0.06,
                                         sub_1_arrow_alpha = 1,
                                         sub_2_arrow_alpha = 1,
                                         y_axis_limits = c(0,.22),
                                         y_axis_breaks = c(0,0.1,0.2),
                                         y_axis_labels = c(0,0.1,0.2),
                                         point_size = 4)

PRAD_density <- ggdraw(PRAD_density) + draw_label(label = "PRAD",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/PRAD_density.png",plot = PRAD_density,base_height = 5/6,base_width = 6.5/4)


# READ ---- 

READ_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "READ",
                                         main_size = common.size,
                                         sub_1 = "APC E1322*",
                                         sub_2 = "APC E1209*",
                                         sub_1_lab_y = .15,
                                         sub_1_lab_x = 2e5,
                                         sub_1_arrow_x =   2613352.63,
                                         sub_1_arrow_y = 0.03,
                                         sub_2_lab_y = .065,
                                         sub_2_lab_x = 2e5,
                                         sub_2_arrow_x =    2513352.63,
                                         sub_2_arrow_y = .02,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),
                                         x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),extender_percent = 0.05,sub_1_arrow_alpha = 1,sub_2_arrow_alpha = 1,
                                         y_axis_limits = c(0,.3),
                                         y_axis_breaks = c(0,0.1,0.2),
                                         y_axis_labels = c(0,0.1,0.2),
                                         point_size = 4)

READ_density <- ggdraw(READ_density) + draw_label(label = "READ",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/READ_density.png",plot = READ_density,base_height = 5/6,base_width = 6.5/4)



# SKCMM ---- 

SKCMM_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                          this_tumor = "SKCMM",
                                          main_size = common.size,
                                          sub_1 = "NRAS Q61R",
                                          sub_2 = "NRAS Q61K",
                                          sub_1_lab_y = .4,
                                          sub_1_lab_x = 2e6,
                                          sub_1_arrow_x =   4707860.13,
                                          sub_1_arrow_y = 0.1,
                                          sub_2_lab_y = .2,
                                          sub_2_lab_x = 1e6,
                                          sub_2_arrow_x =    4438279.31,
                                          sub_2_arrow_y = .08,
                                          x_axis_breaks = c(1,1e4,1e5,1e6,1e7),
                                          x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6),expression(10^7)),
                                          extender_percent = 0.02,
                                          sub_1_arrow_alpha = 1,
                                          sub_2_arrow_alpha = 1,
                                          y_axis_limits = c(0,.65),
                                          y_axis_breaks = c(0,0.2,0.4,0.6),
                                          y_axis_labels = c(0,0.2,0.4,0.6),
                                          point_size = 4)

SKCMM_density <- ggdraw(SKCMM_density) + draw_label(label = "SKCMM",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/SKCMM_density.png",plot = SKCMM_density,base_height = 5/6,base_width = 6.5/4)



# SKCMP ---- 

SKCMP_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                          this_tumor = "SKCMP",
                                          main_size = common.size,
                                          sub_1 = "BRAF V600E",
                                          sub_2 = "KIT K642E",
                                          sub_1_lab_y = .4,
                                          sub_1_lab_x = 1e6,
                                          sub_1_arrow_x =   3492084.99,
                                          sub_1_arrow_y = 0.1,
                                          sub_2_lab_y = .3,
                                          sub_2_lab_x = 1e5,
                                          sub_2_arrow_x =    167735.46,
                                          sub_2_arrow_y = .1,
                                          x_axis_breaks = c(1,1e4,1e5,1e6,1e7),
                                          x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6),expression(10^7)),
                                          extender_percent = 0.02,
                                          sub_1_arrow_alpha = 1,
                                          sub_2_arrow_alpha = 1,
                                          y_axis_limits = c(0,.65),
                                          y_axis_breaks = c(0,0.2,0.4,0.6),
                                          y_axis_labels = c(0,0.2,0.4,0.6),
                                          point_size = 4)

SKCMP_density <- ggdraw(SKCMP_density) + draw_label(label = "SKCMP",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/SKCMP_density.png",plot = SKCMP_density,base_height = 5/6,base_width = 6.5/4)



# STAD ----

STAD_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "STAD",
                                         main_size = common.size,
                                         sub_1 = "RHOA Y42S",
                                         sub_2 = "PIK3CA H1047R",
                                         sub_1_lab_y = .13,
                                         sub_1_lab_x = 5e5,
                                         sub_1_arrow_x =65728.42,
                                         sub_1_arrow_y = 0.05,
                                         sub_2_lab_y = .25,
                                         sub_2_lab_x = 5e4,
                                         sub_2_arrow_x = 27539.21,
                                         sub_2_arrow_y = .06,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),
                                         x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),
                                         y_axis_limits = c(0,.4),
                                         y_axis_breaks = c(0,0.1,0.2,0.3),
                                         y_axis_labels = c(0,0.1,0.2,0.3),
                                         point_size = 4,extender_percent = 0.06)

STAD_density <- ggdraw(STAD_density) + draw_label("STAD",x = .5,y=.9,size = common.size)

save_plot(filename = "figures/STAD_density.png",plot = STAD_density,base_height = 5/6,base_width = 6.5/4)


# THCA ---- 

THCA_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "THCA",
                                         main_size = common.size,
                                         sub_1 = "BRAF V600E",
                                         sub_2 = "NRAS Q61R",
                                         sub_1_lab_y = .07,
                                         sub_1_lab_x = 1e7,
                                         sub_1_arrow_x =   33029682.79,
                                         sub_1_arrow_y = 0.025,
                                         sub_2_lab_y = .05,
                                         sub_2_lab_x = 2e6,
                                         sub_2_arrow_x =    4601059.2,
                                         sub_2_arrow_y = .025,
                                         x_axis_breaks = c(1,1e5,1e6,1e7),
                                         x_axis_labels = c(expression(10^0),expression(10^5),expression(10^6),expression(10^7)),extender_percent = 0.02,
                                         sub_1_arrow_alpha = 1,
                                         sub_2_arrow_alpha = 1,
                                         y_axis_limits = c(0,.13),
                                         y_axis_breaks = c(0,0.1,0.2),
                                         y_axis_labels = c(0,0.1,0.2),
                                         point_size = 4)

THCA_density <- ggdraw(THCA_density) + draw_label(label = "THCA",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/THCA_density.png",plot = THCA_density,base_height = 5/6,base_width = 6.5/4)


# UCEC ---- 

UCEC_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "UCEC",
                                         main_size = common.size,
                                         sub_1 = "PTEN R130G",
                                         sub_2 = "NFE2L2 R34G",
                                         sub_1_lab_y = .4,
                                         sub_1_lab_x = 5e5,
                                         sub_1_arrow_x =   219762.83,
                                         sub_1_arrow_y = 0.09,
                                         sub_2_lab_y = .3,
                                         sub_2_lab_x = 5e4,
                                         sub_2_arrow_x =    149243.04,
                                         sub_2_arrow_y = .09,
                                         x_axis_breaks = c(1,1e4,1e5,1e6,1e7),
                                         x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6),expression(10^7)),
                                         extender_percent = 0.045,
                                         sub_1_arrow_alpha = 1,
                                         sub_2_arrow_alpha = 1,
                                         y_axis_limits = c(0,.65),
                                         y_axis_breaks = c(0,0.2,0.4,0.6),
                                         y_axis_labels = c(0,0.2,0.4,0.6),
                                         point_size = 4)

UCEC_density <- ggdraw(UCEC_density) + draw_label(label = "UCEC",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/UCEC_density.png",plot = UCEC_density,base_height = 5/6,base_width = 6.5/4)


# legend
legend_plot <- ggplot(data = subset(combined_all_data.noNA, tumor_type=="LUAD"), 
                      aes(x=gamma_epistasis)) + 
  geom_density(aes(fill=synonymous),position = "stack") + scale_fill_discrete(name="Substitution\ntype",labels=c("Non-synonymous", "Synonymous")) + theme(legend.title=element_text(size=common.size) , legend.text=element_text(size=common.size))

legend_for_plot <- get_legend(legend_plot)

#### combining density plots  ----

plot_combine <- plot_grid(   UCEC_density,
                             SKCMM_density,
                             COAD_density,
                             SKCMP_density,
                             STAD_density,
                             CESC_density,
                             BRCA_density,
                             BLCA_density,
                             LUAD_density,
                             LUSC_density,
                             HNSC_HPVneg_density,
                             READ_density,
                             LIHC_density,
                             LGG_density,
                             GBM_density,
                             ESCA_density,
                             PRAD_density,
                             PAAD_density,
                             OV_density,
                             KIRC_density,
                             HNSC_HPVpos_density,
                             THCA_density,
                             LAML_density,
                             legend_for_plot,
                             ncol=4,align = 'hv')

combined_with_labs <- ggdraw(plot = plot_grid(NULL,plot_combine,NULL,NULL,nrow = 2,ncol = 2,rel_heights = c(1,0.02),rel_widths = c(0.02,1))) + 
  geom_text(x=.5,y=0.01,label="Selection intensity",size=common.size*(5/14)) + 
  geom_text(x=.01,y=0.5,label="Density",size=common.size*(5/14),angle=90)

save_plot(filename = "figures/combined_density_plot.png",base_width = 6.5,plot = plot_combine,dpi=600,base_height = 10)

save_plot(filename = "figures/combined_density_plot_w_labs.png",base_width = 6.5,plot = combined_with_labs,base_height = 5)




length(which(combined_all_data.noNA$synonymous==T))/nrow(combined_all_data.noNA)

length(which(combined_all_data.noNA$synonymous==F))/nrow(combined_all_data.noNA)
