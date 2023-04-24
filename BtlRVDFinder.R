#Btl RVD Finder
#Morgan Carter 2023 R v 4.2.2
#Takes amino acid fasta files of putative Btl proteins and extracts the RVD sequence. 
#It also identifies other important motifs and names them based on strain and repeat number.


# Packages ----------------------------------------------------------------
setwd("") #Change for your use 
library(seqinr)
library(stringr)


# Inputs and Setup ------------------------------------------------------------------
fasta_file <- read.fasta("Final_All_Btls.fa", as.string = T, set.attributes = F, forceDNAtolower = T)
Btl_info_table = data.frame(matrix(ncol=10, nrow=length(fasta_file))) #Starts data frame for our information to fill in
colnames(Btl_info_table) <- c( "fastaName", "Btl", "Repeat_no", "StartAA", "PseudoRepeat", "FirstRepeat", "LastRepeat", "NLS", "Complete", "RVDs")
Inc_btl_tracker <- c("starting_filler") #tracker for incomplete sequences to be named sequentially

# Outputs -----------------------------------------------------------------
rvd_out_file = "Btl_RVDs.fa"
aa_out_file = "Btl_aa.fa"
Btl_info_file = "Btl_info.csv"
aa_intact_out_file = "Btl_intact_aa.fa" #only going to be Btl proteins that are predicted complete/intact


# Main Function -----------------------------------------------------------
#This extracts all the data pieces from each Btl sequence and names the Btl

for(i in 1:length(fasta_file)){
  RVD_seq <- c()
  Btl_info_table$StartAA[i] <- substr(fasta_file[[i]], 1, 1) #check starting amino acid to see if complete (m)
  if(grepl("i(k|e)k..gg", fasta_file[[i]])==T){ #to find out if an HS/HY pseudo repeat is first and pulls that
    HS_rep<- regmatches(fasta_file[[i]], gregexpr("i(k|e)k..gg", fasta_file[[i]]))
    HS_rep_rvd <- c(str_sub(unname(unlist(HS_rep)), 4,5)) #grab just the "rvd"
    Btl_info_table$PseudoRepeat[i] <- toupper(HS_rep_rvd) #add to the info table
  }else{
    Btl_info_table$PseudoRepeat[i] <- "None"
  }
  if(grepl("(s|c|n|d|g)(y|c|h)d(c|g|y)(g|a)a", fasta_file[[i]])==T){#to find out if a YD/CD repeat is first and pulls that
    YD_rep<- regmatches(fasta_file[[i]], gregexpr("(s|c|n|d|g)(y|c|h)d(c|g|y)(g|a)a", fasta_file[[i]]))
    YD_rep_rvd <-c(str_sub(unname(unlist(YD_rep)), 2,3))
    RVD_seq <-
      c(YD_rep_rvd) #grab just the YD/HD
    Btl_info_table$FirstRepeat[i] <- toupper(YD_rep_rvd)
  }else{
    Btl_info_table$FirstRepeat[i] <- "Other"
  }
  repeat_chunk <- regmatches(fasta_file[[i]], gregexpr("(i|m|v|d)(a|t|h|v).{1,3}(g|s)(g|s)(s|a|t)(p|q|r|w)", fasta_file[[i]]))#pulls the repeat chunks with the RVD
  repeat_chunk_char <- unname(unlist(repeat_chunk)) #convert from list to character vector
  for(j in 1:length(repeat_chunk_char)){ #for each repeat chunk in the list of repeat chunks from this protein
      rep_mod <- repeat_chunk_char[j] #pull repeat chunks one by one
      rep_mod <- str_sub(rep_mod, 4, -5) #takes off the first and last amino acids - just RVD
      n=2-nchar(rep_mod) #looking to see if 2 amino acids are left
      miss_res <- str_c(replicate(n, "*"), collapse = "") #make a vector of the right number of *s
      rep_mod <- paste(rep_mod, miss_res, sep="") #create a string of the RVD with added *s if necessary
      RVD_seq <- c(RVD_seq, rep_mod) #add to vector of RVDs for this protein
  }
  repeat_num <- length(RVD_seq)
  Btl_info_table$Repeat_no[i] <- repeat_num
  if(grepl("n\\*",RVD_seq[repeat_num])==F){ #looking for last RVD (typicall N*)
    Btl_info_table$LastRepeat[i] <- "Other"
  }else{
    Btl_info_table$LastRepeat[i] <- "N*"
  }
  RVD_seq <- str_c(RVD_seq, collapse=" ")  
  Btl_info_table$RVDs[i] <- RVD_seq
  if(grepl(".irk",fasta_file[[i]])==T){ #looking for the NLS (QIRK or RIRK)
    Btl_info_table$NLS[i] <- toupper(unlist(regmatches(fasta_file[[i]], gregexpr(".irk", fasta_file[[i]]))))
  }else{
    Btl_info_table$NLS[i] <- c("None")
  }
  if(str_detect("None",Btl_info_table$NLS[i])==T|(str_detect("None",Btl_info_table$PseudoRepeat[i])==T&str_detect("Other",Btl_info_table$FirstRepeat[i])==T)|str_detect("m", Btl_info_table$StartAA[i])==F){ #if there is no NLS and it doesn't start with a start codon or have either pseudo/first repeat, call it incomplete
    Btl_info_table$Complete[i] <- "No"
    if(str_detect("None",Btl_info_table$PseudoRepeat[i])==T & str_detect("Other",Btl_info_table$FirstRepeat[i])==T){
      Btl_info_table$FirstRepeat[i] <- NA
    }
    if(str_detect("None",Btl_info_table$NLS[i])==T & str_detect("Other",Btl_info_table$LastRepeat[i])==T){
      Btl_info_table$LastRepeat[i] <- NA
    }
  }
  if(grepl("Btl",names(fasta_file[i]))==F){ #if not named with Btl in name then format it to give it a "btl" name
    strain_name <- str_remove(names(fasta_file[i]), "_metscf_[0-9]*_[0-9]") #removing info about the contig from ZygoLife data
    strain_num <- str_sub(strain_name, -2) #last 2 digits in the strain name
    strain_num <- str_remove(strain_num, "_") #if strain number is only 1 digit, HKI strains
    if("No"%in%Btl_info_table$Complete[i]){
      btl_name <- paste("Btl", "A", "_", strain_num, sep="") #this mess is for adding a letter as a stand in for the repeat no in incomplete proteins
      k=1
      while(btl_name%in%Inc_btl_tracker){ #checking to see if an incomplete one is already in that strain and ups the letter
        k=k+1
        btl_name <- paste("Btl", toupper(letters[k]), "_", strain_num, sep="")
      }
      Inc_btl_tracker <- c(Inc_btl_tracker, btl_name)
    }else{
      Btl_info_table$Complete[i] <- "Yes" #For complete Btl proteins, making a new name with the repeat number
      btl_name <- paste("Btl", repeat_num, "_", strain_num, sep="")
    }
    fasta_head <- paste(">", btl_name, " : ", names(fasta_file[i]),sep="") #creates fasta header format with name of protein
  }else{ #just take the name in the file
    btl_name <- names(fasta_file[i])
    fasta_head <- paste(">", btl_name, sep="")
    Btl_info_table$Complete[i] <- "Yes"
  }
  Btl_info_table$fastaName[i] <- gsub(">", "", fasta_head)
  fasta <- paste(fasta_head, toupper(RVD_seq), sep="\n") #creates fasta entry with header and RVD sequence
  write(fasta, rvd_out_file, append=T)
  Btl_info_table$Btl[i] <- btl_name
  aa_fasta <- paste(fasta_head, toupper(fasta_file[[i]]), sep="\n")
  if(Btl_info_table$Complete[i]=="Yes"){
    write(aa_fasta, aa_intact_out_file, append=T)
  }
  write(aa_fasta, aa_out_file, append=T)
}
write.csv(Btl_info_table, file=Btl_info_file, row.names = F) #you will get an warning message about the first repeat replacement length but ignore that