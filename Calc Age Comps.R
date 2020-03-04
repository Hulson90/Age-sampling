#########################
# Script to estimate population #'s at age for GOA/AI species
# This script builds from previously calculated pop'n est's at size

#########################
# Set directories and read in necessary data
path<-getwd()
pathD<-paste0(path,"/Data")
pathR<-paste0(path,"/AC results")

sizepop<-read.csv(paste0(pathD,"/sizepop.csv"))
#specimen<-read.csv(paste0(pathD,"/specimen.csv"))
specimen<-read.csv(paste0(pathD,"/specimen_full.csv"))
specimen<-subset(specimen,is.na(specimen$AGE)==FALSE)
agepop_RACE<-read.csv(paste0(pathD,"/agepop.csv"))

# Get data parameters together
species<-sort(unique(specimen$SPECIES_CODE))
yrs<-sort(unique(specimen$YEAR))
sex<-sort(unique(specimen$SEX))

#########################
# Loop thru species
for(sp in 1:length(species)){
  
  sizepop_sp<-subset(sizepop,sizepop$SPECIES_CODE==species[sp])
  specimen_sp<-subset(specimen,specimen$SPECIES_CODE==species[sp])
  agepop_RACE_sp<-subset(agepop_RACE,agepop_RACE$SPECIES_CODE==species[sp])
  
  # Set up results matrices
  ages<-sort(unique(specimen_sp$AGE))
  AGEPOP_M<-matrix(nrow=length(yrs),ncol=length(ages))
  colnames(AGEPOP_M)<-as.character(ages)
  rownames(AGEPOP_M)<-as.character(yrs)
  AGEPOP_F<-matrix(nrow=length(yrs),ncol=length(ages))
  colnames(AGEPOP_F)<-as.character(ages)
  rownames(AGEPOP_F)<-as.character(yrs)
  AGEPOP_U<-matrix(nrow=length(yrs),ncol=length(ages))
  colnames(AGEPOP_U)<-as.character(ages)
  rownames(AGEPOP_U)<-as.character(yrs)
  AGEPOP_T<-matrix(nrow=length(yrs),ncol=length(ages))
  colnames(AGEPOP_T)<-as.character(ages)
  rownames(AGEPOP_T)<-as.character(yrs)
  
  # Set up test matrix for comparison with RACE output
  RACEmatch<-matrix(nrow=length(yrs),ncol=4)
  colnames(RACEmatch)<-c("Year","Male","Female","Unsexed")
  RACEmatch<-as.data.frame(RACEmatch)
  RACEmatch$Year<-yrs
  
  #########################
  # Loop thru years
  for(y in 1:length(yrs)){
    
    sizepop_sp_y<-subset(sizepop_sp,sizepop_sp$YEAR==yrs[y])
    specimen_sp_y<-subset(specimen_sp,specimen_sp$YEAR==yrs[y])
    agepop_RACE_sp_y<-subset(agepop_RACE_sp,agepop_RACE_sp$SURVEY_YEAR==yrs[y])
    
    #########################
    # Loop thru sex
    for(sx in 1:length(sex)){
      
      specimen_sp_y_sx<-subset(specimen_sp_y,specimen_sp_y$SEX==sex[sx])
      agepop_RACE_sp_y_sx<-subset(agepop_RACE_sp_y,agepop_RACE_sp_y$SEX==sex[sx])
      
      # Remove matrices to wipe clean each loop
      if(sx==1 & exists("pop_age_est_M")==TRUE)
        rm(pop_age_est_M)
      if(sx==2 & exists("pop_age_est_F")==TRUE)
        rm(pop_age_est_F)
      if(sx==3 & exists("pop_age_est_U")==TRUE)
        rm(pop_age_est_U)
      
      # Test if there's specimen data for particular sex
      if(length(specimen_sp_y_sx$SEX) == 0) {
        cat(paste("No specimen data for sex", sx, "\n"))
        next
      }
      
      # If sex unknown and there is specimen data then use all specimen data
      if(sx==3)
        specimen_sp_y_sx<-specimen_sp_y
      
      # If there is no sizecomp data, we are wasting our time
      if(sex[sx] == 1 & sum(sizepop_sp_y$MALES)==0) {
        cat(paste("No sizecomp data for sex", sx, "\n"))
        next
      }
      if(sex[sx] == 2 & sum(sizepop_sp_y$FEMALES)==0) {
        cat(paste("No sizecomp data for sex", sx, "\n"))
        next
      }
      if(sex[sx] == 3 & sum(sizepop_sp_y$UNSEXED)==0) {
        cat(paste("No sizecomp data for sex", sx, "\n"))
        next
      }
      
      # Get vector of possible lengths
      lenlist<-seq(from = min(sizepop_sp_y$LENGTH), to = max(sizepop_sp_y$LENGTH), by = 10)
      
      # Get number of ages by length and age
      age_num <- tapply(specimen_sp_y_sx$AGE, list(specimen_sp_y_sx$LENGTH, specimen_sp_y_sx$AGE), length)
      age_num[is.na(age_num)] <- 0
      
      # Turn these into fractions
      age_frac <- age_num/apply(age_num, 1, sum)
      
      # Find lengths from age data where there is no sizecomp data
      no.lengths <- unique(sort(specimen_sp_y_sx$LENGTH))[is.na(match(unique(sort(specimen_sp_y_sx$LENGTH)), sizepop_sp_y$LENGTH))]
      
      # Find lengths from sizecomp data where there is no age data.
      no.ages <- sizepop_sp_y$LENGTH[is.na(match(sizepop_sp_y$LENGTH, unique(specimen_sp_y_sx$LENGTH)))]
      if(length(no.ages) == 0)
        no.ages <- sizepop_sp_y$LENGTH
      no.age.sizecomp <- sizepop_sp_y[match(no.ages, sizepop_sp_y$LENGTH),]
      
      # Nothing else we can do when there are age data with no sizecomp data, so get rid of these records
      age_frac <- age_frac[is.na(match(as.numeric(dimnames(age_frac)[[1]]), no.lengths)),  ]
      
      # Estimate numbers by age and length
      if(sex[sx] ==1)
        pop_age_est_M <- age_frac * sizepop_sp_y$MALES[match(as.numeric(dimnames(age_frac)[[1]]), as.numeric(sizepop_sp_y$LENGTH), nomatch = 0, incomparables = no.lengths)]
      if(sex[sx] ==2)
        pop_age_est_F <- age_frac * sizepop_sp_y$FEMALES[match(as.numeric(dimnames(age_frac)[[1]]), as.numeric(sizepop_sp_y$LENGTH), nomatch = 0, incomparables = no.lengths)]
      if(sex[sx] ==3)
        pop_age_est_U <- age_frac * sizepop_sp_y$UNSEXED[match(as.numeric(dimnames(age_frac)[[1]]), as.numeric(sizepop_sp_y$LENGTH), nomatch = 0, incomparables = no.lengths)]
      
      #########################
      # End sex loop   
    }
    
    # Now sum up the numbers at age for all lengths and remove any 0s, and check to see if matches with RACE output
    
    # Males
    if(exists("pop_age_est_M")==TRUE){if(length(pop_age_est_M)>0){
      age_est_M <- apply(pop_age_est_M, 2,sum)
      if(length(which(age_est_M==0))>0)
        age_est_M <- age_est_M[-which(age_est_M==0)]
      agepop_RACE_sp_y_M<-subset(agepop_RACE_sp_y,agepop_RACE_sp_y$SEX==sex[1])
      test_M<-matrix(nrow=length(names(age_est_M)),ncol=3)
      colnames(test_M)<-c("Age","Calc","RACE")
      test_M[,1]<-as.numeric(names(age_est_M))
      test_M[,2]<-age_est_M
      for(i in 1:length(test_M[,3])){
        if(length(match(test_M[i,1],agepop_RACE_sp_y_M$AGE))>0)
          test_M[i,3]<-agepop_RACE_sp_y_M$AGEPOP[match(test_M[i,1],agepop_RACE_sp_y_M$AGE)]
      }
      test_M<-as.data.frame(test_M)
      RACEmatch$Male[y]<-max(abs(test_M$Calc-test_M$RACE)/test_M$RACE)
      AGEPOP_M[y,match(as.numeric(names(age_est_M)),ages)]<-age_est_M
      AGEPOP_M[y,is.na(AGEPOP_M[y,])] <- 0
    }}
    
    # Females
    if(exists("pop_age_est_F")==TRUE){if(length(pop_age_est_F)>0){
      age_est_F <- apply(pop_age_est_F, 2,sum)
      if(length(which(age_est_F==0))>0)
        age_est_F <- age_est_F[-which(age_est_F==0)]
      agepop_RACE_sp_y_F<-subset(agepop_RACE_sp_y,agepop_RACE_sp_y$SEX==sex[2])
      test_F<-matrix(nrow=length(names(age_est_F)),ncol=3)
      colnames(test_F)<-c("Age","Calc","RACE")
      test_F[,1]<-as.numeric(names(age_est_F))
      test_F[,2]<-age_est_F
      for(i in 1:length(test_F[,3])){
        if(length(match(test_F[i,1],agepop_RACE_sp_y_F$AGE))>0)
          test_F[i,3]<-agepop_RACE_sp_y_F$AGEPOP[match(test_F[i,1],agepop_RACE_sp_y_F$AGE)]
      }
      test_F<-as.data.frame(test_F)
      RACEmatch$Female[y]<-max(abs(test_F$Calc-test_F$RACE)/test_F$RACE)
      AGEPOP_F[y,match(as.numeric(names(age_est_F)),ages)]<-age_est_F
      AGEPOP_F[y,is.na(AGEPOP_F[y,])] <- 0
    }}
    
    # Unsexed
    if(exists("pop_age_est_U")==TRUE){if(length(pop_age_est_U)>0){
      age_est_U <- apply(pop_age_est_U, 2,sum)
      if(length(which(age_est_U==0))>0)
        age_est_U <- age_est_U[-which(age_est_U==0)]
      agepop_RACE_sp_y_U<-subset(agepop_RACE_sp_y,agepop_RACE_sp_y$SEX==sex[3])
      test_U<-matrix(nrow=length(names(age_est_U)),ncol=3)
      colnames(test_U)<-c("Age","Calc","RACE")
      test_U[,1]<-as.numeric(names(age_est_U))
      test_U[,2]<-age_est_U
      for(i in 1:length(test_U[,3])){
        if(length(match(test_U[i,1],agepop_RACE_sp_y_U$AGE))>0)
          test_U[i,3]<-agepop_RACE_sp_y_U$AGEPOP[match(test_U[i,1],agepop_RACE_sp_y_U$AGE)]
      }
      test_U<-as.data.frame(test_U)
      RACEmatch$Unsexed[y]<-max(abs(test_U$Calc-test_U$RACE)/test_U$RACE)
      AGEPOP_U[y,match(as.numeric(names(age_est_U)),ages)]<-age_est_U
      AGEPOP_U[y,is.na(AGEPOP_U[y,])] <- 0
    }}
    
    
    
    #########################
    # End year loop
  }
  
  # Test if there's disagreement with RACE calcs (larger than 5% for the maximum difference at age - there are some years that are different with RACE that are unexplainable)
  if(length(which(RACEmatch$Male > 0.05)) > 0)
    cat(paste("Disagreement between calcs and RACE for Males of species", species[sp], "\n"))
  if(length(which(RACEmatch$Female > 0.05)) > 0)
    cat(paste("Disagreement between calcs and RACE for Females of species", species[sp], "\n"))
  if(length(which(RACEmatch$Unsexed > 0.05)) > 0)
    cat(paste("Disagreement between calcs and RACE for Unsexed of species", species[sp], "\n"))
  
  # Finalize results matrices
  AGEPOP_T <- matrix(mapply(sum,AGEPOP_M,AGEPOP_F,AGEPOP_U,MoreArgs=list(na.rm=TRUE)),ncol=length(ages))
  colnames(AGEPOP_T)<-as.character(ages)
  rownames(AGEPOP_T)<-as.character(yrs)
  
  # Write results matrices
  write.csv(RACEmatch,paste0(pathR,"/RACEmatch_",species[sp],".csv"))
  write.csv(AGEPOP_M,paste0(pathR,"/AGEPOP_M_",species[sp],".csv"))
  write.csv(AGEPOP_F,paste0(pathR,"/AGEPOP_F_",species[sp],".csv"))
  write.csv(AGEPOP_U,paste0(pathR,"/AGEPOP_U_",species[sp],".csv"))
  write.csv(AGEPOP_T,paste0(pathR,"/AGEPOP_T_",species[sp],".csv"))
  
  #########################
  # End species loop 
}









