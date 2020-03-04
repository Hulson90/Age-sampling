#########################
# Script to estimate population #'s at size and age for GOA/AI species

#########################
# Set directories and read in necessary data
path<-getwd()
pathD<-paste0(path,"/Data")
pathR<-paste0(path,"/SzAC results")

lfreq<-read.csv(paste0(pathD,"/lfreq.csv"))
#lfreq<-read.csv(paste0(pathD,"/lfreq_full.csv"))
sizepop_RACE<-read.csv(paste0(pathD,"/sizepop.csv"))
specimen<-read.csv(paste0(pathD,"/specimen_full.csv"))
specimen<-subset(specimen,is.na(specimen$AGE)==FALSE)
agepop_RACE<-read.csv(paste0(pathD,"/agepop.csv"))
cpue<-read.csv(paste0(pathD,"/CPUE.csv"))
strata<-read.csv(paste0(pathD,"/strata.csv"))

# Get data parameters together
species<-sort(unique(specimen$SPECIES_CODE))
yrs<-sort(unique(specimen$YEAR))
sex<-sort(unique(specimen$SEX))
stratum<-sort(unique(strata$STRATUM))

#########################
# Loop thru species
for(sp in 1:length(species)){
  
  lfreq_sp<-subset(lfreq,lfreq$SPECIES_CODE==species[sp])
  sizepop_RACE_sp<-subset(sizepop_RACE,sizepop_RACE$SPECIES_CODE==species[sp])
  specimen_sp<-subset(specimen,specimen$SPECIES_CODE==species[sp])
  agepop_RACE_sp<-subset(agepop_RACE,agepop_RACE$SPECIES_CODE==species[sp])
  cpue_sp<-subset(cpue,cpue$SPECIES_CODE==species[sp])
  
  # Set up results matrices (length then age comp)
  lengths<-sort(unique(lfreq_sp$LENGTH))
  SZPOP_M<-matrix(nrow=length(yrs),ncol=length(lengths))
  colnames(SZPOP_M)<-as.character(lengths)
  rownames(SZPOP_M)<-as.character(yrs)
  SZPOP_F<-matrix(nrow=length(yrs),ncol=length(lengths))
  colnames(SZPOP_F)<-as.character(lengths)
  rownames(SZPOP_F)<-as.character(yrs)
  SZPOP_U<-matrix(nrow=length(yrs),ncol=length(lengths))
  colnames(SZPOP_U)<-as.character(lengths)
  rownames(SZPOP_U)<-as.character(yrs)
  
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
  
  # Set up test matrix for comparison for size/age pop'n est with RACE output
  RACEmatch_Sz<-matrix(nrow=length(yrs),ncol=4)
  colnames(RACEmatch_Sz)<-c("Year","Male","Female","Unsexed")
  RACEmatch_Sz<-as.data.frame(RACEmatch_Sz)
  RACEmatch_Sz$Year<-yrs
  
  RACEmatch_Ag<-matrix(nrow=length(yrs),ncol=4)
  colnames(RACEmatch_Ag)<-c("Year","Male","Female","Unsexed")
  RACEmatch_Ag<-as.data.frame(RACEmatch_Ag)
  RACEmatch_Ag$Year<-yrs
  
  #########################
  # Loop thru years
  for(y in 1:length(yrs)){
    
    lfreq_sp_y<-subset(lfreq_sp,lfreq_sp$YEAR==yrs[y])
    cpue_sp_y<-subset(cpue_sp,cpue_sp$YEAR==yrs[y])
    sizepop_RACE_sp_y<-subset(sizepop_RACE_sp,sizepop_RACE_sp$YEAR==yrs[y])
    specimen_sp_y<-subset(specimen_sp,specimen_sp$YEAR==yrs[y])
    agepop_RACE_sp_y<-subset(agepop_RACE_sp,agepop_RACE_sp$SURVEY_YEAR==yrs[y])
    
    
    ########################################################################################
    # Estimate pop'n @ size
    ########################################################################################
    SZPOP_M_st<-matrix(nrow=length(stratum),ncol=length(lengths))
    rownames(SZPOP_M_st)<-stratum
    colnames(SZPOP_M_st)<-lengths
    SZPOP_F_st<-matrix(nrow=length(stratum),ncol=length(lengths))
    rownames(SZPOP_F_st)<-stratum
    colnames(SZPOP_F_st)<-lengths
    SZPOP_U_st<-matrix(nrow=length(stratum),ncol=length(lengths))
    rownames(SZPOP_U_st)<-stratum
    colnames(SZPOP_U_st)<-lengths

      #########################
      # Loop thru strata
      for(st in 1:length(stratum)){
        
        # Subset data to strata level
        strata_st<-subset(strata,strata$STRATUM==stratum[st])        
        cpue_sp_y_st<-subset(cpue_sp_y,cpue_sp_y$STRATUM==stratum[st])
        hls_cpue<-unique(cpue_sp_y_st$HAULJOIN)
        lfreq_sp_y_st<-subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls_cpue)

        # Subset data to sex-specific (M=males, F=females, U=unsexed)
        lfreq_sp_y_M_st<-subset(lfreq_sp_y_st,lfreq_sp_y_st$SEX==1)
        lfreq_sp_y_F_st<-subset(lfreq_sp_y_st,lfreq_sp_y_st$SEX==2)
        lfreq_sp_y_U_st<-subset(lfreq_sp_y_st,lfreq_sp_y_st$SEX==3)

        
        # Determine number of hauls of catch with lengths
        count<-length(unique(c(lfreq_sp_y_M_st$HAULJOIN,lfreq_sp_y_F_st$HAULJOIN,lfreq_sp_y_U_st$HAULJOIN)))
        
        # Identify hauls with catch but no lengths
        hls_l<-unique(c(lfreq_sp_y_M_st$HAULJOIN,lfreq_sp_y_F_st$HAULJOIN,lfreq_sp_y_U_st$HAULJOIN))
        hls_c<-cpue_sp_y_st$HAULJOIN[which(is.na(cpue_sp_y_st$CATCHJOIN)==FALSE)]
        hls_nol<-hls_c[which(is.na(match(hls_c,hls_l)==TRUE))]
      
        # Calc pop'n #'s in strata
        st_num<-mean(cpue_sp_y_st$NUMCPUE)*strata_st$AREA
        
        # Calc CPUE ratio among hauls
        cprat<-tapply(cpue_sp_y_st$NUMCPUE,cpue_sp_y_st$HAULJOIN,mean)/sum(cpue_sp_y_st$NUMCPUE)
        
        # Calc Total lengths sampled by haul
        n_st<-tapply(lfreq_sp_y_st$FREQUENCY,lfreq_sp_y_st$HAULJOIN,sum)
        
        # Calc sex-specific numbers at length
        n_h_M<-tapply(lfreq_sp_y_M_st$FREQUENCY,list(lfreq_sp_y_M_st$HAULJOIN,lfreq_sp_y_M_st$LENGTH),sum)
        n_h_M[is.na(n_h_M)] <- 0
        n_h_F<-tapply(lfreq_sp_y_F_st$FREQUENCY,list(lfreq_sp_y_F_st$HAULJOIN,lfreq_sp_y_F_st$LENGTH),sum)
        n_h_F[is.na(n_h_F)] <- 0
        n_h_U<-tapply(lfreq_sp_y_U_st$FREQUENCY,list(lfreq_sp_y_U_st$HAULJOIN,lfreq_sp_y_U_st$LENGTH),sum)
        n_h_U[is.na(n_h_U)] <- 0

        # Sex-specific ratio of total
        ratio_h_M<-n_h_M/as.vector(n_st[match(as.numeric(rownames(n_h_M)),as.numeric(names(n_st)))])
        ratio_h_F<-n_h_F/as.vector(n_st[match(as.numeric(rownames(n_h_F)),as.numeric(names(n_st)))])
        ratio_h_U<-n_h_U/as.vector(n_st[match(as.numeric(rownames(n_h_U)),as.numeric(names(n_st)))])

        if(length(hls_nol)>0){

        # Estimate size comp for hauls with catch that did not sample lengths
        ratio_h_M_unk<-colSums(ratio_h_M)/count
        ratio_h_F_unk<-colSums(ratio_h_F)/count
        ratio_h_U_unk<-colSums(ratio_h_U)/count
        total<-sum(ratio_h_M_unk,ratio_h_F_unk,ratio_h_U_unk)
        ratio_h_M_unk<-ratio_h_M_unk/total
        ratio_h_F_unk<-ratio_h_F_unk/total
        ratio_h_U_unk<-ratio_h_U_unk/total
        
        # Add unkown size com hauls to sex-specific ratio of total
        ratio_h_M_unk_add<-matrix(ratio_h_M_unk,nrow=length(hls_nol),ncol=length(ratio_h_M_unk),byrow=TRUE)
        rownames(ratio_h_M_unk_add)<-hls_nol
        colnames(ratio_h_M_unk_add)<-colnames(ratio_h_M)
        ratio_h_M<-rbind(ratio_h_M,ratio_h_M_unk_add)
        ratio_h_F_unk_add<-matrix(ratio_h_F_unk,nrow=length(hls_nol),ncol=length(ratio_h_F_unk),byrow=TRUE)
        rownames(ratio_h_F_unk_add)<-hls_nol
        colnames(ratio_h_F_unk_add)<-colnames(ratio_h_F)
        ratio_h_F<-rbind(ratio_h_F,ratio_h_F_unk_add)
        ratio_h_U_unk_add<-matrix(ratio_h_U_unk,nrow=length(hls_nol),ncol=length(ratio_h_U_unk),byrow=TRUE)
        rownames(ratio_h_U_unk_add)<-hls_nol
        colnames(ratio_h_U_unk_add)<-colnames(ratio_h_U)
        ratio_h_U<-rbind(ratio_h_U,ratio_h_U_unk_add)

        }
        
        # Put it all together to get numbers-at-sex-at-length by strata, and put it in results matrix
        szpop_M<-round(colSums(ratio_h_M*as.vector(cprat[match(as.numeric(rownames(ratio_h_M)),as.numeric(names(cprat)))])*st_num),digits=0)
        SZPOP_M_st[st,match(as.numeric(names(szpop_M)),lengths)]<-szpop_M
        szpop_F<-round(colSums(ratio_h_F*as.vector(cprat[match(as.numeric(rownames(ratio_h_F)),as.numeric(names(cprat)))])*st_num),digits=0)
        SZPOP_F_st[st,match(as.numeric(names(szpop_F)),lengths)]<-szpop_F
        szpop_U<-round(colSums(ratio_h_U*as.vector(cprat[match(as.numeric(rownames(ratio_h_U)),as.numeric(names(cprat)))])*st_num),digits=0)
        SZPOP_U_st[st,match(as.numeric(names(szpop_U)),lengths)]<-szpop_U

        # End stratum loop
      }

    
    # Now sum up across strata and see if it matches with RACE output
    
    # Males
    SZPOP_M_st[is.na(SZPOP_M_st)] <- 0
    SZPOP_M_y<-colSums(SZPOP_M_st)
    test_M_sz<-matrix(nrow=length(lengths),ncol=3)
    colnames(test_M_sz)<-c("Length","Calc","RACE")
    test_M_sz[,1]<-lengths
    test_M_sz[,2]<-SZPOP_M_y
    test_M_sz[match(sizepop_RACE_sp_y$LENGTH,lengths),3]<-sizepop_RACE_sp_y$MALES
    test_M_sz[is.na(test_M_sz)] <- 0
    test_M_sz<-as.data.frame(test_M_sz)
    RACEmatch_Sz$Male[y]<-max(abs(test_M_sz$Calc[which(test_M_sz$Calc>0)]-test_M_sz$RACE[which(test_M_sz$RACE>0)])/test_M_sz$RACE[which(test_M_sz$RACE>0)])
    SZPOP_M[y,]<-SZPOP_M_y
    
    # Females
    SZPOP_F_st[is.na(SZPOP_F_st)] <- 0
    SZPOP_F_y<-colSums(SZPOP_F_st)
    test_F_sz<-matrix(nrow=length(lengths),ncol=3)
    colnames(test_F_sz)<-c("Length","Calc","RACE")
    test_F_sz[,1]<-lengths
    test_F_sz[,2]<-SZPOP_F_y
    test_F_sz[match(sizepop_RACE_sp_y$LENGTH,lengths),3]<-sizepop_RACE_sp_y$FEMALES
    test_F_sz[is.na(test_F_sz)] <- 0
    test_F_sz<-as.data.frame(test_F_sz)
    RACEmatch_Sz$Female[y]<-max(abs(test_F_sz$Calc[which(test_F_sz$Calc>0)]-test_F_sz$RACE[which(test_F_sz$RACE>0)])/test_F_sz$RACE[which(test_F_sz$RACE>0)])
    SZPOP_F[y,]<-SZPOP_F_y
    
    # Unsexed
    SZPOP_U_st[is.na(SZPOP_U_st)] <- 0
    SZPOP_U_y<-colSums(SZPOP_U_st)
    test_U_sz<-matrix(nrow=length(lengths),ncol=3)
    colnames(test_U_sz)<-c("Length","Calc","RACE")
    test_U_sz[,1]<-lengths
    test_U_sz[,2]<-SZPOP_U_y
    test_U_sz[match(sizepop_RACE_sp_y$LENGTH,lengths),3]<-sizepop_RACE_sp_y$UNSEXED
    test_U_sz[is.na(test_U_sz)] <- 0
    test_U_sz<-as.data.frame(test_U_sz)
    RACEmatch_Sz$Unsexed[y]<-max(abs(test_U_sz$Calc[which(test_U_sz$Calc>0)]-test_U_sz$RACE[which(test_U_sz$RACE>0)])/test_U_sz$RACE[which(test_U_sz$RACE>0)])
    SZPOP_U[y,]<-SZPOP_U_y
    
    # Set Size pop'n ests for age pop'n est script
    sizepop_sp_y<-cbind(lengths,SZPOP_M_y,SZPOP_F_y,SZPOP_U_y)
    colnames(sizepop_sp_y)<-c("LENGTH","MALES","FEMALES","UNSEXED")
    sizepop_sp_y<-as.data.frame(sizepop_sp_y)
    
    
    ########################################################################################
    # Estimate pop'n @ age
    ########################################################################################
    
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
      RACEmatch_Ag$Male[y]<-max(abs(test_M$Calc-test_M$RACE)/test_M$RACE)
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
      RACEmatch_Ag$Female[y]<-max(abs(test_F$Calc-test_F$RACE)/test_F$RACE)
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
      RACEmatch_Ag$Unsexed[y]<-max(abs(test_U$Calc-test_U$RACE)/test_U$RACE)
      AGEPOP_U[y,match(as.numeric(names(age_est_U)),ages)]<-age_est_U
      AGEPOP_U[y,is.na(AGEPOP_U[y,])] <- 0
    }}
    
    
    
    #########################
    # End year loop
  }
  
  
  # Test if there's disagreement with RACE calcs (flag if larger than 5% for the maximum difference at size)
  if(length(which(RACEmatch_Sz$Male > 0.05)) > 0)
    cat(paste("Disagreement between calcs and RACE for Male Size pop'n for species", species[sp], "\n"))
  if(length(which(RACEmatch_Sz$Female > 0.05)) > 0)
    cat(paste("Disagreement between calcs and RACE for Female Size pop'n for species", species[sp], "\n"))
  if(length(which(RACEmatch_Sz$Unsexed > 0.05)) > 0)
    cat(paste("Disagreement between calcs and RACE for Unsexed Size pop'n for species", species[sp], "\n"))
  
  # Test if there's disagreement with RACE calcs (flag if larger than 5% for the maximum difference at age)
  if(length(which(RACEmatch_Ag$Male > 0.05)) > 0)
    cat(paste("Disagreement between calcs and RACE for Male Age pop'n for species", species[sp], "\n"))
  if(length(which(RACEmatch_Ag$Female > 0.05)) > 0)
    cat(paste("Disagreement between calcs and RACE for Female Age pop'n for species", species[sp], "\n"))
  if(length(which(RACEmatch_Ag$Unsexed > 0.05)) > 0)
    cat(paste("Disagreement between calcs and RACE for Unsexed Age pop'n for species", species[sp], "\n"))
  
  # Finalize results matrices
  
  SZPOP_T <- matrix(mapply(sum,SZPOP_M,SZPOP_F,SZPOP_U,MoreArgs=list(na.rm=TRUE)),ncol=length(lengths))
  colnames(SZPOP_T)<-as.character(lengths)
  rownames(SZPOP_T)<-as.character(yrs)
  
  AGEPOP_T <- matrix(mapply(sum,AGEPOP_M,AGEPOP_F,AGEPOP_U,MoreArgs=list(na.rm=TRUE)),ncol=length(ages))
  colnames(AGEPOP_T)<-as.character(ages)
  rownames(AGEPOP_T)<-as.character(yrs)
  
  # Write results matrices
  write.csv(RACEmatch_Sz,paste0(pathR,"/RACEmatch_Sz_",species[sp],".csv"))
  write.csv(SZPOP_M,paste0(pathR,"/SZPOP_M_",species[sp],".csv"))
  write.csv(SZPOP_F,paste0(pathR,"/SZPOP_F_",species[sp],".csv"))
  write.csv(SZPOP_U,paste0(pathR,"/SZPOP_U_",species[sp],".csv"))
  write.csv(SZPOP_T,paste0(pathR,"/SZPOP_T_",species[sp],".csv"))
  
  write.csv(RACEmatch_Ag,paste0(pathR,"/RACEmatch_Ag_",species[sp],".csv"))
  write.csv(AGEPOP_M,paste0(pathR,"/AGEPOP_M_",species[sp],".csv"))
  write.csv(AGEPOP_F,paste0(pathR,"/AGEPOP_F_",species[sp],".csv"))
  write.csv(AGEPOP_U,paste0(pathR,"/AGEPOP_U_",species[sp],".csv"))
  write.csv(AGEPOP_T,paste0(pathR,"/AGEPOP_T_",species[sp],".csv"))
  
  #########################
  # End species loop 
}









