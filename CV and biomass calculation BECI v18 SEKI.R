# BECI CUBIC VOLUME AND BIOMASS EQUATIONS

# Summary: Calculates cubic volumn and biomass for trees, 
# using BECI equations. 

# Workflow: calculate CV for every tree. Calculate biomass for 
# bole (stem), bark, and branches.  

# Several Helper functions are defined below, but to calculate tree
# Biomass and Cubic Volumn information, all you need to do is run
# the function "Biomass" in the "MAIN" section. 
# Several diagnostic equations are also provided, which you can run
# after running the "Biomass" function. 

#################################
####### HELPER FUNTIONS #########
#################################

# Parses data file given DATA, CODEFILE, and DIR.  
# Returns DATA, with a new column "SPP6" where large redwoods are 
# speparated from smaller redwoods because they have different 
# bark biomass equations. SESE and SEGI redwoods are renamed to 
# SESE.l/SEGI.l if they are large (DBH > 100 cm). 
# Also performs unit conversions when necessary, as some 
# equations require metric while others require imperial units.
# DATA file is merged with CODEFILE, so that the returned DATA file 
# contains equation information (equation numbers, density, species
# common name, etc.) for each tree. 

ParseData <- function(data, codefile, dir = "") {
  #set's working directory to dir 
  if(!(dir == "")){setwd(dir)}
  infile <- read_csv(data) 
  codefile <- read_csv(codefile)
  #performs unit conversions where necessary
  if(!("Ht_ft" %in% colnames(infile))) {infile$Ht_ft <- infile$Ht_m/0.3048
  } else if (!("Ht_m" %in% colnames(infile))) {infile$Ht_m<-infile$Ht_ft*0.3048
  }
  if(!("DBH_in" %in% colnames(infile))) {infile$DBH_in <- infile$DBH_cm/2.54
  } else if (!("DBH_cm" %in% colnames(infile))) {infile$DBH_cm<-infile$DBH_in*2.54
  }
  if(!("SPP4" %in% colnames(infile))) {infile$SPP4<-infile$Species
  }
  #Identifies redwoods < 100 cm DBH
  redtest<-infile$SPP4 =="SEGI" | infile$SPP4 == "SESE"          
  dbhtest<-infile$DBH_cm>100                                     
  test<-redtest&dbhtest
  codefile$SPP4 = NULL 
  infile<-mutate(infile, SPP6 = as.factor(if_else(test, paste(SPP4,".l",sep=""),SPP4)))
  #merges code file with datafile
  mergefile<-merge(infile, codefile, by.x = "SPP6",by.y = "SPP6")
  return(mergefile)
}


# Calculates tree cubic volumn using BECI given parsed DATA. 
# Returns DATA with a new column "CVTS_ft" containing the estimated 
# cubic volumn for each tree, in cubic feet.

CalcCVTS<-function(data){
  #Separates trees according to their equation number
  for (i in  c(3, 5, 6, 8, 16, 17, 18, 19, 20, 21, 23, 24, 26, 28, 32, 33, 34, 37, 38, 40, 42, 43)) {
    assign(paste("subset.eq", i, sep = ""), subset(data, data$Eq_CV == i))
  }
  #Helper function that calculates TMP_DBH, BA, and BA_TMP for small trees (trees < 6 inches in DBH)
  SmallBATMP<-function(table) {
    table$TMP_DBH <- 6.0
    table$BA <- table$DBH_in^2 * 0.005454154
    table$BA_TMP <- (6^2)*0.005454154
    return(table)
  } 
  #Small trees (< 6 inch DBH)
  #Calculates cubic volumn for small Dougals firs using equation 3
  
  small3 <- subset(subset.eq3, subset.eq3$DBH_in < 6.0)
  if (dim(small3)[1]>0) {
    small3 = SmallBATMP(small3)
    small3$CF4_TMPcalc <- 0.248569 + 0.0253524 * (small3$Ht_ft/small3$TMP_DBH) - 0.0000560175 * (small3$Ht_ft^2/small3$TMP_DBH)
    small3[,"CF4_TMP"] <- NA
    for (i in 1:nrow(small3)){
      if(small3$CF4_TMPcalc[i] < 0.3) {small3$CF4_TMP[i] = 0.3
      } else if (small3$CF4_TMPcalc[i] > 0.4) {small3$CF4_TMP[i] = 0.4
      } else {small3$CF4_TMP[i] = small3$CF4_TMPcalc[i]
      }
    }
    small3$CV4_TMP <- small3$CF4_TMP * small3$BA_TMP * small3$Ht_ft 
    small3$Tarif_TMPcalc <- (small3$CV4_TMP * 0.912733) / (small3$BA_TMP - 0.087266)
    small3[,"Tarif_TMP"] <- NA
    for (i in 1:nrow(small3)){
      if(small3$Tarif_TMPcalc[i] <= 0.0) {small3$Tarif_TEMP[i] = 0.01
      } else {small3$Tarif_TMP[i] = small3$Tarif_TMPcalc[i]
      }
    }
    small3$TARIFcacl <- small3$Tarif_TMP * (0.5 * (6.0-small3$DBH_in)^2 + (1.0 + 0.063 * (6.0 - small3$DBH_in)^2))
    small3[,"TARIF"] <- NA
    for (i in 1:nrow(small3)){
      if(small3$TARIFcacl[i] <= 0.0) {small3$TARIF[i] = 0.01
      } else {small3$TARIF[i] = small3$TARIFcacl[i]
      }
    }
    small3$TERM <- ((1.033 * (1.0 + 1.382937 * exp(-4.015292 * (small3$DBH_in/10.0)))) * (small3$BA + 0.087266) - 0.174533)
    small3$CVTS_ft <- small3$TARIF * small3$TERM
  } #end if small Eq 3
  #Calculates CV for small Jefferey, Ponderosa and Foothill pines using equation 5
  small5 <- subset(subset.eq5, subset.eq5$DBH_in < 6.0)
  if (dim(small5)[1]>0){
    small5 = SmallBATMP(small5)
    small5$CF4_TMPcalc <- 0.402060 - 0.899914 * (1/small5$TMP_DBH)
    small5[,"CF4_TMP"] <- NA
    for (i in 1:nrow(small5)){
      if(small5$CF4_TMPcalc[i] < 0.3) {small5$CF4_TMP[i] = 0.3
      } else if (small5$CF4_TMPcalc[i] > 0.4) {small5$CF4_TMP[i] = 0.4
      } else {small5$CF4_TMP[i] = small5$CF4_TMPcalc[i]
      }
    }
    small5$CV4_TMP <- small5$CF4_TMP * small5$BA_TMP * small5$Ht_ft 
    small5$Tarif_TMPcalc <- (small5$CV4_TMP * 0.912733) / (small5$BA_TMP - 0.087266)
    small5[,"Tarif_TMP"] <- NA

    for (i in 1:nrow(small5)){
      if(small5$Tarif_TMPcalc[i] <= 0.0) {small5$Tarif_TEMP[i] = 0.01
      } else {small5$Tarif_TMP[i] = small5$Tarif_TMPcalc[i]
      }
    }
    small5$TARIFcacl <- small5$Tarif_TMP * (0.5 * (6.0-small5$DBH_in)^2 + (1.0 + 0.063 * (6.0 - small5$DBH_in)^2))
    small5[,"TARIF"] <- NA
    for (i in 1:nrow(small5)){
      if(small5$TARIFcacl[i] <= 0.0) {small5$TARIF[i] = 0.01
      } else {small5$TARIF[i] = small5$TARIFcacl[i]
      }
    }
    small5$TERM <- ((1.033 * (1.0 + 1.382937 * exp(-4.015292 * (small5$DBH_in/10.0)))) * (small5$BA + 0.087266) - 0.174533)
    small5$CVTS_ft <- small5$TARIF * small5$TERM
  } #end if small Eq 5

  #Calculates CV for small Lodgepole pines using equation 16
  small16 <- subset(subset.eq16, subset.eq16$DBH_in < 6.0)
  if (dim(small16)[1]>0){
    small16 = SmallBATMP(small16)
    small16$CF4_TMPcalc <- 0.422709 - 0.0000612236 * (small16$Ht_ft^2/small16$TMP_DBH)
    small16[,"CF4_TMP"] <- NA
    for (i in 1:nrow(small16)){
      if(small16$CF4_TMPcalc[i] < 0.3) {small16$CF4_TMP[i] = 0.3
      } else if (small16$CF4_TMPcalc[i] > 0.4) {small16$CF4_TMP[i] = 0.4
      } else {small16$CF4_TMP[i] = small16$CF4_TMPcalc[i]
      }
    }
    small16$CV4_TMP <- small16$CF4_TMP * small16$BA_TMP * small16$Ht_ft 
    small16$Tarif_TMPcalc <- (small16$CV4_TMP * 0.912733) / (small16$BA_TMP - 0.087266)
    small16[,"Tarif_TMP"] <- NA
    for (i in 1:nrow(small16)){
      if(small16$Tarif_TMPcalc[i] <= 0.0) {small16$Tarif_TEMP[i] = 0.01
      } else {small16$Tarif_TMP[i] = small16$Tarif_TMPcalc[i]
      }
    }
    small16$TARIFcacl <- small16$Tarif_TMP * (0.5 * (6.0-small16$DBH_in)^2 + (1.0 + 0.063 * (6.0 - small16$DBH_in)^2))
    small16[,"TARIF"] <- NA
    for (i in 1:nrow(small16)){
      if(small16$TARIFcacl[i] <= 0.0) {small16$TARIF[i] = 0.01
      } else {small16$TARIF[i] = small16$TARIFcacl[i]
      }
    }
    small16$TERM <- ((1.033 * (1.0 + 1.382937 * exp(-4.015292 * (small16$DBH_in/10.0)))) * (small16$BA + 0.087266) - 0.174533)
    small16$CVTS_ft <- small16$TARIF * small16$TERM
  } # end of small Eq 16
  #Calculates CV for small California Red Firs and Noble Firs using equation 18
  small18 <- subset(subset.eq18, subset.eq18$DBH_in < 6.0)
  if (dim(small18)[1]>0) {
    small18 = SmallBATMP(small18)
    small18$CF4_TMPcalc <- 0.231237 + 0.028176 * (small18$Ht_ft/small18$TMP_DBH)
    small18[,"CF4_TMP"] <- NA
    for (i in 1:nrow(small18)){
      if(small18$CF4_TMPcalc[i] < 0.3) {small18$CF4_TMP[i] = 0.3
      } else if (small18$CF4_TMPcalc[i] > 0.4) {small18$CF4_TMP[i] = 0.4
      } else {small18$CF4_TMP[i] = small18$CF4_TMPcalc[i]
      }
    }
    small18$CV4_TMP <- small18$CF4_TMP * small18$BA_TMP * small18$Ht_ft 
    small18$Tarif_TMPcalc <- (small18$CV4_TMP * 0.912733) / (small18$BA_TMP - 0.087266)
    small18[,"Tarif_TMP"] <- NA
    for (i in 1:nrow(small18)){
      if(small18$Tarif_TMPcalc[i] <= 0.0) {small18$Tarif_TEMP[i] = 0.01
      } else {small18$Tarif_TMP[i] = small18$Tarif_TMPcalc[i]
      }
    }
    small18$TARIFcacl <- small18$Tarif_TMP * (0.5 * (6.0-small18$DBH_in)^2 + (1.0 + 0.063 * (6.0 - small18$DBH_in)^2))
    small18[,"TARIF"] <- NA
    for (i in 1:nrow(small18)){
      if(small18$TARIFcacl[i] <= 0.0) {small18$TARIF[i] = 0.01
      } else {small18$TARIF[i] = small18$TARIFcacl[i]
      }
    }
    small18$TERM <- ((1.033 * (1.0 + 1.382937 * exp(-4.015292 * (small18$DBH_in/10.0)))) * (small18$BA + 0.087266) - 0.174533)
    small18$CVTS_ft <- small18$TARIF * small18$TERM
  } # end of small Eq 18

  #Calculates CV for small Incense Cedars using equation 19
  small19 <- subset(subset.eq19, subset.eq19$DBH_in < 6.0)
  if (dim(small19)[1]>0){
    small19 = SmallBATMP(small19)
    small19$CF4_TMPcalc <- 0.225786 + 4.44236 * (1/small19$Ht_ft)
    small19[,"CF4_TMP"] <- NA
    for (i in 1:nrow(small19)){
      if(small19$CF4_TMPcalc[i] < 0.27) {small19$CF4_TMP[i] = 0.27
      } else {small19$CF4_TMP[i] = small19$CF4_TMPcalc[i]
      }
    }
    small19$CV4_TMP <- small19$CF4_TMP * small19$BA_TMP * small19$Ht_ft 
    small19$Tarif_TMPcalc <- (small19$CV4_TMP * 0.912733) / (small19$BA_TMP - 0.087266)
    small19[,"Tarif_TMP"] <- NA
    for (i in 1:nrow(small19)){
      if(small19$Tarif_TMPcalc[i] <= 0.0) {small19$Tarif_TEMP[i] = 0.01
      } else {small19$Tarif_TMP[i] = small19$Tarif_TMPcalc[i]
      }
    }
    small19$TARIFcacl <- small19$Tarif_TMP * (0.5 * (6.0-small19$DBH_in)^2 + (1.0 + 0.063 * (6.0 - small19$DBH_in)^2))
    small19[,"TARIF"] <- NA
    for (i in 1:nrow(small19)){
      if(small19$TARIFcacl[i] <= 0.0) {small19$TARIF[i] = 0.01
      } else {small19$TARIF[i] = small19$TARIFcacl[i]
      }
    }
    small19$TERM <- ((1.033 * (1.0 + 1.382937 * exp(-4.015292 * (small19$DBH_in/10.0)))) * (small19$BA + 0.087266) - 0.174533)
    small19$CVTS_ft <- small19$TARIF * small19$TERM
  } # end if small Eq 19
  #Calculates CV for small sugar pines and small western white pines using equation 20
  small20 <- subset(subset.eq20, subset.eq20$DBH_in < 6.0)

  if (dim (small20)[1]>0){
    small20 = SmallBATMP(small20)
    small20$CF4_TMPcalc <- 0.358550 - 0.488134 * (1/small20$TMP_DBH)
    small20[,"CF4_TMP"] <- NA
    for (i in 1:nrow(small20)){
      if(small20$CF4_TMPcalc[i] < 0.3) {small20$CF4_TMP[i] = 0.3
      } else if (small20$CF4_TMPcalc[i] > 0.4) {small20$CF4_TMP[i] = 0.4
      } else {small20$CF4_TMP[i] = small20$CF4_TMPcalc[i]
      }
    }
    small20$CV4_TMP <- small20$CF4_TMP * small20$BA_TMP * small20$Ht_ft 
    small20$Tarif_TMPcalc <- (small20$CV4_TMP * 0.912733) / (small20$BA_TMP - 0.087266)
    small20[,"Tarif_TMP"] <- NA
    for (i in 1:nrow(small20)){
      if(small20$Tarif_TMPcalc[i] <= 0.0) {small20$Tarif_TEMP[i] = 0.01
      } else {small20$Tarif_TMP[i] = small20$Tarif_TMPcalc[i]
      }
    }
    small20$TARIFcacl <- small20$Tarif_TMP * (0.5 * (6.0-small20$DBH_in)^2 + (1.0 + 0.063 * (6.0 - small20$DBH_in)^2))
    small20[,"TARIF"] <- NA
    for (i in 1:nrow(small20)){
      if(small20$TARIFcacl[i] <= 0.0) {small20$TARIF[i] = 0.01
      } else {small20$TARIF[i] = small20$TARIFcacl[i]
      }
    }
    small20$TERM <- ((1.033 * (1.0 + 1.382937 * exp(-4.015292 * (small20$DBH_in/10.0)))) * (small20$BA + 0.087266) - 0.174533)
    small20$CVTS_ft <- small20$TARIF * small20$TERM
  } # end of small Eq 20

  #Calculates CV for small white firs and grand firs using equation 23
  small23 <- subset(subset.eq23, subset.eq23$DBH_in < 6.0)
  if (dim(small23)[1]>0){
    small23 = SmallBATMP(small23)
    small23$CF4_TMPcalc <- 0.299039 + 1.91272 * (1/small23$Ht_ft) + 0.0000367217 * (small23$Ht_ft^2/small23$TMP_DBH)
    small23[,"CF4_TMP"] <- NA
    for (i in 1:nrow(small23)){
      if(small23$CF4_TMPcalc[i] < 0.3) {small23$CF4_TMP[i] = 0.3
      } else if (small23$CF4_TMPcalc[i] > 0.4) {small23$CF4_TMP[i] = 0.4  #
      } else {small23$CF4_TMP[i] = small23$CF4_TMPcalc[i]
      }
    }
    small23$CV4_TMP <- small23$CF4_TMP * small23$BA_TMP * small23$Ht_ft 
    small23$Tarif_TMPcalc <- (small23$CV4_TMP * 0.912733) / (small23$BA_TMP - 0.087266)
    small23[,"Tarif_TMP"] <- NA
    for (i in 1:nrow(small23)){
      if(small23$Tarif_TMPcalc[i] <= 0.0) {small23$Tarif_TEMP[i] = 0.01
      } else {small23$Tarif_TMP[i] = small23$Tarif_TMPcalc[i]
      }
    }
    small23$TARIFcacl <- small23$Tarif_TMP * (0.5 * (6.0-small23$DBH_in)^2 + (1.0 + 0.063 * (6.0 - small23$DBH_in)^2))
    small23[,"TARIF"] <- NA
    for (i in 1:nrow(small23)){
      if(small23$TARIFcacl[i] <= 0.0) {small23$TARIF[i] = 0.01
      } else {small23$TARIF[i] = small23$TARIFcacl[i]
      }
    }
    small23$TERM <- ((1.033 * (1.0 + 1.382937 * exp(-4.015292 * (small23$DBH_in/10.0)))) * (small23$BA + 0.087266) - 0.174533)
    small23$CVTS_ft <- small23$TARIF * small23$TERM
  } #end if small Eq 23 


# Large Trees (>6 inches DBH)
#Calculates CV for large Douglas firs using equation 3
  large3 <- subset(subset.eq3, subset.eq3$DBH_in >= 6.0)
  if (dim(large3)[1]>0){
    large3$BA <- (large3$DBH_in^2)*0.005454154
    large3$CF4calc <- 0.248569 + 0.0253524 * (large3$Ht_ft/large3$DBH_in) - 0.0000560175 * (large3$Ht_ft^2/large3$DBH_in)
    large3[,"CF4"] <- NA
    for (i in 1:nrow(large3)){
      if(large3$CF4calc[i] < 0.3) {large3$CF4[i] = 0.3
      } else if (large3$CF4calc[i] > 0.4) {large3$CF4[i] = 0.4
      } else {large3$CF4[i] = large3$CF4calc[i]
      }
    }
    large3$CV4 <- large3$CF4 * large3$BA * large3$Ht_ft
    large3$Tarifcalc <- (large3$CV4 * 0.912733) / (large3$BA - 0.087266)
    large3[,"TARIF"] <- NA
    for (i in 1:nrow(large3)){
      if(large3$Tarifcalc[i] <= 0.0) {large3$TARIF[i] = 0.01
      } else {large3$TARIF[i] = large3$Tarifcalc[i]
      }
    }
    large3$TERM <- ((1.033 * (1.0 + 1.382937 * exp(-4.015292 * (large3$DBH_in/10.0)))) * (large3$BA + 0.087266) - 0.174533)
    large3$CVTS_ft <- (large3$CV4 * large3$TERM)/ (large3$BA - 0.087266)
  } #end if large Eq 3

#Calculates CV for large Jefferey, Ponderosa, and white pines (eq. 5)
  large5 <- subset(subset.eq5, subset.eq5$DBH_in >= 6.0)
  if (dim(large5)[1]>0){
    large5$BA <- (large5$DBH_in^2)*0.005454154
    large5$CF4calc <- 0.402060 - 0.899914 * (1/large5$DBH_in)
    large5[,"CF4"] <- NA
    for (i in 1:nrow(large5)){
      if(large5$CF4calc[i] < 0.3) {large5$CF4[i] = 0.3
      } else if (large5$CF4calc[i] > 0.4) {large5$CF4[i] = 0.4
      } else {large5$CF4[i] = large5$CF4calc[i]
      }
    }
    large5$CV4 <- large5$CF4 * large5$BA * large5$Ht_ft
    large5$Tarifcalc <- (large5$CV4 * 0.912733) / (large5$BA - 0.087266)
    large5[,"TARIF"] <- NA
    for (i in 1:nrow(large5)){
      if(large5$Tarifcalc[i] <= 0.0) {large5$TARIF[i] = 0.01
      } else {large5$TARIF[i] = large5$Tarifcalc[i]
      }
    }
    large5$TERM <- ((1.033 * (1.0 + 1.382937 * exp(-4.015292 * (large5$DBH_in/10.0)))) * (large5$BA + 0.087266) - 0.174533)
    large5$CVTS_ft <- (large5$CV4 * large5$TERM)/ (large5$BA - 0.087266)
  } #end if large Eq 5

  #Calculates CV for large Lodgepole pines (eq. 16)
  large16 <- subset(subset.eq16, subset.eq16$DBH_in >= 6.0)
  if (dim(large16)[1]>0){
    large16$BA <- (large16$DBH_in^2)*0.005454154
    large16$CF4calc <- 0.422709 - 0.0000612236 * (large16$Ht_ft^2/large16$DBH_in)
    large16[,"CF4"] <- NA
    for (i in 1:nrow(large16)){
      if(large16$CF4calc[i] < 0.3) {large16$CF4[i] = 0.3
      } else if (large16$CF4calc[i] > 0.4) {large16$CF4[i] = 0.4
      } else {large16$CF4[i] = large16$CF4calc[i]
      }
    }
    large16$CV4 <- large16$CF4 * large16$BA * large16$Ht_ft
    large16$Tarifcalc <- (large16$CV4 * 0.912733) / (large16$BA - 0.087266)
    large16[,"TARIF"] <- NA
    for (i in 1:nrow(large16)){
      if(large16$Tarifcalc[i] <= 0.0) {large16$TARIF[i] = 0.01
      } else {large16$TARIF[i] = large16$Tarifcalc[i]
      }
    }
    large16$TERM <- ((1.033 * (1.0 + 1.382937 * exp(-4.015292 * (large16$DBH_in/10.0)))) * (large16$BA + 0.087266) - 0.174533)
    large16$CVTS_ft <- (large16$CV4 * large16$TERM)/ (large16$BA - 0.087266)
  } #end of large Eq 16
  
  #Calculates CV for large California Red Firs and Noble Firs (Eq. 18)
  large18 <- subset(subset.eq18, subset.eq18$DBH_in >= 6.0)
  if (dim(large18)[1]>0){
  
    large18$BA <- (large18$DBH_in^2)*0.005454154
    large18$CF4calc <- 0.231237 + 0.028176 * (large18$Ht_ft/large18$DBH_in)
    large18[,"CF4"] <- NA
    for (i in 1:nrow(large18)){
      if(large18$CF4calc[i] < 0.3) {large18$CF4[i] = 0.3
      } else if (large18$CF4calc[i] > 0.4) {large18$CF4[i] = 0.4
      } else {large18$CF4[i] = large18$CF4calc[i]
      }
    }
    large18$CV4 <- large18$CF4 * large18$BA * large18$Ht_ft
    large18$Tarifcalc <- (large18$CV4 * 0.912733) / (large18$BA - 0.087266)
    large18[,"TARIF"] <- NA
    for (i in 1:nrow(large18)){
      if(large18$Tarifcalc[i] <= 0.0) {large18$TARIF[i] = 0.01
      } else {large18$TARIF[i] = large18$Tarifcalc[i]
      }
    }
    large18$TERM <- ((1.033 * (1.0 + 1.382937 * exp(-4.015292 * (large18$DBH_in/10.0)))) * (large18$BA + 0.087266) - 0.174533)
    large18$CVTS_ft <- (large18$CV4 * large18$TERM) / (large18$BA - 0.087266)
    } #end if large Eq 18
    
    #Calculates CV for large Incense Cedars (eq. 19)
    large19 <- subset(subset.eq19, subset.eq19$DBH_in >= 6.0)
    if (dim(large19)[1]>0){
    large19$BA <- (large19$DBH_in^2)*0.005454154
    large19$CF4calc <- 0.225786 + 4.44236 * (1/large19$Ht_ft)
    large19[,"CF4"] <- NA
    for (i in 1:nrow(large19)){
      if(large19$CF4calc[i] < 0.27) {large19$CF4[i] = 0.27
      } else {large19$CF4[i] = large19$CF4calc[i]
      }
    }
    large19$CV4 <- large19$CF4 * large19$BA * large19$Ht_ft
    large19$Tarifcalc <- (large19$CV4 * 0.912733) / (large19$BA - 0.087266)
    large19[,"TARIF"] <- NA
    for (i in 1:nrow(large19)){
      if(large19$Tarifcalc[i] <= 0.0) {large19$TARIF[i] = 0.01
      } else {large19$TARIF[i] = large19$Tarifcalc[i]
      }
    }
    large19$TERM <- ((1.033 * (1.0 + 1.382937 * exp(-4.015292 * (large19$DBH_in/10.0)))) * (large19$BA + 0.087266) - 0.174533)
    large19$CVTS_ft <- (large19$CV4 * large19$TERM) / (large19$BA - 0.087266)
  } # end if large Eq 19
  
  #Calculates CV for large Sugar pines and Western white pines (eq. 20)
  large20 <- subset(subset.eq20, subset.eq20$DBH_in >= 6.0)
  if (dim(large20)[1]>0){
    large20$BA <- (large20$DBH_in^2)*0.005454154
    large20$CF4calc <- 0.358550 - 0.488134 * (1/large20$DBH_in)
    large20[,"CF4"] <- NA
    for (i in 1:nrow(large20)){
      if(large20$CF4calc[i] < 0.3) {large20$CF4[i] = 0.3
      } else if (large20$CF4calc[i] > 0.4) {large20$CF4[i] = 0.4
      } else {large20$CF4[i] = large20$CF4calc[i]
      }
    }
    large20$CV4 <- large20$CF4 * large20$BA * large20$Ht_ft
    large20$Tarifcalc <- (large20$CV4 * 0.912733) / (large20$BA - 0.087266)
    large20[,"TARIF"] <- NA
    for (i in 1:nrow(large20)){
      if(large20$Tarifcalc[i] <= 0.0) {large20$TARIF[i] = 0.01
      } else {large20$TARIF[i] = large20$Tarifcalc[i]
      }
    }
    large20$TERM <- ((1.033 * (1.0 + 1.382937 * exp(-4.015292 * (large20$DBH_in/10.0)))) * (large20$BA + 0.087266) - 0.174533)
    large20$CVTS_ft <- (large20$CV4 * large20$TERM) / (large20$BA - 0.087266)
  } # end if large Eq 20
  
  #Calculates CV for large White Firs and Grand Firs (eq. 23)
  large23 <- subset(subset.eq23, subset.eq23$DBH_in >= 6.0)
  if (dim(large23)[1]>0){
    large23$BA <- (large23$DBH_in^2)*0.005454154
    large23$CF4calc <- 0.299039 + 1.91272 * (1/large23$Ht_ft) + 0.0000367217 * (large23$Ht_ft^2/large23$DBH_in)
    large23[,"CF4"] <- NA
    for (i in 1:nrow(large23)){
      if(large23$CF4calc[i] < 0.3) {large23$CF4[i] = 0.3
      } else if (large23$CF4calc[i] > 0.4) {large23$CF4[i] = 0.4
      } else {large23$CF4[i] = large23$CF4calc[i]
      }
    }
    
    large23$CV4 <- large23$CF4 * large23$BA * large23$Ht_ft
    large23$Tarifcalc <- (large23$CV4 * 0.912733) / (large23$BA - 0.087266)
    large23[,"TARIF"] <- NA
    for (i in 1:nrow(large23)){
      if(large23$Tarifcalc[i] <= 0.0) {large23$TARIF[i] = 0.01
      } else {large23$TARIF[i] = large23$Tarifcalc[i]
      }
    }
    large23$TERM <- ((1.033 * (1.0 + 1.382937 * exp(-4.015292 * (large23$DBH_in/10.0)))) * (large23$BA + 0.087266) - 0.174533)
    large23$CVTS_ft <- (large23$CV4 * large23$TERM) / (large23$BA - 0.087266)
  } # end if large Eq 23
  
  #For equations with 1 step
  
  #Calculates CV for Western Hemlock using equation 6
  subset.eq6$CVTS_ft <- 10^(-2.72170 + (2.00857*log10(subset.eq6$DBH_in)) + (1.08620*log10(subset.eq6$Ht_ft)) - (0.00568*subset.eq6$DBH_in))
  
  #Calculates CV for California Nutmeg and Pacific Yew using equation 8
  subset.eq8$CVTS_ft <- 10^(-2.464614 + (1.701993*log10(subset.eq8$DBH_in)) + (1.067038*log10(subset.eq8$Ht_ft)))
  
  #Calculates CV for Unknown Conifers and mountain hemlock using equation 17
  subset.eq17$CVTS_ft <- 0.001106485*(subset.eq17$DBH_in)^1.8140497*(subset.eq17$Ht_ft)^1.2744923
  
  #Calculates CV for Western Juniper using equation 21
  subset.eq21$CVTS_ft <- 0.005454154 * (0.30708901 + 0.00086157622 * subset.eq21$Ht_ft - 0.0037255243 * subset.eq21$DBH_in * 
                                          (subset.eq21$Ht_ft/(subset.eq21$Ht_ft-4.5))) * subset.eq21$DBH_in^2 * subset.eq21$Ht_ft *
                                             (subset.eq21$Ht_ft/(subset.eq21$Ht_ft-4.5))^2
  
  #Calculates CV for large and small giant sequoias using equation 24
  subset.eq24$CVTS_ft <- exp(-6.2597 + 1.9967 * log(subset.eq24$DBH_in) + 0.9642 * log(subset.eq24$Ht_ft))
  
  #Calculates CV for Pacific Dogwood, White Alder, and unknown hardwoods, using equation 26
  subset.eq26$CVTS_ft <- 10^(-2.672775 + (1.920617*log10(subset.eq26$DBH_in)) + (1.074024*log10(subset.eq26$Ht_ft)))
  
  #Calculates CV for quaking aspen, using equation 28
  subset.eq28$CVTS_ft <- 10^(-2.635360 + (1.946034*log10(subset.eq28$DBH_in)) + (1.024793*log10(subset.eq28$Ht_ft)))
  
  #Calculates CV for Golden Chinkapin using equation equation 32
  subset.eq32$CVTS_ft <- 0.0120372263*(subset.eq32$DBH_in)^2.02232*(subset.eq32$Ht_ft)^0.68638
  
  #Calculates CV for Bay Laurel using equation 33
  subset.eq33$CVTS_ft <- 0.0057821322*(subset.eq33$DBH_in)^1.94553*(subset.eq33$Ht_ft)^0.88389
  
  #Calculates CV for Tan Oak using equation 34
  subset.eq34$CVTS_ft <- 0.0058870024*(subset.eq34$DBH_in)^1.94165*(subset.eq34$Ht_ft)^0.86562
  
  #Calculates CV for Bigleaf maple using equation 37
  subset.eq37$CVTS_ft <- 0.0101786350*(subset.eq37$DBH_in)^2.22462*(subset.eq37$Ht_ft)^0.575619
  
  #Calculates CV for California Black Oak using equation 38
  subset.eq38$CVTS_ft <- 0.0070538108*(subset.eq38$DBH_in)^1.97437*(subset.eq38$Ht_ft)^0.85034
  
  #Calculates CV for Pacific Madrone and Willows using equation 40
  subset.eq40$CVTS_ft <- 0.0067322665*(subset.eq40$DBH_in)^1.96628*(subset.eq40$Ht_ft)^0.83458
  
  #Calculates CV for Canyon Live Oak using equation 42
  subset.eq42$CVTS_ft <- 0.0097438611*(subset.eq42$DBH_in)^2.20527*(subset.eq42$Ht_ft)^0.61190
  
  #Calculates CV for California Live Oak using equation 42
  subset.eq43$CVTS_ft <- 0.0065261029*(subset.eq43$DBH_in)^2.31958*(subset.eq43$Ht_ft)^0.62528
  
  
  #Re-combine into single data frame
  #set k equal to the number of columns in data
  k<-dim(data)[2]
  data.small <- rbind(small3, small5, small16, small18, small19, small20, small23)
  if (dim(data.small)[1]>0){
    data.small <- cbind(data.small[,1:k],data.small[,"CVTS_ft"])
    colnames(data.small) [k+1] <- "CVTS_ft"
  } # end if data.small
  
  data.large <- rbind(large3, large5, large16, large18, large19, large20, large23)
  if (dim(data.large)[1]>0){
    data.large <- cbind(data.large[,1:k],data.large[,"CVTS_ft"])
    colnames(data.large) [k+1] <- "CVTS_ft"
  } # end if data.large
  
  data.1step<-rbind(subset.eq6,subset.eq8,subset.eq17,subset.eq21,subset.eq24,subset.eq26,subset.eq28,
                        subset.eq32,subset.eq33,subset.eq34,subset.eq37, subset.eq38,subset.eq40, subset.eq42, subset.eq43)
  if (dim(data.1step)[1]>0){
    data.1step<-cbind(data.1step[,1:k], data.1step[,"CVTS_ft"])
    colnames(data.1step) [k+1] <- "CVTS_ft"
  }#end if data.1step
  
  data <- rbind(data.small, data.large, data.1step)

  return (data)
} 


# Calculates Stem Biomass given parsed DATA and cubic volume 
# calculations from CalcCVTS. Returns DATA with two new columns:
# estimated Stem Biomass in both Kg ("Stembiom_kg") and Tons 
# ("Stembiom_tons"). 

CalcStemBiomass<-function (data){
  data$Stembiom_tons <- (data$CVTS_ft * data$Density) / 2000
  data$Stembiom_kg <- data$Stembiom_tons * 907.18474
  return(data)
}


# Calculates Bark Biomass given parsed DATA. Returns DATA 
# with a new column "barkbiom_kg" containing estimated 
# bark biomass in kg. 

CalcBarkBiomass<-function (data){
  
  BarkSet<-unique(data$Eq_bark)
  
   
  for (i in  c(0,1,2,4,5,8,9,10,11,12,13,14,15, 16,17,18,20,21)) {
    assign(paste("bark", i, sep = ""), subset(data, data$Eq_bark == i))
  }
  #If statement includes only branch equations present

  
  #DBH in cm and HT in m
  
  # No separate branch and biomass calcs for hardwoods. Components included in total. (eq. 0)
  if( 0 %in% BarkSet)bark0$barkbiom_kg<-0 
  
  #Bark biomass for White Firs (eq. 1)
  if( 1 %in% BarkSet)bark1$barkbiom_kg <- exp(2.1069 + 2.7271 * log(bark1$DBH_cm)) / 1000
  
  #Bark biomass for Grand Firs (eq. 2)
  if( 2 %in% BarkSet)bark2$barkbiom_kg <- 0.6 + 16.4 * (bark2$DBH_cm/100)^2 * bark2$Ht_m
  
  #Bark biomass for California red fir (eq.4)
  if( 4 %in% BarkSet)bark4$barkbiom_kg <- exp(1.47146 + 2.8421 * log(bark4$DBH_cm)) / 1000
  
  #Bark biomass for Noble firs (eq. 5)
  if( 5 %in% BarkSet)bark5$barkbiom_kg <- exp(2.79189 + 2.4313 * log(bark5$DBH_cm)) / 1000
  
  #Bark biomas for Douglas firs (eq.8)
  if( 8 %in% BarkSet)bark8$barkbiom_kg <- exp(-4.3103 + 2.4300 * log(bark8$DBH_cm))
  
  #Bark biomass for Jefferey pines, Ponderosa pines, and Foothill pines (eq. 9)
  if( 9 %in% BarkSet)bark9$barkbiom_kg <- exp(-3.6263 + 1.34077 * log(bark9$DBH_cm) + 0.8567 * log(bark9$Ht_m))
  
  #Bark biomass for Sugar pines (eq. 10)
  if( 10 %in% BarkSet)bark10$barkbiom_kg <- exp(2.183174 + 2.6610 * log(bark10$DBH_cm)) / 1000
  
  #Bark biomass for Western White Pines (eq. 11)
  if( 11 %in% BarkSet)bark11$barkbiom_kg <- 1.2 + 11.2 * (bark11$DBH_cm/100)^2 * bark11$Ht_m 
  
  #Eq 12
  if( 12 %in% BarkSet)bark12$barkbiom_kg <- exp(-13.3146 + 2.8594 * log(bark12$DBH_cm))*1000
  
  #Eq 13
  if( 13 %in% BarkSet)bark13$barkbiom_kg <- 0.336 + 0.00058 * bark13$DBH_cm^2 * bark13$Ht_m
  
  #Eq 14
  if( 14 %in% BarkSet)bark14$barkbiom_kg <- 3.2 + 9.1 * (bark14$DBH_cm/100)^2 * bark14$Ht_m
  
  #Eq 15
  if( 15 %in% BarkSet)bark15$barkbiom_kg <- exp(-4.371 + 2.259*log(bark15$DBH_cm))
  
  #Eq 16
  if( 16 %in% BarkSet)bark16$barkbiom_kg <- exp(-10.175 + 2.6333 * log(bark16$DBH_cm * pi))
  
  #Eq 17
  if( 17 %in% BarkSet)bark17$barkbiom_kg <- exp(7.189689 + 1.5837 * log(bark17$DBH_cm)) / 1000
  
  #Eq 18
  if( 18 %in% BarkSet)bark18$barkbiom_kg <- 1.3 + 27.6 * (bark18$DBH_cm/100)^2 * bark18$Ht_m
  
  #Eq 20
  if( 20 %in% BarkSet)bark20$barkbiom_kg <- exp(-4.6424 + 2.4617 * log(bark20$DBH_cm))
  
  #Eq 21
  if( 21 %in% BarkSet)bark21$barkbiom_kg <- 0.9 + 27.4 * (bark21$DBH_cm/100)^2 * bark21$Ht_m
  
  data <- rbind(bark0,bark1,bark2,bark4,bark5,bark8,bark9,bark10,bark11,bark12,bark13,bark14,bark15,bark16,bark17,bark18,bark20,bark21)
  
  return (data)
} 


# Calculates Branch Biomass given parsed DATA. Returns DATA 
# with a new column "branchbiom_kg" containing estimated
# bark biomass in kg. 

CalcBranchBiomass<-function (data){
  BranchSet<-unique(data$Eq_branch)
  for (i in  c(0,1,3,6,7,8,9,10,11,12,13,14,16,17)) {
    assign(paste("branch", i, sep = ""), subset(data, data$Eq_branch == i))
  }
  #Eq 0
  if( 0 %in% BranchSet)branch0$branchbiom_kg<-0
  #Eq 1
  if( 1 %in% BranchSet)branch1$branchbiom_kg <- 13.0 +12.4 * (branch1$DBH_cm/100)^2 * branch1$Ht_m
  #Eq3
  if( 3 %in% BranchSet)branch3$branchbiom_kg <- exp(-4.1817 + 2.3324 * log(branch3$DBH_cm))
  #Eq 6
  if( 6 %in% BranchSet)branch6$branchbiom_kg <- exp(-3.6941 + 2.1382 * log(branch6$DBH_cm))
  #Eq 7
  if( 7 %in% BranchSet)branch7$branchbiom_kg <- exp(-4.1068 + 1.5177 * log(branch7$DBH_cm) + 1.0424 * log(branch7$Ht_m))
  #Eq 8
  if( 8 %in% BranchSet)branch8$branchbiom_kg <- exp(-7.637 + 3.3648 * log(branch8$DBH_cm))
  #Eq 9
  if( 9 %in% BranchSet)branch9$branchbiom_kg <- 9.5 + 16.8 * (branch9$DBH_cm/100)^2 * branch9$Ht_m
  #Eq 10
  if( 10 %in% BranchSet)branch10$branchbiom_kg <- 0.199 + 0.00381 * branch10$DBH_cm^2 * branch10$Ht_m 
  #Eq 11
  if( 11 %in% BranchSet)branch11$branchbiom_kg <- 7.8 + 12.3 * (branch11$DBH_cm/100)^2 * branch11$Ht_m
  #Eq 12
  if( 12 %in% BranchSet)branch12$branchbiom_kg <- exp(-4.570+2.271*log(branch12$DBH_cm))
  #Eq 13
  if( 13 %in% BranchSet)branch13$branchbiom_kg <- exp(-7.2775 + 2.3337 * log(branch13$DBH_cm * pi))
  #Eq 14
  if( 14 %in% BranchSet)branch14$branchbiom_kg <- 1.7 + 26.2 * (branch14$DBH_cm/100)^2 * branch14$Ht_m
  #Eq 16
  BF <- (exp(-4.5648 + 2.6232 * log (branch16$DBH_cm)))* (1/(2.7638 + 0.062 * branch16$DBH_cm^1.3364))
  if( 16 %in% BranchSet)branch16$branchbiom_kg <-exp(-4.5648 + 2.6232 * log (branch16$DBH_cm))-BF
  #Eq 17
  if( 17 %in% BranchSet)branch17$branchbiom_kg <- exp(-5.2581 + 2.6045*log(branch17$DBH_cm))
  data <- rbind(branch0,branch1,branch3,branch6,branch7,branch8,branch9,branch10,branch11,branch12,branch13,branch14,branch16,branch17)
  return(data)
} #end of CalcBranchBiomass
 

# Calculates total biomass given parsed DATA which includes 
# stem biomass, bark biomass, and branch biomass. Returns 
# DATA with new columns containing summed biomass in Kg 
# ("Sumbiom_kg") and lbs ("Sumbiom_lbs").

CalcSumBiomass <- function (data){
  data$Sumbiom_kg <- data$Stembiom_kg + data$barkbiom_kg + data$branchbiom_kg
  data$Sumbiom_lbs <- data$Sumbiom_kg*2.2046226

  return(data)
}


# Given DATA, CODEFILE, and DIR, runs all other helper functions 
# and returns DATA file containing all Biomass and Cubic 
# Volume calculations, as well as another column containing 
# estimated Above Ground Living Biomass "AGL.MG.ha" in Mg/ha
# if original data table coontains a column labeled "TPH" 
# containing TPH measurements 

CalcPieces <-function (data, codefile, dir = ""){
  mergefile <- ParseData(data, codefile, dir)
  see1<-CalcCVTS(mergefile)
  see2<-CalcStemBiomass(see1)
  see3<-CalcBarkBiomass(see2)
  see4<-CalcBranchBiomass(see3)
  see5<-CalcSumBiomass(see4)
  ###Convert to AGL per tree in Mg/ha
  if ("TPH" %in% colnames(see5)) {
    see5$AGL.Mg.ha<-see5$Sumbiom_kg/1e3*see5$TPH
  }
  see5
} 



#################################
############# MAIN ##############
#################################

# Input: DATA file, CODEFILE (provided), DIR (optional)

# CODEFILE contains information needed to calculate tree Cubic Volume 
# and Biomass using the BECI equations (equation codes for stem,
# bark, and branches, and wood density value). By default, uses 
# "BECI EQN LUT.csv" which is included in the repository.

# DIR is your working directory path. Leave blank if you are running
# this program from your working directory. 

# OUTPUT is what you would like to name your output files. Tree
# biomass file will be named output_tree, and plot biomass file 
# will be named output_plot. Default is 'XXXX_biomass"

# DATA file must be a CSV file containing at least the following: 
# 1. A column named "Species" or "SPP4" containing the species code 
# for each tree
# 2. A column named "DBH_in" and/or "DBH_cm" containing the trees' 
# DBH measurement in inches and/or cm respectively
# 3. A column named "Ht_m" and/or "Ht_ft" containing the trees' 
# height in meters or feet respectively.

# RETURNS a new table containing original tree information along 
# with equation information and tree biomass data. Writes this
# table into DIR ("XXXX_tree_biomass.csv"). 

# If provided DATA includes a "Plot" column and a "Live" column 
# (binary indicator for tree is living), then function also writes
# a file containing plot biomass data and summary statistics into 
# DIR ("XXX_tree_biomass.csv"). 

# If provided DATA includes a "TPH" column, then AGL in mg/ha calcualtions
# will be included in returned table. 


Biomass<-function(data, codefile = "BECI EQN LUT.csv", dir = "", output = "XXXX_biomass") { 
  require (tidyverse)
  require(readxl)
  # Runs all previously defined biomass and cubic volumn equatins 
  AGL.tree<-CalcPieces(data, codefile, dir)
  # Creates a file containing plot summary statistics if possible
  if(("Plot" %in% colnames(AGL.tree)) && "Live" %in% colnames(AGL.tree)) {
    AGL.plot<-AGL.tree%>%
      group_by (Plot,Live)%>%
      summarise (AGL_Mg_ha = sum(AGL.Mg.ha),
                TreeHt.x = mean (Ht_m),
                TreeHt.sd = sd (Ht_m),
                TreeHt.max = max(Ht_m))
    mean(AGL.plot$AGL_Mg_ha)
    min(AGL.plot$AGL_Mg_ha)
    max(AGL.plot$AGL_Mg_ha)
    # Writes file containing plot summary statistics into directory
    write_csv(AGL.plot, paste(output, "_plot", ".csv", sep = ""))
  }
  # Writes file containing tree biomass information into directory 
  write_csv(AGL.tree, paste(output, "_tree", ".csv", sep = ""))
  return(AGL.tree)
}


data = Biomass("Complicated Test Data.csv", output = "complicated")

#################################
######## DIAGNOSTICS ############
#################################

# Plots Above Ground Living Biomass vs. Total Biomass given DATA,
# which is the table returned by BIOMASS. Needs original table to have
# TPH  column

BIOM_vs_AGL <- function(DATA) {
  return(qplot(Sumbiom_kg,AGL.Mg.ha, data=AGL.tree))
}
AGL.tree = Biomass("SEKI_biomass_data.csv", "BECI EQN LUT.csv")


# Plots Above Ground Living Biomass vs. Total Biomass for 
# each tree species, given DATA, which is the table returned by
# the BIOMASS function. Needs original table to have TPH column

BIOM_AGL_by_Species <- function(DATA) {
  require(lattice)
  xyplot (Sumbiom_kg~AGL.Mg.ha|SPP6, data=AGL.tree)
}





