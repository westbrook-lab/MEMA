library(testthat)
library(MEMA)

##Testing building the M matrix for RUV normalization
#Test two unique rows
data <- data.table::data.table(Barcode = "B1",
                               Well = c("WL1","WL2"),
                               SignalType = "Signal",
                               CellLine = "ACells",
                               Drug = "none",
                               Ligand=c("L1","L2"))
testthat::compare(createRUVM(data),
                  matrix(data = c(1,0,0,1),
                         nrow = 2,
                         dimnames = list(list("B1_WL1_ACells_L1_none","B1_WL2_ACells_L2_none"),
                                         list("ACells_L1_none","ACells_L2_none"))))

#Test one unique row and two repeated rows
data <- data.table::data.table(Barcode = c("B1","B1","B2"),
                               Well = c("WL1","WL2","WL2"),
                               SignalType = "Signal",
                               CellLine = "ACells",
                               Drug = "none",
                               Ligand=c("L1","L2","L2"))
testthat::compare(createRUVM(data),
                  matrix(data = c(1,0,0,0,1,1),
                         nrow = 3,
                         dimnames = list(list("B1_WL1_ACells_L1_none","B1_WL2_ACells_L2_none","B2_WL2_ACells_L2_none"),
                                         list("ACells_L1_none","ACells_L2_none"))))

#Test one unique row and two repeated rows with pipe names
data <- data.table::data.table(Barcode = c("B1","B1","B2"),
                               Well = c("W","W","W"),
                               SignalType = "Signal",
                               CellLine = "ACells",
                               Drug = c("none","none","drugA"),
                               Ligand=c("L1|2","L2","L2"))
testthat::compare(createRUVM(data),
                  matrix(data = c(1,0,0,0,1,0,0,0,1),
                         nrow = 3,
                         dimnames = list(list("B1_W_ACells_L1|2_none","B1_W_ACells_L2_none","B2_W_ACells_L2_drugA"),
                                         list("ACells_L1|2_none","ACells_L2_none","ACells_L2_drugA"))))

