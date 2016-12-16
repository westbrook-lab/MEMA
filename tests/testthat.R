library(testthat)
library(MEMA)

##Testing building the M matrix for RUV normalization
#Test two unique rows
data <- data.table::data.table(BWL = c("B1WL1","B1WL2"),
                               SignalType = "Signal",
                               Ligand=c("L1","L2"))
testthat::compare(createRUVM(data),
                  matrix(data = c(1,0,0,1),
                         nrow = 2,
                         dimnames = list(list("B1WL1","B1WL2"),
                                         list("L1","L2"))))

#Test one unique row and two repeated rows
data <- data.table::data.table(BWL = c("B1WL1","B1WL2","B2WL2"),
                               SignalType = "Signal",
                               Ligand=c("L1","L2","L2"))
testthat::compare(createRUVM(data),
                  matrix(data = c(1,0,0,0,1,1),
                         nrow = 3,
                         dimnames = list(list("B1WL1","B1WL2","B2WL2"),
                                         list("L1","L2"))))

#Test one unique row and two repeated rows with pipe names
data <- data.table::data.table(BWL = c("B1WL1|2","B1WL2","B2WL2"),
                               SignalType = "Signal",
                               Ligand=c("L1|2","L2","L2"))
testthat::compare(createRUVM(data),
                  matrix(data = c(1,0,0,0,1,1),
                         nrow = 3,
                         dimnames = list(list("B1WL1|2","B1WL2","B2WL2"),
                                         list("L1|2","L2"))))



#Test two unique rows
data <- data.table::data.table(BWLD = c("B1WL1D1","B1WL2D1"),
                               SignalType = "Signal",
                               LD=c("L1","L2"))
testthat::compare(createRUVMGeneral(data, unitID = "BWLD", uniqueID = "LD"),
                  matrix(data = c(1,0,0,1),
                         nrow = 2,
                         dimnames = list(list("B1WL1D1","B1WL2D1"),
                                         list("L1","L2"))))

#Test replicates as they would be in a 96 well array with replicates
# and drug+ligand treatments and pipes in names
data <- data.table::data.table(BW = c("B1WA01","B1WA02","B1WB01","B1WB02","B1WC01"),
                               SignalType = "Signal",
                               LD=c("L1D1","L1D2","L1D2","L2|1D1","L1D1"))
testthat::compare(createRUVMGeneral(data, unitID = "BW", uniqueID = "LD"),
                  matrix(data = c(1,0,0,0,1,0,1,1,0,0,0,0,0,1,0),
                         nrow = 5,
                         dimnames = list(list("B1WA01","B1WA02","B1WB01","B1WB02","B1WC01"),
                                         list("L1D1","L1D2","L2|1D1"))))
