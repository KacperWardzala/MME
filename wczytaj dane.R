pedigree <- read.csv2('http://theta.edu.pl/wp-content/uploads/2023/01/pedigree.csv', header = FALSE)
head(pedigree)
laktacje <- read.csv2('http://theta.edu.pl/wp-content/uploads/2023/01/laktacje.csv',stringsAsFactors = FALSE)
head(laktacje)

library(dplyr)

laktacje

laktacje$calvingDate <- as.Date(laktacje$calvingDate, format="%d.%m.%Y")
laktacje$inseminationDate <- as.Date(laktacje$inseminationDate, format="%d.%m.%Y")

ifl <- laktacje %>%
  group_by(cowId) %>%
  summarise(IFL = as.numeric(max(inseminationDate) - min(inseminationDate)))

head(ifl)

icf <- laktacje %>%
  filter(success == 1) %>%
  group_by(cowId) %>%
  summarise(ICF = as.numeric(calvingDate - min(inseminationDate)))
head(icf)

library(data.table)
library(AGHmatrix)

# Wczytanie danych
pedigreeOrg <- fread('http://theta.edu.pl/wp-content/uploads/2023/01/pedigree.csv', header = FALSE)

uSireId = unique(pedigreeOrg$V2)
uDamId = unique(pedigreeOrg$V3)
uCowId = unique(pedigreeOrg$V1)

orginalId = c(uSireId, uDamId, uCowId)
orginalId = as.data.table(orginalId)
colnames(orginalId) = "orginalId"

orginalId = merge(orginalId, pedigreeOrg, by.x = "orginalId",
                  by.y = "V1", all.x = TRUE)

orginalId$V2[is.na(orginalId$V2)] = 0
orginalId$V3[is.na(orginalId$V3)] = 0
orginalId$V4[is.na(orginalId$V4)] = 0
orginalId$V5[is.na(orginalId$V5)] = 0

setorder(orginalId, V4)

orginalId$newId = 1:nrow(orginalId)
newId = orginalId[, c(1, 6)]
colnames(newId) = c("V2", "newSireId")

orginalId = merge(orginalId, newId, by = "V2", all.x = TRUE)
orginalId$newSireId[is.na(orginalId$newSireId)] = 0

colnames(newId) = c("V3", "newDamId")

orginalId = merge(orginalId, newId, by = "V3", all.x = TRUE)
orginalId$newDamId[is.na(orginalId$newDamId)] = 0

orginalId = orginalId[, c(3, 2, 1, 4, 5, 6, 7, 8)]
setorder(orginalId, newId)

A = Amatrix(orginalId[, c(6, 7, 8)])

#laktaccjeee

data = fread('http://theta.edu.pl/wp-content/uploads/2023/01/laktacje.csv')


data$inseminationDate <- as.Date(data$inseminationDate, format = "%d.%m.%Y")

firstInsemination = data[, .(firstIns = min(inseminationDate)),
                         by = cowId]
lastInsemination = data[, .(lastIns = max(inseminationDate)),
                        by = cowId]

IFL = merge(firstInsemination, lastInsemination, by = "cowId",
            all.x = TRUE)
IFL$phenotype = IFL$lastIns - IFL$firstIns

IFL = IFL[, -c(2, 3)]

IFL = merge(IFL, orginalId[, c(1, 6)], by.x = "cowId",
            by.y = "orginalId", all.x = TRUE)

setorder(IFL, newId)



ICF = icf %>%
  arrange(match(icf[[1]], IFL[[1]]))  # sortowanie icf według kolejności IFL w pierwszej kolumnie
ICF[,3] = IFL[,3]



#### fixed effects ####
breed = model.matrix(~ as.factor(orginalId$V5[321:499]) - 1)
birthYear = model.matrix(~ as.factor(orginalId$V4[321:499]) - 1)

X = cbind(breed, birthYear)

#### Z matrix ####
Z = matrix(0, nrow(IFL), nrow(A))
I = diag(nrow(IFL))
Z[1:179, 321:499] = I

Z
X

#G i R brane z treści zadania
G = rbind(c(235.55, -238.38), c(-283.38, 412.24))
G

R = rbind(c(3632.90, -1107.40),c(-1107.40, 8134.30))
R