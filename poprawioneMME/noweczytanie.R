library(dplyr)
library(data.table)
library(AGHmatrix)

# Wczytanie danych z laktacji
laktacje = read.csv2("http://theta.edu.pl/wp-content/uploads/2023/01/laktacje.csv", sep = ';')

# Zmiana formatu na datę
laktacje$inseminationDate = as.Date(laktacje$inseminationDate, format = "%d.%m.%Y")
laktacje$calvingDate = as.Date(laktacje$calvingDate, format = "%d.%m.%Y")
head(laktacje)
dim(laktacje)

# Wczytanie danych z rodowodu
pedigree = read.csv2("http://theta.edu.pl/wp-content/uploads/2023/01/pedigree.csv", sep = ';',header = F)
colnames(pedigree) = c("cow", "sire", "dam", "YOB", "heard")
head(pedigree)
dim(pedigree)

## Wczytanie danych genetycznych

genetyka = read.csv2("http://theta.edu.pl/wp-content/uploads/2024/01/Genotypes.csv", sep = ';', header=F)
head(genetyka)

# Przekształcanie pedigree


# Stworzenie tabeli unikatowych ID (starych)
org_id <- data.table(id = unique(c(pedigree$sire, pedigree$dam, pedigree$cow)))
colnames(org_id) = "organised_id"
org_id = merge(org_id, pedigree, by.x = "organised_id", by.y = "cow", all.x = T)
org_id[is.na(org_id)] <- 0

# Wyznaczenie nowych ID (po prostu liczby)
setorder(org_id, YOB)
org_id$new_id = 1:nrow(org_id)

# Zamiana starych ID na nowe ID
org_id <- org_id %>%
  rowwise() %>%
  mutate(
    dam = if_else(
      dam != "0" | dam %in% organised_id,  # Using OR instead of AND
      as.character(ifelse(dam != "0", org_id$new_id[org_id$organised_id == dam], 'HERE')),
      as.character(dam)
    ),
    sire = if_else(
      sire != "0" | sire %in% organised_id,  # Using OR instead of AND
      as.character(ifelse(sire != "0", org_id$new_id[org_id$organised_id == sire], 'HERE')),
      as.character(sire)
    )
  ) %>%
  ungroup()

# Wybór danych do macierzy spokrewnień
A_data = as.data.table(org_id[,c(6, 2, 3)])
A_data[] <- lapply(A_data, function(x) as.integer(as.character(x)))

head(A_data)
A = Amatrix(A_data)
dim(A)
heatmap(A)


# Przekształcanie genetyki - zmiana nazw id na nowe i posortowanie

new_gen  <- merge(genetyka, org_id, by.x = "V1", by.y = "organised_id", all.x = T)

setorder(new_gen, new_id)
id_num <- new_gen[, 7006]

new_gen <- as.matrix(new_gen[,c(2:7001)])
rownames(new_gen) <- id_num




# Wyliczenie wartości fenotypowych
ins_times = c()
calv_times = c()
cows = c()

for(cow in unique(laktacje$cowId)){
  subset = laktacje[laktacje$cowId == cow,]
  first_ins = subset[1,]$inseminationDate
  last_ins = subset[nrow(subset),]$inseminationDate
  
  
  ins_time = last_ins - first_ins
  ins_times = c(ins_times, ins_time)
  calv_time = subset[nrow(subset),]$calvingDate - first_ins
  calv_times = c(calv_times, calv_time)
  cows = c(cows, cow)
  
}

# IFL ostatnia inseminacja  - pierwsza ins
# ICF calving - pierwsza inseminacja
phenotypes = data.frame(cows, ins_times, calv_times)
colnames(phenotypes) = c("cow", "IFL", "ICF")
merged_phenotypes = merge(phenotypes, org_id,by.x = "cow", by.y = "organised_id",)
merged_phenotypes_lactations = merge(merged_phenotypes, laktacje, by.x = 'cow', by.y = 'cowId')
head(merged_phenotypes_lactations)


dim(merged_phenotypes_lactations)

IFL = phenotypes$IFL
ICF = phenotypes$ICF



#Wyznaczenie maceirzy efektów stałych X i efektow losowych Z

heard = model.matrix(~ as.factor(org_id$heard[321:499]) -1)
birth = model.matrix(~ as.factor(org_id$YOB[321:499]) -1)

# Macierz efektów losowych
X = cbind(heard, birth)

# Macierz efektów stałych
Z = matrix(0, nrow = nrow(phenotypes), ncol = nrow(A))

# Macierz jednostkowa
I = diag(nrow(phenotypes))
Z[1:179, 321:499] = I


