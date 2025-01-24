library(MASS)
# Potrzebne biblioteki
library(AGHmatrix)
library(dplyr)
library(data.table)
library(lme4)
library(sommer)
library(ggplot2)
library(MASS)

# modele mieszane z genetyką jednocechowe

mme5 = function(y, X, Z, A, sigma_a, sigma_e) {
  alpha = sigma_e / sigma_a
  invA = A #ginv(A)
  C = rbind(cbind(t(X)%*%X, t(X)%*%Z),
            cbind(t(Z)%*%X, t(Z)%*%Z+invA*c(alpha)))
  rhs = rbind(t(X)%*%y, t(Z)%*%y)
  invC = ginv(C)
  estimators = invC%*%rhs
  list(C = C, est = estimators)
}




# Obliczenie macierzy G i H

Gm <- Gmatrix(SNPmatrix = new_gen,
             method = "VanRaden",
             missingValue = -9,
             maf = 0.01)

invG <- ginv(Gm)

invA22 <- ginv(A)[321:499, 321:499]


invH <- rbind(
  cbind(matrix(0, nrow =320, ncol =320), matrix(0, nrow =320, ncol =179)),
  cbind(matrix(0, nrow =179, ncol =320), invG-invA22))

dim(invH)

# Bierzemy tutaj invH zamiast a po prostu

efects_IFL_MT = mme5(y = IFL, X = X, Z = Z, A = invH, g11, z11)

C = efects_IFL_MT$C
C22 = C[20:518, 20:518]
invC22_IFL_GEN = ginv(C22)



(r2 = diag(1 - invC22_IFL_GEN*(as.numeric(z11)/as.numeric(g11))))
r2[r2 < 0] = 0
r_IFL_GEN = sqrt(r2)

plot(r_IFL_GEN, ylim = c(0, 1), main = 'Dokładność oszacowania dla IFL z genotypami')


#dla ICF
efects_ICF_MT = mme5(y = ICF, X = X, Z = Z, A = invH, g22, z22)

C = efects_ICF_MT$C
C22 = C[20:518, 20:518]
invC22_ICF_GEN = ginv(C22)



(r2 = diag(1 - invC22_ICF_GEN*(as.numeric(z22)/as.numeric(g22))))
r2[r2 < 0] = 0
r_ICF_GEN = sqrt(r2)

plot(r_ICF_GEN, ylim = c(0, 1), main = 'Dokładność oszacowania dla ICF z genotypami')

cbind(r_IFL_GEN, r_ICF_GEN)

#estymacje
IFL_gen_est = efects_IFL_MT$est[340:518]
ICF_gen_est = efects_ICF_MT$est[340:518]


#standaryzacje
IFL_std = 10*(IFL- mean(IFL)) / sd(IFL) + 100
ICF_std = 10*(ICF- mean(ICF)) / sd(ICF) + 100

IFL_gen_est_std = 10*(IFL_gen_est - mean(IFL_gen_est)) / sd(IFL_gen_est) + 100
ICF_gen_est_std = 10*(ICF_gen_est - mean(ICF_gen_est)) / sd(ICF_gen_est) + 100

# Dopasowanie regresji liniowej dla IFL
model_IFL <- lm(IFL_std ~ IFL_gen_est_std)  # Regresja liniowa
cor_IFL <- cor(IFL_gen_est_std, IFL_std)    # Korelacja

# Wykres z linią trendu dla IFL
plot(IFL_gen_est_std, IFL_std, 
     main = "Obserwowane IFL vs estymowana wartość hodowlana z genotypami",
     #ylim = c(60, 160), 
     xlab = "Estymowana wartość IFL (standaryzowana)", 
     ylab = "Obserwowana wartość IFL (standaryzowana)")
abline(model_IFL, col = "red", lwd = 2)  # Dodanie linii trendu w kolorze czerwonym



# Dopasowanie regresji liniowej dla ICF
model_ICF <- lm(ICF_std ~ ICF_gen_est_std)  # Regresja liniowa
cor_ICF <- cor(ICF_gen_est_std, ICF_std)    # Korelacja

# Wykres z linią trendu dla IFL
plot(ICF_gen_est_std, ICF_std, 
     main = "Obserwowane IFL vs estymowana wartość hodowlana z genotypami",
     #ylim = c(60, 160), 
     xlab = "Estymowana wartość IFL (standaryzowana)", 
     ylab = "Obserwowana wartość IFL (standaryzowana)")
abline(model_ICF, col = "blue", lwd = 2)  # Dodanie linii trendu w kolorze czerwonym


