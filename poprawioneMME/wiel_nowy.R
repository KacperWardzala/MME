y1 = IFL
y2 = ICF



# Tworzenie macierzy wystąpień dla efetków stałych
heard = model.matrix(~ as.factor(org_id$heard[321:499]) -1)
birth = model.matrix(~ as.factor(org_id$YOB[321:499]) -1)
X1= cbind(heard, birth)
X2 = X1
X = matrix(0, nrow = 2 * nrow(X1), ncol = 2 * ncol(X1)) #  TO MOŻNA ROBIĆ TYLKO JEŚLI MACIERZE WYSTĄPIEŃ SĄ IDENTYCZNE DLA OBU CZYNNIKÓW
X[1:nrow(X1), 1:ncol(X1)] = X1
X[(nrow(X1)+1):nrow(X), (ncol(X1)+1):ncol(X)] = X2



# W tym modelu Z jest rozszerzona, ale jako, ze jest dla tych sanych osobników to wystarczy ją rozszerzyć
Z = matrix(0, nrow = 2 * nrow(phenotypes), ncol =  2 * nrow(A)) # TO MOŻNA ROBIĆ TYLKO JEŚLI MACIERZE WYSTĄPIEŃ SĄ IDENTYCZNE DLA OBU CZYNNIKÓW
I = diag(nrow(phenotypes))
Z[1:179, 321:499] = I
Z[180:358, 820:998] = I

y_m = as.matrix(c(y1,y2))


library(MASS)

mme3 = function(y, X, Z, A, G, R) {
  invA = ginv(A)
  invG = ginv(G)
  
  
  R = kronecker(R, diag(179))
  invR = ginv(R)
  
  
  C = rbind(cbind(t(X)%*%invR%*%X, t(X)%*%invR%*%Z),
            cbind(t(Z)%*%invR%*%X, t(Z)%*%invR%*%Z+kronecker(invG, invA)))
  print(dim(t(Z)%*%invR%*%Z+kronecker(invG, invA)))
  rhs = rbind(t(X)%*%invR%*%y, t(Z)%*%invR%*%y)
  invC = ginv(C)
  estimators = invC%*%rhs
  list(C = C, est = estimators)
}



results = mme3(y_m, X, Z, A, G, R)
results


est_IFL = results$est[340:518]
est_ICF = results$est[858:1036]

#standaryzacja

est_IFL_std = 10*(est_IFL-mean(est_IFL))/sd(est_IFL) + 100
est_ICF_std = 10*(est_ICF-mean(est_ICF))/sd(est_ICF) + 100


cbind(est_IFL_std, IFL_std)
cbind(est_ICF_std, ICF_std)

#dokładnośc oceny wielocechowej

C = as.matrix(mme3(y_m, X, Z, A, G, R)$C)
invC = ginv(C)

invC22 = invC[39:nrow(invC), 39:nrow(invC)]

trait1 = diag(invC22)[1:499]
trait2 = diag(invC22)[500:998]

(r2_IFL = (g11-trait1) / g11)
(r2_ICF = (g22-trait2) / g22)
r2_IFL[r2_IFL < 0] = 0
r2_ICF[r2_ICF < 0] = 0

rICF = sqrt(r2_IFL)
rIFL = sqrt(r2_ICF)

cbind(rIFL, rIFL)

rICF_wielo = rIFL
rIFL_wielo = rICF

plot(rIFL_wielo, ylim = c(0, 1), main = 'wartości r dla cechy IFL westymowanej modelem wielocechowym')
plot(rICF_wielo, ylim = c(0, 1), main = 'wartości r dla cechy ICF westymowanej modelem wielocechowym')


# Dopasowanie regresji liniowej dla IFL
model_IFL <- lm(IFL_std ~ est_IFL_std)  # Regresja liniowa
cor_IFL <- cor(est_IFL_std, IFL_std)    # Korelacja

# Wykres z linią trendu dla IFL
plot(est_IFL_std, IFL_std, 
     main = "Obserwowane IFL vs estymowana wartość hodowlana",
     ylim = c(90, 160), 
     xlab = "Estymowana wartość IFL (standaryzowana)", 
     ylab = "Obserwowana wartość IFL (standaryzowana)")
abline(model_IFL, col = "red", lwd = 2)  # Dodanie linii trendu w kolorze czerwonym

# Dopasowanie regresji liniowej dla ICF
model_ICF <- lm(ICF_std ~ est_ICF_std)  # Regresja liniowa
cor_ICF <- cor(est_ICF_std, ICF_std)    # Korelacja

# Wykres z linią trendu dla ICF
plot(est_ICF_std, ICF_std, 
     main = "Obserwowane ICF vs estymowana wartość hodowlana",
     ylim = c(90, 160), 
     xlab = "Estymowana wartość ICF (standaryzowana)", 
     ylab = "Obserwowana wartość ICF (standaryzowana)")
abline(model_ICF, col = "blue", lwd = 2)  # Dodanie linii trendu w kolorze niebieskim

boxplot(list(r_IFL,r_ICF, rIFL_wielo, rICF_wielo), names = c('r IFL', 'r ICF', 'r IFL mt', 'r ICF mt')  ,main = 'Porównanie wartości r')

