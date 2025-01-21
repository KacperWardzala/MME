y1 = as.numeric(y)
y2 = y_ICF



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

y_m = as.matrix(c(y1, y2))


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


est_IFL = results$est[321:499]
est_ICF = results$est[820:998]

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

(r2_IFL = (6277.935-trait1) / 6277.935)
(r2_ICF = (6017.857-trait2) / 6017.857)
r2_IFL[r2_IFL < 0] = 0
r2_ICF[r2_ICF < 0] = 0

rICF = sqrt(r2_IFL)
rIFL = sqrt(r2_ICF)

cbind(rIFL, rIFL)

rICF_wielo = rIFL
rIFL_wielo = rICF

plot(rIFL_wielo, ylim = c(0, 1))
plot(rICF_wielo, ylim = c(0, 1))


plot(IFL_std, est_IFL_std, main = 'Obserwowane IFLvs estymowana wartosc hodowlana')
plot(ICF_std, est_ICF_std, main = 'Obserwowane ICF vs estymowana wartosc hodowlana')



