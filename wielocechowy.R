y1 = as.numeric(y)
y2 = y_ICF

G = rbind(c(235.55, -238.38), c(-283.38, 412.24))
G

R = rbind(c(3632.90, -1107.40),c(-1107.40, 8134.30))
R

Z1 = matrix(0, nrow(IFL), nrow(A))
Z2 = Z1

Z = matrix(0, 2*nrow(IFL), 2*nrow(A))
Z[1:179, 1:499] = Z1
Z[180:358, 500:998] = Z2
Z

breed = model.matrix(~ as.factor(orginalId$V5[321:499]) - 1)
birthYear = model.matrix(~ as.factor(orginalId$V4[321:499]) - 1)

X1 = cbind(breed, birthYear)
X2 = X1

X = matrix(0, 2*179, 2*19)
X[1:179, 1:19] = X1
X[180:358, 20:38] = X2

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

y_m = as.matrix(c(y1, y2))


results = mme3(y_m, X, Z, A, G, R)
results

dim(results$C)

dim(results$est)

est_1 = results$est[5:12]
est_2 = results$est[13:20]

est_1
est_2

#standaryzacja

est_1_new = 10*(est_1-mean(est_1))/sd(est_1) + 100
est_2_new = 10*(est_2-mean(est_2))/sd(est_2) + 100

mean(est_1_new)
sd(est_1_new)

cbind(est_1_new, est_2_new)
cor(est_1_new, est_2_new)

#dokładnośc oceny wielocechowej

C = as.matrix(mme3(y_m, X, Z, A, G, R)$C)
invC = ginv(C)

invC22 = invC[39:1036, 39:1036]

trait1 = diag(invC22)[1:499]
trait2 = diag(invC22)[500:998]

(r2_1 = (235.55-trait1) / 235.55)
(r2_2 = (412.24-trait2) / 412.24)
r2_1[r2_1 < 0] = 0
r2_2[r2_2 < 0] = 0

(r1 = sqrt(r2_1))
(r2 = sqrt(r2_2))

cbind(r1, r2)

r1_wielo = r1
r2_wielo = r2
