library(MASS)

y = IFL$phenotype
y = as.matrix(as.numeric(y))

G = Gmatrix(
  SNPmatrix = t(genetyka),
  method = "VanRaden",
  missingValue = -9,
  maf = 0.01)

Z = I

mme4 = function(y, X, Z1, Z2, A, G, sigma_a, sigma_e, n) {
  sigma_g = sigma_a / n
  alpha1 = sigma_e / sigma_a
  alpha2 = sigma_e / sigma_g
  invA = ginv(A)
  invG = ginv(G)
  C = rbind(cbind(t(X)%*%X, t(X)%*%Z1, t(X)%*%Z2),
            cbind(t(Z1)%*%X, t(Z1)%*%Z1+invA*c(alpha1), t(Z1)%*%Z2),
            cbind(t(Z2)%*%X, t(Z2)%*%Z1, t(Z2)%*%Z2 + invG*c(alpha2)))
  print(dim(t(X)%*%X))
  print(dim(cbind(t(Z1)%*%X)))
  print(dim(t(Z2)%*%X))
  rhs = rbind(t(X)%*%y, t(Z1)%*%y, t(Z2)%*%y)
  invC = ginv(C)
  estimators = invC%*%rhs
  list(C = C, est = estimators)
}

(est = mme4(y, X, Z, Z, A[321:499, 321:499], G, 0.02086693, 0.007582836, 179)$est)


C = as.matrix(mme4(y, X, Z, Z, A[321:499, 321:499], G, 0.02086693, 0.007582836, 179)$C)
invC = ginv(C)

invC22 = invC[19:198, 19:198]
(r2 = diag(1 - invC22*(0.007582836/0.02086693)))
r2[r2 < 0] = 0
(r = sqrt(r2))