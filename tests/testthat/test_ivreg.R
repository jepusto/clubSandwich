context("ivreg objects")
library(AER)
data("CigarettesSW", package = "AER")
CigarettesSW <- within(CigarettesSW, {
  rprice <- price/cpi
  rincome <- income/population/cpi
  tdiff <- (taxs - tax)/cpi
})

obj_un <- ivreg(log(packs) ~ log(rprice) + log(rincome) | log(rincome) + tdiff + I(tax/cpi),
             data = CigarettesSW)
obj <- ivreg(log(packs) ~ log(rprice) + log(rincome) | log(rincome) + tdiff + I(tax/cpi),
            data = CigarettesSW, weights = population)
summary(obj)

X <- model.matrix(obj, component = "regressors")
Z <- model.matrix(obj, component = "instruments")
y <- log(CigarettesSW$packs)
w <- weights(obj)

# check unweighted iv
XZ <- model.matrix(obj_un, component = "projected")
ZtZ_inv <- chol2inv(chol(t(Z) %*% Z))
XZ_check <- Z %*% ZtZ_inv %*% t(Z) %*% X
all.equal(XZ, XZ_check, check.attributes=FALSE)
all.equal(coef(obj_un), lm.fit(XZ, y)$coefficients)
all.equal(bread(obj_un), chol2inv(chol(t(XZ) %*% XZ)) * nobs(obj_un), check.attributes=FALSE)
hii <- diag(X%*% chol2inv(chol(t(XZ) %*% XZ)) %*% t(XZ))
all.equal(hatvalues(obj_un), hii)

# check weighted iv
XZ <- model.matrix(obj, component = "projected")
ZwZ_inv <- chol2inv(chol(t(Z) %*% (w * Z)))
XZ_check <- Z %*% ZwZ_inv %*% t(Z) %*% (w * X)
all.equal(XZ, XZ_check, check.attributes=FALSE)
all.equal(coef(obj), lm.wfit(XZ, y, w)$coefficients)
all.equal(bread(obj), chol2inv(chol(t(XZ) %*% (w * XZ))) * nobs(obj), check.attributes=FALSE)
hii <- diag(X%*% chol2inv(chol(t(XZ) %*% (w * XZ))) %*% t(w * XZ))
all.equal(hatvalues(obj), hii) # does not agree because hatvalues doesn't work with weighting
