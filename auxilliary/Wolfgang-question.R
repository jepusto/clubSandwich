############################################################################

### install/load packages

# load the 'nlme' package
library(nlme)

# install/load the 'remotes' package
if (!require("remotes", quietly=TRUE)) {
  install.packages("remotes")
  library("remotes")
}

# install/load the 'esmpack' package
if (!require("esmpack", quietly=TRUE)) {
  remotes::install_github("wviechtb/esmpack")
  library(esmpack)
}


############################################################################

### data preparation

# download the dataset (if it isn't already downloaded)
if (!file.exists("auxilliary/viechtbauer2025_data_esmda_example.txt")) download.file("https://www.wvbauer.com/lib/exe/fetch.php/suppl:viechtbauer2025_data_esmda_example.txt", destfile="auxilliary/viechtbauer2025_data_esmda_example.txt")

# read in the data
dat <- read.table("auxilliary/viechtbauer2025_data_esmda_example.txt", header=TRUE, sep="\t", na.strings="", as.is=TRUE)

# examine the first six rows of the dataset
head(dat)

# order the dataset by subject id and obs within subject id (it already is, but just in case)
dat <- sort_by(dat, ~ id + obs)

# compute variables posaff and negaff (average of the 3 PA and 6 NA items)
dat$posaff <- combitems(c("mood_cheerf", "mood_relaxed", "mood_satisfi"), data=dat)
dat$negaff <- combitems(c("mood_irritat", "mood_anxious", "mood_down", "mood_guilty", "mood_insecur", "mood_lonely"), data=dat)

# create a lagged version of the posaff variable
dat$lposaff <- lagvar(posaff, id=id, day=day, data=dat)

# create a lagged version of the negaff variable
dat$lnegaff <- lagvar(negaff, id=id, day=day, data=dat)

# create a lagged version of the eventpl variable
dat$leventpl <- lagvar(eventpl, id=id, day=day, data=dat)

# keep only rows with complete data on the relevant variables
dat <- dat[complete.cases(dat[c("posaff", "lposaff", "negaff", "lnegaff", "eventpl", "leventpl")]),]

# restructure the relevant parts of the dataset into a 'very long' format
tmp <- dat[c("id", "obs", "posaff", "negaff", "eventpl", "lposaff", "lnegaff", "leventpl")]
tmp <- reshape(tmp, direction="long", idvar="row", varying=list(c("posaff", "negaff", "eventpl")))
tmp$row <- NULL
tmp <- tmp[order(tmp$id, tmp$obs, tmp$time),]
tmp <- rename("time", "var", tmp)
tmp$outcome <- c("posaff", "negaff", "eventpl")
tmp <- rename("posaff", "y", tmp)
rownames(tmp) <- NULL
head(tmp, 12)

############################################################################
### model fitting

# [...]

# fit a simpler model with correlated random intercepts for each outcome but
# no random slopes (note: convergence takes a minute or two)
print(system.time(
  res2 <- lme(
    y ~ outcome + outcome:lposaff + outcome:lnegaff + outcome:leventpl - 1,
    random = ~ outcome - 1 | id,
    correlation = corSymm(form = ~ var | id/obs), 
    weights = varIdent(form = ~ 1 | outcome),
    data=tmp, 
    control=list(maxIter=5000, msMaxIter=5000, msMaxEval=5000, msVerbose=TRUE)
  )
))

summary(res2)

# obtain the cluster-robust inference results
res3 <- coef_test(res2, vcov="CR2")

vcovCR(res2, type = "CR2")

tmp$obs_id <- as.character(tmp$obs)

res3 <- update(res2, correlation = corSymm(form = ~ var | id/obs_id))
coef_test(res3, vcov = "CR2")


mod_formula <- nlme::getGroupsFormula(res2$modelStruct$corStruct)
grps <- stats::model.frame(mod_formula, data = nlme::getData(res2))
grps1 <- apply(grps, 1, paste, collapse = "/")
head(grps1)

grps2 <- lapply(grps, as.character)
grps2 <- do.call(paste, args = c(grps2, sep = "/"))
head(grps2)


file.remove("auxilliary/viechtbauer2025_data_esmda_example.txt")

