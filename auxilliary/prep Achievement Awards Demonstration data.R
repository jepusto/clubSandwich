library(dplyr)
library(tidyr)
library(foreign)

filenames <- paste0("auxilliary/AngristLavy_AERdata/base",c("99","00","01","02"), ".dta")

possible_units <- c(0,18,20,22,24)

qntl <- function(x, groups = 4) {
  qtl <- quantile(x, (0:groups) / groups)
  if (length(unique(qtl)) == (groups + 1)) {
    print(paste(qtl, collapse = " "))
    cut(x, breaks = qtl, labels = 1:groups, include.lowest = TRUE, right = FALSE, ordered_result = TRUE)
  } else {
    NA
  }
}

AchievementAwardsRCT <- 
  lapply(filenames, read.dta) %>%
  bind_rows() %>%
  mutate(year = factor(year, levels = c(99,0,1,2), labels = 1999:2002),
         school_type = factor(ifelse(semrel, "Religious", ifelse(semarab, "Arab","Secular"))),
         student_id = paste0(year,"-",student_id),
         sex = ifelse(boy==1, "Boy","Girl"),
         attempted = possible_units[1 + att18 + att20 + att22 + att24],
         awarded = possible_units[1 + awr18 + awr20 + awr22 + awr24]) %>%
  select(school_id, school_type, pair, treated, year, student_id, pair, 
         sex, siblings = m_ahim, immigrant = ole5, father_ed = educav, mother_ed = educem,
         Bagrut_status = zakaibag, attempted, awarded , achv_math, achv_english = achv_eng, achv_hebrew = achv_hib,
         lagscore) %>%
  group_by(year, sex) %>%
  mutate(qrtl = qntl(lagscore, groups = 4),
         half = qntl(lagscore, groups = 2)) %>%
  ungroup()

with(AchievementAwardsRCT, table(year, sex))
with(AchievementAwardsRCT, table(year, qrtl, sex))
with(AchievementAwardsRCT, table(year, half, sex))

save(AchievementAwardsRCT, file = "data/AchievementAwardsRCT.RData", compress = "xz")
head(AchievementAwardsRCT)
