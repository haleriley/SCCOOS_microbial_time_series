# AOU_RF_V2 murdered darlings
# 2023-11-02
# RJH


# o2.days <- as.character(strptime(o2$Date, format = "%Y-%m-%d"))
# aou.daily.mean <- tapply(o2[,'aou'], INDEX = o2.days, FUN = mean)
# aou.corrected.daily.mean <- tapply(o2[,'aou.corrected'], INDEX = o2.days, FUN = mean)
# 
# names(aou.daily.mean) <- unique(o2.days)
# names(aou.corrected.daily.mean) <- unique(o2.days)
# 
# o2.daily.state <- rep(NA, length(aou.daily.mean))
# o2.daily.state[which(aou.daily.mean >= 0)] <- 'A'
# o2.daily.state[which(aou.daily.mean < 0)] <- 'H'
# names(o2.daily.state) <- unique(o2.days)


# predictors <- colnames(asv.train)[order(colSums(asv.train), decreasing = T)] # set predictors to asv global edge number

predictors <- colnames(asv.train)[order(colSums(asv.train), decreasing = T)[1:hyper.grid$n.edges[i]]]

predictors <- colnames(asv.train)[order(colSums(asv.train), decreasing = T)[1:selected.params$n.edges]]




