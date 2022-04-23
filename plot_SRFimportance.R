library(ggplot2)
library(reshape2)
library(RColorBrewer)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

var_imp <- read.table("YCU_PJI_preope_SRFimportance.txt")

var_imp_top <- var_imp[head(rev(order(rowSums(var_imp))),10),]
var_imp_melt <- melt(as.matrix(var_imp_top))
var_imp_melt_summary <- summarySE(var_imp_melt, measurevar="value", groupvars=c("Var1"))
var_imp_melt_summary$Var1 <- factor(var_imp_melt_summary$Var1,levels=rev(as.vector(rownames(var_imp_top))))

ggplot(var_imp_melt_summary,aes(x=Var1,y=value)) +
geom_bar(width=0.8,position = position_dodge(width = 0.8),stat="identity") +
geom_errorbar(aes(ymin=value-ci, ymax=value+ci),width=.2,position=position_dodge(.9))+
theme_classic(base_size = 16) + coord_flip()
ggsave("YCU_PJI_preope_SRFimportance.pdf",width=6,height=5)
