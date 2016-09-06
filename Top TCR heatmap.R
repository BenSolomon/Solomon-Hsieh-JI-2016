library(plyr); library(dplyr); library(vegan); library(reshape2); 
library(ggplot2); library(RColorBrewer)

average.reads.f <- function(data, tis) {
  sub <- subset(data, tissue == tis)
  lev <- levels(sub$population)
  tally <- data.frame()
  for (i in 1:length(lev)){
    sub.sub <- subset(sub, population == lev[i])
    df.c <- dcast(sub.sub, trvcdr3 ~ mouseNum, fun.aggregate = sum, value.var = "numReads")
    df.c <- data.frame(df.c[,1], decostand(df.c[-1], "total", 2))
    df.c <- data.frame(df.c[,1], apply(df.c[-1], 1, mean), rep(lev[i], length(df.c[,1])))
    tally <- rbind(tally, df.c)
  }
  names(tally) <- c("trvcdr3", "frequency", "population")
  return(tally)
}

top.tcr.f <- function(data, number = 25) {
  data$population <- factor(data$population)
  groups <- levels(data$population)
  tcrs <- vector()
  for (i in 1:length(groups)){
    x <- subset(data, population == groups[i])
    x <- x[order(x$frequency, decreasing = T), ]
    tcrs <- c(tcrs, as.character(x$trvcdr3[1:number]))
  }
  return(tcrs)
}


df <- read.csv("TCRsample.csv", header = T)

df.1 <- average.reads.f(df, "mLN") #Generates average population frequencies across mice for given tissue

top <- top.tcr.f(df.1, 25) #Finds top 25 TCRs for each population

df.s <- df.1[df.1[,1] %in% top, ] #Limits TCR data set to only the top 25 found above
names(df.s) <- c("tcr", "frequency", "population")
df.s$tcr <- factor(df.s$tcr)
df.s$tcr <- factor(df.s$tcr, levels = levels(df.s$tcr)[match(top, levels(df.s$tcr))]) #Reorders TCR levels to be in order of "top 25" group

blank <- expand.grid(tcr = levels(df.s$tcr), population = levels(df.s$population))

df.s <- merge(df.s, blank, all = T)
df.s$frequency[is.na(df.s$frequency)] <- 0

s <- summary(df.s$frequency[df.s$frequency > 0])
s <- quantile(df.s$frequency[df.s$frequency > 0], probs = seq(0,0.9,0.125))

df.s$groups <- cut(df.s$frequency, breaks = c(0, s, 1), right = F) 
#Breaks TCR frequencies into frequency groups since distribution is skewed

my.cols <- c("white", brewer.pal(8, "Reds"))

ggplot(data = df.s, aes(x = tcr, y = population, fill = groups)) + 
    theme_classic()+
    geom_tile(color = "grey50") +
    scale_fill_manual(values = my.cols) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0)) +
    theme(panel.grid=element_blank(), panel.border=element_blank())
