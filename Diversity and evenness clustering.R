library(reshape2);library(ggplot2);library(vegan);library(gplots);library(RColorBrewer)

##Functions
numbreaks.f <- function(x, char){
  b <- nchar(gsub(char, "", x))
  a <- nchar(x)
  return(a-b+1)
}

namesplit.f <- function(x, char) {
  x <- as.character(x)
  c <- numbreaks.f(x[1], char)
  r <- length(x)
  tally <- as.data.frame(matrix(nrow = r, ncol = c))
  for (i in 1:length(x)){
    for (j in 1:numbreaks.f(x[1], char)) {
      tally[i,j] <- strsplit(x[i], char)[[1]][j]
    }
  }
  return(tally)
}

df <- read.csv("TCRsample.csv", header = T)
df.c <- dcast(df, trvcdr3 ~ mouseNum + tissue + population, fun.aggregate = sum, value.var = "numReads")
df.c <- as.data.frame(t(df.c[,-1]))

len <- 50
scale <- seq(0,2, length.out = len)
df.div <- renyi(df.c, scales = scale, hill = F)
df.even <- df.div
for (i in 2:len){
  df.even <- cbind(df.even, df.even[,i]/df.even[,1])
}
df.even <- df.even[,(len+1):ncol(df.even)]
names(df.even) <- 1:ncol(df.even)

namesplit.f(rownames(df.even), "_")[,3]
color.key <- c("CD44hi" = brewer.pal(4, "Set1")[1], 
               "Foxp3+" = brewer.pal(4, "Set1")[2],
               "RORgt+" = brewer.pal(4, "Set1")[3],
               "Foxp3+RORgt+" = brewer.pal(4, "Set1")[4])
col.color <- color.key[namesplit.f(rownames(df.even), "_")[,3]]

hmcol <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(1000)

#Diversity clustering
heatmap.2(as.matrix(dist(df.div)),
          distfun = function(x) as.dist(x),
          col = hmcol,
          trace = "none",
          margins = c(15,15),
          ColSideColors = col.color)
legend("topright",
       legend = names(color.key),
       col = color.key,
       lty= 1, lwd = 5,
       border=FALSE, bty="n", y.intersp = 0.8, cex=0.8)

#Evennes clustering
heatmap.2(cor(t(df.even)), 
          trace = "none", 
          col = hmcol,
          margins = c(15,15),
          ColSideColors = col.color)
legend("topright",
       legend = names(color.key),
       col = color.key,
       lty= 1, lwd = 5,
       border=FALSE, bty="n", y.intersp = 0.8, cex=0.8)

