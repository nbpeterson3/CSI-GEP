library(jaccard)
library(dplyr)
library(igraph)
library(vroom)
library(elbow)
library(ggplot2)

df = read.table("benchmarking/inmf/output/median_kl.csv", header = TRUE, sep = ",")
names(df) = c("median_KL", "K")

#Run a loess regression on JL data to smooth the curve
sl <- df %>% ggplot(aes(K, median_KL)) + geom_point() + geom_smooth(method = "loess", span = 0.3, method.args = list(degree = 1))
gbuild <- ggplot_build(sl)
elbow.df <- data.frame(X = gbuild$data[[2]]$x, Y = gbuild$data[[2]]$y)
#Find inflection point
ipoint <- elbow(data = elbow.df, plot = F)
#Since the data was run on smoothed data, retreive the closest rank to the x intercept calculated by elbow package
optimum.rank <- df$K[which.min(abs(df$K - ipoint$X_selected))]

print(optimum.rank)
