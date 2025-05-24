#plot for sequential z-test
library(ggplot2)

samplesize = 100
true_mean = 0.3

time = seq(1, samplesize,1)
samples = rnorm(samplesize, mean = true_mean, sd = 1)
alpha = 0.05
quantile =qnorm(1 - alpha, mean = 0, sd = 1)

seq_test <- function(n){
  numerator = sum(samples[1:n]) - sqrt(samplesize) * quantile
  difference = sqrt(samplesize - n)
  return(pnorm(numerator/difference))
}
Y = numeric(samplesize - 1)
for(i in 1:(samplesize)){
  Y[i]=seq_test(i)
}
#use ggplot to plot the sequential test
df = data.frame(time = time, Y = Y)
ggplot(df, aes(x = time, y = Y)) +
  geom_line(color = 'black') +
  labs(title = 'Sequential Z-Test', x = 'Time', y = 'Sequential Test Value') +
  theme_minimal()


