#plot for sequential z-test for testing mu<=0 against mu>0
library(ggplot2)

#parameters
samplesize = 500 #number of samples
true_mean = 0.3 #true mean of the samples
alpha = 0.05 #significance level

set.seed(9) #for reproducibility

###function to compute the sequential z-test
seq_test <- function(n){
  numerator = sum(samples[1:n]) - sqrt(samplesize) * quantile
  difference = sqrt(samplesize - n)
  return(pnorm(numerator/difference))
}


time = seq(1, samplesize,1)
###generate samples from a normal distribution with the true mean
samples = rnorm(samplesize, mean = true_mean, sd = 1)

###compute the quantile for the sequential z-test
quantile =qnorm(1 - alpha, mean = 0, sd = 1)

###compute the sequential z-test for all samples
Y = numeric(samplesize) #vector to store the sequential test values
for(i in 1:(samplesize)){
  Y[i]=seq_test(i)
}
#use ggplot to plot the sequential test
df = data.frame(time = time, Y = Y)
ggplot(df, aes(x = time, y = Y)) +
  geom_line(color = 'black') +
  labs(title = 'Sequential z-Test', x = 'Time', y = 'Sequential Test Value') +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size = 12))

#find the first time point where the sequential test value exceeds 0.9999
which(Y >= 0.9999)[1]
