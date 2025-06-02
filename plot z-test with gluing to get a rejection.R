#plot the time we stop for z-tests with gluing to get a decisions
#i.e. If the seq. z-test reaches the threshold, we define a new test based on then next 10 data points to get a decision
#for testing mü <= 0 against mü > 0

library(ggplot2)
library(patchwork)
library(latex2exp)

set.seed(9) #for reproducibility

###normal sequential test defined from z-test at step n of level alpha (as in Chapter 3)
seq_test <- function(samples, n, sample_size, alpha){
  #compute the upper quantile of the z-test
  quantile = qnorm(1 - alpha, mean = 0, sd = 1)
  
  numerator = sum(samples[1:n]) - sqrt(sample_size) * quantile 
  difference = sqrt(sample_size - n)
  return(pnorm(numerator/difference))
}

###compute the step where we get a rejection for the normal seq. test, or set to sample_size + penalty if we never reject
stop_time_normal <- function(samples, sample_size, alpha){
  #compute the normal sequential test for each step
  Y = numeric(sample_size)
  for(i in 1:(sample_size)){
    Y[i] = seq_test(samples, i, sample_size, alpha)
  }
  #If we never reject, return sample_size + penalty to indicate a non-rejection
  if(length(which(Y >= threshold)) == 0) {
    return (sample_size + penalty)
  }
  #else find the index such that the sequential test is >= our boundary
  else{
    return(which(Y >= threshold)[1])
  }
}

###compute the step where the glued seq. test rejects, or set to sample_size + penalty if we never reject
stop_time_glued <- function(samples, sample_size, alpha, boundary){
  Y = numeric(sample_size)
  #compute the normal sequential test for each step
  for(i in 1:(sample_size)){
    Y[i] = seq_test(samples, i, sample_size, alpha)
  }
  
  #find the index such that the value of the sequential test is >= boundary
  if(length(which(Y >= boundary)) == 0) {
    return (sample_size + penalty)
  }
  else{
    stop_time <- which(Y >= boundary)[1]
  }
  #if there are less than size_new data points left, return the value of the normal seq. test
  if(stop_time + size_new > sample_size && Y[sample_size] == 1) {
    return (sample_size)
  }
  if(stop_time + size_new > sample_size && Y[sample_size] == 0) {
    return (sample_size + penalty)
  }
  
  ##if there are at least 10 data points left, we define a new test based on the next 10 data points
  alpha = Y[stop_time] #save the alpha of the sequential test at the step where we reach the boundary
  quantile = qnorm(1 - Y[stop_time], mean = 0, sd = 1) #quantile for z-test with new alpha
  
  #compute the new test based on the next 10 data points
  newtest = 1/sqrt(size_new) * sum(samples[(stop_time + 1):(stop_time + size_new)])
  #reject the null hypothesis if the new test is greater than the quantile
  if(newtest > quantile){
    return(stop_time + size_new)
  }
  #else we do not reject
  else{
    return(sample_size + penalty)
  }
}


#parameters
boundary = 0.9 #If the seq. z-test reaches the boundary, we define a new test
threshold = 0.9999 #If the seq. test is greater than this threshold, we reject the null hypothesis
penalty = 10 #penalty for not rejecting the null hypothesis, usefull to plot the results
alpha = 0.05 #starting alpha of z-test (and sequential test)
sample_size = 500 #size of one sample
rep = 1000 #number of repetitions for each value of mü
size_new = 15 #size of the new test,i.e. how many data points we consider to get a decision after we stoped

#plot the stop times
x = c(0, 0.1, 0.2, 0.3, 0.4) #values of mü (min, max, frequency)


#define matrix to save stop times for both the normal seq. tests and the glued seq. tests 
stoping_times_normal = matrix(0, nrow = length(x), ncol = rep)
stoping_times_glued = matrix(0, nrow = length(x), ncol = rep)

for (i in x){
  for(j in 1:rep){
    #get samples from the true distribution
    samples = rnorm(sample_size, mean = i, sd = 1)
    
    #compute stop times for the normal seq. test and the glued seq. test
    stoping_times_normal[which(x == i), j] = stop_time_normal(samples, sample_size, alpha)
    stoping_times_glued[which(x == i), j] = stop_time_glued(samples, sample_size, alpha, boundary)
    
  }
}

###make Histogram-plots for different values of mü
p <- vector("list", length(x)) #list to store the plots of size length(x)
average_normal = numeric(length(x)) #vector to store the average decision time for the normal seq. test
average_glued = numeric(length(x)) #vector to store the average decision time for the glued seq. test
nonrej_normal = numeric(length(x)) #vector to store the number of non-rejections for the normal seq. test
nonrej_glued = numeric(length(x)) #vector to store the number of non-rejections for the glued seq. test
for(i in 1:length(x)){
  #get the values for the current mü
  values_glued = stoping_times_glued[i,]
  values_normal = stoping_times_normal[i,]
  
  #compute the average decision time for the normal seq. test and the glued seq. test
  average_normal[i] = mean(values_normal[values_normal < sample_size])
  average_glued[i] = mean(values_glued[values_glued < sample_size])
  
  #save the amount of non rejectionss for the normal seq. test and the glued seq. test
  nonrej_normal[i] = length(which(values_normal == sample_size + penalty))
  nonrej_glued[i] = length(which(values_glued == sample_size + penalty))
  
  ###plot the histogram of the stop times for the normal seq. test and the glued seq. test
  # Combine both data sets into one with a 'Method' column
  df <- rbind(
    data.frame(x = values_normal, Method = "Normal seq. test"),
    data.frame(x = values_glued, Method = "Glued seq. test")
  )
  
  # Plot in a histogram
  p[[i]] <- ggplot(df, aes(x = x, fill = Method)) +
    geom_histogram(binwidth = 5, alpha = 0.5, position = "identity") +
    labs(title = TeX(paste("$\\mu =", x[i], "$")),subtitle = paste0("For rejection: Avg (Normal): ", round(average_normal[i]),
                                                                    ", Avg (Glued): ", round(average_glued[i])), x = "Decision time", y = "Count") +
    theme_minimal() +
    theme(
      plot.subtitle = element_text(
        size = 10,         # change to desired size
        face = "italic",   # optional: italic, bold, etc.
        color = "gray37" # optional: color
      )) + 
    theme(
      plot.title = element_text(margin = margin(b = 0)),       # reduce bottom margin of title
      plot.subtitle = element_text(margin = margin(t = 0))     # reduce top margin of subtitle
    ) +
    scale_fill_manual(values = c("Normal seq. test" = "red", "Glued seq. test" = "blue")) +
    guides(fill = guide_legend(title = NULL)) +
    coord_cartesian(xlim = c(-5, sample_size + penalty + 10))
}

#combine all plots such that it has only one legend
combined_plot <- wrap_plots(p, ncol = 1, guides = "collect") & theme(legend.position = "bottom")
plot(combined_plot)

#Average decision time if the tests do reject
average_normal
average_glued

#percentage of non-rejection/rep
nonrej_normal/rep * 100
nonrej_glued/rep * 100
