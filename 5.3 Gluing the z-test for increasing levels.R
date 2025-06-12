#Section 5.3 gluing the z-test for increasing levels
#plot the time we stop for gluing the seq. z-test if we define a new test each time we reach 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 and then conclude with a tradional z-test based on 15 data points 
#i.e. we compute the seq. z-test for the data, if we reach a new desired level, we define a new one and so on
#for testing mü <= 0 against mü > 0

set.seed(9) #for reproducibility

###load libraries
library(ggplot2)
library(patchwork)
library(latex2exp)

###set parameters
rep = 2000 #how often we take samples
sample_size = 500 #size of one sample
alpha = 0.05 #starting level of z-test (and sequential test)
penalty = 10 #penalization for the stop time if we don't reject the null, used to make plots interpretable
threshold = 0.9999 #value, for which we reject
x = c(0, 0.1, 0.2, 0.3, 0.4) #values of mü
k = 1 #set to 2, to see how many of the non-rejections of the glued test, happen after 1/4 of the data is observed (does will appear at value 520 in the plots)

###define functions
##compute the sequential test from the z-test (as in Chapter 3) and compute the time we stop
#sequential test defined from the z-tests for data points start:end
#N = number of data points for traditional z-test, alpha = level of the test
seq_test <- function(samples, start, end, N, alpha){
  #compute the quantile for the level alpha sequential test
  quantile = qnorm(1 - alpha, mean = 0, sd = 1)
  
  numerator = sum(samples[start:end]) - sqrt(N)*quantile
  denominator = sqrt(N - (end-start + 1))
  return(pnorm(numerator/denominator))
}

#define function that computes the normal seq. z-test and returns the first index where the seq. test is over the threshold
norm_seq_test <- function(samples, comp_size, alpha){
  #define a vector to store the values of the seq.test and fill it with the sequential test values
  seq_tests = c(n = comp_size)
  for (i in 1:comp_size){
    seq_tests[i] = seq_test(samples, 1, i, comp_size, alpha)
  }
  
  #If there is no value in seq_tests that is >= threshold, return the comp_size + penalization to indicate a non-rejection
  if (length(which(seq_tests >= threshold)) == 0){
    return(comp_size + penalty)
  }
  #else return the index of the first value in seq_tests that is >= threshold, so leads to a rejection
  else{
    return(which(seq_tests >= threshold)[1])
  }
}

##compute the time we stop for the glued test
##i.e. we define the seq.test for the z-test based on comp_size/4 data points and stop as soon as it hits 0.2.
##Then we define a new seq.test for the next 1/2 of the data points and stop as soon as it hits 0.4.
##Then we define a new seq.test for the rest of the data points and stop as soon as it hits 0.6.
#We repeat this until we reach 0.9 and then we take the next 15 data points and compute the seq.test for all of them

#function that takes samples, computes the seq. test and gives back the first index where the seq. test is over the threshold
glued_seq_test <- function(samples, comp_size, alpha){
  seq_tests = c(n = comp_size)
  
  #compute sequential test using the first 1/2 of the data points
  for(i in 1:(ceiling(comp_size/4))){
    #compute the sequential test for all steps
    seq_tests[i] = seq_test(samples, 1, i, ceiling(comp_size/4), alpha)
  }
  #If we do not reach 0.2, return the comp_size + penalization to indicate a non-rejection
  if(length(which(seq_tests >= 0.2)) == 0) {
    return(comp_size + k*penalty)
  }
  #else, stopping_time = index of first value that is >= 0.2
  stopping_time = which(seq_tests >= 0.2)[1]
  #if we stop at the end and thus reject, return the last index of the current test
  if(stopping_time + 1 > ceiling(comp_size/4)){
    return(stopping_time)
  }
  #compute the sequential test of level >=0.2 using the first 1/2 of the data, starting from where we stopped with the first test
  for(i in (stopping_time + 1):(ceiling(comp_size/2))){
    #compute the sequential test for all steps
    seq_tests[i] = seq_test(samples, stopping_time + 1, i, (ceiling(comp_size/2)) - (stopping_time + 1) + 1, seq_tests[stopping_time])
  }
  #If we do not reach 0.5, return the comp_size + penalization to indicate a non-rejection
  if(length(which(seq_tests >= 0.5)) == 0) {
    return(comp_size + penalty)
  }
  #else, stopping_time = index of first value that is >= 0.5
  stopping_time = (which(seq_tests >= 0.5)[1])
  #if we stop at the end, return the last index of the current test
  if(stopping_time + 1 > ceiling(comp_size/2)){
    return(stopping_time)
  }
  
  #compute the sequential test of level >=0.5 using the first 3/4 of the data, starting from where we stopped with the previous test
  for(i in (stopping_time + 1):(ceiling(3*comp_size/4))){
    #compute the sequential test for all steps
    seq_tests[i] = seq_test(samples, stopping_time + 1, i, (ceiling(3*comp_size/4)) - (stopping_time + 1) + 1, seq_tests[stopping_time])
  }
  #If we do not reach 0.7, return the comp_size + penalization to indicate a non-rejection
  if(length(which(seq_tests >= 0.7)) == 0) {
    return(comp_size + penalty)
  }
  #else, stopping_time = index of first value that is >= 0.7
  stopping_time = (which(seq_tests >= 0.7)[1])
  #if we stop at the end, return the last index of the current test
  if(stopping_time + 1 > ceiling(3*comp_size/4)){
    return(stopping_time)
  }
  
  #For the rest of the data, compute the sequential test for level = value of seq. test at stopping time >= 0.7
  for(i in (stopping_time + 1):(comp_size)){
    #compute the sequential test for all steps
    seq_tests[i] = seq_test(samples, stopping_time + 1, i, comp_size - (stopping_time + 1) + 1, seq_tests[stopping_time])
  }
  #If we do not reach 0.9, return the comp_size + penalization to indicate a non-rejection
  if(length(which(seq_tests >= 0.9)) == 0) {
    return(comp_size + penalty)
  }
  #else, stopping_time = index of first value that is >= 0.9
  stopping_time = (which(seq_tests >= 0.9)[1])
  
  #if we stop at the end, return the last index of the current test
  if(stopping_time + 1 > comp_size){
    return(stopping_time)
  }
  
  #if there are no 15 data points left, return the last value of the current test
  if((stopping_time + 15) >= comp_size){
    if(seq_tests[comp_size] == 1){
      return(comp_size)
    }
    return(comp_size + penalty)
  }
  
  ##if there are at least 15 data points left, we define a new test based on the next 15 data points
  alpha = seq_tests[stopping_time] #save the alpha of the sequential test at the step where we reach the boundary
  quantile = qnorm(1 - alpha, mean = 0, sd = 1) #quantile for z-test with new alpha
  
  #compute the new test based on the next 15 data points
  newtest = 1/sqrt(15) * sum(samples[(stopping_time + 1):(stopping_time + 15)])
  #reject the null hypothesis if the new test is greater than the quantile
  if(newtest > quantile){
    return(stopping_time + 15)
  }
  #else we do not reject
  else{
    return(comp_size + penalty)
  }
}

###simulate
##compute the stop times for both tests for several samples
#define matrix to save stop times for both the normal seq. tests and the glued seq. tests 
stoping_times_normal = matrix(0, nrow = length(x), ncol = rep)
stoping_times_glued = matrix(0, nrow = length(x), ncol = rep)

for (i in x){
  for(j in 1:rep){
    #get samples from the true distribution
    samples = rnorm(sample_size, mean = i, sd = 1)
    
    #compute stop times for the normal seq. test and the glued seq. test
    stoping_times_normal[which(x == i), j] = norm_seq_test(samples, sample_size, alpha)
    stoping_times_glued[which(x == i), j] = glued_seq_test(samples, sample_size, alpha)
  }
}

###plot
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
  
  #save the amount of non rejections for the normal seq. test and the glued seq. test
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

