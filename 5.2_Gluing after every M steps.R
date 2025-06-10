#plot the time we stop for the seq. z-test if we glue it after every M steps
#i.e. computing the sequential z-test based on step_size data points,
#where the sequential z-test is defined on the z-test with step_size + 1 values
#plot it for rep samples for different values of mü
#for testing mü <= 0 against mü > 0

set.seed(9) #for reproducibility

###load libraries
library(ggplot2)
library(patchwork)
library(latex2exp)


###set parameters
rep = 2000 #how often we take samples
sample_size = 500 #size of one sample
step_size = 50 #after every step_size data points, we define a new test (define s.t. sample_size = step_size * m for some m)
alpha = 0.05 #starting level of z-test (and sequential test)
penalty = 10 #penalization for the stop time if we don't reject the null, used to make plots interpretable
threshold = 0.9999 #value, for which we reject
x = c(0, 0.1, 0.2, 0.3, 0.4) #values of mü


###define functions
##compute the sequential test from the z-test (as in Chapter 3) and compute the time we stop
#sequential test defined from the z-tests for data points start:end
#N = number of data points for traditional z-test
seq_test <- function(samples, start, end, N, alpha){
  #compute the quantile for the level alpha sequential test
  quantile = qnorm(1 - alpha, mean = 0, sd = 1)
  
  #take a set of step_size data points and define z-test if we have step_size + 1 data points
  numerator = sum(samples[start:end]) - sqrt(N)*quantile
  denominator = sqrt(N - (end-start + 1))
  return(pnorm(numerator/denominator))
}

#define function that computes the normal seq. z-test and returns the first index where the seq. test is over the threshold
norm_seq_test <- function(samples, comp_size, alpha){
  #compute the sequential test for all steps
  #define a vector to store the values of the seq.test
  #and fill it with the sequential test values
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

##sequential test with gluing: define seq. z-test for step_size + 1 data points and define a new test after step_size points
##i.e. we define the seq. test for a sample of size step_size + 1, but stop one step before its finished
##then we can use the value of the seq. test after step_size data points as a new level

#function that takes samples, computes the seq. test and gives back the first index where the seq. test is over the threshold
glued_seq_test <- function(samples, comp_size, step_size, alpha){
  #compute quantile for the first z-test
  quantile = qnorm(1 - alpha, mean = 0, sd = 1)
  
  #define a vector to store the values of the seq.test
  seq_test_glued = c(n = comp_size)
  
  #loop for a new seq.test after every step_size
  for(j in 1:(comp_size/step_size)){
    #loop for next step_size data points, to compute the values of the current sequential test
    #if we are at the last step_size, we need to compute the seq. test for the z-test based on step_size data points, s.t. it is either 0 or 1
    if(j == (comp_size/step_size)){
      for(i in 1:step_size){
        seq_test_glued[(j-1)*step_size + i] = seq_test(samples, (j-1)*step_size + 1, (j-1)*step_size + i, step_size, alpha)
      }
    }
    #else we compute the seq. test for the z-test based on step_size + 1 data points
    else{
      for(i in 1:step_size){
        #compute the sequential test for each step_size data points
        seq_test_glued[(j-1)*step_size + i] = seq_test(samples, (j-1)*step_size + 1, (j-1)*step_size + i, step_size + 1, alpha)
      }
    }
    #use the last value of the seq. test as the new alpha for the next seq. test
    alpha = seq_test_glued[(j*step_size)]
  }

  #If there is no value in seq_tests that is >= threshold, return the comp_size + penalization to indicate a non-rejection
  if(length(which(seq_test_glued >= threshold)) == 0){
    return(comp_size + penalty)
  }
  #return the index of the first value where the glued seq. test is >= threshold
  else{
    return(which(seq_test_glued >= threshold)[1])
  }
} 

###simulate
#define matrix to save stop times for both the normal seq. tests and the glued seq. tests 
stoping_times_normal = matrix(0, nrow = length(x), ncol = rep)
stoping_times_glued = matrix(0, nrow = length(x), ncol = rep)

for (i in x){
  for(j in 1:rep){
    #get samples from the true distribution
    samples = rnorm(sample_size, mean = i, sd = 1)
    
    #compute stop times for the normal seq. test and the glued seq. test
    stoping_times_normal[which(x == i), j] = norm_seq_test(samples, sample_size, alpha)
    stoping_times_glued[which(x == i), j] = glued_seq_test(samples, sample_size, step_size, alpha)
  }
}
###plot
##make Histogram-plots for different values of mü
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
