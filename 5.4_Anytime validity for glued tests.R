#5.4 Anytime Validity of the Simulations
#plot the rejection rates of the glued sequential tests from Sections 5.1, 5.2 and 5.3

set.seed(9) #for reproducibility

###load libraries
library(ggplot2)

###set parameters
alphas = seq(0.01, 0.1, 0.005) #alpha values for the tests
rep = 2000 #number of repetitions for each alpha
sample_size = 500#size of each of the samples
threshold = 0.9999 #threshold for the sequential test
step_size = 50 #size of the step for test 2
size_new = 15 #size ot the new z-test in test1


###define functions
#sequential test function used for the glued test in section 5.1
seq_test_1 <- function(samples, n, sample_size, alpha){
  #compute the upper quantile of the z-test
  quantile = qnorm(1 - alpha, mean = 0, sd = 1)
  
  numerator = sum(samples[1:n]) - sqrt(sample_size) * quantile 
  difference = sqrt(sample_size - n)
  return(pnorm(numerator/difference))
}

#sequential test function used for the glued tests in sections 5.2 and 5.3
seq_test_23 <- function(samples, start, end, N, alpha){
  #compute the quantile for the level alpha sequential test
  quantile = qnorm(1 - alpha, mean = 0, sd = 1)
  
  #take a set of step_size data points and define z-test if we have step_size + 1 data points
  numerator = sum(samples[start:end]) - sqrt(N)*quantile
  denominator = sqrt(N - (end-start + 1))
  return(pnorm(numerator/denominator))
}


##compute all three glued sequential tests and return if they reject(1) or not (0)
#compute the glued sequential test from section 5.1
test_1 <- function(samples, sample_size, alpha, boundary){
  Y = numeric(sample_size)
  #compute the normal sequential test for each step
  for(i in 1:(sample_size)){
    Y[i] = seq_test_1(samples, i, sample_size, alpha)
  }
  
  #find the index such that the value of the sequential test is >= boundary
  if(length(which(Y >= boundary)) == 0) {
    return (0)
  }
  else{
    stop_time <- which(Y >= boundary)[1]
  }
  #if there are less than size_new data points left, return the value of the normal seq. test
  if(stop_time + size_new > sample_size && Y[sample_size] == 1) {
    return (Y[sample_size])
  }
  if(stop_time + size_new > sample_size && Y[sample_size] == 0) {
    return (0)
  }
  
  ##if there are at least 10 data points left, we define a new test based on the next 10 data points
  alpha = Y[stop_time] #save the alpha of the sequential test at the step where we reach the boundary
  quantile = qnorm(1 - Y[stop_time], mean = 0, sd = 1) #quantile for z-test with new alpha
  
  #compute the new test based on the next 10 data points
  newtest = 1/sqrt(size_new) * sum(samples[(stop_time + 1):(stop_time + size_new)])
  #reject the null hypothesis if the new test is greater than the quantile
  if(newtest > quantile){
    return(1)
  }
  #else we do not reject
  else{
    return(0)
  }
}

#compute the glued sequential test from section 5.2
test_2 <- function(samples, comp_size, step_size, alpha){
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
        seq_test_glued[(j-1)*step_size + i] = seq_test_23(samples, (j-1)*step_size + 1, (j-1)*step_size + i, step_size, alpha)
      }
    }
    #else we compute the seq. test for the z-test based on step_size + 1 data points
    else{
      for(i in 1:step_size){
        #compute the sequential test for each step_size data points
        seq_test_glued[(j-1)*step_size + i] = seq_test_23(samples, (j-1)*step_size + 1, (j-1)*step_size + i, step_size + 1, alpha)
      }
    }
    #use the last value of the seq. test as the new alpha for the next seq. test
    alpha = seq_test_glued[(j*step_size)]
  }
  
  #If there is no value in seq_tests that is >= threshold, return the comp_size + penalization to indicate a non-rejection
  if(length(which(seq_test_glued >= threshold)) == 0){
    return(0)
  }
  #return the index of the first value where the glued seq. test is >= threshold
  else{
    return(1)
  }
}

#compute the glued sequential test from section 5.3
test_3 <- function(samples, comp_size, alpha){
  seq_tests = c(n = comp_size)
  
  #compute sequential test using the first 1/2 of the data points
  for(i in 1:(ceiling(comp_size/4))){
    #compute the sequential test for all steps
    seq_tests[i] = seq_test(samples, 1, i, ceiling(comp_size/4), alpha)
  }
  #If we do not reach 0.2, return the comp_size + penalization to indicate a non-rejection
  if(length(which(seq_tests >= 0.2)) == 0) {
    return(0)
  }
  #else, stopping_time = index of first value that is >= 0.2
  stopping_time = which(seq_tests >= 0.2)[1]
  #if we stop at the end and thus reject, return the last index of the current test
  if(stopping_time + 1 > ceiling(comp_size/4)){
    return(1)
  }
  #compute the sequential test of level >=0.2 using the first 1/2 of the data, starting from where we stopped with the first test
  for(i in (stopping_time + 1):(ceiling(comp_size/2))){
    #compute the sequential test for all steps
    seq_tests[i] = seq_test(samples, stopping_time + 1, i, (ceiling(comp_size/2)) - (stopping_time + 1) + 1, seq_tests[stopping_time])
  }
  #If we do not reach 0.5, return the comp_size + penalization to indicate a non-rejection
  if(length(which(seq_tests >= 0.5)) == 0) {
    return(0)
  }
  #else, stopping_time = index of first value that is >= 0.5
  stopping_time = (which(seq_tests >= 0.5)[1])
  #if we stop at the end, return the last index of the current test
  if(stopping_time + 1 > ceiling(comp_size/2)){
    return(1)
  }
  
  #compute the sequential test of level >=0.5 using the first 3/4 of the data, starting from where we stopped with the previous test
  for(i in (stopping_time + 1):(ceiling(3*comp_size/4))){
    #compute the sequential test for all steps
    seq_tests[i] = seq_test(samples, stopping_time + 1, i, (ceiling(3*comp_size/4)) - (stopping_time + 1) + 1, seq_tests[stopping_time])
  }
  #If we do not reach 0.7, return the comp_size + penalization to indicate a non-rejection
  if(length(which(seq_tests >= 0.7)) == 0) {
    return(0)
  }
  #else, stopping_time = index of first value that is >= 0.7
  stopping_time = (which(seq_tests >= 0.7)[1])
  #if we stop at the end, return the last index of the current test
  if(stopping_time + 1 > ceiling(3*comp_size/4)){
    return(1)
  }
  
  #For the rest of the data, compute the sequential test for level = value of seq. test at stopping time >= 0.7
  for(i in (stopping_time + 1):(comp_size)){
    #compute the sequential test for all steps
    seq_tests[i] = seq_test(samples, stopping_time + 1, i, comp_size - (stopping_time + 1) + 1, seq_tests[stopping_time])
  }
  #If we do not reach 0.9, return the comp_size + penalization to indicate a non-rejection
  if(length(which(seq_tests >= 0.9)) == 0) {
    return(0)
  }
  #else, stopping_time = index of first value that is >= 0.9
  stopping_time = (which(seq_tests >= 0.9)[1])
  
  #if we stop at the end, return the last index of the current test
  if(stopping_time + 1 > comp_size){
    return(1)
  }
  
  #if there are no 15 data points left, return the last value of the current test
  if((stopping_time + 15) >= comp_size){
    if(seq_tests[comp_size] == 1){
      return(1)
    }
    return(0)
  }
  
  ##if there are at least 15 data points left, we define a new test based on the next 15 data points
  alpha = seq_tests[stopping_time] #save the alpha of the sequential test at the step where we reach the boundary
  quantile = qnorm(1 - alpha, mean = 0, sd = 1) #quantile for z-test with new alpha
  
  #compute the new test based on the next 15 data points
  newtest = 1/sqrt(15) * sum(samples[(stopping_time + 1):(stopping_time + 15)])
  #reject the null hypothesis if the new test is greater than the quantile
  if(newtest > quantile){
    return(1)
  }
  #else we do not reject
  else{
    return(0)
  }
}

###simulate
#Initialize matrices
rejections_1 <- matrix(0, nrow = length(alphas), ncol = rep)
rejections_2 <- matrix(0, nrow = length(alphas), ncol = rep)
rejections_3 <- matrix(0, nrow = length(alphas), ncol = rep)

#Loop through repetitions and alphas
for (i in 1:rep){
  samples<- rnorm(sample_size, mean = 0, sd = 1)
  for(j in seq_along(alphas)){
    alpha <- alphas[j]
    rejections_1[j, i] <- test_1(samples, sample_size, alpha, threshold)
    rejections_2[j, i] <- test_2(samples, sample_size, step_size, alpha)
    rejections_3[j, i] <- test_3(samples, sample_size, alpha)    
  }
}

#compute empirical Type I errors
rejection_rate_1 = rowSums(rejections_1)/rep
rejection_rate_2 = rowSums(rejections_2)/rep
rejection_rate_3 = rowSums(rejections_3)/rep

#prepare data frame for ggplot
rejection_rate_data = data.frame(
  alpha = rep(alphas, 3), # We now have 3 sets of data
  rejection_rate = c(rejection_rate_1, rejection_rate_2, rejection_rate_3),
  test_type = rep(c("5.1 gluing for early rejection", "5.2 gluing every 50 steps", "5.3 gluing for increasing levels"), each = length(alphas))
)

#plot the rejection rates
ggplot(rejection_rate_data, aes(x = alpha, y = rejection_rate, color = test_type)) +
  geom_line() +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  labs(title = "Type I Error Plot",
       x = "Level",
       y = "Rejection Rate") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "orange", "darkgreen")) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
