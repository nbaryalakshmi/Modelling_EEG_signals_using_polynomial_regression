#Determine a suitable mathematical model in explaining the relationship 
#between the audio signal y and the two brain signals x1 and x2
#X.csv file contains EEG signals.
#'x1' is prefrontal cortex and 'x2' is auditory cortex.
data_eeg_x = read.csv("X.csv")
data_eeg_x = data.matrix(data_eeg_x)
cat("No. of rows and columns in brain signal data:",nrow(data_eeg_x),",",ncol(data_eeg_x))
#y.csv file contains sound signal 'y'.
data_sound_y = read.csv("y.csv")
data_sound_y = data.matrix(data_sound_y)
cat("No. of rows and columns in sound signal data:",nrow(data_sound_y),",",ncol(data_sound_y))
#time.csv contains sampling time of all three signals in 'seconds'.
#2 minutes of signal with sampling frequency of 20Hz.
data_time = read.csv("time.csv")
data_time = data.matrix(data_time)

######################### Model 1 y=θ1*x1^3+θ2*x2^5+θbias+ε ######################### 
#Task 2.1#Estimate model parameters
matrix_of_ones = matrix(1, nrow(data_eeg_x), 1)
model1_X = cbind(data_eeg_x[,"x1"]^3, data_eeg_x[,"x2"]^5,matrix_of_ones)
model1_theta_hat = solve(t(model1_X) %*% model1_X) %*% t(model1_X) %*% data_sound_y[,"y"]
cat("Model 1 parameter values for θ1, θ2 and θbias:", model1_theta_hat)
#Task 2.2#RSS
model1_y_hat = model1_X %*% model1_theta_hat
model1_error = data_sound_y[,"y"] - model1_y_hat
model1_rss = sum(model1_error^2)
cat("Model 1 residual sum of squared errors:",model1_rss)
#Task 2.3#log-likelihood
model1_num_samples=nrow(data_eeg_x)
model1_resdl_varnc=model1_rss/(model1_num_samples-1)
model1_log_liklhd=-(model1_num_samples/2)*log(2*pi)-(model1_num_samples/2)*log(model1_resdl_varnc)-(1/(2*model1_resdl_varnc))*model1_rss
cat("Model 1 log-likelihood: ",model1_log_liklhd)
#Task 2.4#Aic & Bic
no_estmtd_params=3
model1_aic=(2*no_estmtd_params)-(2*model1_log_liklhd)
model1_bic=(no_estmtd_params*log(model1_num_samples))-(2*model1_log_liklhd)
cat("Model 1 Akaike information criterion: ",model1_aic)
cat("Model 1 Bayesian information criterion: ",model1_bic)
#Task 2.5
#Q-Q plot
qqnorm(model1_error, main = "Q-Q Plot for Model 1 residual",
       xlab = "Theoretical Quantiles", ylab = "Empirical Quantiles",
       plot.it = TRUE, datax = FALSE,col="blue")
qqline(model1_error,datax = FALSE, distribution = qnorm,col="red")
#histogram
hist(model1_error,breaks=10,col="blue",xlab ="Model 1 residual",
     main="Histogram plot for Model 1 residual")
#Shapiro-Wilk test
shapiro.test(model1_error)
#Kolmogorov-Smirnov test
ks.test(model1_error,"pnorm")

###################################################
########################## Model 2#y=θ1*x1^4 + θ2*x2^2 + θbias + ε ######################### 
#Task 2.1#Estimate model parameters
model2_X = cbind(data_eeg_x[,"x1"]^4, data_eeg_x[,"x2"]^2, matrix_of_ones)
model2_theta_hat = solve(t(model2_X) %*% model2_X) %*% t(model2_X) %*% data_sound_y[,"y"]
cat("Model 2 parameter values for θ1, θ2 and θbias:", model2_theta_hat)
#Task 2.2#RSS
model2_y_hat = model2_X %*% model2_theta_hat
model2_error = data_sound_y[,"y"] - model2_y_hat
model2_rss = sum(model2_error^2)
cat("Model 2 residual sum of squared errors:",model2_rss)
#Task 2.3#log-likelihood
model2_num_samples=nrow(data_eeg_x)
model2_resdl_varnc=model2_rss/(model2_num_samples-1)
model2_log_liklhd=-(model2_num_samples/2)*log(2*pi)-(model2_num_samples/2)*log(model2_resdl_varnc)-(1/(2*model2_resdl_varnc))*model2_rss
cat("Model 2 log-likelihood: ",model2_log_liklhd)
#Task 2.4#Aic & Bic
no_estmtd_params=3
model2_aic=(2*no_estmtd_params)-(2*model2_log_liklhd)
model2_bic=(no_estmtd_params*log(model2_num_samples))-(2*model2_log_liklhd)
cat("Model 2 Akaike information criterion: ",model2_aic)
cat("Model 2 Bayesian information criterion: ",model2_bic)
#Task 2.5
#Q-Q plot
qqnorm(model2_error, main = "Q-Q Plot for Model 2 residual",
       xlab = "Theoretical Quantiles", ylab = "Empirical Quantiles",
       plot.it = TRUE, datax = FALSE,col="blue")
qqline(model2_error,datax = FALSE, distribution = qnorm,col="red")
#histogram
hist(model2_error,breaks=10,col="blue",xlab ="Model 2 residual",
     main="Histogram plot for Model 2 residual")
#Shapiro-Wilk test
shapiro.test(model2_error)
#Kolmogorov-Smirnov test
ks.test(model2_error,"pnorm")

###################################################
########################## Model 3 #y=θ1*x1^3 + θ2*x2 + θ3*x1 + θbias + ε ######################### 
#Task 2.1#Estimate model parameters
model3_X = cbind(data_eeg_x[,"x1"]^3, data_eeg_x[,"x2"], data_eeg_x[,"x1"], matrix_of_ones)
model3_theta_hat = solve(t(model3_X) %*% model3_X) %*% t(model3_X) %*% data_sound_y[,"y"]
cat("Model 3 parameter values for θ1, θ2, θ3 and θbias:", model3_theta_hat)
#Task 2.2#RSS
model3_y_hat = model3_X %*% model3_theta_hat
model3_error = data_sound_y[,"y"] - model3_y_hat
model3_rss = sum(model3_error^2)
cat("Model 3 residual sum of squared errors:",model3_rss)
#Task 2.3#log-likelihood
model3_num_samples=nrow(data_eeg_x)
model3_resdl_varnc=model3_rss/(model3_num_samples-1)
model3_log_liklhd=-(model3_num_samples/2)*log(2*pi)-(model3_num_samples/2)*log(model3_resdl_varnc)-(1/(2*model3_resdl_varnc))*model3_rss
cat("Model 3 log-likelihood: ",model3_log_liklhd)
#Task 2.4#Aic & Bic
no_estmtd_params=4
model3_aic=(2*no_estmtd_params)-(2*model3_log_liklhd)
model3_bic=(no_estmtd_params*log(model3_num_samples))-(2*model3_log_liklhd)
cat("Model 3 Akaike information criterion: ",model3_aic)
cat("Model 3 Bayesian information criterion: ",model3_bic)
#Task 2.5
#Q-Q plot
qqnorm(model3_error, main = "Q-Q Plot for Model 3 residual",
       xlab = "Theoretical Quantiles", ylab = "Empirical Quantiles",
       plot.it = TRUE, datax = FALSE,col="blue")
qqline(model3_error,datax = FALSE, distribution = qnorm,col="red")
#histogram
hist(model3_error,breaks=10,col="blue",xlab ="Model 3 residual",
     main="Histogram plot for Model 3 residual")
#Shapiro-Wilk test
shapiro.test(model3_error)
#Kolmogorov-Smirnov test
ks.test(model3_error,"pnorm")

###################################################
#########################  Model 4 #y=θ1*x1 + θ2*x1^2 + θ3*x1^3 + θ4*x2^3 + θbias + ε ######################### 
#Task 2.1#Estimate model parameters
model4_X = cbind(data_eeg_x[,"x1"], data_eeg_x[,"x1"]^2, data_eeg_x[,"x1"]^3, data_eeg_x[,"x2"]^3, matrix_of_ones)
model4_theta_hat = solve(t(model4_X) %*% model4_X) %*% t(model4_X) %*% data_sound_y[,"y"]
cat("Model 4 parameter values for θ1, θ2, θ3, θ4 and θbias:", model4_theta_hat)
#Task 2.2#RSS
model4_y_hat = model4_X %*% model4_theta_hat
model4_error = data_sound_y[,"y"] - model4_y_hat
model4_rss = sum(model4_error^2)
cat("Model 4 residual sum of squared errors:",model4_rss)
#Task 2.3#log-likelihood
model4_num_samples=nrow(data_eeg_x)
model4_resdl_varnc=model4_rss/(model4_num_samples-1)
model4_log_liklhd=-(model4_num_samples/2)*log(2*pi)-(model4_num_samples/2)*log(model4_resdl_varnc)-(1/(2*model4_resdl_varnc))*model4_rss
cat("Model 4 log-likelihood: ",model4_log_liklhd)
#Task 2.4#Aic & Bic
no_estmtd_params=5
model4_aic=(2*no_estmtd_params)-(2*model4_log_liklhd)
model4_bic=(no_estmtd_params*log(model4_num_samples))-(2*model4_log_liklhd)
cat("Model 4 Akaike information criterion: ",model4_aic)
cat("Model 4 Bayesian information criterion: ",model4_bic)
#Task 2.5
#Q-Q plot
qqnorm(model4_error, main = "Q-Q Plot for Model 4 residual",
       xlab = "Theoretical Quantiles", ylab = "Empirical Quantiles",
       plot.it = TRUE, datax = FALSE,col="blue")
qqline(model4_error,datax = FALSE, distribution = qnorm,col="red")
#histogram
hist(model4_error,breaks=10,col="blue",xlab ="Model 4 residual",
     main="Histogram plot for Model 4 residual")
#Shapiro-Wilk test
shapiro.test(model4_error)
#Kolmogorov-Smirnov test
ks.test(model4_error,"pnorm")

###################################################
######################### Model 5 #y=θ1*x1^3 + θ2*x1^4 + θ3*x2 + θbias + ε ######################### 
#Task 2.1#Estimate model parameters
model5_X = cbind(data_eeg_x[,"x1"]^3, data_eeg_x[,"x1"]^4, data_eeg_x[,"x2"], matrix_of_ones)
model5_theta_hat = solve(t(model5_X) %*% model5_X) %*% t(model5_X) %*% data_sound_y[,"y"]
cat("Model 5 parameter values for θ1, θ2, θ3 and θbias:", model5_theta_hat)
#Task 2.2#RSS
model5_y_hat = model5_X %*% model5_theta_hat
model5_error = data_sound_y[,"y"] - model5_y_hat
model5_rss = sum(model5_error^2)
cat("Model 5 residual sum of squared errors:",model5_rss)
#Task 2.3#log-likelihood
model5_num_samples=nrow(data_eeg_x)
model5_resdl_varnc=model5_rss/(model5_num_samples-1)
model5_log_liklhd=-(model5_num_samples/2)*log(2*pi)-(model5_num_samples/2)*log(model5_resdl_varnc)-(1/(2*model5_resdl_varnc))*model5_rss
cat("Model 5 log-likelihood: ",model5_log_liklhd)
#Task 2.4#Aic & Bic
no_estmtd_params=4
model5_aic=(2*no_estmtd_params)-(2*model5_log_liklhd)
model5_bic=(no_estmtd_params*log(model5_num_samples))-(2*model5_log_liklhd)
cat("Model 5 Akaike information criterion: ",model5_aic)
cat("Model 5 Bayesian information criterion: ",model5_bic)
#Task 2.5
#Q-Q plot
qqnorm(model5_error, main = "Q-Q Plot for Model 5 residual",
       xlab = "Theoretical Quantiles", ylab = "Empirical Quantiles",
       plot.it = TRUE, datax = FALSE,col="blue")
qqline(model5_error,datax = FALSE, distribution = qnorm,col="red")
#histogram
hist(model5_error,breaks=10,col="blue",xlab ="Model 5 residual",
     main="Histogram plot for Model 5 residual")
#Shapiro-Wilk test
shapiro.test(model5_error)
#Kolmogorov-Smirnov test
ks.test(model5_error,"pnorm")

###################################################
#Task 2.6#select best regression model
cat("Model 3 is the best performing model among the five given models based on the results.")

#Task 2.7#train-test-split
#install.packages("dplyr")
library(dplyr)
set.seed(1)
#create id column to split by id
matrix_id=matrix(1:nrow(data_eeg_x), nrow(data_eeg_x), 1)
colnames(matrix_id)="Id"
data_eeg_sound=cbind(matrix_id,data_eeg_x,data_sound_y)
data_eeg_sound=data.frame(data_eeg_sound)
#get 70% data to train
data_train=data_eeg_sound%>%dplyr::sample_frac(0.70)
#get 30% data to test
data_test=dplyr::anti_join(data_eeg_sound,data_train,by="Id")
cat("No. of rows and columns in original dataset: ",dim(data_eeg_sound))
cat("No. of rows and columns in train dataset: ",dim(data_train))
cat("No. of rows and columns in test dataset: ",dim(data_test))

#Task 2.7.1#Estimate model parameters using the training dataset
#Best model#Model 3 #y=θ1*x1^3 + θ2*x2 + θ3*x1 + θbias + ε
matrix_of_ones = matrix(1, nrow(data_train), 1)
modelBest_X = cbind(data_train[,"x1"]^3, data_train[,"x2"], data_train[,"x1"], matrix_of_ones)
modelBest_theta_hat = solve(t(modelBest_X) %*% modelBest_X) %*% t(modelBest_X) %*% data_train[,"y"]
cat("Best Model (train data) parameter values for θ1, θ2, θ3 and θbias:", modelBest_theta_hat)
#Task 2.7.2#Compute the model's prediction on the testing data
matrix_of_ones_test = matrix(1, nrow(data_test), 1)
modelBest_X_test = cbind(data_test[,"x1"]^3, data_test[,"x2"], data_test[,"x1"], matrix_of_ones_test)
modelBest_y_hat = modelBest_X_test %*% modelBest_theta_hat
modelBest_error = data_test[,"y"] - modelBest_y_hat
modelBest_rss = sum(modelBest_error^2)
cat("Best Model (Test data) residual sum of squared errors:",modelBest_rss)
cat("Predicted Values:", modelBest_y_hat)
#Task 2.7.3#Compute 95% confidence interval, plot them with error bars 
#together with model prediction and testing data samples.
#Va(yhat)=sigma^2*x(X^T*X)^-1*x^T
modeltest_num_samples=nrow(data_test)
modeltest_resdl_varnc=modelBest_rss/(modeltest_num_samples-1)#variance
X_trnsp_X_invrs=solve(t(modelBest_X_test)%*%modelBest_X_test)
var_yhat=matrix(1,modeltest_num_samples,1)
for(i in 1:modeltest_num_samples)
{
  x_i=matrix(modelBest_X_test[i,],1,4)
  #print(x_i)
  var_yhat[i,1]=modeltest_resdl_varnc*x_i%*%X_trnsp_X_invrs%*%t(x_i)
}
#conf_intrvl=yhat(+/-)1.96*sqrt(Var(yhat))
conf_intrvl=1.96*(sqrt(var_yhat))
cf_upper=modelBest_y_hat+conf_intrvl
cf_lower=modelBest_y_hat-conf_intrvl
#plot confidence interval
#split time data
data_time_id=cbind(matrix_id,data_time)
data_time_id=data.frame(data_time_id)
data_time_train=data_time_id%>%dplyr::sample_frac(0.70)
data_time_test=dplyr::anti_join(data_time_id,data_time_train,by="Id")
#all points
plot(x=data_time_test[,"time"],modelBest_y_hat,col="orange",
     main="Confidence interval (all data)",
     xlab="Sampling Time in seconds",ylab="Predicted values")#predicted
#points(x=data_time_test[,"time"],data_test[,"y"],col="green")#actual
segments(data_time_test[,"time"],cf_lower,data_time_test[,"time"],cf_upper,col = "blue",lwd=2)#ci
legend(0,95,legend=c("Predicted values","Confidence Interval"),fill=c("orange","blue"))
#first 50 points
rng=1:50
plot(x=data_time_test[rng,"time"],modelBest_y_hat[rng,],col="orange",lwd=2,
     main="Confidence interval (first 50 points)",
     xlab="Sampling Time in seconds",ylab="Predicted values")#predicted
points(x=data_time_test[rng,"time"],data_test[rng,"y"],col="green",lwd=2)#actual
segments(data_time_test[rng,"time"],cf_lower[rng],data_time_test[rng,"time"],cf_upper[rng],
         col = "blue",lwd=2)#ci
legend(0,95,legend=c("Predicted values","Confidence Interval","Actual values"),
       fill=c("orange","blue","green"))

#install.packages("plotrix")
#library("plotrix")# Load library plotrix
#rng=1:25
#plotCI(x=data_test[rng,"x1"],y=modelBest_y_hat[rng,],li=cf_lower[rng],ui=cf_upper[rng],col="orange")#predicted and ci
#points(data_test[rng,"x1"],data_test[rng,"y"],col="blue")#test samples
#
#plotCI(x=data_test[rng,"x2"],y=data_test[rng,"y"],li=cf_lower[rng],ui=cf_upper[rng],col="blue")
#points(data_test[rng,"x2"],modelBest_y_hat[rng],col="red")

#Task 3
#compute posterior distribution of best model using rejection abc
#Model 3: y = 2.715713*𝑥1^3+-3.15135*𝑥2+4.18139*𝑥1+-6.6514+ε
#The two parameters with largest absolute values are 2.715713 and 4.18139
set.seed(1)
prior_x1cube=runif(n=1500,min=2.715713-1,max=2.715713+1)
prior_x1=runif(n=1500,min=4.18139-1,max=4.18139+1)
#perform rejection abc
selctd_x1cube=matrix(1,modeltest_num_samples,1)
selctd_x1=matrix(1,modeltest_num_samples,1)
res_all=matrix(1,1500,3)
colnames(res_all)=c("x1cube","x1","varnc")

for(i in 1:1500)
{
  theta_abc=matrix(c(prior_x1cube[i],-3.15135,prior_x1[i],-6.6514),4,1)
  y_hat_abc = modelBest_X_test %*% theta_abc
  model_abc_error = data_test[,"y"] - y_hat_abc
  model_abc_rss = sum(model_abc_error^2)
  model_abc_varnc=model_abc_rss/(modeltest_num_samples-1)#variance
  res_all[i,]=c(prior_x1cube[i],prior_x1[i],model_abc_varnc)
}
c=0
for(i in 1:1500)
{
  if(res_all[i,"varnc"]<=modeltest_resdl_varnc+2.6)
  {
    c=c+1
    selctd_x1cube[c,1]=res_all[i,"x1cube"]
    selctd_x1[c,1]=res_all[i,"x1"]
  }
  if(c>=modeltest_num_samples)
    break
}
#print("c")
#print(c)

#plot joint and marginal posterior distribution 
#selctd_x1cube#2.715713*𝑥1^3
dt_mean=mean(selctd_x1cube)
dt_sd=sd(selctd_x1cube)
dt_densty1=dnorm(selctd_x1cube,mean=dt_mean,sd=dt_sd)
plot(selctd_x1cube,dt_densty1,lty=1,xlab="First Parameter",ylab="Distribution and Density",col="blue",
     main="Marginal Posterior Distribution (First Parameter)")
dt_dstrbtn1=pnorm(selctd_x1cube,mean=dt_mean,sd=dt_sd)
points(selctd_x1cube,dt_dstrbtn1,lty=1,col="green")
legend(2,1.2,legend=c("Density","Distribution"),
       fill=c("blue","green"))
hist(selctd_x1cube,breaks=10,xlab="First Parameter",ylab="Frequency",
     main="Histogram for marginal Posterior Distribution (First Parameter)")

#selctd_x1[c,1]#4.18139*𝑥1
dt_mean=mean(selctd_x1)
dt_sd=sd(selctd_x1)
dt_densty2=dnorm(selctd_x1,mean=dt_mean,sd=dt_sd)
dt_dstrbtn2=pnorm(selctd_x1,mean=dt_mean,sd=dt_sd)
plot(selctd_x1,dt_dstrbtn2,lty=1,xlab="Second Parameter",ylab="Distribution and Density",col="green",
     main="Marginal Posterior Distribution (Second Parameter)")
points(selctd_x1,dt_densty2,lty=1,col="blue")
legend(2,1.2,legend=c("Density","Distribution"),
       fill=c("blue","green"))
hist(selctd_x1,breaks=10,xlab="Second Parameter",ylab="Frequency",
     main="Histogram for marginal Posterior Distribution (Second Parameter)")

#Joint Posterior
plot(selctd_x1cube, selctd_x1,lty=1,xlab="First Parameter",ylab="Second Parameter",col="blue",
     main="Joint Posterior Distribution plot")
########## 
#install.packages("LaplacesDemon")
library(LaplacesDemon)
joint.density.plot(selctd_x1cube, selctd_x1,
                   Title="Joint Density Plot for first and second parameter",
                   contour=TRUE, color=FALSE)


