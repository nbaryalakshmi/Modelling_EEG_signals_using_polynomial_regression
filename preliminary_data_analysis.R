################ Read data ################ 
#X.csv file contains EEG signals.
#'x1' is prefrontal cortex and 'x2' is auditory cortex.
data_eeg = read.csv("X.csv")
data_eeg = data.matrix(data_eeg)
#Rename columns
colnames(data_eeg)=c("prefrontal","auditory")
#check if null values present
cat("Sum of null values in EEG signal data:",sum(is.na(data_eeg)))
print("First five rows in EEG signals file:")
print(data_eeg[1:5,])

#y.csv file contains sound signal 'y'.
data_sound = read.csv("y.csv")
data_sound = data.matrix(data_sound)
#Rename columns
colnames(data_sound)=c("sound_signal")
#check if null values present
cat("Sum of null values in Sound signal data:",sum(is.na(data_sound)))
print("First five rows in Sound signal file:")
print(data_sound[1:5,])

#time.csv contains sampling time of all three signals in 'seconds'.
#2 minutes of signal with sampling frequency of 20Hz.
data_time = read.csv("time.csv")
data_time = data.matrix(data_time)
cat("Sum of null values in Time data:",sum(is.na(data_time)))
print("First five rows in time data file:")
print(data_time[1:5,])

#All signals are subject to additive noise(assume independent
#and identically distributed Gaussian with zero mean)
#with unknown variance due to distortions during recording.

################ Preliminary data analysis ################ 
#Time series plots (audio and EEG signals)
#Prefrontal
plot(data_time,data_eeg[,"prefrontal"], type="l",
     xlab="Time in seconds (Sampling frequency of 20Hz)",
     ylab="EEG signal (Prefrontal)",
     main="Time series plot for EEG signal (Prefrontal)")
#Auditory
plot(data_time,data_eeg[,"auditory"], type="l",
     xlab="Time in seconds (Sampling frequency of 20Hz)",
     ylab="EEG signal (Auditory)",
     main="Time series plot for EEG signal (Auditory)")
#Sound
plot(data_time,data_sound, type="l",
     xlab="Time in seconds (Sampling frequency of 20Hz)",
     ylab="Sound signal",
     main="Time series plot for Sound signal")

#Distribution for each signal
func_density=function(prm_data, x_label, y_label)
{
  dt_mean=mean(prm_data)
  dt_sd=sd(prm_data)
  dt_densty=dnorm(prm_data,mean=dt_mean,sd=dt_sd)
  plot(prm_data,dt_densty,lty=1,xlab=x_label,ylab=y_label,
       main="Standard Normal Density")
}
func_distribution=function(prm_data, x_label, y_label)
{
  dt_mean=mean(prm_data)
  dt_sd=sd(prm_data)
  dt_dstrbtn=pnorm(prm_data,mean=dt_mean,sd=dt_sd)
  plot(prm_data,dt_dstrbtn,lty=1,xlab=x_label,ylab=y_label,
       main="Standard Normal Distribution")
}
#Prefrontal Distribution
func_density(data_eeg[,"prefrontal"],"Prefrontal","Density")
func_distribution(data_eeg[,"prefrontal"],"Prefrontal","Distribution")

#Auditory Distribution
func_density(data_eeg[,"auditory"],"Auditory","Density")
func_distribution(data_eeg[,"auditory"],"Auditory","Distribution")

#Sound Distribution
func_density(data_sound,"Sound","Density")
func_distribution(data_sound,"Sound","Distribution")

#Central Tendencies, Scale, Skewness
fun_central_tend=function(prm_data,x_name)
{
  dt_mean=mean(prm_data)
  dt_median=median(prm_data)
  hist(prm_data,breaks=10,xlab=x_name,main=paste("Histogram for",x_name))
  abline(v = dt_mean, lwd = 5 , col="cyan")
  abline(v = dt_median, lwd = 3 , col = "orange")
  dt_mode = unique(prm_data)
  dt_mode = dt_mode[which.max(tabulate(match(prm_data, dt_mode)))]
  abline(v = dt_mode, lwd = 3 , col="blue")
  cat("Mean: ",dt_mean)
  cat("\nMedian: ",dt_median)
  cat("\nMode: ",dt_mode)
}
fun_scale=function(prm_data)
{
  dt_max = max(prm_data)
  dt_min = min(prm_data)
  dt_range = dt_max - dt_min
  dt_variance = var(prm_data)
  dt_std = sd(prm_data)
  cat("Sample Maximum: ",dt_max,"\nSample Minimum: ",dt_min,"\nSample Range: ",dt_range)
  cat("Sample Variance: ",dt_variance,"\nSample Standard Deviation: ",dt_std)
}
fun_skewness=function(prm_data)
{
  n = length(prm_data)
  dt_skewness = (sqrt(n) * sum( (prm_data - mean(prm_data))^3 ) ) / (sqrt( sum( (prm_data - mean(prm_data))^2 ) ))^3
  cat("Skewness: %s",dt_skewness)
}

#Prefrontal
fun_central_tend(data_eeg[,"prefrontal"],"Prefrontal")
fun_scale(data_eeg[,"prefrontal"])
fun_skewness(data_eeg[,"prefrontal"])

#Auditory
fun_central_tend(data_eeg[,"auditory"],"Auditory")
fun_scale(data_eeg[,"auditory"])
fun_skewness(data_eeg[,"auditory"])

#Sound
fun_central_tend(data_sound,"Sound")
fun_scale(data_sound)
fun_skewness(data_sound)

#Box plot
boxplot(data_eeg)
#outliers for prefrontal
boxplot.stats(data_eeg[,"prefrontal"])$out
#outliers for auditory
boxplot.stats(data_eeg[,"auditory"])$out
#boxplot for sound
boxplot(data_sound)
boxplot.stats(data_sound)$out

#Correlation and scatter plots
#scatter plot prefrontal~sound
plot(data_eeg[,"prefrontal"],data_sound, 
     ylab="Sound signal",xlab="EEG signal (Prefrontal)",
     main="Correlation and scatter plot for EEG signal (Prefrontal) and Sound signal")
#regression line lm(y~x,data=data)
abline(lm(data_sound~data_eeg[,"prefrontal"]),col="red", lwd=2)
cat("Correlation between EEG signal (Prefrontal) and Sound signal",cor(data_eeg[,"prefrontal"],data_sound))

#scatter plot Auditory~sound
plot(data_eeg[,"auditory"],data_sound,
     ylab="Sound signal",xlab="EEG signal (Auditory)",
     main="Correlation and scatter plot for EEG signal (Auditory) and Sound signal")
#regression line lm(y~x,data=data)
abline(lm(data_sound~data_eeg[,"auditory"]),col="blue", lwd=2)
cat("Correlation between EEG signal (Auditory) and Sound signal",cor(data_eeg[,"auditory"],data_sound))



