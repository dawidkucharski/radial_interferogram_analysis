#----------------------new with ImageJ---------------------------------------------------------------------------------------------------------------------------------------------------
library(pracma)
library("car")
library(tikzDevice)
csv_data <- read.csv(file ="/Users/DawidKucharski/Library/Mobile Documents/com~apple~CloudDocs/Dokumenty/Twyman-Green_experiment/Image_processing/#1/Values_corrupted.csv", header = TRUE, sep = ",",stringsAsFactors = TRUE)

localmin <- function(x) {rlex <- rle(x<0)
rlex$values <- seq_along(rlex$values) *rlex$values
ix <- inverse.rle(rlex)
ixt <- tapply(seq_along(ix)[ix>0], ix[ix>0], function(i,xi) i[which.min(xi[i])],xi=x)
unname(ixt)
}
plot(csv_data$Normalized_Integrated_Intensity,type="l")
movingAverage <- function(x, n=1, centered=FALSE) {
  
  if (centered) {
    before <- floor  ((n-1)/2)
    after  <- ceiling((n-1)/2)
  } else {
    before <- n-1
    after  <- 0
  }
  
  # Track the sum and count of number of non-NA items
  s     <- rep(0, length(x))
  count <- rep(0, length(x))
  
  # Add the centered data 
  new <- x
  # Add to count list wherever there isn't a 
  count <- count + !is.na(new)
  # Now replace NA_s with 0_s and add to total
  new[is.na(new)] <- 0
  s <- s + new
  
  # Add the data from before
  i <- 1
  while (i <= before) {
    # This is the vector with offset values to add
    new   <- c(rep(NA, i), x[1:(length(x)-i)])
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # Add the data from after
  i <- 1
  while (i <= after) {
    # This is the vector with offset values to add
    new   <- c(x[(i+1):length(x)], rep(NA, i))
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # return sum divided by count
  s/count
}
i <- 2
while(i <=24000) {                  # Start while-loop
  csv_data[ , i] <- movingAverage(csv_data[ , i], 1, TRUE)
  #plot(csv_data[,i],type="l")
  fit1<-smooth.spline(csv_data[ , i], w = NULL, df=20)
  #pred2<-predict(fit1, deriv = 2)
  #fit2<-smooth.spline(pred2, w = NULL, df=100)
  #plot(pred2,type = 'l')
  #plot(fit1,col='blue',lwd=2,lty=1, type="l", ylab=expression(I), xlab= ("Pixels"))
  #grid()
  
  find_peaks <- function (x, m = 3){
    shape <- diff(sign(diff(x, na.pad = FALSE)))
    pks <- sapply(which(shape < 0), FUN = function(i){
      z <- i - m + 1
      z <- ifelse(z > 0, z, 1)
      w <- i + m + 1
      w <- ifelse(w < length(x), w, length(x))
      if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
    })
    pks <- unlist(pks)
    pks
  }
  peaks<-find_peaks(fit1$y, m = 10)
 # peaks<-peaks[(peaks >= 50) & (peaks < 350)]
  #bad <- c(1, 350)
  #peaks<-subset(peaks, !(peaks %in% bad))
  #print(peaks)
  D <- data.frame(x = seq(1,length(peaks),1), y=c(2*peaks))
  reg = lm(D$y^2~D$x)
  plot(D$y^2~D$x)
  abline(lm(D$y^2~D$x))
  ox<-round(abs(round(coefficients(reg)[1],2))/abs(round(coefficients(reg)[2],2)),2)
  print(ox)
  write.table(ox, file="/Users/DawidKucharski/Library/Mobile Documents/com~apple~CloudDocs/Dokumenty/Twyman-Green_experiment/Image_processing/#1/ox_corrupted.csv", append=TRUE,sep=",",col.names=FALSE,row.names=FALSE)
  
  i <- i + 1
}
ox = read.csv("/Users/DawidKucharski/Library/Mobile Documents/com~apple~CloudDocs/Dokumenty/Twyman-Green_experiment/Image_processing/#1/ox.csv", header = FALSE)
#plot(csv_data$slice_2,type="l")
ox[]<-apply(ox,2,function(x) ifelse(x>1,abs(ox$V1-1),x))
epsilon<-1-ox

ox_corrupted = read.csv("/Users/DawidKucharski/Library/Mobile Documents/com~apple~CloudDocs/Dokumenty/Twyman-Green_experiment/Image_processing/#1/ox_corrupted.csv", header = FALSE)
#plot(csv_data$slice_2,type="l")
ox_corrupted[]<-apply(ox_corrupted,2,function(x) ifelse(x>1,abs(ox_corrupted$V1-1),x))
epsilon_corrupted<-1-ox_corrupted

df1<-data.frame(epsilon$V1,epsilon_corrupted$V1)
df1$Absolute_Distance<-abs(df1$epsilon.V1-df1$epsilon_corrupted.V1)
df1

plot(epsilon$V1, xlim=c(1500,2500),type="b")
abline(h = 0.02, col = "darkgreen") 
abline(h = 0.96, col = "darkgreen") 
# Function for amplitude-unwrapping the signal between 0 and 1
amplitude_unwrap <- function(signal) {
  # Initialize the unwrapped signal with the input signal
  unwrapped_signal <- signal
  # Detect and correct jumps between 0 and 1 in amplitude
  for (i in 2:length(signal)) {
    diff_val <- signal[i] - signal[i - 1]
    if (abs(diff_val) > 0.5) {
      # Unwrap the jump
      unwrapped_signal[i:length(signal)] <- unwrapped_signal[i:length(signal)] - sign(diff_val)
    }}
  return(unwrapped_signal)
}
#Example signal (replace this with your actual signal)
input_signal <- epsilon$V1
# Perform amplitude unwrapping
unwrapped_signal <- amplitude_unwrap(input_signal)
# Plot both wrapped and unwrapped signals in the same plot
plot(c(1:length(input_signal)), input_signal, type = "l", col = "black",ylim = c(0,105),
     xlab = "Sample Index", ylab = "Signal", main = "Wrapped and Unwrapped Signals")
lines(c(1:length(unwrapped_signal)), unwrapped_signal, col = "red")
#legend("topright", legend = c("Wrapped Signal", "Unwrapped Signal"), col = c("blue", "red"), lty = 1)

x<-seq(from = 0, to = (length(epsilon$V1)-1), by = 1)
A<- ((max(unwrapped_signal)-min(unwrapped_signal))/2)
C<-((max(unwrapped_signal)+min(unwrapped_signal))/2)
raw.fft = fft(unwrapped_signal)
truncated.fft = raw.fft[seq(1, length(unwrapped_signal)/2 - 1)]
truncated.fft[1] = 0
omega = which.max(abs(truncated.fft)) * 2 * pi / length(unwrapped_signal)
res<- nls(unwrapped_signal ~ A * cos(omega * x + phi) + C, data = data.frame(x, unwrapped_signal),
          start=list(A = A, omega = 0.0017, phi=pi, C=C), trace = TRUE)
y_est<-predict(res,x)
lines(c(1:length(unwrapped_signal)), y_est, col = "blue")
plot(x=c(1:length(unwrapped_signal)), y=(unwrapped_signal-y_est), type="l",lwd=2, ylab=(expression(P*" ["*mu*"W]")), xlab= (expression(t*" ["*min*"]")),col="green")
abs(max(abs(unwrapped_signal-y_est))-min(abs(unwrapped_signal-y_est)))



#Roughness------------------------------------
# Example profile data (replace with your own)
profile <- abs(unwrapped_signal-y_est)  # Replace with your actual profile data

# Function to calculate the mean of the profile
calculate_mean <- function(profile_data) {
  mean(profile_data)
}

# Function for baseline detrending (remove mean)
detrend_profile <- function(profile_data) {
  mean_value <- calculate_mean(profile_data)
  detrended_profile <- profile_data - mean_value
  return(detrended_profile)
}

# Spline smoothing function
spline_smoothing <- function(profile_data) {
  smoothed_profile <- smooth.spline(seq_along(profile_data), profile_data)
  return(predict(smoothed_profile))
}

# Function to calculate Ra (Average Roughness)
calculate_Ra <- function(profile_data) {
  sum(abs(profile_data)) / length(profile_data)
}

# Function to calculate Rz (Average Maximum Height)
calculate_Rz <- function(profile_data) {
  max_val <- max(profile_data)
  min_val <- min(profile_data)
  return(max_val + abs(min_val))
}

# Detrend the profile
detrended_profile <- detrend_profile(profile)

# Apply spline smoothing to the detrended profile
smoothed_profile <- spline_smoothing(detrended_profile)

# Assuming 'your_data' contains your dataset

# Perform linear regression to obtain the linear trend
linear_model <- lm(profile ~ seq_along(profile))

# Extract the slope and intercept of the linear trend
linear_slope <- coef(linear_model)[2]
linear_intercept <- coef(linear_model)[1]

# Detrend the data by subtracting the linear trend
detrended_linear <- profile - (linear_slope * seq_along(profile) + linear_intercept)

# Perform polynomial regression (change degree as needed)
#poly_degree <- 3
#nonlinear_model <- lm(detrended_linear ~ poly(seq_along(detrended_linear), poly_degree))

# Extract the fitted values (non-linear trend)
#fitted_values <- fitted(nonlinear_model)

# Detrend the data by subtracting the non-linear trend
#detrended_nonlinear <- detrended_linear - fitted_values

# Apply spline smoothing to the detrended profile
smoothed_profile <- spline_smoothing(detrended_linear)

# Calculate Ra and Rz for the smoothed profile
Ra <- calculate_Ra(smoothed_profile$y)
Rz <- calculate_Rz(smoothed_profile$y)

# Plot the raw profile, detrended profile, and smoothed profile with calculated values
plot(266*detrended_linear, type = "l",
     ylab = "", xlab = "Index",cex.main=1.6,cex.axis=1.6,cex.lab=1.6, lwd=1.5,col=grey(.5))
title(ylab=(expression(d*" ["*n*"m]")), mgp=c(2.3,1,0),cex.lab=2)

lines(266*smoothed_profile$y, col="black", lwd=3,lty=1)  # Smoothed profile

# Add text for roughness parameters on the plot
text(10000, max(266*detrended_linear), paste0("Ra:", round(Ra*266,2),"nm"), col = "black", cex = 0.8)
text(10000, max(266*detrended_linear) - 30, paste0("Rz:", round(Rz*266,2),"nm"), col = "black", cex = 0.8)


