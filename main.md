#========================================================
# 0 # START
# -Geopolitics-of-Renewable-Energy-time-varying-interactions-between-geopolitical-risk-and-renewable
#Code for the paper

#========================================================
# 1 # Cubic Spline Interpolation
#========================================================
library(tseries)
library(zoo)

#Converting quarterly renewable energy production for the Netherlands (NL_REP) data into monthly using cubic spline interpolation
#download data
date<- read_excel(file.choose())

#prepare the data
date$Date<- as.Date.character(date$Date)
date$REP<- as.numeric(date$NL_REP)
names(date)<- c("DATE", "VALUE")
head(date) 

#Add the missing months
DateSeq <- seq(date$DATE[1],tail(date$DATE,1),by="1 month")

# interpolate: na.approx linearly interpolates NA values in a time series, na.spline will use cubic spline interpolation.
gerMonthly <- data.frame(DATE=DateSeq, Interp.Value=spline(date, method="natural", xout=DateSeq)$y)
?spline
RE_cons<- merge(date, gerMonthly, by='DATE', all.y = T)

#ploting the converted values 
ts.plot(RE_cons$Interp.Value)

#export data as an excel file
write_xlsx(RE_cons,"/Users/macbook.12/Desktop/RUG 1 year/Research/R code\\REPmonthly.xlsx")

#========================================================
# 2 # Seasonal Adjustment
#========================================================
##R code for Deseasonalizing a series
#Deseasonalizing, then using monthly changes

#Upload the data
data<-read_excel(file.choose())
head(data)
data1<- data$US_REP #Renewables Energy production in the US
data2<- data$NL_REP #Renewables Energy production in the NL

#make a time series Renewable energy production (US_REP & NL_REP)
tsdata1<-ts(data1,start=c(1999,3),frequency = 12)
tsdata2<-ts(data2,start=c(1999,3),frequency = 12)
ts.plot(tsdata1)
ts.plot(tsdata2)

#Now decompose using loess smothing
#note the final product has more than just new time series
tsdata11<-stl(tsdata1,s.window = "periodic") #Renewables US
head(tsdata11$time.series)
plot(tsdata11)

tsdata22<-stl(tsdata2,s.window = "periodic") #Renewables NL
head(tsdata22$time.series)
plot(tsdata22)

#Append deseasonalised series to new
US_REPds<-tsdata1-tsdata11$time.series[,1]
plot(US_REPds)
dataREP_USnew<-cbind(tsdata1,US_REPds)
colnames(dataREP_USnew)[1]<-"US_REP"
head(dataREP_USnew)

NL_REPds<-tsdata1-tsdata22$time.series[,1]
plot(NL_REPds)
dataREP_NLnew<-cbind(tsdata1,NL_REPds)
colnames(dataREP_USnew)[1]<-"NL_REP"
head(dataREP_NLnew)

#plot, separately and together
plot(dataREP_USnew)
plot(dataREP_USnew[,1],xlab="",ylab="",lwd=4,col="dark grey",main="Renewables Energy Production US")
par(new=TRUE)
plot(dataREP_USnew[,2],xlab="",ylab="",lwd=2,col="black",lty=2)
legend("bottomright",lty = c(1,2),col=c("dark grey","black"),legend = c("Not adjusted","Seasonally adjusted"))

plot(dataREP_NLnew)
plot(dataREP_NLnew[,1],xlab="",ylab="",lwd=4,col="dark grey",main="Renewables Energy Production NL")
par(new=TRUE)
plot(dataREP_NLnew[,2],xlab="",ylab="",lwd=2,col="black",lty=2)
legend("bottomright",lty = c(1,2),col=c("dark grey","black"),legend = c("Not adjusted","Seasonally adjusted"))

#export as an excel file
adjusteddata<- cbind(dataREP_USnew, dataREP_NLnew)
adjusteddata<- as.data.frame(adjusteddata)
write_xlsx(adjusteddata,"/Users/macbook.12/Desktop/RUG 1 year/Research/Data\\file.xlsx")

#========================================================
# 3 # Vector Error Correction Model and Cointegration
#========================================================
library(vars)  # vec2var
library(tsDyn) # VECM
library(urca)
library(vars)
library(mFilter)
library(tseries)
library(TSstudio)
library(foREPast)
library(tidyverse)
library(xml2)
library(readxl)
library(writexl)

#Loading the Dataset & checking for stationarity
data <- read_excel(file.choose())
head(data)
str(data)
rep <- ts(data$NL_REPds, start = c(1999,3), frequency = 12) #Renewable energy production in NL, seasonally adjusted 
elc<- ts(data$NL_ELEC, start = c(1999,3), frequency = 12) # Real electrisity prices
cli <- ts(data$NL_CLI, start = c(1999,3), frequency = 12) #Composite Leading Indicator, as a proxy for GDP for NL
gas <- ts(data$EU_NGASPreal, start = c(1999,3), frequency = 12) #Real Crude gas Price for EU, transformed into real price with the use of the Dutch CPI, seasonally adjusted 
risk <- ts(data$RISK, start = c(1999,3), frequency = 12) #Risk Index averaged for 28 biggest natural gas exporting countries for the years 1999-2020

pp.test(rep) #p-value = 0.01
pp.test(cli) #p-value = 0.2267
pp.test(elc) #p-value = 0.4863
pp.test(gas) #p-value = 0.4028
pp.test(risk) #p-value = 0.051
#Results:CLI, ELC, GAS and RISK are non-stationary

#Cheking wheather we can still run our VAR model in levels given the non-stationarity 
# (resoning: when variables are non stationary, a VAR model in levels is not appropriate since it is a spurious regression which is a non-interpretable regression. However, although variables are non stationary but when cointegrations exist, a VAR model in levels can be estimated which has a long-term interpretation.)

# level data 
data(data)
lev <- cbind(data$NL_REPds, data$EU_NGASPreal, data$NL_ELEC, data$RISK , data$NL_CLI)
nr_lev <- nrow(lev)

# the sample period
yq <- expand.grid(1:12, 1999:2020)[1:nr_lev,]
yq<- yq[ -c(1,2), ]
colnames(yq) <- c("m", "yyyy")
rownames(yq) <- NULL
lev<- as.data.frame(lev)

# Cointegration Test
#if r = 0, estimate VAR in differences → period.
#if r > 0, apply the cajorls() function to the ca.jo() output to get estimated parameters of the VECM model in differences

coint_ca.jo <- ca.jo(lev, ecdet = "none", type  = "eigen", K = 7, spec = "longrun", season = 12, dumvar = NULL)
summary(coint_ca.jo) # r=0, as test statistics for r = 0 is 25.99 < 38.78(10pct), 33.32(5pct), 30.84(1pct)
coint_ca.jo
?ca.jo
#Conclusion: In the case of no cointegration (r=0), since all variables are non-stationary in level, the above VECM model reduces to a VAR model with growth variables. We astimate VAR in diffrences.  
#========================================================
# 4 # The SVAR model 
#========================================================
#Loading the Dataset
data <- read_excel(file.choose())
head(data)
str(data)

#========================================================
#Creating the Time Series Objectives for NL
rep <- ts(log(data$NL_REPds), start = c(1999,3), frequency = 12) #Renewable energy production in NL, seasonally adjusted 
elc<- ts(log(data$NL_ELEC), start = c(1999,3), frequency = 12) #Electricity Price 
cli <- ts(log(data$NL_CLI), start = c(1999,3), frequency = 12) #Composite Leading Indicator, as a proxy for GDP for NL
gas <- ts(log(data$EU_NGASPreal), start = c(1999,3), frequency = 12) #Real Crude gas Price for EU, transformed into real price with the use of the Dutch CPI, seasonally adjusted 
risk <- ts(log(data$RISK), start = c(1999,3), frequency = 12) #Risk Index averaged for 28 biggest natural gas exporting countries for the years 1999-2020

svardata1<- cbind(risk, gas, elc, cli, rep)
colnames(svardata1) <- cbind("Risk", "Natural gas", "Electricity", "Cli", "Renewables")
head(svardata1)
plot(svardata1, main = "Time series for the SVAR Model (level)")

#========================================================
#Checking for stationarity
#Phillips-Perron Unit Root Test of level variables
pp.test(rep) #p-value = 0.01
pp.test(cli) #p-value = 0.2267
pp.test(gas) #p-value = 0.4028
pp.test(risk) #p-value = 0.051
pp.test(elc) #p-value = 0.4863
#Results:Cli, ELC, GAS and RISK are non-stationary

#The Augmented Dickey-Fuller test of level variables
lt.adf.lv.none<- list(
  REP = ur.df(rep, type = "none", selectlags = c("Fixed")),
  Cli = ur.df(cli, type = "none", selectlags = c("Fixed")),
  GAS = ur.df(gas, type = "none", selectlags = c("Fixed")),
  ELC = ur.df(elc, type = "none", selectlags = c("Fixed")),
  RISK = ur.df(risk, type = "none", selectlags = c("Fixed")))

print(lt.adf.lv.none)
summary(lt.adf.lv.none$REP) #t statistics = 0.497, p-value = 0.62
summary(lt.adf.lv.none$Cli) #t statistics = -0.470, p-value = 0.639
summary(lt.adf.lv.none$GAS) #t statistics = -0.6443, p-value = 0.52  
summary(lt.adf.lv.none$ELC) #t statistics = 0.499, p-value = 0.618  
summary(lt.adf.lv.none$RISK) #t statistics = 0.499, p-value = 0.618
#Results: non-stationary

lt.adf.lv.drift<- list(
  REP = ur.df(REP, type = "drift", selectlags = c("BIC")),
  Cli = ur.df(Cli, type = "drift", selectlags = c("BIC")),
  GAS = ur.df(gas, type = "drift", selectlags = c("BIC")),
  RISK = ur.df(risk, type = "drift", selectlags = c("BIC")))
print(lt.adf.lv.drift)

lt.adf.lv.trend<- list(
  REP = ur.df(REP, type = "trend", selectlags = c("BIC")),
  Cli = ur.df(Cli, type = "trend", selectlags = c("BIC")),
  GAS = ur.df(gas, type = "trend", selectlags = c("BIC")),
  RISK = ur.df(risk, type = "trend", selectlags = c("BIC")))
print(lt.adf.lv.trend)

#Transforming variables into the first differenced term 
rep<- diff(rep, 12)
cli<- diff(cli, 12)
gas<- diff(gas, 12)
risk<- diff(risk, 12)
elc<- diff(elc, 12)
svardata2<- cbind(risk, gas, elc, cli, rep)
colnames(svardata2) <- cbind("Risk", "Gas", "Electricity", "Cli", "Renewables")
head(svardata2)

rep<- diff(rep)
cli<- diff(cli)
gas<- diff(gas)
risk<- diff(risk)
elc<- diff(elc)
svardata2<- cbind(risk, gas, elc, cli, rep)
colnames(svardata2) <- cbind("Risk", "Gas", "Electricity", "Cli", "Renewables")
head(svardata2)

#Checking for stationarity
#Phillips-Perron Unit Root Test of variables in the first differenced term 
pp.test(rep) #p-value = 0.01
pp.test(cli) #p-value = 0.0201 
pp.test(gas) #p-value = 0.0153
pp.test(risk) #p-value = 0.01
pp.test(elc) #p-value = 0.000
#Results: stationary

#The Augmented Dickey-Fuller test of variables in the first differenced term 
lt.adf.lv.none<- list(
  REP = ur.df(rep, type = "none", selectlags = c("Fixed")),
  Cli = ur.df(cli, type = "none", selectlags = c("Fixed")),
  GAS = ur.df(gas, type = "none", selectlags = c("Fixed")),
  ELC = ur.df(elc, type = "none", selectlags = c("Fixed")),
  RISK = ur.df(risk, type = "none", selectlags = c("Fixed")))

print(lt.adf.lv.none)
summary(lt.adf.lv.none$REP) #t statistics = -5.550, p-value = 7.39e-08  
summary(lt.adf.lv.none$Cli) #t statistics = -3.719, p-value = 0.000248  
summary(lt.adf.lv.none$GAS) #t statistics = -3.485, p-value = 0.000582
summary(lt.adf.lv.none$RISK) #t statistics = -3.926, p-value = 0.000112
summary(lt.adf.lv.none$ELC) #t statistics = -3.926, p-value = 0.000112
#Results: stationary

lt.adf.lv.drift<- list(
  REP = ur.df(REP, type = "drift", selectlags = c("BIC")),
  Cli = ur.df(Cli, type = "drift", selectlags = c("BIC")),
  GAS = ur.df(gas, type = "drift", selectlags = c("BIC")),
  RISK = ur.df(risk, type = "drift", selectlags = c("BIC")))
print(lt.adf.lv.drift)

lt.adf.lv.trend<- list(
  REP = ur.df(REP, type = "trend", selectlags = c("BIC")),
  Cli = ur.df(Cli, type = "trend", selectlags = c("BIC")),
  GAS = ur.df(gas, type = "trend", selectlags = c("BIC")),
  RISK = ur.df(risk, type = "trend", selectlags = c("BIC")))
print(lt.adf.lv.trend)

#========================================================
#Correlation Matrix
library(Hmisc)
res2 <- rcorr(as.matrix(svardata2))
res2

library(corrplot)
# Insignificant correlation are crossed
corrplot(res2$r, type="upper", order="hclust", 
         p.mat = res2$P, sig.level = 0.01, insig = "blank")

#========================================================
#Time Series Plots
plot<- plot.ts(svardata1, main = "NL - Level variables used in the analysis")
plot<- plot.ts(svardata2, main = "Time series for the SVAR Model (log difference)")
ts_plot(rep, title = "Renewable Energy Production", Xtitle = "Time", Ytitle = "Renewables Consumption")
ts_plot(Cli, title = "Composite Leading Indicator", Xtitle = "Time", Ytitle = "Cli")
ts_plot(gas, title = "Real Crude gas Price", Xtitle = "Time", Ytitle = "Real gas price")
ts_plot(risk, title = "Risk Index", Xtitle = "Time", Ytitle = "Risk")

#Descriptive statistic.
summary(svardata2)
library(pastecs)
res <- stat.desc(svardata2)
round(res, 2)
?stat.desc

# calculate skewness & kurtosis in r
library(moments)
skewness(svardata2)
kurtosis(svardata2)

#========================================================
#Buidling the VAR Model & selection of lag
head(svardata2)
info.var <- VARselect(svardata2, lag.max = 30, type = "trend", season = 12)
info.var$selection
info.var$criteria
?VAR
#========================================================
#Estimate the reduced-form VAR
var.est1 <- VAR(svardata2, p = 9, type = "both", season = 12)
summary(var.est1)

new<- (BQ(var.est1))
plot(irf(var.est1, n.ahead = 36, boot = TRUE, cumulative = TRUE))
plot(irf(new, n.ahead = 36, cumulative = TRUE))
#========================================================
#Determine the persistence of the model
# autocorrelation and  partial autocorrelation functions
par(mfrow=c(2,2))
acfrep<- acf(rep, main = "ACF for Renewable Energy Production")
pacfrep<- pacf(rep, main = "PACF for Renewable Energy Production")
acfcli<- acf(cli, main = "ACF for Composite Leading Indicator")
pacfcli<- pacf(cli, main = "PACF for Composite Leading Indicator")

par(mfrow=c(2,2))
acfgas<- acf(gas, main = "ACF for Real Natural Gas Price")
pacfgas<- pacf(gas, main = "PACF fpr Real Natural Gas Price")
acfrisk<- acf(risk, main = "ACF for Political Risk Index")
pacfrisk<- pacf(risk, main = "PACF for Political Risk Index")

par(mfrow=c(1,2))
acfelc<- acf(elc, main = "ACF for Electricity Price")
pacfelc<- pacf(elc, main = "PACF for Electricity Price")


#Serial Correlation 
Serial1<- serial.test(var.est1, lags.pt = 24, type = "PT.asymptotic")
Serial1 #p-value =  > 0,05, meaning that there is no serial correlation being observed, 
#the model passes the test

#Heteroskedasticity
Arch1<- arch.test(var.est1, lags.multi = 12, multivariate.only = TRUE)
Arch1 #p-value = 0,000 < 0,05, meaning that the model suffers from heteroskedastisity

#Normal distribution of the residuals 
Norm1<- normality.test(var.est1, multivariate.only = TRUE)
Norm1 #the model passed the skewness test, JB-Test and Kurtosis tests

#Testing for structural breaks in the Residuals
Stability1<- stability(var.est1, type = "OLS-CUSUM") #the system is stable
plot(Stability1)
#stable 

#========================================================
#Granger Causality
download.packages("grangers")
library(grangers)

GrangerREP<- causality(var.est1, cause = "Renewables")
print(GrangerREP)
#p-value = 0,000 < 0,05, so we can reject Granger causality H0: Renewables do not Granger-cause Cli Gas Risk; 
#p-value = 0,009 < 0,05, so we can reject H0: No instantaneous causality between: Renewables and Risk Gas Cli
GrangerCli<- causality(var.est1, cause = "Cli")
GrangerCli 
#p-value = 0,000 < 0,05, so we can reject Granger causality H0: Cli do not Granger-cause Renewables Gas Risk
#p-value = 0,009 < 0,05, so we can reject H0: No instantaneous causality between: Cli and Risk Gas Renewables
GrangerGas<- causality(var.est1, cause = "Gas")
GrangerGas 
#p-value = 0,059 > 0,05, so we cannot reject Granger causality H0: Gas do not Granger-cause Renewables Cli Risk
#p-value = 0,984 > 0,05, so we cannot reject H0: No instantaneous causality between: Gas and Risk Cli Renewables
GrangerRISK<- causality(var.est1, cause = "Risk")
GrangerRISK 
#p-value = 0,842 > 0,05, so we cannot reject Granger causality H0: Risk do not Granger-cause Gas Cli Renewables
#p-value = 0,901 > 0,05, so we cannot reject H0: No instantaneous causality between: Risk and Gas Cli Renewables
GrangerELEC<- causality(var.est1, cause = "Electrisity")
GrangerELEC
#we can still check the Impulse Response Functions 
#In principle, the two estimations methodology has different estimations motive. The Granger Causality test is use to verify causal diREPtion between two variables only. 
#For instance, Whether X must leads for Y to follow and vice versa. 
#While Impulse Response Function (VIF) test the decomposition of shocks impact on a variable at a given period. 
#That is, what is the effect X shock on Y? Then, VIF is appropriate methodology to use.

grangertest(risk, rep, order = 48)
grangertest(risk, gas, order = 48)
grangertest(risk, elc, order = 48)
grangertest(risk, cli, order = 48)

grangertest(gas, rep, order = 48)
grangertest(gas, risk, order = 48)
grangertest(gas, elc, order = 48)
grangertest(gas, cli, order = 48)

grangertest(elc, rep, order = 48)
grangertest(elc, risk, order = 48)
grangertest(elc, gas, order = 48)
grangertest(elc, cli, order = 48)

grangertest(rep, elc, order = 48)
grangertest(rep, risk, order = 48)
grangertest(rep, gas, order = 48)
grangertest(rep, cli, order = 48)

grangertest(cli, elc, order = 48)
grangertest(cli, risk, order = 48)
grangertest(cli, gas, order = 6)
grangertest(cli, rep, order = 48)

#========================================================
#Variance Decomposition
FEVDvar1<- fevd(var.est1, n.ahead = 20)
plot(FEVDvar1)

#========================================================
#VAR FoREPasting 
ForcastVAR1<- predict(var.est1, n.ahead = 6, ci = 0,95)
fanchart(ForcastVAR1, names = "Renewables")
fanchart(ForcastVAR1, names = "Cli")
fanchart(ForcastVAR1, names = "Gas")
fanchart(ForcastVAR1, names = "Risk")

#========================================================
#========================================================
#SVAR Model
#Setting the Restrictions for the SVAR model
#Contemporaneous coefficients
summary(var.est1)

a.mat <- diag(5)
a.mat[2, 1] <- NA
a.mat[3, 1] <- NA
a.mat[4, 1] <- NA
a.mat[5, 1] <- NA
a.mat[3, 2] <- NA
a.mat[4, 2] <- NA
a.mat[5, 2] <- NA
a.mat[4, 3] <- NA
a.mat[5, 3] <- NA
a.mat[5, 4] <- NA

print(a.mat)

b.mat <- diag(5)
print(b.mat)
#========================================================
#Buidling the Model
svar.one <- SVAR(var.est1, estmethod = "direct", Amat = a.mat)
svar.one
summary (svar.one)

#========================================================
#Impulse response functions Plot
library(devtools)

#create function to extract the irf to plot them with ggplot later
extract_varirf <- function(...){
  
  varirf_object <- list(...) #list one or more varirf input objects
  
  get_vec_length <- function(list_item){nrow(list_item[[1]][[1]])}
  
  if (!("varirf" %in% mapply(class, varirf_object))){
    stop("this function only accepts 'varirf' class objects")
  }
  
  if (length(unique(mapply(class, varirf_object)))!=1){
    stop("all input items must be 'varirf' class objects")
  }    
  if (length(unique(mapply(get_vec_length, varirf_object)))!=1){
    stop("all irf vectors must have the same length")   
  }  
  
  period <- as.data.frame(0:(nrow(varirf_object[[1]][[1]][[1]])-1)) 
  names(period) <- "period"
  
  for (l in 1:length(varirf_object)){
    for (i in 1:3){
      for (j in 1:dim(varirf_object[[l]][[i]][[1]])[2]){
        for (k in 1:length(varirf_object[[l]][[1]])){
          temp_colname <- paste(names(varirf_object[[l]][i]), #vector type (irf, lower, or upper)
                                names(varirf_object[[l]][[i]])[k], #impulse name
                                colnames(varirf_object[[l]][[i]][[k]])[j], #response name
                                sep = "_")
          
          temp <- as.data.frame(varirf_object[[l]][[i]][[k]][, j]) #extracts the vector
          
          names(temp) <- temp_colname #add the column name (vectortype_impulse_reponse)
          period <- cbind(period, temp) 
        }
        
      }
    }
  }
  names(period) <- tolower(names(period))
  return(period)
}
irf_all<- irf(svar.one, n.ahead = 24, ortho = TRUE, boot = TRUE, cumulative = TRUE, ci = 0.9)

#extract irf 
multiple_varirf <- extract_varirf(irf_all)
head(multiple_varirf)

#plot irf with ggplot 
names(multiple_varirf)

#Risk
gas_risk <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_gas_risk, ymin=lower_gas_risk, ymax=upper_gas_risk)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("real natural gas - political risk index")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
gas_risk

cli_risk <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_cli_risk, ymin=lower_cli_risk, ymax=upper_cli_risk)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("cli - political risk index")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
cli_risk

electricity_risk <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_electricity_risk, ymin=lower_electricity_risk, ymax=upper_electricity_risk)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("electrisity price - political risk index")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
electricity_risk

renewables_risk <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_renewables_risk, ymin=lower_renewables_risk, ymax=upper_renewables_risk)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("renewable energy gerenarion - political risk index")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
renewables_risk

risk_risk <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_risk_risk, ymin=lower_risk_risk, ymax=upper_risk_risk)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("political risk index- political risk index")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
risk_risk

library(patchwork)
risk_irf <-  gas_risk / electricity_risk / renewables_risk / risk_risk / cli_risk
risk_irf

#Natural gas price
gas_gas <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_gas_gas, ymin=lower_gas_gas, ymax=upper_gas_gas)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("real natural gas - real natural gas")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
gas_gas

cli_gas <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_cli_gas, ymin=lower_cli_gas, ymax=upper_cli_gas)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("cli - real natural gas")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
cli_gas

electricity_gas <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_electricity_gas, ymin=lower_electricity_gas, ymax=upper_electricity_gas)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("electrisity price - real natural gas")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
electricity_gas

renewables_gas <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_renewables_gas, ymin=lower_renewables_gas, ymax=upper_renewables_gas)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("renewable energy gerenarion - real natural gas")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
renewables_gas

risk_gas <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_risk_gas, ymin=lower_risk_gas, ymax=upper_risk_gas)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("political risk index- real natural gas")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
risk_gas

gas_irf <-  gas_gas / electricity_gas / renewables_gas / risk_gas / cli_gas
gas_irf

#Electricity price
gas_electricity <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_gas_electricity, ymin=lower_gas_electricity, ymax=upper_gas_electricity)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("real natural gas - electrisity price")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
electricity_gas

cli_electricity <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_cli_electricity, ymin=lower_cli_electricity, ymax=upper_cli_electricity)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("cli - electrisity price")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
cli_electricity

electricity_electricity <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_electricity_electricity, ymin=lower_electricity_electricity, ymax=upper_electricity_electricity)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("electrisity price - electrisity price")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
electricity_electricity

renewables_electricity <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_renewables_electricity, ymin=lower_renewables_electricity, ymax=upper_renewables_electricity)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("renewable energy gerenarion - electrisity price")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
renewables_electricity

risk_electricity <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_risk_electricity, ymin=lower_risk_electricity, ymax=upper_risk_electricity)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("political risk index- electrisity price")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
risk_electricity

electricity_irf <-  gas_electricity / electricity_electricity / renewables_electricity / risk_electricity / cli_electricity
electricity_irf

#Renewable energy generation 
gas_renewables <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_gas_renewables, ymin=lower_gas_renewables, ymax=upper_gas_renewables)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("real natural gas - renewable energy gerenarion")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
renewables_gas

cli_renewables <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_cli_renewables, ymin=lower_cli_renewables, ymax=upper_cli_renewables)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("cli - renewable energy gerenarion")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
cli_renewables

electricity_renewables <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_electricity_renewables, ymin=lower_electricity_renewables, ymax=upper_electricity_renewables)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("electrisity price - renewable energy gerenarion")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
electricity_renewables

renewables_renewables <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_renewables_renewables, ymin=lower_renewables_renewables, ymax=upper_renewables_renewables)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("renewable energy gerenarion - renewable energy gerenarion")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
renewables_renewables

risk_renewables <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_risk_renewables, ymin=lower_risk_renewables, ymax=upper_risk_renewables)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("political risk index - renewable energy gerenarion")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
risk_renewables

renewables_irf <-  gas_renewables / electricity_renewables / renewables_renewables / risk_renewables / cli_renewables
renewables_irf

#CLI
gas_cli <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_gas_cli, ymin=lower_gas_cli, ymax=upper_gas_cli)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("real natural gas - cli")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
cli_gas

cli_cli <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_cli_cli, ymin=lower_cli_cli, ymax=upper_cli_cli)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("cli - cli")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
cli_cli

electricity_cli <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_electricity_cli, ymin=lower_electricity_cli, ymax=upper_electricity_cli)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("electrisity price - cli")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
electricity_cli

renewables_cli <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_renewables_cli, ymin=lower_renewables_cli, ymax=upper_renewables_cli)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("renewable energy gerenarion - cli")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
renewables_cli

risk_cli <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_risk_cli, ymin=lower_risk_cli, ymax=upper_risk_cli)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("political risk index - cli")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
risk_cli

cli_irf <-  gas_cli / electricity_cli / renewables_cli / risk_cli / cli_cli
cli_irf

#plot all irf together
gas_irf
electricity_irf
renewables_irf
risk_irf
cli_irf

#========================================================
#SVAR Variance Decomposition
FEVDsvar1<- fevd(svar.one, n.ahead = 60)
FEVDsvar1
plot(FEVDsvar1)


#========================================================
#Risk of Only the direct importers for the Dutch natural gas market (Norway, Russia, United Kingdom)
#Additional SVAR model for regional analysis of the risk factor
rep <- ts(log(data$NL_REPds), start = c(1999,3), frequency = 12) #Renewable energy production in NL, seasonally adjusted 
cli <- ts(log(data$NL_CLI), start = c(1999,3), frequency = 12) #Composite Leading Indicator, as a proxy for GDP for NL
elc<- ts(log(data$NL_ELEC), start = c(1999,3), frequency = 12) #Electrisity price
gas <- ts(log(data$EU_NGASPnom), start = c(1999,3), frequency = 12) #Real Crude gas Price for EU, transformed into real price with the use of the Dutch CPI, seasonally adjusted 
risk <- ts(log(data$RISK_NLimporters), start = c(1999,3), frequency = 12) #Risk Index for only the direct importers for the Dutch natural gas market (Norway, Russia, United Kingdom) 

svardata11<- cbind(risk, gas, elc, cli, rep)
colnames(svardata11) <- cbind("Risk", "Natural gas", "Electricity", "Cli", "Renewables")
head(svardata11)

#========================================================
#Checking for stationarity
#Phillips-Perron Unit Root Test of level variables
pp.test(rep) #p-value = 0.01
pp.test(cli) #p-value = 0.2267
pp.test(gas) #p-value = 0.4028
pp.test(risk) #p-value = 0.051
pp.test(elc) #p-value = 0.4863
#Results:Cli, ELC, GAS and RISK are non-stationary

#The Augmented Dickey-Fuller test of level variables
lt.adf.lv.none<- list(
  REP = ur.df(rep, type = "none", selectlags = c("Fixed")),
  Cli = ur.df(cli, type = "none", selectlags = c("Fixed")),
  GAS = ur.df(gas, type = "none", selectlags = c("Fixed")),
  ELC = ur.df(elc, type = "none", selectlags = c("Fixed")),
  RISK = ur.df(risk, type = "none", selectlags = c("Fixed")))

print(lt.adf.lv.none)
summary(lt.adf.lv.none$REP) #t statistics = 0.497, p-value = 0.62
summary(lt.adf.lv.none$Cli) #t statistics = -0.470, p-value = 0.639
summary(lt.adf.lv.none$GAS) #t statistics = -0.6443, p-value = 0.52  
summary(lt.adf.lv.none$ELC) #t statistics = 0.499, p-value = 0.618  
summary(lt.adf.lv.none$RISK) #t statistics = 0.499, p-value = 0.618
#Results: non-stationary

#Transforming variables into the first differenced term 
rep<- diff(rep)
cli<- diff(cli)
gas<- diff(gas)
risk<- diff(risk)
elc<- diff(elc)
svardata22<- cbind(risk, gas, elc, cli, rep)
colnames(svardata22) <- cbind("Risk", "Gas", "Electricity", "Cli", "Renewables")
head(svardata22)

#Checking for stationarity
#Phillips-Perron Unit Root Test of variables in the first differenced term 
pp.test(rep) #p-value = 0.01
pp.test(cli) #p-value = 0.0201 
pp.test(gas) #p-value = 0.0153
pp.test(risk) #p-value = 0.01
pp.test(elc) #p-value = 0.2065!!!
#Results: stationary

#The Augmented Dickey-Fuller test of variables in the first differenced term 
lt.adf.lv.none<- list(
  REP = ur.df(rep, type = "none", selectlags = c("Fixed")),
  Cli = ur.df(cli, type = "none", selectlags = c("Fixed")),
  GAS = ur.df(gas, type = "none", selectlags = c("Fixed")),
  ELC = ur.df(elc, type = "none", selectlags = c("Fixed")),
  RISK = ur.df(risk, type = "none", selectlags = c("Fixed")))

print(lt.adf.lv.none)
summary(lt.adf.lv.none$REP) #t statistics = -5.550, p-value = 7.39e-08  
summary(lt.adf.lv.none$Cli) #t statistics = -3.719, p-value = 0.000248  
summary(lt.adf.lv.none$GAS) #t statistics = -3.485, p-value = 0.000582
summary(lt.adf.lv.none$RISK) #t statistics = -3.926, p-value = 0.000112
summary(lt.adf.lv.none$ELC) #t statistics = -3.926, p-value = 0.000112
#Results: stationary
#========================================================
#Time Series Plots
plot<- plot.ts(risk, main = "Political risk index for Norway, Russia & United Kingdom")

#========================================================
#Buidling the VAR Model & selection of lag
head(svardata22)
info.var <- VARselect(svardata22, lag.max = 30, type = "both")
info.var$selection
info.var$criteria

#========================================================
#Estimate the reduced-form VAR
var.est2 <- VAR(svardata22, p = 9, type = "both", season = 12)
summary(var.est2)

plot(irf(var.est2, n.ahead = 24, boot = TRUE, cumulative = TRUE, ci = 0.9))

#========================================================
#Determine the persistence of the model
acf(rec, main = "ACF for Renewable Energy Consumption")
pacf(rec, main = "PACF for Renewable Energy Consumption")

acf(cli, main = "ACF for Composite Leading Indicator")
pacf(cli, main = "PACF for Composite Leading Indicator")

acf(gas, main = "ACF for Real Crude gas Price")
pacf(gas, main = "PACF for Real Crude gas Price")

acf(risk.mena, main = "Risk Index for MENA region")
pacf(risk.mena, main = "Risk Index for MENA region")

acf(risk.eca, main = "Risk Index for ECA region")
pacf(risk.eca, main = "Risk Index for ECA region")

acf(risk.amer, main = "Risk Index for AMER region")
pacf(risk.amer, main = "Risk Index for AMER region")

acf(risk.other, main = "Risk Index for Other region")
pacf(risk.other, main = "Risk Index for Other region")

#Serial Correlation 
Serial2<- serial.test(var.est2, lags.pt = 12, type = "PT.asymptotic")
Serial2 #p-value = 0,008 < 0,05, meaning that some serial correlation is observed

#Heteroskedasticity
Arch2<- arch.test(var.est2, lags.multi = 12, multivariate.only = TRUE)
Arch2 #p-value = 1 > 0,05, meaning that the model does not suffer from heteroskedastisity
#the model passes the test

#Normal distribution of the residuals 
Norm2<- normality.test(var.est2, multivariate.only = TRUE)
Norm2 #the model passed the skewness test, JB-Test and Kurtosis tests

#Testing for structural breaks in the Residuals
Stability2<- stability(var.est2, type = "OLS-CUSUM") #the system is stable
plot(Stability2)
#stable 

#========================================================
#Granger Causality
GrangerRISK<- causality(var.est2, cause = "Risk")
GrangerRISK
#p-value = 0,000 < 0,05, so we can reject Granger causality H0: Renewables do not Granger-cause CLI Gas Risk; 
#p-value = 0,006 < 0,05, so we can reject H0: No instantaneous causality between: Renewables and Risk Gas Cli

#select order = 1, 3, 6, 12, 24, 48
grangertest(risk, rep, order = 48)
grangertest(risk, gas, order = 48)
grangertest(risk, elc, order = 48)
grangertest(risk, cli, order = 48)

#========================================================
#Impulse Response Functions 
#Risk
rent<- irf(var.est2, n.ahead = 36, boot = TRUE, cumulative = TRUE)
plot(rent)

#========================================================
#========================================================

#SVAR Model
#Setting the Restrictions for the SVAR model
#Contemporaneous coefficients
summary(var.est2)

a.mat <- diag(5)
a.mat[2, 1] <- NA
a.mat[3, 1] <- NA
a.mat[4, 1] <- NA
a.mat[5, 1] <- NA
a.mat[3, 2] <- NA
a.mat[4, 2] <- NA
a.mat[5, 2] <- NA
a.mat[4, 3] <- NA
a.mat[5, 3] <- NA
a.mat[5, 4] <- NA

print(a.mat)

#========================================================
#Buidling the Model
svar.two <- SVAR(var.est2, estmethod = "direct", Amat = a.mat)
svar.two
summary (svar.two)

#========================================================
#Impulse response functions
irf_all<- irf(svar.two, n.ahead = 24, ortho = TRUE, boot = TRUE, cumulative = TRUE, ci = 0.9)

#extract irf 
multiple_varirf <- extract_varirf(irf_all)
head(multiple_varirf)

#plot irf with ggplot 
names(multiple_varirf)

#Risk
risk_risk <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_risk_risk, ymin=lower_risk_risk, ymax=upper_risk_risk)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("political risk index- political risk index")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
risk_risk

risk_gas <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_risk_gas, ymin=lower_risk_gas, ymax=upper_risk_gas)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("political risk index- real natural gas")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
risk_gas

risk_electricity <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_risk_electricity, ymin=lower_risk_electricity, ymax=upper_risk_electricity)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("political risk index- electrisity price")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
risk_electricity

risk_renewables <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_risk_renewables, ymin=lower_risk_renewables, ymax=upper_risk_renewables)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("political risk index - renewable energy gerenarion")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
risk_renewables

risk_cli <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_risk_cli, ymin=lower_risk_cli, ymax=upper_risk_cli)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("political risk index - cli")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
risk_cli

library(patchwork)
risk_irf <-  risk_risk / risk_gas / risk_electricity / risk_renewables / risk_cli
risk_irf

#========================================================
#SVAR Variance Decomposition
FEVDsvar2<- fevd(var.est2, n.ahead = 38)
FEVDsvar2
plot(FEVDsvar2)


#========================================================
#Risk index for only Russia 
#Additional SVAR model for regional analysis of the risk factor
rep <- ts(log(data$NL_REPds), start = c(1999,3), frequency = 12) #Renewable energy production in NL, seasonally adjusted 
cli <- ts(log(data$NL_CLI), start = c(1999,3), frequency = 12) #Composite Leading Indicator, as a proxy for GDP for NL
elc<- ts(log(data$NL_ELEC), start = c(1999,3), frequency = 12) #Electrisity price
gas <- ts(log(data$EU_NGASPnom), start = c(1999,3), frequency = 12) #Real Crude gas Price for EU, transformed into real price with the use of the Dutch CPI, seasonally adjusted 
risk <- ts(log(data$Russia), start = c(1999,3), frequency = 12) #Risk Index for only the direct importers for the Dutch natural gas market (Norway, Russia, United Kingdom) 

svardata11<- cbind(risk, gas, elc, cli, rep)
colnames(svardata11) <- cbind("Risk", "Natural gas", "Electricity", "Cli", "Renewables")
head(svardata11)

#========================================================
#Checking for stationarity
#Phillips-Perron Unit Root Test of level variables
pp.test(rep) #p-value = 0.01
pp.test(cli) #p-value = 0.2267
pp.test(gas) #p-value = 0.4028
pp.test(risk) #p-value = 0.051
pp.test(elc) #p-value = 0.4863
#Results:Cli, ELC, GAS and RISK are non-stationary

#The Augmented Dickey-Fuller test of level variables
lt.adf.lv.none<- list(
  REP = ur.df(rep, type = "none", selectlags = c("Fixed")),
  Cli = ur.df(cli, type = "none", selectlags = c("Fixed")),
  GAS = ur.df(gas, type = "none", selectlags = c("Fixed")),
  ELC = ur.df(elc, type = "none", selectlags = c("Fixed")),
  RISK = ur.df(risk, type = "none", selectlags = c("Fixed")))

print(lt.adf.lv.none)
summary(lt.adf.lv.none$REP) #t statistics = 0.497, p-value = 0.62
summary(lt.adf.lv.none$Cli) #t statistics = -0.470, p-value = 0.639
summary(lt.adf.lv.none$GAS) #t statistics = -0.6443, p-value = 0.52  
summary(lt.adf.lv.none$ELC) #t statistics = 0.499, p-value = 0.618  
summary(lt.adf.lv.none$RISK) #t statistics = 0.499, p-value = 0.618
#Results: non-stationary

#Transforming variables into the first differenced term 
rep<- diff(rep)
cli<- diff(cli)
gas<- diff(gas)
risk<- diff(risk)
elc<- diff(elc)
svardata22<- cbind(risk, gas, elc, cli, rep)
colnames(svardata22) <- cbind("Risk", "Gas", "Electricity", "Cli", "Renewables")
head(svardata22)

#Checking for stationarity
#Phillips-Perron Unit Root Test of variables in the first differenced term 
pp.test(rep) #p-value = 0.01
pp.test(cli) #p-value = 0.0201 
pp.test(gas) #p-value = 0.0153
pp.test(risk) #p-value = 0.01
pp.test(elc) #p-value = 0.2065!!!
#Results: stationary

#The Augmented Dickey-Fuller test of variables in the first differenced term 
lt.adf.lv.none<- list(
  REP = ur.df(rep, type = "none", selectlags = c("Fixed")),
  Cli = ur.df(cli, type = "none", selectlags = c("Fixed")),
  GAS = ur.df(gas, type = "none", selectlags = c("Fixed")),
  ELC = ur.df(elc, type = "none", selectlags = c("Fixed")),
  RISK = ur.df(risk, type = "none", selectlags = c("Fixed")))

print(lt.adf.lv.none)
summary(lt.adf.lv.none$REP) #t statistics = -5.550, p-value = 7.39e-08  
summary(lt.adf.lv.none$Cli) #t statistics = -3.719, p-value = 0.000248  
summary(lt.adf.lv.none$GAS) #t statistics = -3.485, p-value = 0.000582
summary(lt.adf.lv.none$RISK) #t statistics = -3.926, p-value = 0.000112
summary(lt.adf.lv.none$ELC) #t statistics = -3.926, p-value = 0.000112
#Results: stationary
#========================================================
#Time Series Plots
plot<- plot.ts(risk, main = "Politica risk index for Russia")

#========================================================
#Buidling the VAR Model & selection of lag
head(svardata22)
info.var <- VARselect(svardata22, lag.max = 30, type = "both")
info.var$selection
info.var$criteria

#========================================================
#Estimate the reduced-form VAR
var.est2 <- VAR(svardata22, p = 9, type = "both", season = 12)
summary(var.est2)

plot(irf(var.est2, n.ahead = 24, boot = TRUE, cumulative = TRUE))

#========================================================
#Determine the persistence of the model
acf(rec, main = "ACF for Renewable Energy Consumption")
pacf(rec, main = "PACF for Renewable Energy Consumption")

acf(cli, main = "ACF for Composite Leading Indicator")
pacf(cli, main = "PACF for Composite Leading Indicator")

acf(gas, main = "ACF for Real Crude gas Price")
pacf(gas, main = "PACF for Real Crude gas Price")

acf(risk.mena, main = "Risk Index for MENA region")
pacf(risk.mena, main = "Risk Index for MENA region")

acf(risk.eca, main = "Risk Index for ECA region")
pacf(risk.eca, main = "Risk Index for ECA region")

acf(risk.amer, main = "Risk Index for AMER region")
pacf(risk.amer, main = "Risk Index for AMER region")

acf(risk.other, main = "Risk Index for Other region")
pacf(risk.other, main = "Risk Index for Other region")

#Serial Correlation 
Serial2<- serial.test(var.est2, lags.pt = 12, type = "PT.asymptotic")
Serial2 #p-value = 0,008 < 0,05, meaning that some serial correlation is observed

#Heteroskedasticity
Arch2<- arch.test(var.est2, lags.multi = 12, multivariate.only = TRUE)
Arch2 #p-value = 1 > 0,05, meaning that the model does not suffer from heteroskedastisity
#the model passes the test

#Normal distribution of the residuals 
Norm2<- normality.test(var.est2, multivariate.only = TRUE)
Norm2 #the model passed the skewness test, JB-Test and Kurtosis tests

#Testing for structural breaks in the Residuals
Stability2<- stability(var.est2, type = "OLS-CUSUM") #the system is stable
plot(Stability2)
#stable 

#========================================================
#Granger Causality
GrangerRISK<- causality(var.est2, cause = "Risk")
GrangerRISK
#p-value = 0,000 < 0,05, so we can reject Granger causality H0: Renewables do not Granger-cause CLI Gas Risk; 
#p-value = 0,006 < 0,05, so we can reject H0: No instantaneous causality between: Renewables and Risk Gas Cli

#select order = 1, 3, 6, 12, 24, 48
grangertest(risk, rep, order = 48)
grangertest(risk, gas, order = 48)
grangertest(risk, elc, order = 48)
grangertest(risk, cli, order = 48)

#========================================================
#Impulse Response Functions 
#Risk
rent<- irf(var.est2, n.ahead = 36, boot = TRUE, cumulative = TRUE)
plot(rent)

#========================================================
#========================================================

#SVAR Model
#Setting the Restrictions for the SVAR model
#Contemporaneous coefficients
summary(var.est2)

a.mat <- diag(5)
a.mat[2, 1] <- NA
a.mat[3, 1] <- NA
a.mat[4, 1] <- NA
a.mat[5, 1] <- NA
a.mat[3, 2] <- NA
a.mat[4, 2] <- NA
a.mat[5, 2] <- NA
a.mat[4, 3] <- NA
a.mat[5, 3] <- NA
a.mat[5, 4] <- NA

print(a.mat)

#========================================================
#Buidling the Model
svar.two <- SVAR(var.est2, estmethod = "direct", Amat = a.mat)
svar.two
summary (svar.two)

#========================================================
#Impulse response functions
irf_all<- irf(svar.two, n.ahead = 24, ortho = TRUE, boot = TRUE, cumulative = TRUE, ci = 0.9)
plot(irf(svar.two, n.ahead = 24, ortho = TRUE, boot = TRUE, cumulative = TRUE, ci = 0.9))

#extract irf 
multiple_varirf <- extract_varirf(irf_all)
head(multiple_varirf)

#plot irf with ggplot 
names(multiple_varirf)

#Risk
risk_risk <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_risk_risk, ymin=lower_risk_risk, ymax=upper_risk_risk)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("political risk index- political risk index")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
risk_risk

risk_gas <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_risk_gas, ymin=lower_risk_gas, ymax=upper_risk_gas)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("political risk index- real natural gas")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
risk_gas

risk_electricity <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_risk_electricity, ymin=lower_risk_electricity, ymax=upper_risk_electricity)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("political risk index- electrisity price")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
risk_electricity

risk_renewables <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_risk_renewables, ymin=lower_risk_renewables, ymax=upper_risk_renewables)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("political risk index - renewable energy gerenarion")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
risk_renewables

risk_cli <- multiple_varirf %>% 
  ggplot(aes(x=period, y=irf_risk_cli, ymin=lower_risk_cli, ymax=upper_risk_cli)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("political risk index - cli")+
  ylab(NULL)+
  xlab(NULL) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))
risk_cli

library(patchwork)
risk_irf <-  risk_risk / risk_gas / risk_electricity / risk_renewables / risk_cli
risk_irf

#========================================================
#SVAR Variance Decomposition
FEVDsvar2<- fevd(var.est2, n.ahead = 38)
FEVDsvar2
plot(FEVDsvar2)

#========================================================
#END
#========================================================
