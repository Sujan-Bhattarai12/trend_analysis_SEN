# Name ---------------Sujan Bhattarai
# Project Name---------use logistic regression to calculate probability for precipitaion
                     #phase 
# Date................. May 4, 2023
#::::::::::::::::::::::::::::::::::::::::::::

## to install CRHMr package for updated R, you first need these packages
#install.packages(c("operators", "topmodel", "DEoptim", "XML"))

### the installation of CRHMr package requires download of Ecohydrology package
#install.packages("devtools")
#library(devtools)
#install_github("CentreForHydrology/CRHMr")

#load all required packages that will be used for analysis
library(CRHMr)
library(lubridate)
library(ggplot2)
library(dplyr)
library(trend)
library(changepoint)
library(patchwork)
library(caret)
library(mblm)
library(zoo)
library(pracma)
library(car)
library(gridExtra)
library(xts)

#read WFDEI files
trhup.file <- file.path("C:/Users/sujan/OneDrive/Desktop/R_codes/precipitation_phase/Langtang_WFDEI_teauQsiQlip_1979_2018.obs")
trhup <- readObsFile(trhup.file, timezone = 'MST', quiet = F)

## see the correlation between them 
# cor_matrix <- cor(trhup[, -1])
# corrplot(cor_matrix, method = "circle", type = "upper", order = "hclust")

#wind is measured at 5m above the ground. Multiply it by 0.83 to obtain wind at 2m. 
#Geonor height is at 2m above the ground

trhup$u.1 <- trhup$u.1  
trhup <- changeEatoRH(trhup)
trhup <- subset(trhup, select = c(datetime, t.1, rh.1, u.1, p.1))

#now read precip file
# p <- readObsFile(p.file, timezone = 'MST', quiet = F)
# 
# trhup <- assembleObs(trhu,p, quiet = F)
# rain snow phase determination 
trhup.phase <- phaseCorrect(trhup, Tcol = 1, RHcol = 2, Ucol = 3, Pcol = 4,
                            shield = 0, quiet = F) 

writeObsFile(trhup.phase , file.path("C:/Users/sujan/OneDrive/Desktop/R_codes/precipitation_phase",
                                     'rainsnow_1979_2016.obs'),
             comment = 'Precipitation Phase - harder and pomeroy 2013')

#make new columns for months and days, add seasons to make it easy for  analysis and visualizations 
data <- trhup.phase %>% 
  mutate(phaseCorrectRainRatio= round(phaseCorrectRainRatio, digits=2),
         year= year(datetime),
         month= month(datetime, label = TRUE),
         day= day(datetime),
         period= if_else(hour(datetime)>5 & hour(datetime)< 19, "Day", "Night"),
         season = ifelse(
           month <= "Feb"  | month == "Dec", "winter",
           ifelse( month == "Mar" | month == "Apr" | month == "May", "pre_monsoon",
                   ifelse(month == "Oct" | month == "Nov", "post_monsoon", "monsoon")))) %>% 
  filter(year != 2019 & p.1 != 0)

clone_copy <- data

##### summarized yearly data for exponential moving average
moving_average <- data %>% 
  group_by(year) %>% 
  summarize(temp= mean(t.1), rh= mean(rh.1), wind=mean(u.1), precip=sum(p.1),
            rain_ratio= mean(phaseCorrectRainRatio)) %>% ungroup()

moving_average$temp_ema        <- movavg(moving_average$temp,       n=5, type='e')
moving_average$rh_ema          <- movavg(moving_average$rh,         n=5, type='e')
moving_average$wind_ema        <- movavg(moving_average$wind,       n=5, type='e')
moving_average$precip_ema      <- movavg(moving_average$precip,     n=5, type='e')
moving_average$rain_ratio_ema  <- movavg(moving_average$rain_ratio, n=5, type='e')

###smoother precipitation curve as it shows cyclicity
acf(moving_average$precip_ema)
moving_average$smooth_precip <- moving_average$precip - moving_average$precip_ema
acf(moving_average$smooth_precip)

###trend test 
mk.test(moving_average$temp_ema)
sens.slope(moving_average$temp_ema)
pettitt.test(moving_average$temp_ema)

mk.test(moving_average$rh_ema)
sens.slope(moving_average$rh_ema)
pettitt.test(moving_average$rh_ema)

mk.test(moving_average$wind_ema)
sens.slope(moving_average$wind_ema)
pettitt.test(moving_average$wind_ema)

mk.test(moving_average$smooth_precip)
sens.slope(moving_average$smooth_precip)
pettitt.test(moving_average$smooth_precip)
acf(moving_average$precip_ema)

mk.test(moving_average$rain_ratio_ema)
sens.slope(moving_average$rain_ratio_ema)
pettitt.test(moving_average$rain_ratio_ema)


#######plot graphs for all variables
mkmodel_temp     <- mblm(temp_ema ~ year, data=moving_average)
mkmodel_wind     <- mblm(wind_ema ~ year, data=moving_average)
mkmodel_rh       <- mblm(rh_ema ~ year, data=moving_average)
mkmodel_precip   <- mblm(precip_ema ~ year, data=moving_average)
mkmodel_smoother <- mblm(smooth_precip ~ year, data=moving_average)

###legend
colors <- c("variable_color"= "blue", "weight_color"="lightblue", 'abline_color'= 'red')

#prepare graphs
make_graph <- function(type, variable, weight, xlab, ylab, mkmodel, legend.position="FALSE") {
  plot <- moving_average %>%
    ggplot(aes(year, {{variable}})) +
    xlab(xlab) +
    ylab(ylab) +
    geom_line(data = moving_average, aes(year, {{weight}}, color = 'col2')) + 
    geom_point(data= moving_average, aes(year, {{weight}}, color= "col2"))+
    geom_abline(intercept = coef(mkmodel)[1], slope = coef(mkmodel)[2] ) +
    theme_bw() +
    scale_x_continuous(breaks = seq(1979, 2019, by = 10)) +
    theme(axis.title = element_text(size = 12)) +
    theme(legend.position = legend.position) +
    scale_color_manual(values = c('col2' = '#E69F00'))
  
  if (identical(type, geom_bar)) {
    plot <- plot + geom_bar(stat = "identity", fill = "#0072B2", alpha=0.5)
  } else {
    plot <- plot + type(stat = "identity", aes(color = "col1")) +
      scale_color_manual(values = c('col1' = '#0072B2', 'col2' = '#E69F00'))
  }
  plot <- plot +
    geom_abline(intercept = coef(mkmodel)[1], slope = coef(mkmodel)[2], color = "red")
  return(plot)
}

png("trend_test.tiff", height = 5, width = 7, units = "in", res = 300)
grid.arrange(make_graph( type = geom_line, variable= temp,    weight = temp_ema,   xlab="", ylab = "Temperature °C",  mkmodel= mkmodel_temp),
             make_graph( type = geom_line, variable= wind,    weight = wind_ema,   xlab="", ylab = expression("Wind (ms"^-1*")"), mkmodel= mkmodel_wind),
             make_graph( type = geom_step, variable= rh,     weight = rh_ema,     xlab="", ylab =  "Relative Humidity (%)", mkmodel= mkmodel_rh),
             make_graph( type = geom_bar,  variable =smooth_precip, weight = precip_ema, xlab="", ylab = expression("Precipitation (mm)"), mkmodel=mkmodel_smoother))
dev.off()


### Do MK trend test, Sen's slope in all ema columns
calculate_mk_test <- function(data) {
  start <- which(colnames(data) == "temp_ema")
  end   <- which(colnames(data) == "rain_ratio_ema")
  mk_sen_results <- list()
  
  for (i in start:end) {
    column <- data[[i]]
    mk <- mk.test(column)
    sen <- sens.slope(column)
    petite <- pettitt.test(column)
    mk_sen_results[[colnames(data)[i]]] <- list(mk, sen, petite )
  }
  return(mk_sen_results)
} 
calculate_mk_test(moving_average)

# par(mar = c(1, 1, 1, 1)) # this is used to set the graphical plot window: not a part of climate analysis
# cpt.meanvar(annual_mean$annual_average, method="AMOC", minseglen = 10) # this test identifies the change point using AMOC
# the code below this point test the same thing for annual daytime and nighttime
# filtering data only for days(6am and 6pm) requred in model fitting in next step

annual_Day <- data %>% 
  group_by(year) %>%
  filter( period=="Day") %>% 
  summarize(rain_ratio= mean(phaseCorrectRainRatio)) %>% 
  ungroup() %>% 
  mutate(day_rain_ratio_ema = movavg( rain_ratio, n=5, type='e')) 

mk.test(annual_Day$day_rain_ratio_ema)
sens.slope(annual_Day$day_rain_ratio_ema)
pettitt.test(annual_Day$day_rain_ratio_ema)
x <- cpt.meanvar(annual_Day$day_rain_ratio_ema, method="AMOC", minseglen = 10)# identifies changepoint
param.est(x)

#### for night period
annual_Night <- data %>% 
  group_by(year) %>%
  filter( period=="Night" ) %>% 
  summarize(rain_ratio= mean(phaseCorrectRainRatio)) %>% 
  ungroup() %>% 
  mutate(night_rain_ratio_ema = movavg(rain_ratio, n=5, type='e'))

mk.test(annual_Night$night_rain_ratio_ema) 
sens.slope(annual_Night$night_rain_ratio_ema)
pettitt.test(annual_Night$night_rain_ratio_ema)
x <- cpt.meanvar(annual_Night$night_rain_ratio_ema, method="AMOC", minseglen = 10)# identifies changepoint
param.est(x)

# ggplot does not have readymade function for sens slope.. so manual model is created for day and night
# so that it can be used in visualization part
mkmodel1 <- mblm(day_rain_ratio_ema ~ year, data=annual_Day) #for day 
mkmodel2 <- mblm(night_rain_ratio_ema ~ year, data=annual_Night) #for night

# this data is required for plotting for day and night averages
annual_mean <- data %>%
  group_by(year, period) %>%
  filter(p.1 != 0 ) %>%
  summarize(annual_average= mean(phaseCorrectRainRatio))

# plotting trend#for annual day and night
ggplot(annual_mean, aes(year, annual_average, color= period))+
  geom_line()+
  theme_bw()+
  geom_abline(intercept = coef(mkmodel1)[1], slope = coef(mkmodel1)[2], color = '#0072B2')+
  geom_abline(intercept = coef(mkmodel2)[1], slope = coef(mkmodel2)[2], color = '#0072B2')+
  ylab("Annual rainfall ratio")+
  ylim(0.25, 1)+
  scale_x_continuous(breaks = seq(1979, 2019, by = 10))+
  theme(axis.title = element_text(size = 10))+
  geom_line(data=annual_Day, aes(year, day_rain_ratio_ema), color="deeppink")+
  #annotate("text", x = 1999, y = 0.25, label = "sen's slope for diurnal period= 0.0014, p value= 0.02", size= 4)+
  geom_line(data =annual_Night, aes(year, night_rain_ratio_ema), color="deeppink")+  #454545
  #annotate("text", x = 1999, y = 0.3, label = "sen's slope for night period= 0.0008, p value= 0.2", size= 4)+
  scale_colour_manual(name= "", values=colors)+
  theme(legend.position="top")+
  scale_color_manual(values = c('col1' = '#0072B2', 'col2' = '#E69F00'))


ggsave("day-and-night.tiff", height=5, width=4) # saves the previously created ggplot graphics in high quality

# for august and for all other months
#caluculating rainfall ratio trend for aug
# for other months, need to change month value to different one

monthwise_test <- function(x, y) {
  filtered_data <- data %>% 
    filter(month == x &  period==y) %>% 
    group_by(year, month, period) %>% 
    summarize(mean = mean(phaseCorrectRainRatio)) %>% 
    ungroup() %>% 
    mutate(mean_ema = movavg(mean, n=5, type='e')) 
  return(list(filtered_data, mk.test(filtered_data$mean_ema), sens.slope(filtered_data$mean_ema),
              pettitt.test(filtered_data$mean_ema), param.est(cpt.meanvar(filtered_data$mean_ema, method="AMOC", minseglen = 10)) ))
}

monthwise_test("Jan", "Day")
monthwise_test("Feb", "Day")
monthwise_test("Mar", "Day")
monthwise_test("Apr", "Day")
monthwise_test("May", "Day")
monthwise_test("Jun", "Day")
monthwise_test("Jul", "Day")
monthwise_test("Aug", "Day")
monthwise_test("Sep", "Day")
monthwise_test("Oct", "Day")
monthwise_test("Nov", "Day")
monthwise_test("Dec", "Day")

monthwise_test("Jan", "Night")
monthwise_test("Feb", "Night")
monthwise_test("Mar", "Night")
monthwise_test("Apr", "Night")
monthwise_test("May", "Night")
monthwise_test("Jun", "Night")
monthwise_test("Jul", "Night")
monthwise_test("Aug", "Night")
monthwise_test("Sep", "Night")
monthwise_test("Oct", "Night")
monthwise_test("Nov", "Night")
monthwise_test("Dec", "Night")

# par(mar = c(1, 1, 1, 1))
# logistic regression for temperature boundary point
# phase is estimated based on rain ratio.. with greater than 0.5 as 1, and less than 0.5 as 0
#then creating phase column as factor, which is enforced to logistic model

#### data for logistic maping
summary_log <- data %>% 
  group_by(year, month, day) %>% 
  summarize(temp= mean(t.1), 
            rh= mean(rh.1),
            wind=mean(u.1), 
            pp=mean(p.1),
            rain_ratio= mean(phaseCorrectRainRatio),
            season= unique(season)) %>% 
  mutate(phase=as.factor(ifelse(rain_ratio > 0.5, 1, 0))) 

#summary_log$temp_1 <-  rnorm(9546, 0, 0.01) + summary_log$temp
##define the model

log_regression <- glm(phase ~  temp + rh + wind, family = "binomial", data= summary_log)
summary(log_regression)
vif(log_regression)

# plotting logistic graph between temp and rainfall fraction 
ggplot(summary_log, aes(temp, rain_ratio))+
  geom_point(aes(color= season), alpha=0.5)+
  theme_bw()+ ylab("rainfall ratio")+ 
  scale_y_continuous(breaks = seq(0.00, 1.00, by = 0.10))+
  geom_jitter(width = 0.01, alpha=0.1)+
  geom_smooth(method= 'glm',  method.args=list(family="binomial"))+
  theme(axis.title = element_text(size = 14))+ xlab("Temperature °C")+ ylab("Logistic probability")+
  theme(legend.text=element_text(size= 12))

ggsave("logistic.png", width = 7, height= 7, type="cairo", dpi= 300)

## visual representation of rain-ratio with rall variables
visual_relation <- function(variable){
  summary_log %>% 
    ggplot(aes({{variable}}, rain_ratio))+
    geom_point(aes(color= season))+
    geom_smooth(method= "glm",  method.args=list(family="binomial"))+
    theme_bw()+
    geom_jitter( alpha=0.2)
}

visual_relation(rh)
visual_relation(wind)
visual_relation(temp)
visual_relation(pp)

# creating graph for sample year 2000

datetime <- seq(
  from = as.Date("1979/01/01"),
  to = as.Date("2018/12/31"),
  by = "days")

data <- trhup.phase %>% 
  mutate(phaseCorrectRainRatio= round(phaseCorrectRainRatio, digits=2),
         year= year(datetime),
         month= month(datetime, label = TRUE),
         day= day(datetime),
         period= if_else(hour(datetime)>5 & hour(datetime)< 19, "Day", "Night"),
         season = ifelse( month <= "Feb"  | month == "Dec", "winter",
                          ifelse( month == "Mar" | month == "Apr" | month == "May", "pre_monsoon",
                                  ifelse(month == "Oct" | month == "Nov", "post_monsoon", "monsoon")))) %>% 
  filter(year != 2019)

time_series <- data %>% 
  group_by(year, month, day) %>% 
  summarize(rain=sum(phaseCorrectRain), snow= sum(phaseCorrectSnow), temp= mean(t.1))

ts2 <- xts(time_series[4:5], order.by = datetime)
ts2_2000 <- ts2["2000"]

dev.new()
png( "rain and snow.tiff", height = 4, width = 10, units = "in", res = 300)
plot( ts2$rain, lwd = 0.5, xlim=c(1979, 2018), col = "#454545", main = NA,
      yaxis.right = F, ylab = "Precipitation(mm/day)", grid.col=NA, cex=0.7)
lines(ts2$snow, col = 4, lwd = 0.5)
addLegend(legend.loc = "topright", legend.names = c("rainfall", "snowfall"), col = c("#454545", 4),
          lty = 1, lwd = 1.5)

axis(side = 4, at = pretty(ts2_2000$phaseCorrectSnow.sum), col = 4, lwd = 1.6)
dev.off()

############## plot for the sample year

dev.new()
png("rain and snow 2000.tiff", height = 4, width = 10, units = "in", res = 300)

plot(ts2_2000$rain, lwd = 1, xlim(1979, 2018), col = "#454545",
     main = NA, yaxis.right = F, ylab = "Precipitation(mm/day)", cex=0.7,
     grid.col=NA)
lines(ts2_2000$snow, col = 4,lwd = 1)
addLegend(
  legend.loc = "topright",
  legend.names = c("rainfall", "snowfall"),
  col = c("#454545", 4),
  lty = 1,
  lwd = 1.5)
dev.off()

# month_wise grid plot with 40 lines
lines_graph_maker <- function(calculate, column, ylabel, legend){
  monthly_data <- data %>%  
    group_by(year, month) %>% 
    summarize(mean = calculate({{column}}))
  
  monthly_data$year <- factor(monthly_data$year)
  color_palette <- colorRampPalette(c("#0165fc", "yellow", "red"))(length(unique(monthly_data$year)))
  
  ggplot(monthly_data, aes(month, mean, group = year, color = year)) +
    geom_line(size = ifelse(between(as.integer(monthly_data$year), 1, 1) |
                              between(as.integer(monthly_data$year), 40, 40), 0.7, 0.2),
              alpha = 1) + 
    scale_color_manual(values = color_palette) +
    # geom_line(data = monthly_data[1, ], aes(group = year, color = year), size = 2) +
    # geom_line(data = monthly_data[length(monthly_data$year), ], aes(group = year, color = year), size = 2) +
    theme_bw() +
    theme(axis.title = element_text(size = 10)) +
    ylab(ylabel) +
    xlab("") +
    theme(legend.position = legend, legend.text = element_text(size = 12))
}

plot1  <- lines_graph_maker(max,  t.1,  "Maximum Daily Temperature °C", legend="none")
plot2  <- lines_graph_maker(min,  t.1,  "Minimum Daily Temperature °C", legend="none")
plot3  <- lines_graph_maker(mean, rh.1,  "Relative Humidity (%)", legend="none")
plot4  <- lines_graph_maker(mean, u.1,  expression("Wind (ms"^-1*")"), legend="none")
plot5  <- lines_graph_maker(sum,  p.1,  "Precipitation (mm)", legend="none")

png("4plots_combined.tiff", height = 8, width = 8, units = "in", res = 300)
grid.arrange(plot1, plot2, plot3, plot4, plot5)
dev.off()

png("oneplotwithlegend.tiff", height = 7, width = 7, units = "in", res = 300)
lines_graph_maker(sum,  p.1,  "Precipitation (m/s)", legend="right")
dev.off()


####seasonal mathematical measure 

seasonal <- function(x, y) {
  filtered_data <- clone_copy  %>% 
    filter(season == x &  period==y) %>% 
    group_by(year, season, period) %>% 
    summarize(mean = mean(phaseCorrectRainRatio)) %>% 
    ungroup() %>% 
    mutate(mean_ema = movavg(mean, n=5, type='e')) 
  return(list(filtered_data, mk.test(filtered_data$mean_ema), sens.slope(filtered_data$mean_ema),
              pettitt.test(filtered_data$mean_ema), cpt.meanvar(filtered_data$mean_ema, method="AMOC", minseglen = 10)))
}

seasonal("monsoon", "Day")
seasonal("monsoon", "Night")
seasonal('pre_monsoon', "Day")
seasonal('pre_monsoon', "Night")
seasonal('post_monsoon', "Day")
seasonal('post_monsoon', "Night")
seasonal('winter', "Day")
seasonal('winter', "Night")

######rain to snow precipitation ration
volume <- data %>%
  group_by(year) %>% 
  summarize(rain = mean(phaseCorrectRain),
            snow = mean(phaseCorrectSnow)) %>% 
  mutate(ratio=rain/snow) %>% 
  summarize(mean = mean(ratio),
            lower_ci = quantile(ratio, probs = 0.025),
            upper_ci = quantile(ratio, probs = 0.975))


