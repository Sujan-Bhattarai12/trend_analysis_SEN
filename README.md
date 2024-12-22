### Using Multivariate Logistic Regression on Climate Data to Detect Climate Change
##### Highlights of the Findings
The Langtang region receives 57-78% of its annual total precipitation as rain and 22-43% as snow. Most precipitation as rain occurs during June, July, and August.
The rain-snow transition temperature is 1.78 °C.
The rainfall fraction shows an increasing trend annually: daytime (+0.17% per year), nighttime (+0.10% per year), and across all timestamps except during winter nights.
Significant changepoints in rainfall fraction were observed between 1990 and 2000.
##### Abstract
Studies worldwide report changes in rainfall fraction in mountain regions due to climate change, often reflecting a shift from snow to rain precipitation. However, for Nepal’s Hindu Kush Himalaya region, the trend in precipitation phase changes is not well documented. This study focuses on the Langtang region, analyzing changes in snow and rain precipitation over 40 years (1979-2018) using bias-corrected WFDEI climate reanalysis data.

Advanced statistical and physical models, including psychrometric energy models, were used to partition precipitation into rain and snow. Statistical techniques like Exponential Weighted Average, Mann-Kendall's test, Sen's slope estimator, and changepoint analysis (Pettit's method and AMOC) were employed to identify trends and changepoints.

The results reveal an increasing trend in rainfall fraction for most timestamps except winter nights, with significant changepoints between 1990 and 2000. Post-monsoon months showed the highest annual increase, while winter nights exhibited minimal change. The rain-snow transition temperature was identified as 1.78 °C.

##### Methodology and Tools
###### Climate Data Collection
WFDEI climate reanalysis datasets were used for their high accuracy in estimating precipitation in the HKH region. These datasets, covering 1979-2018 at a 0.5°x0.5° resolution, were bias-corrected using quantile mapping with in-situ data from the Kyanjing meteorological station.

###### Psychrometric Energy Balance Method
This method integrates temperature, humidity, and wind data to estimate precipitation phase. It considers energy fluxes between falling hydrometeors and the atmosphere, incorporating latent heat transfer and other thermodynamic processes. The rainfall fraction, representing the ratio of rainfall to total precipitation, was calculated using this method.

The study implemented the psychrometric energy balance model using the CRHMr package in R. The phaseCorrect function processed hourly time series data for temperature, wind, relative humidity, and precipitation to partition precipitation into rain or snow.

###### Data Analysis
Mann-Kendall's Trend Test
The Mann-Kendall test was used to assess the monotonic trends in precipitation data. It evaluates the association of data over time by analyzing concordant and discordant pairs. A positive test statistic indicates an increasing trend, while a negative one suggests a decreasing trend.

###### Sen's Slope Estimator
Sen's slope estimator quantifies the trend magnitude in the data. It calculates the median slope of all pairwise differences in time-series observations, with a positive slope indicating an upward trend.

###### Changepoint Detection
Changepoint analysis identifies significant changes in the statistical properties of time-series data. Two methods were employed:

AMOC (At Most One Changepoint): Identifies a single changepoint in the data using mean and variance analysis.
Pettit's Method: Detects changes in data distribution.

##### Logistic Regression Mapping
Logistic regression models the relationship between temperature and the probability of rain. By fitting temperature data into a logistic function, the study estimated the rain-snow transition temperature as 1.78 °C, representing the point where rainfall and snowfall probabilities are equal.
