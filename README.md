
# Using Multivariate logistic regression on climate data to detect climate change
# Project Background

## Highlights

- Langtang region receives 57-78% of annual total precipitation as rain and (22-43%) as snow. June, July, and August receive most precipitation as rain.
- Rain-snow transition temperature occurs at 1.78 °C.
- Rainfall fraction has an increasing trend: annual daytime (+0.17% per year), annual nighttime (0.10% per year), and for all timestamps except winter nights.
- Changepoints in rainfall fraction occurred during 1990-2000.

## Abstract

Mountain studies worldwide document changes in rainfall fraction due to climate change. Most mountain systems have an increasing trend in rainfall fraction due to shifting snow precipitation to rain. In regards to Nepal, which occupies an 800km long belt of the Hindu Kush Himalaya, precipitation phase trend is not well known. This study thus conducts a phase study on the Langtang region and investigates the snow and rain precipitation to changing climate.

The study uses 40 years (1979-2018) of bias-corrected WFDEI climate reanalysis data and employs a physically based psychrometric energy model to partition precipitation into rain and snow. The study identifies points of statistically significant changes in trends using Exponential Weighted Average, Mann Kendall's test, Sen’s slope estimator, and changepoint analysis. Changepoint is identified using two methods: Pettit's method and AMOC (At Most One Changepoint).

In addition, the study identifies transient temperature for rain-snow using logistic mapping. Based on the results, rainfall fraction had an increasing trend for all timestamps (annual, seasonal, monthly) other than winter nights, with changepoints between 1990-2000. Post-monsoon months (October and November) had the greatest annual increasing rate for daytime at 0.34% and 0.27% respectively. Winter months had the least changes with negligible rising trend especially for the nights. Similarly, the transient temperature occurred at 1.78 °C.

## Methodology and Tools

### Climate Data Collection

The study used the WFDEI climate reanalysis datasets (G. P. Weedon et al., 2011; Graham P. Weedon et al., 2014), for its high accuracy in climate estimates in the HKH region (Khan et al., 2017). WFDEI also has much higher accuracy in precipitation estimates than other reanalysis datasets (Bhattarai et al., 2020; Dahri et al., 2021; L. Li et al., 2013). WFDEI provides the global (land only) meteorological datasets from 1979 to 2018, at 0.5°x 0.5° latitude-longitude grid at hourly resolution. For this study, the required data were archived from Pradhananga et al., (2021). The reanalysis data were bias-corrected using the quantile mapping in R and in-situ data from Kyanjing meteorological station (location in figure 1) (D. Pradhananga et al., 2021).

### Psychrometric Energy Balance Method

The psychrometric energy balance method is a physically based method that integrates temperature, humidity, and wind for estimating precipitation phase. The method is based on a principle that several microphysics phenomena determine the phase of a falling hydrometeor. For instance, the energy fluxes occur between the saturated surface of a falling hydrometeor and unsaturated air in the atmosphere (Stewart, 1992); thus precipitation phase estimation must incorporate latent heat transfer along with other physical processes. The method accounts for such thermodynamics, making it a robust method compared to the temperature-only method for phase estimation. Equation (I) represents the mathematical model for estimating the net temperature of falling hydrometeor under the psychrometric method (Harder & Pomeroy, 2013):

\[ T_i = T_a + \frac{D}{\lambda_t} L (\rho_{Ta} - \rho_{sat}(T_i)) \]

The equation calculates the net temperature for a falling hydrometeor. Based on this, the fraction of rainfall ratio in reference to total precipitation is determined. The rain ratio values range from 0 to 1, where higher values indicate higher rainfall fraction compared to total precipitation (Harder & Pomeroy, 2013):

\[ f_r(T_i) = \frac{\sum_{T_i} \text{rainfall(mm)}}{\sum_{T_i} \text{rainfall(mm)} + \sum_{T_i} \text{snowfall(mm)}} \]

In the study, the psychrometric energy model was enforced using the `CRHMr` package (Shook, 2021) in R (version 4.1.2). The `phaseCorrect` function inside the `CRHMr` package takes in hourly time series data for temperature, wind, relative humidity, and precipitation as inputs and partitions the precipitation into snow or rain based on net hydrometeor temperature.

### Data Analysis

#### Mann-Kendall's (MK) Trend Test

Mann-Kendall trend test checks for the monotonic association of a single variable (e.g., Y) over a time series. The test works for a distribution-free dataset, assuming the absence of seasonality/autocorrelation in the dataset. The test checks for concordance and discordance pairs in time series to calculate the test statistics. Higher values for the test statistic indicate a greater monotonic trend and vice versa. 

The test hypothesis is as follows:

- Null hypothesis (H0): monotonic trend does not exist in time series
- Alternative hypothesis (H1): trend exists (linear, cyclic)

The test was calculated for exponential weighted average values using the `mk.test()` function inside the `trend` package in R at a 5% confidence interval.

#### Sen's Slope Estimator

Sen’s slope estimator quantifies the magnitude based on observational data. The Sen’s slope calculates the median slope, which is the measure of central tendency for non-parametric data:

\[ \beta = \text{Median} \left( \frac{x_i - x_j}{i - j} \right), \forall \, j < i \text{ where } 1 < j < i < n \]

The positive \(\beta\) represents the increasing trend in time-series observations and vice versa. In the study, the Sen’s slope was calculated using the `sens.slope` function under the `trend` package in R.

#### Changepoint Detection

The changepoint analysis identifies the point of statistical properties change in a time-series data. It is based on an assumption if a changepoint exists at point “t” between Y1...Yn, Y1 to Yt differs in statistical measures (mean and variance) with Yt to Yn. In this study, a single changepoint was identified using the AMOC (at most one changepoint) method and Pettit's method. The AMOC method was enforced using the `cpt.meanvar()` function inside the `changepoint` package in R. Pettit's method also assumes variation in distribution of data in time series in changepoint.

#### Logistic Regression Mapping

Logistic regression models the mathematic relationships between explanatory variables and binary response variable. The logistic regression scales the binary values to log odds, which is a fractional measure of event occurring to event not occurring. The log odds can be used to predict probabilities using the model coefficients:

\[ P(\text{rain}) = \frac{1}{1 + e^{(\alpha + \beta T_s)}} \]

Logistic regression is often used in estimating phase point change, where the log odds are fitted into the probability equation. In this study, the logistic model is fitted for temperature and rainfall ratio, and the half-point temperature for the Langtang basin is identified.

## References

1. Bhattarai, et al., 2020
2. Dahri, et al., 2021
3. Harder & Pomeroy, 2013
4. Khan, et al., 2017
5. Kendall, 1975
6. Killick, et al., 2012
7. Pradhananga, et al., 2021
8. Sen, 1974
9. Shook, 2021
10. Stewart, 1992
11. Weedon, et al., 2011; 2014
12. Wambui, et al., 2015
13. Wijngaard, et al., 2003
14. Ogungbenro & Morakinyo, 2014
15. Shahid, et al., 2018
16. Wang, et al., 2019
