anatambea introduction
================
Joseph Hicks
2025-01-08

anatambea is a tool for applying Dust, Odin, and MCState to fit a
semi-stochastic mechanistic disease transmission model to continuously
collected prevalence data using particle Markov chain Monte Carlo
(pMCMC). Specifically, this package has been developed to fit a
semi-stochastic version of ‘malariasimulation’ to monthly malaria
prevalence among pregnant women at first antenatal care (ANC) visit.
This vignette walks through a simple example and describes various
options and settings that can be applied.

## Understanding

In the deterministic version of ‘malariasimulation,’ seasonal variation
in malaria transmission intensity is controlled by a Fourier series
whose coefficients have been previously estimated from historical
rainfall data. Different Fourier series equations have been estimated
for each administrative area (level 1) across 80 countries, but as these
are smoothed averages of rainfall, year-on-year differences in rainfall
and therefore malaria transmission intensity are not captured resulting
in non-realistic, rigidly cyclical model estimates. While this provides
a decent seasonal profile for fitting cross-sectional data collected
every few years (as occurs with Demographic Health Surveys and Malaria
Indicator Surveys), this method is not flexible enough for fitting the
model to monthly prevalence data, as collected at ANC1.

To provide more flexibility, we have replaced this seasonal mechanism
with a stochastic process, in which the emergence rate of adult
mosquitoes from pupae varies in a random-walk. Each month, the change in
mosquito emergence rate is determined by a normally distributed random
variable, *δ*, and a volatility constant, *σ*. To ensure *β* values
remains realistic, the upper limit of the random walk was restricted to
*β_max*. The pMCMC then filters likely trajectories of the emergence
rate by comparing observed prevalence with that estimated by the model,
given the emergence rate trend. The output provides a posterior
distribution of likely model trajectories. Because the relationship
between mosquito emergence rate and infection prevalence is defined by a
mechanistic model, we can extract trends in other relevant indicators,
such as the entomological inoculation rate (EIR) or clinical incidence,
which are difficult to observe directly.

## Input Data

In this simple example, we will fit our model to simulated monthly
prevalence data in children under 5 years old. Data for fitting should
be a data frame with at least three columns:

1.  `month`: Month (here formatted as ‘yearmon’, provided by the package
    ‘zoo’)
2.  `tested`:The number of individuals tested
3.  `positive`: Of those tested, the number of positive test results

Like so,

``` r
data_slim <- dplyr::select(anatambea::data_sim,month,positive,tested)
data_slim
#>       month positive tested
#> 1  2015.000       57    181
#> 2  2015.083       59    234
#> 3  2015.167       54    172
#> 4  2015.250       68    231
#> 5  2015.333       82    274
#> 6  2015.417       68    230
#> 7  2015.500       48    185
#> 8  2015.583       59    239
#> 9  2015.667       60    167
#> 10 2015.750       94    187
#> 11 2015.833       97    217
#> 12 2015.917       44    129
#> 13 2016.000       92    271
#> 14 2016.083       95    231
#> 15 2016.167       78    200
#> 16 2016.250       63    146
#> 17 2016.333      114    255
#> 18 2016.417       89    207
#> 19 2016.500       56    138
#> 20 2016.583       96    254
#> 21 2016.667       76    248
#> 22 2016.750       51    162
#> 23 2016.833       63    278
#> 24 2016.917       53    189
#> 25 2017.000       54    175
#> 26 2017.083       34    186
#> 27 2017.167       52    216
#> 28 2017.250       46    266
#> 29 2017.333       68    250
#> 30 2017.417       84    304
#> 31 2017.500       60    134
#> 32 2017.583      137    284
#> 33 2017.667       82    195
#> 34 2017.750       82    209
#> 35 2017.833       92    249
#> 36 2017.917       99    259
#> 37 2018.000      117    298
#> 38 2018.083       81    184
#> 39 2018.167       74    223
#> 40 2018.250       74    236
#> 41 2018.333       70    212
#> 42 2018.417       75    230
#> 43 2018.500       51    178
#> 44 2018.583       62    257
#> 45 2018.667       41    187
#> 46 2018.750       52    266
#> 47 2018.833       45    214
#> 48 2018.917       10     73
#> 49 2019.000       51    219
#> 50 2019.083       80    244
#> 51 2019.167       86    251
#> 52 2019.250       43    180
#> 53 2019.333       66    210
#> 54 2019.417       52    199
#> 55 2019.500       50    209
#> 56 2019.583       34    174
#> 57 2019.667       46    300
#> 58 2019.750       31    241
#> 59 2019.833       35    265
#> 60 2019.917       43    223
#> 61 2020.000       33    266
```

## Running a pMCMC fitting

The `run_pmcmc` function is the central tool to fit our model to
observed data. It formats the provided data set, sets up required
functions and parameters, and finally runs the pMCMC. Broadly, the flow
of actions is as follows:

1.  Process provided data set, specifically formatting necessary time
    variables.
2.  Format and declare parameters. These include both constant or known
    parameter values needed for malariasimulation as well as pMCMC
    parameters that will be fit.
3.  Set particle filter and pMCMC settings.
4.  Run the pMCMC.
5.  Format output.

``` r
result <- anatambea::run_pmcmc(data_raw=data_slim,
                       n_particles=50,
                       proposal_matrix = matrix(1),
                       target_prev = 0.4,
                       target_prev_group = 'u5',
                       max_param=125,
                       prop_treated = 0.4,
                       n_steps = 10,
                       seed = 1L,
                       start_pf_time = 30*12,
                       comparison = 'u5',
                       initial = 'informed',
                       state_check = 0)
#> ── R CMD INSTALL ───────────────────────────────────────────────────────────────
#>   ─  installing *source* package 'odin.model.stripped.seasonal9f89af4b' ...
#>      ** using staged installation
#>      ** libs
#>      using C compiler: 'gcc.exe (GCC) 13.2.0'
#>      gcc  -I"C:/PROGRA~1/R/R-44~1.1/include" -DNDEBUG     -I"C:/rtools44/x86_64-w64-mingw32.static.posix/include"     -O2 -Wall -gdwarf-2 -mfpmath=sse -msse2 -mstackrealign  -UNDEBUG -Wall -pedantic -g -O0 -c odin.c -o odin.o
#>      odin.c: In function 'odin_model_stripped_seasonal_initial_conditions':
#>      odin.c:2286:10: warning: unused variable 't' [-Wunused-variable]
#>     2286 |   double t = scalar_real(t_ptr, "t");
#>          |          ^
#>      gcc  -I"C:/PROGRA~1/R/R-44~1.1/include" -DNDEBUG     -I"C:/rtools44/x86_64-w64-mingw32.static.posix/include"     -O2 -Wall -gdwarf-2 -mfpmath=sse -msse2 -mstackrealign  -UNDEBUG -Wall -pedantic -g -O0 -c registration.c -o registration.o
#>      gcc -shared -static-libgcc -o odin.model.stripped.seasonal9f89af4b.dll tmp.def odin.o registration.o -LC:/rtools44/x86_64-w64-mingw32.static.posix/lib/x64 -LC:/rtools44/x86_64-w64-mingw32.static.posix/lib -LC:/PROGRA~1/R/R-44~1.1/bin/x64 -lR
#>      installing to C:/Users/jthicks/AppData/Local/Temp/RtmpMHAlEk/devtools_install_45d41b5d2d3a/00LOCK-file45d46f565810/00new/odin.model.stripped.seasonal9f89af4b/libs/x64
#>   ─  DONE (odin.model.stripped.seasonal9f89af4b)
#> ── R CMD INSTALL ───────────────────────────────────────────────────────────────
#>   ─  installing *source* package 'dusta76758ba' ...
#>      ** using staged installation
#>      ** libs
#>      using C++ compiler: 'G__~1.EXE (GCC) 13.2.0'
#>      g++ -std=gnu++17  -I"C:/PROGRA~1/R/R-44~1.1/include" -DNDEBUG  -I'C:/Users/jthicks/AppData/Local/R/win-library/4.4/cpp11/include'   -I"C:/rtools44/x86_64-w64-mingw32.static.posix/include"  -I"C:/Users/jthicks/AppData/Local/R/win-library/4.4/dust/include" -DHAVE_INLINE -fopenmp    -O2 -Wall  -mfpmath=sse -msse2 -mstackrealign  -Wall -pedantic -c cpp11.cpp -o cpp11.o
#>      g++ -std=gnu++17  -I"C:/PROGRA~1/R/R-44~1.1/include" -DNDEBUG  -I'C:/Users/jthicks/AppData/Local/R/win-library/4.4/cpp11/include'   -I"C:/rtools44/x86_64-w64-mingw32.static.posix/include"  -I"C:/Users/jthicks/AppData/Local/R/win-library/4.4/dust/include" -DHAVE_INLINE -fopenmp    -O2 -Wall  -mfpmath=sse -msse2 -mstackrealign  -Wall -pedantic -c dust.cpp -o dust.o
#>      g++ -std=gnu++17 -shared -s -static-libgcc -o dusta76758ba.dll tmp.def cpp11.o dust.o -fopenmp -LC:/rtools44/x86_64-w64-mingw32.static.posix/lib/x64 -LC:/rtools44/x86_64-w64-mingw32.static.posix/lib -LC:/PROGRA~1/R/R-44~1.1/bin/x64 -lR
#>      installing to C:/Users/jthicks/AppData/Local/Temp/RtmpMHAlEk/devtools_install_45d4660117a3/00LOCK-file45d43b9d4210/00new/dusta76758ba/libs/x64
#>   ─  DONE (dusta76758ba)
#> 
#> Optimizing initial EIR based on target prevalence.
#> Initial EIR set to 29.3.
#> Running chain 1 / 1
#> Finished 10 steps in 55 secs
#> Time difference of 60.92886 secs
```

## Visualising output

The figure sizes have been customised so that you can easily put two
images side-by-side.

``` r
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
measure_names <- c('prev_05','clininc_all','EIR','betaa')
months <- data_slim$month

timelength <- length(months)
num_months <- length(result$history[1,1,])
start_month_index <- (num_months-timelength+1)

  mcmc <- result$mcmc[,]
  ar <- 1 - coda::rejectionRate(coda::as.mcmc(mcmc))
  ess <- coda::effectiveSize(coda::as.mcmc(mcmc))
  pars_list <- names(mcmc)

history.dfs <- lapply(measure_names,function(x){
  as.data.frame(t(result$history[x,,start_month_index:num_months]))%>%
        dplyr::mutate(month=months)%>%
    reshape2::melt(id='month')%>%
    dplyr::group_by(month)%>%
    dplyr::summarise(median=median(value),
              upper=quantile(value,probs=0.975),
              lower=quantile(value,probs=0.025))%>%
    dplyr::mutate(measure=x)
})
names(history.dfs) <- measure_names

library('ggplot2')

ggplot(data=history.dfs[['prev_05']])+
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper),fill='lightblue')+
  geom_line(aes(x=month,y=median),color='darkblue')+
  geom_point(data=data_slim,aes(x=month,y=positive/tested),fill='black')+
  theme_minimal()

ggplot(data=history.dfs[['clininc_all']])+
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper),fill='lightblue')+
  geom_line(aes(x=month,y=median),color='darkblue')+
  facet_grid(.~measure)+
  theme_minimal()

ggplot(data=history.dfs[['EIR']])+
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper),fill='lightblue')+
  geom_line(aes(x=month,y=median),color='darkblue')+
  facet_grid(.~measure)+
  scale_y_log10()+
  theme_minimal()

ggplot(data=history.dfs[['betaa']])+
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper),fill='lightblue')+
  geom_line(aes(x=month,y=median),color='darkblue')+
  facet_grid(.~measure)+
  scale_y_log10()+
  theme_minimal()
```

![](C:/Users/jthicks/AppData/Local/Temp/RtmpSu1GzM/preview-1e687169345.dir/mamasante_intro_files/figure-gfm/unnamed-chunk-1-1.png)![](C:/Users/jthicks/AppData/Local/Temp/RtmpSu1GzM/preview-1e687169345.dir/mamasante_intro_files/figure-gfm/unnamed-chunk-1-2.png)![](C:/Users/jthicks/AppData/Local/Temp/RtmpSu1GzM/preview-1e687169345.dir/mamasante_intro_files/figure-gfm/unnamed-chunk-1-3.png)![](C:/Users/jthicks/AppData/Local/Temp/RtmpSu1GzM/preview-1e687169345.dir/mamasante_intro_files/figure-gfm/unnamed-chunk-1-4.png)

You can enable figure captions by `fig_caption: yes` in YAML:

    output:
      rmarkdown::html_vignette:
        fig_caption: yes

Then you can use the chunk option `fig.cap = "Your figure caption."` in
**knitr**.

## More Examples

You can write math expressions, e.g. $Y = X\beta + \epsilon$,
footnotes[^1], and tables, e.g. using `knitr::kable()`.

|                   |  mpg | cyl |  disp |  hp | drat |    wt |  qsec |  vs |  am | gear | carb |
|:------------------|-----:|----:|------:|----:|-----:|------:|------:|----:|----:|-----:|-----:|
| Mazda RX4         | 21.0 |   6 | 160.0 | 110 | 3.90 | 2.620 | 16.46 |   0 |   1 |    4 |    4 |
| Mazda RX4 Wag     | 21.0 |   6 | 160.0 | 110 | 3.90 | 2.875 | 17.02 |   0 |   1 |    4 |    4 |
| Datsun 710        | 22.8 |   4 | 108.0 |  93 | 3.85 | 2.320 | 18.61 |   1 |   1 |    4 |    1 |
| Hornet 4 Drive    | 21.4 |   6 | 258.0 | 110 | 3.08 | 3.215 | 19.44 |   1 |   0 |    3 |    1 |
| Hornet Sportabout | 18.7 |   8 | 360.0 | 175 | 3.15 | 3.440 | 17.02 |   0 |   0 |    3 |    2 |
| Valiant           | 18.1 |   6 | 225.0 | 105 | 2.76 | 3.460 | 20.22 |   1 |   0 |    3 |    1 |
| Duster 360        | 14.3 |   8 | 360.0 | 245 | 3.21 | 3.570 | 15.84 |   0 |   0 |    3 |    4 |
| Merc 240D         | 24.4 |   4 | 146.7 |  62 | 3.69 | 3.190 | 20.00 |   1 |   0 |    4 |    2 |
| Merc 230          | 22.8 |   4 | 140.8 |  95 | 3.92 | 3.150 | 22.90 |   1 |   0 |    4 |    2 |
| Merc 280          | 19.2 |   6 | 167.6 | 123 | 3.92 | 3.440 | 18.30 |   1 |   0 |    4 |    4 |

Also a quote using `>`:

> “He who gives up \[code\] safety for \[code\] speed deserves neither.”
> ([via](https://twitter.com/hadleywickham/status/504368538874703872))

[^1]: A footnote here.
