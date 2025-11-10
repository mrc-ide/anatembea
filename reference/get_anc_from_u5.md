# Return ANC prevalence given prevalence in \<5 year old children

`get_anc_from_u5` Return ANC prevalence given prevalence in \<5 year old
children

## Usage

``` r
get_anc_from_u5(prev_u5, avg_prev)
```

## Arguments

- prev_u5:

  Prevalence in under 5 year old children

- avg_prev:

  Average level of prevalence in children under 5 years old to set level
  of immunity within population

## Examples

``` r
get_anc_from_u5(prev_u5=0.3,avg_prev=0.2)
#> $prev_preg_pg
#> [1] 0.3584845
#> 
#> $prev_preg_sg
#> [1] 0.2536578
#> 
#> $prev_preg_mg
#> [1] 0.2067636
#> 
#> $prev_preg_all
#> [1] 0.2400854
#> 
```
