# match admin region

`admin_match` Matches the user input admin unit and country with data

## Usage

``` r
admin_match(admin_unit = NULL, country = NULL, admin_units_seasonal)
```

## Arguments

- admin_unit:

  Character for admin region. Some fuzzy logic will be used to match. If
  not provided then no seasonality is introduced. Default = NULL

- country:

  Character for country within which admin unit is in. Default = NULL

- admin_units_seasonal:

  Dataframe of seasonality data for country and admin unit
