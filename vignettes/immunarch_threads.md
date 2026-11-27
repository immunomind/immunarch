# Limit the number of threads used by immunarch


## What is happening under the hood?

`immunarch` uses `immundata`, which uses `duckplyr`, which runs queries
in **DuckDB**. It is quite a journey from the data to your plots.

DuckDB can use many CPU cores by default to run faster, but sometimes
you want to limit this.

## Why limit threads?

- On shared machines (servers, CI), many threads can slow down other
  users.
- In tutorials, fewer threads make examples more predictable and easier
  to reproduce.
- Lower threads = lower CPU load (but slower queries).

## How to limit the number of threads

Put this near the start of your tutorial/script:

``` r
# Limit the number of CPU threads used by DuckDB in this R session
duckplyr::db_exec("SET threads TO 1")
```

Change later (example: use 4 threads) or reset to default:

``` r
duckplyr::db_exec("SET threads TO 4")   # use 4 threads
# or, if supported in your environment:
duckplyr::db_exec("RESET threads")      # back to DuckDB default
```

## References

- `duckplyr` - `dplyr` powered by DuckDB: https://duckplyr.tidyverse.org
  (CRAN: duckplyr).

- `db_exec` - https://duckplyr.tidyverse.org/reference/db_exec.html

- **DuckDB** configuration docs - thread/memory settings and more:
  https://duckdb.org/docs/stable/configuration/overview.html
