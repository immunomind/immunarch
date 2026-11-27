FROM r-base

# Use the sid suite for consistent versions across system dependencies
RUN apt-get update && apt-get install -y -t sid --no-install-recommends \
    build-essential file libcurl4-openssl-dev libcairo2-dev libxml2-dev libssl-dev \
    libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    libnlopt-dev \
    && rm -rf /var/lib/apt/lists/*

# Install the package installer and alluvial plotting support
RUN R -e "install.packages(c('remotes', 'ggalluvial'))"

# Use Posit binaries for DuckDB on supported architectures
RUN R -e 'install.packages("DBI"); \
    repos <- getOption("repos"); \
    arch <- R.version$arch; \
    if (arch %in% c("x86_64", "aarch64")) { \
      r_series <- paste(R.version$major, sub("[.].*$", "", R.version$minor), sep = "."); \
      repos <- sprintf("https://p3m.dev/cran/latest/bin/linux/manylinux_2_28-%s/%s", arch, r_series); \
    }; \
    install.packages("duckdb", repos = repos, dependencies = FALSE); \
    stopifnot(requireNamespace("duckdb", quietly = TRUE))'

# Copy source files to the image
COPY . /immunarch-src/

# Install Immunarch from source
RUN R -e 'remotes::install_local("/immunarch-src", dependencies = NA, upgrade = "never")'

# Delete Immunarch source from the image
RUN rm -rf /immunarch-src

# Check that immunarch loads
RUN R -e 'library("immunarch")'
