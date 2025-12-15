FROM --platform=linux/amd64 rocker/verse:4.4.1

# Install system dependencies required for R packages
# libgsl-dev: for gsl/energy
# cmake: for xgboost (sometimes needed)
# libglpk-dev: for optimization packages
# libxt-dev: for some plotting
RUN apt-get update && apt-get install -y \
    libgsl-dev \
    cmake \
    libglpk-dev \
    libxt-dev \
    && rm -rf /var/lib/apt/lists/*

# Install renv
ENV RENV_VERSION 1.0.7
RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"

# Create app directory
WORKDIR /home/rstudio/ml4time2event

# Copy lockfile and DESCRIPTION
COPY renv.lock DESCRIPTION ./

# Use Posit Package Manager for Linux binaries to speed up installation
# This overrides the lockfile cran mirror with a binary-capable one for Ubuntu Jammy (which rocker/verse:4.4.1 is based on)
ENV RENV_CONFIG_REPOS_OVERRIDE https://packagemanager.posit.co/cran/__linux__/jammy/latest

# Restore dependencies
# We set RENV_PATHS_LIBRARY to a location inside the image so it doesn't conflict with bind mounts if any
ENV RENV_PATHS_LIBRARY /home/rstudio/renv/library
RUN R -e "renv::restore(prompt = FALSE)"

# Copy the rest of the package source
COPY . .

# Install the package
# We use remotes::install_local which should respect the installed dependencies
RUN R -e "devtools::install('.', dependencies = FALSE, upgrade = 'never', build_vignettes = TRUE)"

# Expose port for RStudio if we switch to that, but default is just R or RStudio
# rocker/verse starts RStudio by default on port 8787 if no command is given?
# Actually rocker/verse inherits from rocker/rstudio.
# So if we just let it run, it will start RStudio Server.
# User wants "dockerized version... as its working form".
# RStudio is a great working form.

EXPOSE 8787

# The default command of rocker/verse is `/init` which starts S6 overlay -> RStudio.
# We can keep that.
