# Dockerized ml4time2event

This directory contains the configuration to run `ml4time2event` in a Docker container.

## Prerequisites

- Docker
- Docker Compose

## Quick Start

1. **Build and Run**:
   Run the following command to build the image and start the container:

   ```bash
   docker-compose up --build
   ```

2. **Access RStudio**:
   - Open your browser and go to [http://localhost:8787](http://localhost:8787).
   - **Username**: `rstudio`
   - **Password**: `ml4time2event` (configured in `docker-compose.yml`)

## Directory Structure

- The current directory is mounted to `/home/rstudio/ml4time2event` inside the container.
- Changes made to files in this directory are reflected inside the container instantly.
- **Note**: The package is installed into the container's system library or `renv` library. If you modify code and want to reinstall the package *within* the container, run:
  ```r
  devtools::install()
  ```
  inside the RStudio console.

## Dependencies

The Docker image uses `renv` to restore dependencies from `renv.lock`. The `renv` library folder inside the container is stored separately from your local `renv/library` to avoid conflicts between operating systems (macOS vs Linux).

If you add new dependencies:
1. Run `renv::snapshot()` inside the container.
2. The `renv.lock` file will be updated on your host machine (since the directory is mounted).
3. Rebuild the image if you want the dependencies baked in: `docker-compose build`.

## Troubleshooting

- **Permissions**: If you encounter permission issues with files created inside the container, check the `ROOT=true` environment variable or user ID mapping.
- **xgboost version**: If `renv::restore()` fails on `xgboost` due to version mismatch or compilation errors, verify `renv.lock` or try installing it manually in the container.
