# Crop-specific SIF from overlapping OCO-2 target-mode soundings

Analysis code for:

> Jones, T. S., Dooley, J. F., Goebel, R., Moldavchuk, L., Ali, A., & Chen, J. *Resolving Crop-Specific Solar-Induced Fluorescence using Overlapping OCO-2 Target-Mode Soundings.* Submitted to *Geophysical Research Letters*.

The scripts in `code/` invert overlapping OCO-2 target-mode soundings over the Lamont, Oklahoma target onto a section-road grid and onto crop classes from the USDA Cropland Data Layer. They then compare the results with SIFnet, GOSIF and HLS NDVI, and test the response to drought.

Input data are not included. The sources are listed in the paper's Data Availability section. The downscaled and SIFnet data will be added here before publication.

## Requirements

R 4.3.1 with these packages: data.table, exactextractr, ggplot2, httr, jsonlite, lubridate, Matrix, ncdf4, patchwork, reshape2, sf, systemfonts, terra. The figure scripts use the Barlow font (Google Fonts), registered from a `fonts/` folder by `fig_fonts.R`.

## Running it

Every script sets its project folder near the top:

```r
proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
```

Some scripts also use this path directly in their `source()` calls. Put the project in that location, or find and replace the path in all scripts.

`code/run_all.sh` runs the main analysis in order, from subsetting the OCO-2 files through the inversions, comparisons and drought analysis. It does not run the figure scripts or the supporting analyses listed below.

### Downloads and credentials

No credentials are stored in these scripts. To download data you need:

- `subset_oco2_to_extent.R`: a NASA Earthdata login in `~/.netrc`. Version 11r files were read from a private storage bucket, so `bucket_dir` is a placeholder (`gs://YOUR-BUCKET/...`) in this copy. That is the only change from the code as run. Get the files from NASA GES DISC and point `bucket_dir` at your own copy.
- `hls_ndvi_lamont.R`: a NASA Earthdata login in `~/.netrc`.
- `fetch_crop_progress.R`: a USDA NASS Quick Stats API key in `~/.nass_api_key`.

## Figures and tables

| Item | Script | Notes |
|---|---|---|
| Figure 1 | `make_fig1_v3.R` | Uses outputs of `make_fig1_patches.R` and `crop_patch_sizes.R` |
| Figure 2 | `make_fig2.R` | Uses panels built by `make_fig1.R` |
| Figure 3 | `make_fig3.R` | Uses outputs of `bootstrap_climatology.R` |
| Figure 4 | `make_fig4.R` | Uses outputs of the comparison scripts |
| Figure S1 | `make_plots.R` | Crop rotation, from `crop_area_rotation.R` |
| Figure S2 | `make_fig_s1_vza.R` | VZA dependence, from `vza_check.R` |
| Figure S3 | `effective_resolution.R` | |
| Figure S4 | `make_plots.R` | |
| Figure S5 | `drought_analysis.R` | |
| Figure S6 | `make_fig5.R` | The output is named `fig5_manuscript` |
| Figure S7 | `drought_by_stage_usda.R` | |
| Table of target counts | `table_target_counts.R` | |

Some script names and comments use figure numbers from earlier drafts. This table has the final numbers.

## Scripts

**Main pipeline, in `run_all.sh` order:** `subset_oco2_to_extent.R`, `make_sub_road_crop.R`, `downscale_sif_inversion.R`, `add_crop_to_sif.R`, `downscale_sif_by_crop.R`, `downscale_sif_by_crop_cdl30.R`, `compare_sifnet.R`, `compare_gosif.R`, `drought_analysis.R`, `hls_ndvi_lamont.R`, `ndvi_by_crop.R`, `compare_ndvi.R`, `crop_area_rotation.R`, `make_plots.R`

**Shared helpers:** `crop_classes.R`, `vza_correction.R`, `crop_progress_calendar.R`, `active_season_bands.R`, `strip_colors.R`, `fig_fonts.R`

**Supporting analyses:**

- Inversion setup and sensitivity: `derive_corr_length.R`, `effective_resolution.R`, `gaussian_footprint_test.R`, `vza_check.R`, `compare_vza_correction.R`, `contrast_amplitude.R`, `table_target_counts.R`
- Crop mosaic: `crop_patch_sizes.R`
- Uncertainty: `bootstrap_climatology.R`, `bootstrap_drought.R`
- Drought and growth stage: `dsci_histogram.R`, `drought_by_stage.R`, `drought_by_stage_ndvi.R`, `drought_by_stage_usda.R`, `drought_stage_ndvi_on_targets.R`, `plot_stage_contrast.R`, `crop_condition_sif.R`, `fetch_crop_progress.R`

**Figures:** `make_fig1.R`, `make_fig1_patches.R`, `make_fig1_v3.R`, `make_fig2.R`, `make_fig3.R`, `make_fig4.R`, `make_fig5.R`, `make_fig_s1_vza.R`
