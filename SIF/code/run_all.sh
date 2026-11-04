#!/bin/sh
# Rerun the whole analysis in order, after new OCO-2 data or CDL years.
# Each step's log goes to logs/<step>.log; stops at the first failure.
#   sh run_all.sh            # from the project folder
#
# Order matters: the crop columns are added to the per-cell GeoPackage after
# the per-cell inversion rewrites it, and the 740 nm crop inversion feeds the
# SIFNet comparison; the drought step writes the drought index compare_ndvi.R uses.
set -e
cd "$(dirname "$0")"
mkdir -p logs

step() {
  name=$1; shift
  echo "$(date +%H:%M) $name"
  Rscript "$@" > "logs/$name.log" 2>&1 || { echo "  FAILED - see logs/$name.log"; exit 1; }
}

step subset_oco2           subset_oco2_to_extent.R    # new OCO-2 days only (resumes)
step sub_road_crop         make_sub_road_crop.R       # crop maps for new CDL years
step inversion_cells       downscale_sif_inversion.R
step add_crop              add_crop_to_sif.R
step inversion_crop_grid   downscale_sif_by_crop.R
step inversion_crop_757    downscale_sif_by_crop_cdl30.R 757
step inversion_crop_740    downscale_sif_by_crop_cdl30.R 740
step compare_sifnet        compare_sifnet.R
step compare_gosif         compare_gosif.R            # after compare_sifnet (joint summary)
step drought               drought_analysis.R
step hls_ndvi              hls_ndvi_lamont.R          # new HLS scenes only (resumes)
step ndvi_by_crop          ndvi_by_crop.R
step compare_ndvi          compare_ndvi.R
step crop_area_rotation    crop_area_rotation.R
step plots                 make_plots.R
echo "$(date +%H:%M) done"
