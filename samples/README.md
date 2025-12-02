This folder contains example sample lists for the outlier summary pipeline.

## CSV Format

**Required columns:**
- `sample`: numeric or string sample identifier (e.g., `25393`, `23003`)

**Optional columns:**
- `category` or `Cell line`: Cell line name (e.g., `PEO1`, `HCT116`)
- `feature1`: First feature annotation (e.g., `HRD`, `HRP`)
- `feature2`: Second feature annotation (e.g., `Cisplatin sensitive`)

## Example Files

- `samples_summary_outliers.csv` - Basic sample list
- `ALLSAMPLES_metadata.csv` - Full sample list with metadata

## Usage

Run the pipeline scripts with a sample CSV:

```bash
# Step 1: Generate outlier summaries
Rscript scripts/01_generate_outlier_summaries.R samples/my_samples.csv

# Step 2: Combine summaries with metadata
Rscript scripts/02_combine_outlier_summaries.R samples/my_samples.csv /output/path /data/path 100

# Step 3: Visualize results
Rscript scripts/03_visualize_outlier_summary.R /output/path/outlier_summary_table_meta_combined.csv "" /output/path
```

See `scripts/README.md` for detailed pipeline documentation.

## Notes

- Scripts expect sample objects at `obj_base/SLX-<sample>_<bin_size>.rds`
- Adjust `obj_base` and `out_base` in scripts if your paths differ
