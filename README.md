# sc_analysis

Single-cell copy number analysis pipeline for processing scAbsolute and scUnique results.

## Quick Start

### Outlier Summary Pipeline

Generate and visualize outlier summaries from scAbsolute single-cell data.

```bash
# Step 1: Generate outlier summaries for each sample
Rscript scripts/01_generate_outlier_summaries.R \
    samples/ALLSAMPLES_sample.csv \
    /Volumes/Fl/sc_analysis/scAboslute-obj \
    /Volumes/Fl/sc_analysis/post-scAbsolute \
    100

# Step 2: Combine all summaries into one CSV with metadata
Rscript scripts/02_combine_outlier_summaries.R \
    samples/ALLSAMPLES_sample.csv \
    /Volumes/Fl/sc_analysis/all_samples_23july2024/sample_test \
    /Volumes/Fl/sc_analysis/post-scAbsolute \
    100

# Step 3: Generate visualization plots
Rscript scripts/03_visualize_summary.R \
    /Volumes/Fl/sc_analysis/all_samples_23july2024/sample_test/outlier_summary_table_meta_combined.csv \
    /Volumes/Fl/sc_analysis/all_samples_23july2024/sample_test \
    feature1 \
    Sample
```

See [scripts/README.md](scripts/README.md) for detailed pipeline documentation.

---

## Legacy Workflow

For each result from scAbsolute:
1. `Rscript scripts/01_generate_outlier_summaries.R samples/ALLSAMPLES_23355_6.csv`
2. GET_summary_of_dropouts
3. PLOT_dropouts

Then, for the results from scUnique...
