# Post-scUnique visualization pipeline

## What it does

For one completed scUnique sample, this pipeline produces:

1. A final-copy-number heatmap ordered by the scUnique evolutionary tree.
2. A histogram and sorted per-cell plot of validated private events.
3. Per-cell event counts and a concise run summary.

The tree receives a wide dedicated margin. Cell labels are shortened by
automatically removing their shared prefix; full names remain in the CSV.

## Run

From the repository root:

```bash
./run_post_scunique.sh \
    <scunique_result_dir> \
    [output_dir] \
    [prefix]
```

If `output_dir` is omitted, results go to
`<scunique_result_dir>/post_scunique/`. If the input directory contains exactly
one `*.finalCN.RDS`, the prefix is inferred automatically.

## Required inputs

The result directory must contain matching files:

- `<prefix>.finalCN.RDS`
- `<prefix>.tree.RDS`
- `<prefix>.df_pass_post.RDS`

## Meaning of `freq == 1`

scUnique calculates `freq` by grouping validated events by chromosome, start,
and end, then counting how many cells carry each locus. Therefore `freq == 1`
means the event is private to one cell within the analyzed cohort. A single cell
can carry zero, one, or many different private events.

## Outputs

- `<prefix>_tree_cn_heatmap_labeled.pdf` and `.png`.
- `<prefix>_freq1_unique_events_distribution.pdf` and `.png`.
- `<prefix>_freq1_unique_events_per_cell.csv` with full and shortened names.
- `<prefix>_post_scunique_summary.csv` with cohort/event totals and the removed
  label prefix.

The PDF heatmap is the preferred version for inspecting cell labels and tree
branches at full resolution.
