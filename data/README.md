# data/

Datasets shipped with this repository for the real-data experiments in the
paper.

## superconductivity.csv

Critical-temperature data for superconducting materials, used in the
high-dimensional real-data experiment (Section 5.3 of the paper).

| Field | Value |
|---|---|
| Source | UCI Machine Learning Repository, dataset ID 464 |
| URL | <https://archive.ics.uci.edu/dataset/464/superconductivty+data> |
| Citation | Hamidieh, K. (2018). *A data-driven statistical model for predicting the critical temperature of a superconductor.* Computational Materials Science, 154, 346–354. |
| License | CC BY 4.0 |
| Rows | 21,263 |
| Columns | 82 (81 features + `critical_temp` target) |
| Preprocessing | As shipped by UCI; no rows dropped, no columns transformed. The runners normalize features at use time via `apply(X, 2, normalize)`. |

The target column (last position, `critical_temp`) is treated as `Y` in the
real-data experiments; the remaining 81 columns are the conditioning vector
`X`.

## diamonds

The Diamonds dataset is loaded at runtime from the **ggplot2** R package
(`ggplot2::diamonds`, n = 53,940). It is therefore not stored in this
directory and is acquired automatically when `Rscript setup.R` installs
`ggplot2`.
