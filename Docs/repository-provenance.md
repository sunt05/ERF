# Repository provenance and contribution policy

## Current status

`sunt05/ERF` is Ting Sun's standalone experimental ERF repository. On
2026-08-05 it left the GitHub fork network rooted at `erf-model/ERF`. GitHub
reports the repository as a non-fork with no parent or source repository. The
detachment is permanent: this repository cannot be reattached to the original
fork network.

The default branch is `development`. Local tooling is configured so that:

- `origin` is `https://github.com/sunt05/ERF.git` and is the default target
  for branches, pushes, pull requests, issues, releases, and other writes;
- `upstream` fetches from `https://github.com/erf-model/ERF.git` for
  comparison and selective synchronization;
- `upstream` has no usable push URL in this checkout; and
- GitHub CLI resolves the repository to `sunt05/ERF` by default.

This repository and its results are not official ERF work and do not imply
upstream review, approval, or endorsement.

## Preserved experimental patch series

The following work is intentionally retained in this repository:

- `sunt05/ERF#1`: SLUCM as a mosaic urban tile, merged into
  `development`;
- `sunt05/ERF#2`: file-driven surface flux forcing and GPU/boundary
  robustness, superseding the closed upstream draft
  `erf-model/ERF#3530`, integrated into `development` on 2026-08-05;
- `sunt05/ERF#3`: repository contribution guidance, merged into
  `development`; and
- `sunt05/ERF#4`: feature-based CI routing, integrated into `development` on
  2026-08-05.

These patch series are maintained for the experimental repository as coherent
local capabilities. Their presence here is not a commitment to submit them
upstream in their current form.

## Possible future upstream contribution

Any future upstream contribution is a separate, explicit decision. The
expected workflow is:

1. review the experimental changes against the then-current
   `erf-model/ERF:development`;
2. split or reduce the work into an upstream-appropriate patch series;
3. validate that series independently of fork-only features;
4. create a separate contribution fork, for example
   `sunt05/ERF-contrib`, because this standalone repository cannot rejoin the
   fork network; and
5. open an upstream pull request only after Ting explicitly authorizes that
   particular contribution.

Until then, all development and pull requests remain in `sunt05/ERF`.
