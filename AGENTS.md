# Repository Agent Instructions

These instructions apply to the entire ERF repository checkout.

## Repository and contribution destination

- Treat `https://github.com/sunt05/ERF.git` (`origin`) as the default repository for all branches, pushes, pull requests, issues, comments, reviews, and other GitHub write operations.
- Treat `https://github.com/erf-model/ERF.git` (`upstream`) as read-only by default. It may be fetched from and compared against, but agents must not create or modify upstream pull requests, issues, comments, reviews, releases, or branches unless Ting explicitly authorizes that exact upstream action in the current conversation.
- A request to create, prepare, or update a pull request means a pull request in `sunt05/ERF`, normally targeting its `development` branch. Do not infer permission to open an upstream pull request.
- Before every GitHub write operation, verify the repository target and use an explicit repository selector such as `--repo sunt05/ERF` where the tool supports it.
- If work may eventually be contributed upstream, first keep it in `sunt05/ERF`. Ask Ting before taking any outward-facing action in `erf-model/ERF`.

## Existing nested instructions

Instruction files under `Submodules/` apply only within their respective submodule trees and do not override this repository-level contribution destination policy.
