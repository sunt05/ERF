# CI routing

ERF's workflows used to fire on every push. A one-line edit to a design note
under `Source/LandSurfaceModel/Noah-MP/dev/` launched the full CUDA, HIP, SYCL,
Windows, Windows-MPI, macOS and Linux matrices — twelve workflows, none of which
could have been affected by the change.

Routing now answers a narrower question: what does *this* change actually need?

## The three pieces

- **`.github/path-filters.yml`** — the categories. Which paths count as core
  solver code, as the Noah-MP feature, as build system, as prose. This is the
  only place patterns are written down.
- **`.github/workflows/detect-changes.yml`** — the policy. A reusable workflow
  that runs `dorny/paths-filter` over those categories and turns them into
  per-workflow `run-*` decisions.
- **`scripts/ci/test-detect-changes.sh`** — the check. Runs that same workflow
  locally under `act`, so routing can be verified before pushing.

Every heavy workflow gains two lines and no policy of its own:

```yaml
jobs:
  changes:
    uses: ./.github/workflows/detect-changes.yml

  <existing-job>:
    needs: changes
    if: needs.changes.outputs.run-<area> == 'true'
```

Keeping the condition a single output reference is deliberate, and enforced:
`detect-changes.yml` exports only `run-*` answers, never the raw categories
behind them. There is deliberately no `needs.changes.outputs.core` to compose
against, because the moment a workflow starts writing its own
`core == 'true' || build == 'true'` the policy has been copied and the copies
drift. A new routing need gets a new `run-*` output. The raw categories remain
visible in each run's routing summary, which is where you want them when you are
working out *why* something was skipped.

## Two tiers

**Feature tier** — every push. Style checks, the baseline Linux GCC build, and
whichever feature workflow matches the paths touched. This is the fast signal
that a branch is not broken.

**Integration tier** — pull requests, pushes to `development`, and manual
dispatch. The full portability matrix: CUDA, HIP, SYCL, Windows, Windows-MPI,
macOS, the ERF CI matrix, and a Doxygen docs rebuild.

A feature workflow also runs at integration tier for *any* code change, not just
changes to its own files, because a core solver edit can break a feature build
without touching a single file the feature owns.

What that means in practice:

- Editing a Noah-MP design note: nothing builds. Docs only.
- Editing Noah-MP sources on a branch: style, Linux GCC, Noah-MP.
- Editing the core solver on a branch: style, Linux GCC.
- Opening the pull request: everything.

## Forcing a run

Every gated workflow accepts `workflow_dispatch`. A manual dispatch is treated
as "the human means it" and runs regardless of what changed — otherwise
dispatching a build on `development`, where there is no merge base to diff
against, would compare against `HEAD~1` and quietly skip.

## Checking routing before you push

```bash
./scripts/ci/test-detect-changes.sh                 # HEAD vs origin/development
./scripts/ci/test-detect-changes.sh <ref-or-sha>    # HEAD vs something else
```

Needs `act` and a running Docker. It executes the real workflow against the real
filter file, so it cannot drift from CI the way a reimplementation would.

## Two things worth knowing

**Feature branches are compared against their merge base with `development`, not
against the previous commit.** Comparing against the previous commit is cheaper,
but `cancel-in-progress` means a follow-up push can cancel the build that would
have covered the earlier source changes; if that follow-up then skips because it
only touched a README, the branch ends up with no green build at all. The cost is
that a docs-only push to a branch carrying source changes rebuilds.

**Never use negation patterns in `path-filters.yml`.** `dorny/paths-filter`
combines patterns with OR, so `!foo.cpp` matches every file that is not
`foo.cpp` and silently turns its category into a catch-all
([dorny/paths-filter#184](https://github.com/dorny/paths-filter/issues/184)).
Where prose lives inside a source tree, match code extensions rather than the
directory — that excludes the prose without a negation.
