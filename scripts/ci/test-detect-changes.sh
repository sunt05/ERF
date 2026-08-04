#!/usr/bin/env bash
# test-detect-changes.sh -- run ERF's CI path routing locally, without pushing.
#
# Runs the real .github/workflows/detect-changes.yml -- the same
# dorny/paths-filter action, reading the same .github/path-filters.yml -- inside
# Docker via nektos/act, so what you see here is what GitHub Actions decides.
#
# Usage:
#   ./scripts/ci/test-detect-changes.sh                 # HEAD vs origin/development
#   ./scripts/ci/test-detect-changes.sh main            # HEAD vs a named ref
#   ./scripts/ci/test-detect-changes.sh abc1234         # HEAD vs a specific commit
#
# Prerequisites:
#   - act     (brew install act)
#   - Docker running
#
# Note on tiers: a workflow_dispatch event counts as the integration tier, so a
# local run reports what a pull request would do. To see what a plain
# feature-branch push would do, read the "Areas touched" block and apply the
# feature-tier rules -- only Style, Linux GCC and the matching feature workflow
# run there.

set -euo pipefail

BASE="${1:-origin/development}"

if ! command -v act >/dev/null 2>&1; then
  echo "act is not installed. brew install act" >&2
  exit 1
fi

# Resolve to a SHA so the container never has to fetch from the remote: there
# are no credentials inside Docker.
if ! BASE_SHA="$(git rev-parse --verify "${BASE}^{commit}" 2>/dev/null)"; then
  echo "Cannot resolve '${BASE}' to a commit. Try: git fetch origin" >&2
  exit 1
fi

echo "Comparing HEAD ($(git rev-parse --short HEAD)) against ${BASE} (${BASE_SHA:0:7})"
echo "---"

ACT_ARGS=(
  workflow_dispatch
  -W .github/workflows/detect-changes.yml
  --input "base=${BASE_SHA}"
  --detect-event
  --container-architecture linux/amd64
  --pull=false
  # node:20 is enough: dorny/paths-filter is a JS action, and the remaining
  # steps are bash and git, both present in that image. The full ubuntu-22.04
  # runner toolchain is not needed to evaluate routing.
  --platform ubuntu-22.04=node:20
  --platform ubuntu-latest=node:20
)

# Git worktrees: .git is a FILE pointing at a directory outside this tree, which
# the container cannot see unless it is mounted at the same absolute path.
GIT_COMMON_DIR_RAW="$(git rev-parse --git-common-dir)"
if [[ "${GIT_COMMON_DIR_RAW}" != ".git" ]]; then
  GIT_COMMON_DIR="$(cd "${GIT_COMMON_DIR_RAW}" && pwd)"
  echo "Worktree detected, mounting: ${GIT_COMMON_DIR}"
  # --bind uses a host bind-mount rather than a copy, so the external .git
  # directory stays reachable.
  ACT_ARGS+=(--bind)
  ACT_ARGS+=(--container-options "-v ${GIT_COMMON_DIR}:${GIT_COMMON_DIR}")
fi

act "${ACT_ARGS[@]}"
