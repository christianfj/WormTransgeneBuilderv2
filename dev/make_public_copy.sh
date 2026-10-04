#!/usr/bin/env bash
# Makes a clean copy of the project for the public repository: only the files tracked in Git, minus the
# private working notes, and WITHOUT the .git directory (so the new repository starts with a fresh history).
#
# Usage (from the project root, with everything committed):
#   bash dev/make_public_copy.sh ../WormTransgeneBuilderv2-public
set -euo pipefail

target="${1:?give the folder to create, for example ../WormTransgeneBuilderv2-public}"
[ -e "$target" ] && { echo "$target already exists; remove it or choose another name." >&2; exit 1; }
[ -z "$(git status --porcelain)" ] || { echo "Uncommitted changes: commit them first." >&2; exit 1; }

# Private working notes and generated-template files that do not belong in the public repository.
exclude='^(CLAUDE\.md|deployment_notes\.txt|docs/ROADMAP\.md|\.here|\.rscignore|dev/01_start\.R|dev/02_dev\.R|dev/03_deploy\.R)$'

mkdir -p "$target"
git ls-files | grep -Ev "$exclude" | rsync -a --files-from=- . "$target"/
echo "Copied $(git ls-files | grep -Evc "$exclude") files to $target"

# Search the copy for things that must not be published.
echo "Searching for secrets and private details..."
patterns='BEGIN [A-Z ]*PRIVATE KEY|\.pem|amazonaws|ec2-[0-9]|/Users/|ghp_[A-Za-z0-9]{20,}|AKIA[0-9A-Z]{16}|password *=|api[_-]?key *='
if grep -rIn -E "$patterns" "$target" --exclude-dir=renv --exclude=renv.lock --exclude=.dockerignore --exclude=.gitignore \
     --exclude=test-package_config.R --exclude=make_public_copy.sh; then
  echo "Possible private details found above. Fix them in the project, commit, and run this again." >&2
  rm -rf "$target"; exit 1
fi
for f in .Rhistory .Renviron .env .DS_Store; do
  if find "$target" -name "$f" | grep -q .; then echo "Found $f in the copy." >&2; rm -rf "$target"; exit 1; fi
done
echo "No private details found."
