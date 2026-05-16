# Bioconda recipe staging area

This directory holds a draft bioconda recipe for BA3. Files here are not
consumed by the BA3 build itself; they are a staging area for submission
to the [bioconda-recipes](https://github.com/bioconda/bioconda-recipes)
repository.

## Layout

```
packaging/bioconda/
  ba3/
    meta.yaml      # package metadata, deps, version, hash
    build.sh      # compile + install script run by conda-build
```

## Submission flow

1. Fork `bioconda/bioconda-recipes` on GitHub.
2. Clone the fork, create a branch.
3. Copy this directory to the fork:
   ```
   cp -r packaging/bioconda/ba3 <bioconda-recipes>/recipes/
   ```
4. (Optional) Test the build locally before opening a PR:
   ```
   conda install -n base -c conda-forge conda-build mamba boa
   cd <bioconda-recipes>
   conda build recipes/ba3 -c bioconda -c conda-forge
   ```
5. Commit, push, open a PR against `bioconda/bioconda-recipes`. The
   bioconda CI will build the recipe on Linux x86_64, Linux aarch64, and
   macOS (Intel + ARM). A reviewer will merge once CI is green.

After the first merge, future version bumps usually only require updating
`version` and `sha256` in `meta.yaml`. Bioconda's bot
(`@BiocondaBot please bump-version`) can often do this automatically when
you tag a new GitHub release.

## When releasing a new BA3 version

1. Tag and publish the GitHub release (existing workflow).
2. Update this recipe:
   - bump `version` in `meta.yaml`
   - update `sha256` (compute with
     `curl -sL https://github.com/brannala/BA3/releases/download/vX.Y.Z/BA3-X.Y.Z.tar.gz | shasum -a 256`)
   - reset `build.number` to `0`
3. Open a PR in your `bioconda-recipes` fork with the bumped recipe.

## Notes on the recipe

- **Compiler**: bioconda's `{{ compiler('cxx') }}` provides `gxx` on Linux
  and `clangxx` on macOS. The Makefile pins `g++`, so `build.sh`
  overrides `CC` with `$CXX`.
- **Includes/libs**: the Makefile hardcodes `/usr/local`. `build.sh`
  overrides `INCLFLAGS` and `LIBFLAGS` to point at `$PREFIX`, where
  conda installs `gsl` and `htslib`.
- **htslib dependency**: this is what brings the static-linking pain on
  generic Linux. Bioconda handles it for us — `htslib` is a first-class
  bioconda package and its runtime deps are pinned automatically.
- **License**: declared as `AGPL-3.0-or-later` to match `include/BA3.h`.
  Note that the README in the repo currently says "GPL"; the actual
  source headers say AGPL. If you want to harmonize this someday it
  would be worth a separate cleanup.
