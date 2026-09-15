# Shared scientific data

This is the repository's only collaborator-facing binary-data namespace. NumPy,
HDF5, and pickle files below this directory are managed by Git LFS; generated
or private data elsewhere remain ignored.

## Lifecycle

- `final-inference/draft/` is mutable development data. Updates must include a
  regenerated `manifest.json` in the same commit.
- `final-inference/vN/` is an immutable paper/release snapshot. Never overwrite
  a frozen version; publish `vN+1` if a correction is necessary.
- A Git commit is the revision identifier during development. Do not create
  `v1` until the analysis, figure, and review are frozen.

## Update the draft

After rebuilding and validating local intermediates:

```bash
python scripts/publish_shared_data.py draft --force-draft
nix-shell -p git-lfs --run 'git lfs ls-files'
git add shared_data/final-inference/draft
```

Review `manifest.json` and ensure every binary is an LFS pointer before commit.
Unchanged content reuses its existing LFS object; changing a large array creates
another full LFS object, so avoid regenerating files unless their scientific
content changed.

## Freeze a paper release

```bash
python scripts/publish_shared_data.py v1
git add shared_data/final-inference/v1
git commit -m 'Freeze final inference dataset v1'
```

The publisher refuses to overwrite an existing frozen `vN` directory.

## Clone without downloading data

```bash
GIT_LFS_SKIP_SMUDGE=1 git clone <repository-url>
```

Fetch only the development inference dataset later:

```bash
git lfs install
git lfs pull --include='shared_data/final-inference/draft/**'
```
