# immunarch documentation deployment

The default branch publishes altdoc at `immunarch.com/`. The `dev` branch publishes altdoc at `immunarch.com/dev/`.
Both builds read the package version from `DESCRIPTION`.

The root build copies a frozen pkgdown site to `/archive/0.9.1/`.
It adds redirects from the old reference and article URLs.
The redirects preserve query strings and page fragments when JavaScript runs.

## First deployment

1. In **Settings → Pages**, set **Source** to **GitHub Actions** and keep the custom domain `immunarch.com`.
2. In the `github-pages` environment, allow the default branch and `dev` to deploy if branch restrictions apply.
3. Push the current package and documentation changes to the default branch. Check the root version and `/reference/repLoad.html#arguments`.
4. Push the same changes to `dev`. Check the preview at `/dev/`.

The workflow stages both sites on `gh-pages`. Each deployment uploads the full tree to GitHub Pages.
The root build retains `dev/`. The dev build updates only `dev/`.

## Frozen archive

`immunarch-0.9.1-pkgdown.tar.gz` is a one-time render of tag `0.9.1`, commit `0b6544a101e9f42ee1bf6a5f14cbd234147e18ad`.
It contains all 18 articles, 82 reference topics, and their computed figures. Help aliases link to their original topics.
The render used R 4.5.3, pkgdown 2.0.7, Clustal W 2.1, and PHYLIP 3.697.
The snapshot uses the archive URL and links its sources to tag `0.9.1`.
The workflow extracts the same files on each root build. It never runs pkgdown.

SHA-256 `d61684c64a42684105cea7b1f54a8c11e348f0892264e043c3164e2f64f26776`
