"""Assemble the root Pages build with the frozen 0.9.1 site and old URLs."""

import html
import json
from pathlib import Path
import re
import sys
import tarfile
from urllib.parse import quote


ARCHIVE_PATH = Path("archive/0.9.1")
SNAPSHOT = Path(__file__).with_name("immunarch-0.9.1-pkgdown.tar.gz")
DOMAIN = "https://immunarch.com"


def redirect_page(destination):
    url = "/" + quote(destination.as_posix(), safe="/")
    escaped_url = html.escape(url, quote=True)
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Documentation moved</title>
  <link rel="canonical" href="{DOMAIN}{escaped_url}">
  <script>
    window.location.replace({json.dumps(url)} + window.location.search + window.location.hash);
  </script>
  <noscript><meta http-equiv="refresh" content="0; url={escaped_url}"></noscript>
</head>
<body><a href="{escaped_url}">Open the documentation</a></body>
</html>
"""


def redirects(site, archive):
    routes = {}
    for section in ("reference", "articles"):
        for page in sorted((archive / section).rglob("*.html")):
            old_path = page.relative_to(archive)
            routes[old_path] = ARCHIVE_PATH / old_path

    original_articles = archive / "articles/web_only"
    for page in sorted(original_articles.rglob("*.html")):
        later_path = Path("articles/web_only_v0") / page.relative_to(original_articles)
        routes[later_path] = page.relative_to(site)
    routes[Path("articles/web_only_v0/introduction.html")] = ARCHIVE_PATH / "articles/v1_introduction.html"
    routes[Path("articles/web_only_v0/data.html")] = ARCHIVE_PATH / "articles/v2_data.html"

    for page in sorted((site / "vignettes").glob("*.html")):
        routes.setdefault(Path("articles") / page.name, page.relative_to(site))

    # altdoc renders one page per Rd topic. pkgdown also served each alias.
    for page in sorted((site / "man").glob("*.html")):
        routes.setdefault(Path("reference") / page.name, page.relative_to(site))
    for rd in sorted(Path(__file__).resolve().parents[2].joinpath("man").glob("*.Rd")):
        topic = Path("man") / f"{rd.stem}.html"
        if (site / topic).is_file():
            for alias in re.findall(r"^\\alias\{([A-Za-z0-9_.-]+)\}", rd.read_text(encoding="utf-8"), re.M):
                routes.setdefault(Path("reference") / f"{alias}.html", topic)

    for old, aliases, current in (
        ("airr_stats", ("chains", "genes", "lengths"), "airr_desc"),
        ("airr_public", ("intersection", "jaccard"), "repsim"),
    ):
        topic = Path("man") / f"{current}.html"
        for name in (old, *(f"{old}_{alias}" for alias in aliases)):
            routes.setdefault(Path("reference") / f"{name}.html", topic)

    for name in ("im_common_args", "register_immunarch_method"):
        routes.setdefault(Path("reference") / f"{name}.html", Path("reference/index.html"))

    return routes


def assemble(site):
    if not (site / "index.html").is_file():
        raise ValueError(f"Build the altdoc site before assembling {site}.")
    archive = site / ARCHIVE_PATH
    if archive.exists():
        raise ValueError(f"Archive directory already exists in {site}.")
    with tarfile.open(SNAPSHOT, "r:gz") as frozen:
        frozen.extractall(archive, filter="data")

    for expected in ("index.html", "reference/repLoad.html", "articles/v1_introduction.html"):
        if not (archive / expected).is_file():
            raise ValueError(f"The 0.9.1 snapshot is missing {expected}.")

    routes = redirects(site, archive)
    # Altdoc owns the current reference index when it publishes that URL.
    if (site / "reference/index.html").is_file():
        routes.pop(Path("reference/index.html"), None)

    for old_path, destination in routes.items():
        if not (site / destination).is_file():
            raise ValueError(f"Missing redirect target for {old_path}: {destination}")
        if (site / old_path).exists():
            raise ValueError(f"Redirect would overwrite an existing page: {old_path}")

    for old_path, destination in routes.items():
        page = site / old_path
        page.parent.mkdir(parents=True, exist_ok=True)
        page.write_text(redirect_page(destination), encoding="utf-8")

    (site / "CNAME").write_text("immunarch.com\n", encoding="utf-8")
    (site / ".nojekyll").touch()
    print(f"Added the frozen 0.9.1 archive and {len(routes)} redirects.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python3 .github/pages/prepare_pages.py <altdoc-output-dir>")
    assemble(Path(sys.argv[1]).resolve())
