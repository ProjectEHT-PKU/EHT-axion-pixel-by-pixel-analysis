"""Documentation structure: reachability from README, relative links and anchors, back links."""
from pathlib import Path
import re
from urllib.parse import unquote

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs" / "en"
INDEX = DOCS / "index.md"
README = ROOT / "README.md"
PAGES = sorted(DOCS.glob("*.md"))

_FENCE = re.compile(r"^(```|~~~)")
_LINK = re.compile(r"(?<!!)\[[^\]]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")
_HEADING = re.compile(r"^(#{1,6})\s+(.*?)\s*#*\s*$")
_EXTERNAL = re.compile(r"^[a-z][a-z0-9+.-]*:", re.IGNORECASE)


def _prose_lines(path):
    """Lines outside fenced code blocks."""
    inside = False
    for line in path.read_text(encoding="utf-8").splitlines():
        if _FENCE.match(line.strip()):
            inside = not inside
            continue
        if not inside:
            yield line


def _slug(text):
    """GitHub-style heading anchor."""
    text = re.sub(r"`([^`]*)`", r"\1", text)
    text = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", text)
    text = text.strip().lower()
    text = re.sub(r"[^\w\- ]", "", text)
    return text.replace(" ", "-")


def anchors(path):
    result, counts = set(), {}
    for line in _prose_lines(path):
        match = _HEADING.match(line)
        if match:
            slug = _slug(match[2])
            n = counts.get(slug, 0)
            counts[slug] = n + 1
            result.add(slug if n == 0 else f"{slug}-{n}")
    return result


def links(path):
    """Relative link targets of a Markdown file (inline code spans removed)."""
    found = []
    for line in _prose_lines(path):
        for target in _LINK.findall(re.sub(r"`[^`]*`", "", line)):
            if not _EXTERNAL.match(target):
                found.append(unquote(target))
    return found


def resolve(source, target):
    file_part, _, anchor = target.partition("#")
    destination = source if not file_part else (source.parent / file_part).resolve()
    return destination, anchor


def test_documentation_pages_exist():
    assert INDEX in PAGES and len(PAGES) > 1
    assert all(re.fullmatch(r"[a-z0-9]+(-[a-z0-9]+)*\.md", p.name) for p in PAGES)
    assert not any(re.match(r"\d", p.name) for p in PAGES), "no numeric prefixes"


@pytest.mark.parametrize("path", [README, *PAGES], ids=lambda p: str(p.relative_to(ROOT)))
def test_relative_links_and_anchors_resolve(path):
    for target in links(path):
        destination, anchor = resolve(path, target)
        assert destination.exists(), f"{path.name}: broken link {target}"
        if anchor:
            assert destination.suffix == ".md", f"{path.name}: anchor on non-Markdown {target}"
            assert anchor in anchors(destination), f"{path.name}: missing anchor {target}"


def test_every_page_is_reachable_from_readme():
    seen, pending = set(), [README]
    while pending:
        page = pending.pop()
        if page in seen:
            continue
        seen.add(page)
        for target in links(page):
            destination, _ = resolve(page, target)
            if destination.suffix == ".md" and destination.is_file():
                pending.append(destination)
    assert INDEX in seen
    assert not [p.name for p in PAGES if p not in seen], "orphan documentation pages"


@pytest.mark.parametrize("path", PAGES, ids=lambda p: p.name)
def test_page_structure_and_back_link(path):
    headings = [line for line in _prose_lines(path) if _HEADING.match(line)]
    assert headings and headings[0].startswith("# "), "page must start with its H1 title"
    assert sum(h.startswith("# ") for h in headings) == 1, "exactly one H1 title"
    assert headings[1] == "## Purpose", "the first section is Purpose"
    assert headings[-1] == "## Related pages", "the last section is Related pages"
    levels = [len(_HEADING.match(h)[1]) for h in headings]
    assert all(b <= a + 1 for a, b in zip(levels, levels[1:])), "no skipped heading levels"
    if path != INDEX:
        assert any(resolve(path, t)[0] == INDEX for t in links(path)), "missing link back to index.md"
