"""Every element id the page script writes to must exist in the HTML."""
import re
from pathlib import Path

HTML = (Path(__file__).resolve().parent.parent / "docs" / "index.html").read_text()


def test_script_targets_exist():
    ids = set(re.findall(r'id="([^"]+)"', HTML))
    used = set(re.findall(r"getElementById\('([^']+)'\)", HTML)) | set(re.findall(r"\$\('([^']+)'\)", HTML))
    dynamic = {"india-svg", "india-legend", "india-info"}   # created by initIndia
    assert not [u for u in used if u not in ids and u not in dynamic]


def test_tags_balanced():
    for tag in ("p", "ul", "li", "section", "div"):
        assert len(re.findall(rf"<{tag}[ >]", HTML)) == len(re.findall(rf"</{tag}>", HTML)), tag
