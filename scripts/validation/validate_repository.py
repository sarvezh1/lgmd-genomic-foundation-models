"""Validate the public scientific repository without network access."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TEXT_SUFFIXES = {".py", ".md", ".json", ".tsv", ".txt", ".cff", ".yaml", ".yml"}
TOP_LEVEL_FILES = {
    ".gitattributes", ".gitignore", "CITATION.cff", "LICENSE", "METHODS.md",
    "README.md", "REPRODUCIBILITY.md", "THIRD_PARTY_NOTICES.md",
    "requirements.txt",
}
TOP_LEVEL_DIRECTORIES = {
    "data", "figures", "licenses", "metadata", "schemas", "scripts", "tables",
    "tests",
}
IGNORED_DIRECTORY_PARTS = {".git", "verification_outputs", "__pycache__"}
IGNORED_BYTECODE_SUFFIXES = {".pyc", ".pyo"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def repository_files(root: Path):
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        relative = path.relative_to(root)
        lower_parts = {part.lower() for part in relative.parts}
        if lower_parts & IGNORED_DIRECTORY_PARTS:
            continue
        if path.suffix.lower() in IGNORED_BYTECODE_SUFFIXES:
            continue
        if path.name.lower() == "desktop.ini":
            continue
        yield path, relative


def validate(root: Path = ROOT) -> tuple[list[str], dict[str, object]]:
    root = root.resolve()
    errors: list[str] = []
    checks: dict[str, object] = {}

    entries = {
        path.name
        for path in root.iterdir()
        if path.name.lower() != "desktop.ini"
        and path.name.lower() not in IGNORED_DIRECTORY_PARTS
        and path.suffix.lower() not in IGNORED_BYTECODE_SUFFIXES
    }
    unexpected = sorted(entries - TOP_LEVEL_FILES - TOP_LEVEL_DIRECTORIES)
    missing = sorted(TOP_LEVEL_FILES - entries)
    if unexpected:
        errors.append(f"unexpected top-level entries: {unexpected}")
    if missing:
        errors.append(f"missing top-level files: {missing}")

    gitignore_rules = {
        line.strip()
        for line in (root / ".gitignore").read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }
    required_cache_rules = {"__pycache__/", "*.py[cod]"}
    missing_cache_rules = sorted(required_cache_rules - gitignore_rules)
    if missing_cache_rules:
        errors.append(f"missing Python cache ignore rules: {missing_cache_rules}")

    manifest = json.loads((root / "metadata/evidence_manifest.json").read_text(encoding="utf-8"))
    evidence = list(manifest["files"].values()) + list(manifest["gene_heldout"]["files"].values())
    bad_evidence = []
    for item in evidence:
        path = root / item["path"]
        if not path.is_file() or sha256(path) != item["sha256"]:
            bad_evidence.append(item["path"])
    if bad_evidence:
        errors.append(f"evidence hash failures: {bad_evidence}")

    canonical = read_tsv(root / "data/canonical_variants.tsv")
    configurations = read_tsv(root / "data/configurations.tsv")
    primary = read_tsv(root / "data/predictions_primary.tsv")
    sensitivity = read_tsv(root / "data/predictions_sensitivity.tsv")
    gene_predictions = read_tsv(root / "data/gene_heldout_predictions.tsv")
    counts = {
        "variants": len(canonical),
        "train": sum(row["split"] == "train" for row in canonical),
        "test": sum(row["split"] == "test" for row in canonical),
        "primary_configurations": sum(row["primary"] == "true" for row in configurations),
        "sensitivity_configurations": sum(row["sensitivity"] == "true" for row in configurations),
        "fixed_predictions": len(primary) + len(sensitivity),
        "gene_heldout_predictions": len(gene_predictions),
    }
    expected = {
        "variants": 6960, "train": 4866, "test": 2094,
        "primary_configurations": 26, "sensitivity_configurations": 2,
        "fixed_predictions": 58632, "gene_heldout_predictions": 16752,
    }
    if counts != expected:
        errors.append(f"scientific count mismatch: {counts}")
    if {row["locus_id"] for row in canonical if row["split"] == "train"} & {
        row["locus_id"] for row in canonical if row["split"] == "test"
    }:
        errors.append("canonical train/test locus overlap")
    checks["scientific_counts"] = counts

    source_rows = read_tsv(root / "metadata/source_to_result.tsv")
    if len(source_rows) != 452 or any(row["verification_status"] == "FAIL" for row in source_rows):
        errors.append("source-to-result traceability is incomplete")
    checks["source_to_result"] = len(source_rows)

    assets = read_tsv(root / "metadata/approved_asset_checksums.tsv")
    bad_assets = []
    for row in assets:
        path = root / row["path"]
        if (
            not path.is_file()
            or path.stat().st_size != int(row["bytes"])
            or sha256(path) != row["sha256"]
        ):
            bad_assets.append(row["path"])
    if len(assets) != 14 or bad_assets:
        errors.append(f"approved asset failure: rows={len(assets)}, bad={bad_assets}")
    checks["approved_assets"] = len(assets)

    allowed_markdown = {
        "README.md", "METHODS.md", "REPRODUCIBILITY.md", "THIRD_PARTY_NOTICES.md"
    }
    markdown = {
        relative.as_posix()
        for path, relative in repository_files(root)
        if path.suffix.lower() == ".md"
    }
    if markdown != allowed_markdown:
        errors.append(f"unexpected Markdown set: {sorted(markdown)}")

    forbidden_paths = [
        "g:" + "\\" + "my drive", "c:" + "\\" + "users",
        "github" + "_repo", "public_release_" + "candidate",
        "approved_" + "candidate2", "results/" + "phase",
        "data/" + "v1/", "mani" + "fests/", "re" + "ports/", "repro" + "duced/",
    ]
    forbidden_names = {
        "desktop.ini", ".env", "credentials.json", "token.json", "id_rsa",
    }
    forbidden_parts = {".pytest_cache", ".mypy_cache", ".ruff_cache"}
    forbidden_suffixes = {".ipynb", ".bak", ".old", ".orig", ".tmp", ".zip", ".7z", ".tar"}
    files = []
    for path, relative in repository_files(root):
        files.append(path)
        lower_parts = {part.lower() for part in relative.parts}
        if path.name.lower() in forbidden_names:
            errors.append(f"prohibited file: {relative.as_posix()}")
        if lower_parts & forbidden_parts:
            errors.append(f"cache path: {relative.as_posix()}")
        if path.suffix.lower() in forbidden_suffixes:
            errors.append(f"prohibited suffix: {relative.as_posix()}")
        if path.suffix.lower() in TEXT_SUFFIXES:
            text = path.read_text(encoding="utf-8-sig", errors="replace").lower()
            for marker in forbidden_paths:
                if marker in text:
                    errors.append(f"internal or obsolete path: {relative.as_posix()}")

    checks["files"] = len(files)
    checks["bytes"] = sum(path.stat().st_size for path in files)
    checks["markdown_files"] = len(markdown)
    checks["status"] = "PASS" if not errors else "FAIL"
    return sorted(set(errors)), checks


def main() -> int:
    errors, checks = validate()
    print(json.dumps({"errors": errors, "checks": checks}, indent=2))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
