#!/usr/bin/env python3
"""
Validate the public code repository before manuscript-associated release.

This script performs lightweight repository checks only. It does not require
the original imaging, transcriptomic, or clinical datasets.
"""

from __future__ import annotations

import json
import py_compile
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

REQUIRED_FILES = [
    "README.md",
    "CITATION.cff",
    ".zenodo.json",
    "CODE_AVAILABILITY.md",
    "LICENSE",
    "MANIFEST.tsv",
]

REQUIRED_DIRS = [
    "R",
    "python",
    "docs",
    "data",
    "environment",
    "tools",
    "config",
]

TEXT_EXTENSIONS = {
    ".R",
    ".r",
    ".py",
    ".md",
    ".txt",
    ".yml",
    ".yaml",
    ".json",
    ".cff",
    ".tsv",
    ".csv",
    ".ini",
    ".cfg",
}

SKIP_DIRS = {
    ".git",
    ".github",
    "__pycache__",
    ".pytest_cache",
    ".Rproj.user",
}

PLACEHOLDERS = [
    "YOUR_USERNAME",
    "OWNER/REPOSITORY",
    "FAMILY_NAME",
    "GIVEN_NAMES",
    "AFFILIATION",
    "0000-0000-0000-0000",
    "ZENODO_VERSION_DOI",
]

FORBIDDEN_SUFFIXES = {
    ".RData",
    ".RDataTmp",
    ".Rhistory",
    ".pyc",
    ".pyo",
}

FORBIDDEN_NAME_PARTS = [
    "test",
    "tmp",
    "temp",
    "backup",
    "old",
    "debug",
]

LOCAL_PATH_PATTERNS = [
    "C:" + "/Users/",
    "C:" + "\\Users\\",
    "D:" + "/",
    "D:" + "\\",
    "/Users/",
    "/home/",
]


def fail(messages: list[str]) -> None:
    print("Repository validation failed:")
    for message in messages:
        print(f" - {message}")
    sys.exit(1)


def iter_repository_files() -> list[Path]:
    files: list[Path] = []
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(ROOT)
        if any(part in SKIP_DIRS for part in rel.parts):
            continue
        files.append(rel)
    return files


def read_text(path: Path) -> str | None:
    try:
        return path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return None


def has_cjk(text: str) -> bool:
    return any(
        "\u4e00" <= ch <= "\u9fff"
        or "\u3400" <= ch <= "\u4dbf"
        or "\uf900" <= ch <= "\ufaff"
        for ch in text
    )


def check_required_items(errors: list[str]) -> None:
    for item in REQUIRED_FILES:
        if not (ROOT / item).is_file():
            errors.append(f"Missing required file: {item}")

    for item in REQUIRED_DIRS:
        if not (ROOT / item).is_dir():
            errors.append(f"Missing required directory: {item}")


def check_json_metadata(errors: list[str]) -> None:
    zenodo_path = ROOT / ".zenodo.json"

    if not zenodo_path.is_file():
        return

    try:
        metadata = json.loads(zenodo_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        errors.append(f".zenodo.json is not valid JSON: {exc}")
        return

    for key in ["title", "description", "upload_type", "license", "creators"]:
        if key not in metadata:
            errors.append(f".zenodo.json is missing key: {key}")

    creators = metadata.get("creators", [])
    if not isinstance(creators, list) or len(creators) == 0:
        errors.append(".zenodo.json must contain at least one creator")


def check_placeholders(errors: list[str]) -> None:
    targets = [
        ROOT / "README.md",
        ROOT / "CITATION.cff",
        ROOT / ".zenodo.json",
        ROOT / "CODE_AVAILABILITY.md",
    ]

    for path in targets:
        if not path.is_file():
            continue
        text = read_text(path)
        if text is None:
            continue
        for placeholder in PLACEHOLDERS:
            if placeholder in text:
                errors.append(f"Placeholder remains in {path.name}: {placeholder}")


def check_cjk_characters(errors: list[str]) -> None:
    for rel in iter_repository_files():
        if rel.suffix not in TEXT_EXTENSIONS:
            continue
        text = read_text(ROOT / rel)
        if text is None:
            continue
        if has_cjk(text):
            errors.append(f"CJK characters detected in text/code file: {rel}")


def check_local_paths(errors: list[str]) -> None:
    for rel in iter_repository_files():
        if rel == Path("tools/validate_repository.py"):
            continue
        if rel.suffix not in TEXT_EXTENSIONS:
            continue
        text = read_text(ROOT / rel)
        if text is None:
            continue
        for pattern in LOCAL_PATH_PATTERNS:
            if pattern in text:
                errors.append(f"Local absolute path detected in {rel}: {pattern}")
                break


def check_forbidden_files(errors: list[str]) -> None:
    for rel in iter_repository_files():
        if rel.suffix in FORBIDDEN_SUFFIXES:
            errors.append(f"Forbidden generated file detected: {rel}")
            continue

        lowered = rel.name.lower()
        if rel.suffix in {".R", ".r", ".py"}:
            for token in FORBIDDEN_NAME_PARTS:
                if token in lowered:
                    errors.append(f"Potential temporary script name detected: {rel}")
                    break


def check_python_syntax(errors: list[str]) -> None:
    for rel in iter_repository_files():
        if rel.suffix != ".py":
            continue
        try:
            py_compile.compile(str(ROOT / rel), doraise=True)
        except py_compile.PyCompileError as exc:
            errors.append(f"Python syntax error in {rel}: {exc.msg}")


def main() -> None:
    errors: list[str] = []

    check_required_items(errors)
    check_json_metadata(errors)
    check_placeholders(errors)
    check_cjk_characters(errors)
    check_local_paths(errors)
    check_forbidden_files(errors)
    check_python_syntax(errors)

    if errors:
        fail(errors)

    print("Repository validation passed.")
    print(f"Checked repository root: {ROOT}")
    print(f"Checked files: {len(iter_repository_files())}")


if __name__ == "__main__":
    main()
