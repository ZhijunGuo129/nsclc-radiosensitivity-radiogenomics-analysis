from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SELF = Path(__file__).resolve()

CJK = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff]")
NON_ASCII_FILENAME = re.compile(r"[^\x20-\x7E]")
WINDOWS_ABSOLUTE_PATH = re.compile(r"(?<!https:)(?<!http:)\b[A-Za-z]:[\\/](?:Users|Documents|Desktop|Downloads|data|project|work|analysis)[\\/]", re.I)
UNIX_HOME_PATH = re.compile(r"(?<!https:)(?<!http:)(?:/Users/[^/\s]+|/home/[^/\s]+)(?:/[^\s\"']*)?")

TEXT_NAMES = {"LICENSE", ".gitignore", ".gitattributes"}
TEXT_SUFFIXES = {
    ".r", ".py", ".md", ".txt", ".cff", ".json", ".yml", ".yaml",
    ".example", ".gitignore", ".rproj", ".tsv"
}
EXCLUDED_SUFFIXES = {
    ".dcm", ".nrrd", ".nii", ".nii.gz", ".rds", ".rdata", ".rhistory",
    ".xlsx", ".xls", ".csv", ".tsv.gz", ".pickle", ".pkl", ".joblib",
    ".h5", ".hdf5", ".zip", ".7z", ".tar", ".gz"
}
EXCLUDED_DIR_NAMES = {
    "__pycache__", ".pytest_cache", ".mypy_cache", ".ruff_cache", ".idea",
    ".vscode", "renv", "packrat", "cache", "tmp", "temp", "logs"
}
SUSPICIOUS_FILE_STEMS = re.compile(
    r"(?:^|[_\-.])(test|testing|tmp|temp|backup|bak|copy|draft|old|obsolete|deprecated|scratch|debug)(?:$|[_\-.])",
    re.I,
)

REQUIRED_TOP_LEVEL = {
    "README.md",
    "LICENSE",
    "CHANGELOG.md",
    "CODE_AVAILABILITY.md",
    "CITATION.cff.template",
    ".zenodo.json.template",
    ".gitignore",
    ".gitattributes",
    "R",
    "python",
    "environment",
    "config",
    "data",
    "docs",
    "tools",
}

BASE_R_PACKAGES = {
    "base", "compiler", "datasets", "graphics", "grDevices", "grid", "methods",
    "parallel", "splines", "stats", "stats4", "tcltk", "tools", "utils"
}
SLICER_ONLY_PYTHON_MODULES = {"slicer", "vtk", "DICOMLib"}
STANDARD_LIBRARY_MODULES = {
    "__future__", "ast", "collections", "csv", "importlib", "json", "math", "os",
    "pathlib", "platform", "re", "shutil", "sys", "traceback", "typing"
}
PYTHON_REQUIREMENT_IMPORT_MAP = {
    "SimpleITK": "simpleitk",
    "numpy": "numpy",
    "pandas": "pandas",
    "pydicom": "pydicom",
    "radiomics": "pyradiomics",
    "openpyxl": "openpyxl",
}


def is_text_file(path: Path) -> bool:
    return path.name in TEXT_NAMES or path.suffix.lower() in TEXT_SUFFIXES


def logical_suffix(path: Path) -> str:
    name = path.name.lower()
    for suffix in (".nii.gz", ".tsv.gz"):
        if name.endswith(suffix):
            return suffix
    return path.suffix.lower()


def strip_r_comments_and_strings(text: str) -> str:
    out: list[str] = []
    quote: str | None = None
    escaped = False
    i = 0
    while i < len(text):
        char = text[i]
        if quote is not None:
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == quote:
                quote = None
            out.append(" ")
            i += 1
            continue
        if char in {'"', "'", "`"}:
            quote = char
            out.append(" ")
            i += 1
            continue
        if char == "#":
            newline = text.find("\n", i)
            if newline == -1:
                out.extend(" " * (len(text) - i))
                break
            out.extend(" " * (newline - i))
            i = newline
            continue
        out.append(char)
        i += 1
    return "".join(out)


def check_r_delimiters(text: str, filename: str) -> list[str]:
    errors: list[str] = []
    stack: list[tuple[str, int]] = []
    pairs = {")": "(", "]": "[", "}": "{"}
    line = 1
    i = 0
    quote: str | None = None
    escaped = False

    while i < len(text):
        char = text[i]
        if char == "\n":
            line += 1
        if quote is not None:
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == quote:
                quote = None
            i += 1
            continue
        if char in {'"', "'", "`"}:
            quote = char
            i += 1
            continue
        if char == "#":
            newline = text.find("\n", i)
            if newline == -1:
                break
            i = newline
            continue
        if char in "([{":
            stack.append((char, line))
        elif char in ")]}" :
            if not stack or stack[-1][0] != pairs[char]:
                errors.append(f"R delimiter mismatch in {filename} at line {line}: {char}")
                return errors
            stack.pop()
        i += 1

    if quote is not None:
        errors.append(f"Unterminated R string in {filename}")
    if stack:
        opener, opener_line = stack[-1]
        errors.append(f"Unclosed R delimiter in {filename}: {opener} from line {opener_line}")
    return errors


def parse_r_dependencies(text: str) -> set[str]:
    packages = set(re.findall(r"\b([A-Za-z][A-Za-z0-9.]*)::", text))
    packages.update(
        re.findall(
            r"(?:library|require|requireNamespace)\(\s*['\"]?([A-Za-z][A-Za-z0-9.]*)",
            text,
        )
    )
    return packages - BASE_R_PACKAGES - {"BiocManager"}


def parse_installed_r_packages(text: str) -> set[str]:
    packages: set[str] = set()
    for block_name in ("cran_packages", "bioconductor_packages"):
        match = re.search(
            rf"{block_name}\s*<-\s*c\((.*?)\)",
            text,
            flags=re.S,
        )
        if match:
            packages.update(re.findall(r"['\"]([^'\"]+)['\"]", match.group(1)))
    return packages


def parse_python_imports(path: Path, text: str) -> tuple[set[str], list[str]]:
    errors: list[str] = []
    modules: set[str] = set()
    try:
        tree = ast.parse(text, filename=path.as_posix())
    except SyntaxError as exc:
        return modules, [f"Python syntax error in {path}: {exc}"]
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            modules.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            modules.add(node.module.split(".")[0])
    return modules, errors


def parse_requirements(path: Path) -> set[str]:
    packages: set[str] = set()
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line or line.startswith("-"):
            continue
        name = re.split(r"[<>=!~\[]", line, maxsplit=1)[0].strip().lower()
        if name:
            packages.add(name)
    return packages


def pipeline_script_paths(text: str) -> set[str]:
    return set(re.findall(r"`((?:R|python)/[^`]+\.(?:R|py))`", text))


def main() -> int:
    errors: list[str] = []
    warnings: list[str] = []

    missing_top = sorted(name for name in REQUIRED_TOP_LEVEL if not (ROOT / name).exists())
    if missing_top:
        errors.append("Missing top-level repository items: " + ", ".join(missing_top))

    r_dependencies: set[str] = set()
    python_imports: set[str] = set()

    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or ".git" in path.parts or path.resolve() == SELF:
            continue

        rel = path.relative_to(ROOT).as_posix()
        lower_parts = {part.lower() for part in path.relative_to(ROOT).parts[:-1]}

        if NON_ASCII_FILENAME.search(rel) or CJK.search(rel):
            errors.append(f"Non-ASCII or CJK filename: {rel}")
        if lower_parts & EXCLUDED_DIR_NAMES:
            errors.append(f"Excluded cache, temporary, or local directory: {rel}")
        if SUSPICIOUS_FILE_STEMS.search(path.stem):
            errors.append(f"Suspicious test, draft, backup, or obsolete filename: {rel}")
        if logical_suffix(path) in EXCLUDED_SUFFIXES:
            errors.append(f"Excluded data, archive, or analysis-output file: {rel}")

        if not is_text_file(path):
            continue

        try:
            text = path.read_text(encoding="utf-8", errors="strict")
        except UnicodeDecodeError as exc:
            errors.append(f"Invalid UTF-8 text file {rel}: {exc}")
            continue

        if CJK.search(text):
            errors.append(f"CJK text detected: {rel}")
        # Templates and documentation intentionally contain abstract placeholder paths.
        if path.suffix.lower() in {".r", ".py"}:
            if WINDOWS_ABSOLUTE_PATH.search(text) or UNIX_HOME_PATH.search(text):
                errors.append(f"Machine-specific absolute path detected: {rel}")

        if path.suffix.lower() == ".py":
            modules, syntax_errors = parse_python_imports(path, text)
            python_imports.update(modules)
            errors.extend(syntax_errors)
        elif path.suffix.lower() == ".r":
            errors.extend(check_r_delimiters(text, rel))
            r_dependencies.update(parse_r_dependencies(text))
            if "install.packages(" in strip_r_comments_and_strings(text) and rel != "environment/install_r_dependencies.R":
                errors.append(f"Package installation call outside the dependency installer: {rel}")

    pipeline_file = ROOT / "docs" / "pipeline.md"
    if pipeline_file.exists():
        for script in sorted(pipeline_script_paths(pipeline_file.read_text(encoding="utf-8"))):
            if not (ROOT / script).exists():
                errors.append(f"Pipeline references a missing script: {script}")

    installer_file = ROOT / "environment" / "install_r_dependencies.R"
    if installer_file.exists():
        installed_r = parse_installed_r_packages(installer_file.read_text(encoding="utf-8"))
        undeclared_r = sorted(r_dependencies - installed_r)
        if undeclared_r:
            errors.append("R dependencies missing from installer: " + ", ".join(undeclared_r))

    requirements_file = ROOT / "environment" / "requirements.txt"
    if requirements_file.exists():
        requirements = parse_requirements(requirements_file)
        third_party_imports = python_imports - STANDARD_LIBRARY_MODULES - SLICER_ONLY_PYTHON_MODULES
        missing_requirements = sorted(
            module
            for module in third_party_imports
            if PYTHON_REQUIREMENT_IMPORT_MAP.get(module, module.lower()) not in requirements
        )
        if missing_requirements:
            errors.append(
                "Python imports missing from requirements.txt: " + ", ".join(missing_requirements)
            )

    rrs_permutation = ROOT / "R" / "06_sensitivity_analyses" / "04_permute_rrs_model.R"
    if rrs_permutation.exists():
        text = rrs_permutation.read_text(encoding="utf-8")
        if re.search(r"abs\s*\(\s*perm_results\$Spearman", text):
            errors.append("Obsolete absolute-correlation two-sided RRS permutation formula detected.")
        required_tokens = {
            "empirical_p_spearman_less",
            "2 * min(empirical_p_spearman_greater, empirical_p_spearman_less)",
        }
        for token in required_tokens:
            if token not in text:
                errors.append(f"Required doubled-tail RRS permutation calculation is missing: {token}")

    placeholder_files = [ROOT / "CITATION.cff.template", ROOT / ".zenodo.json.template"]
    for path in placeholder_files:
        if path.exists() and any(token in path.read_text(encoding="utf-8") for token in ("FAMILY_NAME", "GIVEN_NAMES", "AFFILIATION", "OWNER/REPOSITORY")):
            warnings.append(f"Release metadata still require author completion: {path.name}")

    if errors:
        print("Repository validation failed:")
        for error in errors:
            print(f"- {error}")
        if warnings:
            print("\nWarnings:")
            for warning in warnings:
                print(f"- {warning}")
        return 1

    print("Repository validation passed.")
    print(f"R dependencies detected: {len(r_dependencies)}")
    print(f"Python top-level imports detected: {len(python_imports)}")
    if warnings:
        print("Warnings:")
        for warning in warnings:
            print(f"- {warning}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
