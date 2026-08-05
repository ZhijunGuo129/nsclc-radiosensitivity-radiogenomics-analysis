from importlib import metadata
from pathlib import Path
import platform
import sys

packages = ["numpy", "openpyxl", "pandas", "pydicom", "pyradiomics", "SimpleITK"]
lines = [
    f"Python: {sys.version}",
    f"Platform: {platform.platform()}",
    "",
    "Installed packages:",
]
for package in packages:
    try:
        lines.append(f"{package}=={metadata.version(package)}")
    except metadata.PackageNotFoundError:
        lines.append(f"{package}: not installed")

out = Path("environment") / "captured" / "python_environment.txt"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text("\n".join(lines) + "\n", encoding="utf-8")
