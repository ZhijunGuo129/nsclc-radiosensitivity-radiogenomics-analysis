# -*- coding: utf-8 -*-

# 01_inventory_dicom_series.py
#
# Inventory Lung1 DICOM CT and segmentation series.
#
# Paths are configured through environment variables.

import os
import re
import traceback
from collections import defaultdict

import pandas as pd

try:
    import pydicom
except ImportError:
    raise ImportError(
        "pydicom is not installed. Install pydicom with: python -m pip install pydicom"
    )


def required_env_directory(name):
    value = os.environ.get(name, "").strip()
    if not value:
        raise RuntimeError(f"Environment variable {name} is not set.")
    path = os.path.abspath(os.path.expanduser(value))
    if not os.path.isdir(path):
        raise FileNotFoundError(f"Directory specified by {name} does not exist: {path}")
    return path

############################################################
# 1. Paths
############################################################

lung1_root = required_env_directory("LUNG1_DATA_DIR")
dicom_dir = os.path.join(lung1_root, "TCIA_DICOM")
clinical_csv = os.path.join(lung1_root, "clinical", "Lung1_clinical_clean_initial.csv")
metadata_dir = os.path.join(lung1_root, "metadata")

os.makedirs(metadata_dir, exist_ok=True)

series_inventory_csv = os.path.join(
    metadata_dir,
    "lung1_ct_seg_series_inventory.csv"
)

patient_availability_csv = os.path.join(
    metadata_dir,
    "lung1_patient_ct_seg_availability.csv"
)

candidate_pairs_csv = os.path.join(
    metadata_dir,
    "lung1_ct_seg_candidate_pairs.csv"
)

xlsx_out = os.path.join(
    metadata_dir,
    "lung1_ct_seg_inventory_audit.xlsx"
)

log_csv = os.path.join(
    metadata_dir,
    "lung1_dicom_series_read_log.csv"
)

if not os.path.isdir(dicom_dir):
    raise FileNotFoundError(f"DICOM directory not found: {dicom_dir}")

if not os.path.exists(clinical_csv):
    raise FileNotFoundError(f"Clinical CSV not found: {clinical_csv}")

############################################################
# 2. Helper functions
############################################################

def norm_path(p):
    return os.path.abspath(p).replace("\\", "/")

def safe_get(ds, tag_name, default=""):
    try:
        v = getattr(ds, tag_name, default)
        if v is None:
            return default
        return str(v)
    except Exception:
        return default

def safe_float_get(ds, tag_name, default=None):
    try:
        v = getattr(ds, tag_name, default)
        if v is None:
            return default
        if isinstance(v, (list, tuple)):
            return float(v[0])
        return float(v)
    except Exception:
        return default

def safe_spacing(ds):
    try:
        v = getattr(ds, "PixelSpacing", "")
        if v is None:
            return ""
        return "\\".join([str(x) for x in v])
    except Exception:
        return ""

def extract_patient_id_from_path(path):
    m = re.search(r"LUNG1-\d{3}", path.upper())
    if m:
        return m.group(0)
    return ""

def read_dicom_header(first_file):
    ds = pydicom.dcmread(
        first_file,
        stop_before_pixels=True,
        force=True
    )

    row = {
        "PatientID_tag": safe_get(ds, "PatientID"),
        "StudyInstanceUID": safe_get(ds, "StudyInstanceUID"),
        "SeriesInstanceUID": safe_get(ds, "SeriesInstanceUID"),
        "Modality": safe_get(ds, "Modality"),
        "SeriesDescription": safe_get(ds, "SeriesDescription"),
        "StudyDescription": safe_get(ds, "StudyDescription"),
        "StudyDate": safe_get(ds, "StudyDate"),
        "SeriesDate": safe_get(ds, "SeriesDate"),
        "SeriesNumber": safe_get(ds, "SeriesNumber"),
        "Manufacturer": safe_get(ds, "Manufacturer"),
        "SliceThickness": safe_get(ds, "SliceThickness"),
        "PixelSpacing": safe_spacing(ds),
        "Rows": safe_get(ds, "Rows"),
        "Columns": safe_get(ds, "Columns")
    }

    return row

############################################################
# 3. Load clinical IDs
############################################################

clinical_df = pd.read_csv(clinical_csv)

if "patient_id" in clinical_df.columns:
    clinical_df["patient_id"] = clinical_df["patient_id"].astype(str).str.upper().str.strip()
elif "PatientID" in clinical_df.columns:
    clinical_df["patient_id"] = clinical_df["PatientID"].astype(str).str.upper().str.strip()
else:
    raise ValueError("Cannot find patient_id or PatientID in clinical CSV.")

clinical_ids = sorted(clinical_df["patient_id"].dropna().unique().tolist())

print("\n===== Lung1 clinical IDs =====")
print("Clinical patients:", len(clinical_ids))
print("First 10:", clinical_ids[:10])

############################################################
# 4. Scan DICOM files and group by series directory
############################################################

print("\n===== Scanning DICOM files =====")
print("DICOM dir:", dicom_dir)

series_files = defaultdict(list)

n_all_files = 0
n_dcm_files = 0

for root, dirs, files in os.walk(dicom_dir):
    for fn in files:
        n_all_files += 1

        if fn.lower().endswith(".dcm"):
            full_path = os.path.join(root, fn)
            pid = extract_patient_id_from_path(full_path)

            if pid:
                n_dcm_files += 1
                series_dir = norm_path(root)
                series_files[series_dir].append(norm_path(full_path))

print("All files scanned:", n_all_files)
print("DICOM files with LUNG1 ID:", n_dcm_files)
print("Detected series directories:", len(series_files))

############################################################
# 5. Build series-level inventory
############################################################

inventory_rows = []
log_rows = []

for i, (series_dir, files) in enumerate(sorted(series_files.items()), start=1):
    if i % 50 == 0 or i == 1:
        print(f"Reading series {i}/{len(series_files)}")

    files_sorted = sorted(files)
    first_file = files_sorted[0]

    pid_from_path = extract_patient_id_from_path(series_dir)

    log = {
        "series_dir": series_dir,
        "first_file": first_file,
        "status": "FAILED",
        "message": ""
    }

    try:
        header = read_dicom_header(first_file)

        pid_tag = header.get("PatientID_tag", "")
        pid_final = pid_tag.upper().strip() if pid_tag else pid_from_path

        row = {
            "patient_id": pid_final,
            "patient_id_from_path": pid_from_path,
            "series_dir": series_dir,
            "first_file": first_file,
            "n_files": len(files_sorted),
            "file_size_first": os.path.getsize(first_file)
        }

        row.update(header)

        inventory_rows.append(row)

        log["status"] = "SUCCESS"
        log["message"] = "OK"

    except Exception as e:
        log["status"] = "FAILED"
        log["message"] = str(e)
        traceback.print_exc()

        row = {
            "patient_id": pid_from_path,
            "patient_id_from_path": pid_from_path,
            "series_dir": series_dir,
            "first_file": first_file,
            "n_files": len(files_sorted),
            "file_size_first": os.path.getsize(first_file),
            "PatientID_tag": "",
            "StudyInstanceUID": "",
            "SeriesInstanceUID": "",
            "Modality": "",
            "SeriesDescription": "",
            "StudyDescription": "",
            "StudyDate": "",
            "SeriesDate": "",
            "SeriesNumber": "",
            "Manufacturer": "",
            "SliceThickness": "",
            "PixelSpacing": "",
            "Rows": "",
            "Columns": ""
        }

        inventory_rows.append(row)

    log_rows.append(log)

inventory_df = pd.DataFrame(inventory_rows)
log_df = pd.DataFrame(log_rows)

inventory_df["patient_id"] = inventory_df["patient_id"].astype(str).str.upper().str.strip()
inventory_df["patient_id_from_path"] = inventory_df["patient_id_from_path"].astype(str).str.upper().str.strip()
inventory_df["Modality"] = inventory_df["Modality"].astype(str).str.upper().str.strip()
inventory_df["SeriesDescription"] = inventory_df["SeriesDescription"].astype(str)

print("\n===== Series inventory summary =====")
print("Inventory shape:", inventory_df.shape)

print("\nRead status table:")
print(log_df["status"].value_counts(dropna=False))

print("\nModality table:")
print(inventory_df["Modality"].value_counts(dropna=False))

print("\nSeries count by first 10 patients:")
print(inventory_df.groupby("patient_id").size().head(10))

############################################################
# 6. Patient-level CT/SEG availability
############################################################

availability_rows = []

for pid in clinical_ids:
    sub = inventory_df[inventory_df["patient_id"] == pid].copy()

    n_series = len(sub)
    n_ct = int((sub["Modality"] == "CT").sum())
    n_seg = int((sub["Modality"] == "SEG").sum())
    n_rtstruct = int((sub["Modality"] == "RTSTRUCT").sum())

    availability_rows.append({
        "patient_id": pid,
        "n_series": n_series,
        "n_CT_series": n_ct,
        "n_SEG_series": n_seg,
        "n_RTSTRUCT_series": n_rtstruct,
        "has_CT": n_ct > 0,
        "has_SEG": n_seg > 0,
        "has_RTSTRUCT": n_rtstruct > 0,
        "has_CT_SEG": (n_ct > 0 and n_seg > 0)
    })

availability_df = pd.DataFrame(availability_rows)

print("\n===== Patient CT/SEG availability =====")
print("Clinical patients:", len(availability_df))
print("Patients with CT:", int(availability_df["has_CT"].sum()))
print("Patients with SEG:", int(availability_df["has_SEG"].sum()))
print("Patients with CT+SEG:", int(availability_df["has_CT_SEG"].sum()))

print("\nAvailability status table:")
print(availability_df[["has_CT", "has_SEG", "has_CT_SEG"]].value_counts(dropna=False))

print("\nPatients without CT+SEG:")
print(availability_df.loc[~availability_df["has_CT_SEG"], "patient_id"].tolist())

############################################################
# 7. Create candidate CT-SEG pairs
############################################################

pair_rows = []

for pid in clinical_ids:
    sub = inventory_df[inventory_df["patient_id"] == pid].copy()

    ct_sub = sub[sub["Modality"] == "CT"].copy()
    seg_sub = sub[sub["Modality"] == "SEG"].copy()

    if len(ct_sub) == 0 and len(seg_sub) == 0:
        pair_rows.append({
            "patient_id": pid,
            "pair_status": "NO_CT_NO_SEG",
            "pairing_rule": "",
            "ct_series_dir": "",
            "ct_first_file": "",
            "ct_series_uid": "",
            "ct_study_uid": "",
            "ct_n_files": 0,
            "ct_series_description": "",
            "seg_file_path": "",
            "seg_series_dir": "",
            "seg_series_uid": "",
            "seg_study_uid": "",
            "seg_n_files": 0,
            "seg_series_description": "",
            "message": "No CT or SEG series found"
        })
        continue

    if len(ct_sub) == 0:
        pair_rows.append({
            "patient_id": pid,
            "pair_status": "NO_CT",
            "pairing_rule": "",
            "ct_series_dir": "",
            "ct_first_file": "",
            "ct_series_uid": "",
            "ct_study_uid": "",
            "ct_n_files": 0,
            "ct_series_description": "",
            "seg_file_path": seg_sub.iloc[0]["first_file"] if len(seg_sub) > 0 else "",
            "seg_series_dir": seg_sub.iloc[0]["series_dir"] if len(seg_sub) > 0 else "",
            "seg_series_uid": seg_sub.iloc[0]["SeriesInstanceUID"] if len(seg_sub) > 0 else "",
            "seg_study_uid": seg_sub.iloc[0]["StudyInstanceUID"] if len(seg_sub) > 0 else "",
            "seg_n_files": int(seg_sub.iloc[0]["n_files"]) if len(seg_sub) > 0 else 0,
            "seg_series_description": seg_sub.iloc[0]["SeriesDescription"] if len(seg_sub) > 0 else "",
            "message": "SEG found but no CT series found"
        })
        continue

    if len(seg_sub) == 0:
        ct_pick = ct_sub.sort_values("n_files", ascending=False).iloc[0]

        pair_rows.append({
            "patient_id": pid,
            "pair_status": "NO_SEG",
            "pairing_rule": "largest_CT_only_no_SEG",
            "ct_series_dir": ct_pick["series_dir"],
            "ct_first_file": ct_pick["first_file"],
            "ct_series_uid": ct_pick["SeriesInstanceUID"],
            "ct_study_uid": ct_pick["StudyInstanceUID"],
            "ct_n_files": int(ct_pick["n_files"]),
            "ct_series_description": ct_pick["SeriesDescription"],
            "seg_file_path": "",
            "seg_series_dir": "",
            "seg_series_uid": "",
            "seg_study_uid": "",
            "seg_n_files": 0,
            "seg_series_description": "",
            "message": "CT found but no SEG series found"
        })
        continue

    # Usually one SEG per patient. Pick SEG with largest file size if multiple.
    seg_sub = seg_sub.sort_values(["n_files", "file_size_first"], ascending=False)
    seg_pick = seg_sub.iloc[0]

    seg_study_uid = seg_pick["StudyInstanceUID"]

    ct_same_study = ct_sub[ct_sub["StudyInstanceUID"] == seg_study_uid].copy()

    if len(ct_same_study) > 0:
        ct_pick = ct_same_study.sort_values("n_files", ascending=False).iloc[0]
        rule = "largest_CT_in_same_study_as_selected_SEG"
    else:
        ct_pick = ct_sub.sort_values("n_files", ascending=False).iloc[0]
        rule = "largest_CT_any_study_SEG_study_not_matched"

    pair_rows.append({
        "patient_id": pid,
        "pair_status": "CANDIDATE_PAIR_CREATED",
        "pairing_rule": rule,
        "ct_series_dir": ct_pick["series_dir"],
        "ct_first_file": ct_pick["first_file"],
        "ct_series_uid": ct_pick["SeriesInstanceUID"],
        "ct_study_uid": ct_pick["StudyInstanceUID"],
        "ct_n_files": int(ct_pick["n_files"]),
        "ct_series_description": ct_pick["SeriesDescription"],
        "seg_file_path": seg_pick["first_file"],
        "seg_series_dir": seg_pick["series_dir"],
        "seg_series_uid": seg_pick["SeriesInstanceUID"],
        "seg_study_uid": seg_pick["StudyInstanceUID"],
        "seg_n_files": int(seg_pick["n_files"]),
        "seg_series_description": seg_pick["SeriesDescription"],
        "message": "OK"
    })

pairs_df = pd.DataFrame(pair_rows)

print("\n===== Candidate CT-SEG pairs =====")
print("Pair table shape:", pairs_df.shape)

print("\nPair status table:")
print(pairs_df["pair_status"].value_counts(dropna=False))

print("\nPairing rule table:")
print(pairs_df["pairing_rule"].value_counts(dropna=False))

print("\nPatients without candidate pair:")
print(pairs_df.loc[pairs_df["pair_status"] != "CANDIDATE_PAIR_CREATED", ["patient_id", "pair_status", "message"]])

############################################################
# 8. Save outputs
############################################################

inventory_df.to_csv(series_inventory_csv, index=False, encoding="utf-8-sig")
availability_df.to_csv(patient_availability_csv, index=False, encoding="utf-8-sig")
pairs_df.to_csv(candidate_pairs_csv, index=False, encoding="utf-8-sig")
log_df.to_csv(log_csv, index=False, encoding="utf-8-sig")

with pd.ExcelWriter(xlsx_out, engine="openpyxl") as writer:
    inventory_df.to_excel(writer, sheet_name="series_inventory", index=False)
    availability_df.to_excel(writer, sheet_name="patient_availability", index=False)
    pairs_df.to_excel(writer, sheet_name="candidate_pairs", index=False)
    log_df.to_excel(writer, sheet_name="read_log", index=False)

print("\n===== Lung1 CT-SEG series inventory completed =====")
print("Main output files:")
print(series_inventory_csv)
print(patient_availability_csv)
print(candidate_pairs_csv)
print(xlsx_out)

print("\nKey output files:")
print("1. Modality table")
print("2. Patient CT/SEG availability")
print("3. Pair status table")
print("4. Pairing rule table")
print("5. Patients without candidate pair")
