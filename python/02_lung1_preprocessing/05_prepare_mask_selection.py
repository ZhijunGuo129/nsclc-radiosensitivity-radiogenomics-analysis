# -*- coding: utf-8 -*-

# 05_prepare_mask_selection.py
#
# Copy validated masks and create a targeted re-export list.
#
# Paths are configured through environment variables.

import os
import shutil
import pandas as pd


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
fast_qc_csv = os.path.join(
    lung1_root,
    "metadata",
    "lung1_primary_tumor_mask_qc.csv"
)

pair_csv = os.path.join(
    lung1_root,
    "metadata",
    "lung1_ct_seg_candidate_pairs.csv"
)

fast_root = os.path.join(
    lung1_root,
    "NRRD_tumor_only_fast"
)

final_root = os.path.join(
    lung1_root,
    "NRRD_tumor_only_final"
)

metadata_dir = os.path.join(
    lung1_root,
    "metadata"
)

target_csv = os.path.join(
    metadata_dir,
    "lung1_selected_tumor_mask_targets.csv"
)

final_manifest_csv = os.path.join(
    metadata_dir,
    "lung1_tumor_mask_manifest.csv"
)

summary_csv = os.path.join(
    metadata_dir,
    "lung1_tumor_mask_preparation_summary.csv"
)

os.makedirs(final_root, exist_ok=True)
os.makedirs(metadata_dir, exist_ok=True)

############################################################
# 2. Helper functions
############################################################

def get_first_existing_value(row, candidate_cols, default=""):
    """
    After merge, duplicated columns may become xxx_pair / xxx_fast.
    This function safely retrieves the first available non-empty value.
    """
    for cc in candidate_cols:
        if cc in row.index:
            val = row[cc]
            if pd.notna(val):
                val_str = str(val)
                if val_str.lower() != "nan":
                    return val_str
    return default

def safe_float(x, default=0.0):
    try:
        if pd.isna(x):
            return default
        return float(x)
    except Exception:
        return default

def path_exists_safe(p):
    try:
        return os.path.exists(str(p))
    except Exception:
        return False

############################################################
# 3. Load tables
############################################################

if not os.path.exists(pair_csv):
    raise FileNotFoundError("Pair CSV not found: " + pair_csv)

if not os.path.exists(fast_qc_csv):
    raise FileNotFoundError(
        "FAST QC CSV not found: " + fast_qc_csv +
        "\nThe expected intermediate file was not found. "
        "Check whether this CSV exists in the metadata folder."
    )

pairs = pd.read_csv(pair_csv)
fast_qc = pd.read_csv(fast_qc_csv)

pairs["patient_id"] = pairs["patient_id"].astype(str).str.upper().str.strip()
fast_qc["patient_id"] = fast_qc["patient_id"].astype(str).str.upper().str.strip()

pairs = pairs[pairs["pair_status"] == "CANDIDATE_PAIR_CREATED"].copy()
pairs = pairs.sort_values("patient_id").reset_index(drop=True)

fast_qc = fast_qc.sort_values("patient_id").reset_index(drop=True)

print("\n===== Loaded input tables =====")
print("Candidate CT-SEG pairs:", len(pairs))
print("FAST QC rows:", len(fast_qc))

print("\nPair columns:")
print(list(pairs.columns))

print("\nFAST QC columns:")
print(list(fast_qc.columns))

############################################################
# 4. Merge pair table and FAST QC
############################################################

dat = pairs.merge(
    fast_qc,
    on="patient_id",
    how="left",
    suffixes=("_pair", "_fast")
)

print("\nMerged table shape:")
print(dat.shape)

print("\nMerged columns:")
print(list(dat.columns))

############################################################
# 5. Build final mask manifest and Slicer target list
############################################################

manifest_rows = []
target_rows = []

for _, row in dat.iterrows():

    pid = str(row["patient_id"]).upper().strip()

    fast_status = get_first_existing_value(
        row,
        ["status", "status_fast"],
        default="MISSING_FAST_QC"
    )

    fast_message = get_first_existing_value(
        row,
        ["message_fast", "message", "message_pair"],
        default=""
    )

    fast_mask_path = get_first_existing_value(
        row,
        ["tumor_mask_path", "tumor_mask_path_fast"],
        default=""
    )

    ct_series_dir = get_first_existing_value(
        row,
        ["ct_series_dir", "ct_series_dir_pair"],
        default=""
    )

    seg_series_dir = get_first_existing_value(
        row,
        ["seg_series_dir", "seg_series_dir_pair"],
        default=""
    )

    seg_file_path = get_first_existing_value(
        row,
        ["seg_file_path_pair", "seg_file_path", "seg_file_path_fast"],
        default=""
    )

    tumor_nonzero = safe_float(
        row["tumor_nonzero_voxels"] if "tumor_nonzero_voxels" in row.index else 0,
        default=0
    )

    tumor_volume_cm3 = get_first_existing_value(
        row,
        ["tumor_volume_cm3"],
        default=""
    )

    final_dir = os.path.join(final_root, pid)
    os.makedirs(final_dir, exist_ok=True)

    final_mask_path = os.path.join(
        final_dir,
        f"{pid}_tumor_mask.nrrd"
    )

    use_fast = (
        fast_status == "PASS"
        and fast_mask_path != ""
        and os.path.exists(fast_mask_path)
        and tumor_nonzero > 0
    )

    if use_fast:

        shutil.copy2(fast_mask_path, final_mask_path)

        manifest_rows.append({
            "patient_id": pid,
            "final_status": "FAST_PASS_COPIED",
            "final_mask_path": final_mask_path,
            "source_mask_path": fast_mask_path,
            "needs_slicer": False,
            "fast_status": fast_status,
            "fast_message": fast_message,
            "tumor_nonzero_voxels": tumor_nonzero,
            "tumor_volume_cm3": tumor_volume_cm3
        })

    else:

        # Remove an unreliable final mask before targeted re-export in 3D Slicer.
        if os.path.exists(final_mask_path):
            os.remove(final_mask_path)

        target_rows.append({
            "patient_id": pid,
            "ct_series_dir": ct_series_dir,
            "seg_series_dir": seg_series_dir,
            "seg_file_path": seg_file_path,
            "final_mask_path": final_mask_path,
            "fast_status": fast_status,
            "fast_message": fast_message,
            "fast_tumor_nonzero_voxels": tumor_nonzero,
            "fast_tumor_volume_cm3": tumor_volume_cm3
        })

        manifest_rows.append({
            "patient_id": pid,
            "final_status": "SLICER_TARGET_PENDING",
            "final_mask_path": final_mask_path,
            "source_mask_path": "",
            "needs_slicer": True,
            "fast_status": fast_status,
            "fast_message": fast_message,
            "tumor_nonzero_voxels": tumor_nonzero,
            "tumor_volume_cm3": tumor_volume_cm3
        })

manifest_df = pd.DataFrame(manifest_rows)
target_df = pd.DataFrame(target_rows)

############################################################
# 6. Save outputs
############################################################

manifest_df.to_csv(
    final_manifest_csv,
    index=False,
    encoding="utf-8-sig"
)

target_df.to_csv(
    target_csv,
    index=False,
    encoding="utf-8-sig"
)

summary_df = pd.DataFrame({
    "item": [
        "candidate_pairs",
        "fast_qc_rows",
        "fast_PASS_copied_to_final",
        "needs_Slicer_targeted_export",
        "final_mask_files_existing_now",
        "target_csv",
        "final_mask_root"
    ],
    "value": [
        len(dat),
        len(fast_qc),
        int((manifest_df["final_status"] == "FAST_PASS_COPIED").sum()),
        int((manifest_df["needs_slicer"] == True).sum()),
        int(manifest_df["final_mask_path"].apply(path_exists_safe).sum()),
        target_csv,
        final_root
    ]
})

summary_df.to_csv(
    summary_csv,
    index=False,
    encoding="utf-8-sig"
)

############################################################
# 7. Print summary
############################################################

print("\n===== Prepare selected tumor masks =====")
print(summary_df)

print("\nFinal manifest status table:")
print(manifest_df["final_status"].value_counts(dropna=False))

print("\nTarget fast_status table:")
if len(target_df) > 0:
    print(target_df["fast_status"].value_counts(dropna=False))
else:
    print("No Slicer target cases.")

print("\nFirst 30 Slicer target patients:")
if len(target_df) > 0:
    print(target_df["patient_id"].head(30).tolist())
else:
    print([])

print("\nCheck target key columns:")
if len(target_df) > 0:
    print(target_df[[
        "patient_id",
        "ct_series_dir",
        "seg_series_dir",
        "final_mask_path",
        "fast_status"
    ]].head(10))

print("\nFiles saved:")
print(final_manifest_csv)
print(target_csv)
print(summary_csv)

print("\nDONE.")
