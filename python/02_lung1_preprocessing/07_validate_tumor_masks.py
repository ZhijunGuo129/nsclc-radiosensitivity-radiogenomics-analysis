# -*- coding: utf-8 -*-

# 07_validate_tumor_masks.py
#
# Audit the final Lung1 primary-tumor masks.
#
# Paths are configured through environment variables.

import os
import numpy as np
import pandas as pd
import SimpleITK as sitk
import traceback


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
pair_csv = os.path.join(
    lung1_root,
    "metadata",
    "lung1_ct_seg_candidate_pairs.csv"
)

final_mask_root = os.path.join(
    lung1_root,
    "NRRD_tumor_only_final"
)

image_root = os.path.join(
    lung1_root,
    "NRRD_full"
)

metadata_dir = os.path.join(
    lung1_root,
    "metadata"
)

out_csv = os.path.join(
    metadata_dir,
    "lung1_final_tumor_mask_qc.csv"
)

summary_csv = os.path.join(
    metadata_dir,
    "lung1_final_tumor_mask_qc_summary.csv"
)

xlsx_out = os.path.join(
    metadata_dir,
    "lung1_final_tumor_mask_qc.xlsx"
)

############################################################
# 2. Load candidate patients
############################################################

pairs = pd.read_csv(pair_csv)
pairs["patient_id"] = pairs["patient_id"].astype(str).str.upper().str.strip()
pairs = pairs[pairs["pair_status"] == "CANDIDATE_PAIR_CREATED"].copy()
pairs = pairs.sort_values("patient_id").reset_index(drop=True)

print("\n===== Lung1 final tumor-mask quality control =====")
print("Candidate pairs:", len(pairs))
print("Image root:", image_root)
print("Final tumor mask root:", final_mask_root)

############################################################
# 3. QC loop
############################################################

rows = []

for i, row in pairs.iterrows():

    pid = row["patient_id"]

    if (i + 1) % 25 == 0 or i == 0:
        print(f"QC {i+1}/{len(pairs)}: {pid}")

    image_path = os.path.join(
        image_root,
        pid,
        f"{pid}_image.nrrd"
    )

    mask_path = os.path.join(
        final_mask_root,
        pid,
        f"{pid}_tumor_mask.nrrd"
    )

    one = {
        "patient_id": pid,
        "status": "FAILED",
        "message": "",
        "image_path": image_path,
        "mask_path": mask_path,
        "image_exists": os.path.exists(image_path),
        "mask_exists": os.path.exists(mask_path),
        "image_file_size": os.path.getsize(image_path) if os.path.exists(image_path) else None,
        "mask_file_size": os.path.getsize(mask_path) if os.path.exists(mask_path) else None,
        "image_size": "",
        "mask_size": "",
        "image_spacing": "",
        "mask_spacing": "",
        "same_size": "",
        "same_spacing": "",
        "mask_min": "",
        "mask_max": "",
        "mask_unique_labels": "",
        "n_nonzero_labels": "",
        "mask_nonzero_voxels": "",
        "tumor_volume_cm3": "",
        "bbox_vox_x": "",
        "bbox_vox_y": "",
        "bbox_vox_z": "",
        "bbox_mm_x": "",
        "bbox_mm_y": "",
        "bbox_mm_z": ""
    }

    try:
        if not os.path.exists(image_path):
            raise FileNotFoundError("Image not found: " + image_path)

        if not os.path.exists(mask_path):
            raise FileNotFoundError("Final tumor mask not found: " + mask_path)

        image = sitk.ReadImage(image_path)
        mask = sitk.ReadImage(mask_path)

        image_size = image.GetSize()
        mask_size = mask.GetSize()

        image_spacing = image.GetSpacing()
        mask_spacing = mask.GetSpacing()

        same_size = image_size == mask_size

        same_spacing = all([
            abs(image_spacing[j] - mask_spacing[j]) < 1e-6
            for j in range(3)
        ])

        arr = sitk.GetArrayFromImage(mask)

        labels = np.unique(arr)
        labels_int = [int(x) for x in labels.tolist()]
        labels_nonzero = [x for x in labels_int if x != 0]

        nonzero = int(np.count_nonzero(arr))

        voxel_volume_mm3 = float(mask_spacing[0] * mask_spacing[1] * mask_spacing[2])
        tumor_volume_cm3 = nonzero * voxel_volume_mm3 / 1000.0

        if nonzero > 0:
            coords = np.argwhere(arr != 0)

            zmin, ymin, xmin = coords.min(axis=0)
            zmax, ymax, xmax = coords.max(axis=0)

            bbox_vox_x = int(xmax - xmin + 1)
            bbox_vox_y = int(ymax - ymin + 1)
            bbox_vox_z = int(zmax - zmin + 1)

            bbox_mm_x = bbox_vox_x * mask_spacing[0]
            bbox_mm_y = bbox_vox_y * mask_spacing[1]
            bbox_mm_z = bbox_vox_z * mask_spacing[2]
        else:
            bbox_vox_x = bbox_vox_y = bbox_vox_z = 0
            bbox_mm_x = bbox_mm_y = bbox_mm_z = 0

        one.update({
            "image_size": str(image_size),
            "mask_size": str(mask_size),
            "image_spacing": str(image_spacing),
            "mask_spacing": str(mask_spacing),
            "same_size": same_size,
            "same_spacing": same_spacing,
            "mask_min": int(np.min(arr)),
            "mask_max": int(np.max(arr)),
            "mask_unique_labels": ";".join([str(x) for x in labels_int]),
            "n_nonzero_labels": len(labels_nonzero),
            "mask_nonzero_voxels": nonzero,
            "tumor_volume_cm3": tumor_volume_cm3,
            "bbox_vox_x": bbox_vox_x,
            "bbox_vox_y": bbox_vox_y,
            "bbox_vox_z": bbox_vox_z,
            "bbox_mm_x": bbox_mm_x,
            "bbox_mm_y": bbox_mm_y,
            "bbox_mm_z": bbox_mm_z
        })

        if nonzero <= 0:
            one["status"] = "FAILED"
            one["message"] = "Empty tumor mask"
        elif not same_size:
            one["status"] = "WARNING"
            one["message"] = "Image and mask size mismatch"
        elif not same_spacing:
            one["status"] = "WARNING"
            one["message"] = "Image and mask spacing mismatch"
        elif len(labels_nonzero) != 1:
            one["status"] = "WARNING"
            one["message"] = "Mask is not single-label binary"
        elif labels_nonzero[0] != 1:
            one["status"] = "WARNING"
            one["message"] = "Nonzero mask label is not 1"
        elif tumor_volume_cm3 <= 0:
            one["status"] = "FAILED"
            one["message"] = "Tumor volume <= 0"
        else:
            one["status"] = "PASS"
            one["message"] = "OK"

    except Exception as e:
        one["status"] = "FAILED"
        one["message"] = str(e)
        traceback.print_exc()

    rows.append(one)

    if (i + 1) % 25 == 0:
        pd.DataFrame(rows).to_csv(out_csv, index=False, encoding="utf-8-sig")

qc_df = pd.DataFrame(rows)

############################################################
# 4. Summary
############################################################

vol = pd.to_numeric(qc_df["tumor_volume_cm3"], errors="coerce")
nonzero = pd.to_numeric(qc_df["mask_nonzero_voxels"], errors="coerce")

summary_rows = []

def add_summary(item, value):
    summary_rows.append({"item": item, "value": value})

add_summary("candidate_pairs", len(pairs))
add_summary("final_mask_files_existing", int(qc_df["mask_exists"].sum()))
add_summary("image_files_existing", int(qc_df["image_exists"].sum()))
add_summary("PASS", int((qc_df["status"] == "PASS").sum()))
add_summary("WARNING", int((qc_df["status"] == "WARNING").sum()))
add_summary("FAILED", int((qc_df["status"] == "FAILED").sum()))
add_summary("same_size_count", int((qc_df["same_size"] == True).sum()))
add_summary("same_spacing_count", int((qc_df["same_spacing"] == True).sum()))
add_summary("empty_mask_count", int((nonzero <= 0).sum()))

add_summary("tumor_volume_cm3_min", float(vol.min()))
add_summary("tumor_volume_cm3_q25", float(vol.quantile(0.25)))
add_summary("tumor_volume_cm3_median", float(vol.median()))
add_summary("tumor_volume_cm3_mean", float(vol.mean()))
add_summary("tumor_volume_cm3_q75", float(vol.quantile(0.75)))
add_summary("tumor_volume_cm3_max", float(vol.max()))

summary_df = pd.DataFrame(summary_rows)

print("\n===== Final tumor mask QC status table =====")
print(qc_df["status"].value_counts(dropna=False))

print("\n===== Mask unique label table =====")
print(qc_df["mask_unique_labels"].value_counts(dropna=False).head(20))

print("\n===== Nonzero label count table =====")
print(qc_df["n_nonzero_labels"].value_counts(dropna=False).sort_index())

print("\n===== Tumor volume cm3 summary =====")
print(vol.describe())

print("\n===== Suspicious cases =====")
suspicious = qc_df[qc_df["status"] != "PASS"].copy()
print(suspicious[[
    "patient_id",
    "status",
    "message",
    "mask_unique_labels",
    "mask_nonzero_voxels",
    "tumor_volume_cm3"
]].head(100))

############################################################
# 5. Save outputs
############################################################

qc_df.to_csv(out_csv, index=False, encoding="utf-8-sig")
summary_df.to_csv(summary_csv, index=False, encoding="utf-8-sig")

with pd.ExcelWriter(xlsx_out, engine="openpyxl") as writer:
    summary_df.to_excel(writer, sheet_name="summary", index=False)
    qc_df.to_excel(writer, sheet_name="case_QC", index=False)

    label_tab = qc_df["mask_unique_labels"].value_counts(dropna=False).reset_index()
    label_tab.columns = ["mask_unique_labels", "n_cases"]
    label_tab.to_excel(writer, sheet_name="label_table", index=False)

    suspicious.to_excel(writer, sheet_name="suspicious_cases", index=False)

print("\n===== Lung1 final tumor-mask quality control completed =====")
print("Main output files:")
print(out_csv)
print(summary_csv)
print(xlsx_out)

print("\nKey output files:")
print("1. Final tumor mask QC status table")
print("2. Mask unique label table")
print("3. Nonzero label count table")
print("4. Tumor volume cm3 summary")
print("5. Suspicious cases")
