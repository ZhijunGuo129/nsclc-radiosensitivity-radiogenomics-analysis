# -*- coding: utf-8 -*-

# 03_validate_nrrd_masks.py
#
# Audit converted NRRD images, masks, and DICOM SEG metadata.
#
# Paths are configured through environment variables.

import os
import traceback
import numpy as np
import pandas as pd
import SimpleITK as sitk
import pydicom


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

conv_log_csv = os.path.join(
    lung1_root,
    "metadata",
    "lung1_dicomseg_conversion_log.csv"
)

nrrd_root = os.path.join(
    lung1_root,
    "NRRD_full"
)

metadata_dir = os.path.join(
    lung1_root,
    "metadata"
)

out_csv = os.path.join(
    metadata_dir,
    "lung1_nrrd_mask_qc_and_seg_metadata.csv"
)

summary_csv = os.path.join(
    metadata_dir,
    "lung1_nrrd_mask_qc_summary.csv"
)

xlsx_out = os.path.join(
    metadata_dir,
    "lung1_nrrd_mask_qc_audit.xlsx"
)

os.makedirs(metadata_dir, exist_ok=True)

############################################################
# 2. Helper functions
############################################################

def parse_seg_metadata(seg_file):
    out = {
        "dicom_seg_n_segments": None,
        "dicom_seg_numbers": "",
        "dicom_seg_labels": "",
        "dicom_seg_descriptions": "",
        "dicom_seg_property_types": "",
        "dicom_seg_property_categories": ""
    }

    try:
        ds = pydicom.dcmread(seg_file, stop_before_pixels=True, force=True)

        if not hasattr(ds, "SegmentSequence"):
            return out

        seg_seq = ds.SegmentSequence

        numbers = []
        labels = []
        descriptions = []
        property_types = []
        property_categories = []

        for seg in seg_seq:
            numbers.append(str(getattr(seg, "SegmentNumber", "")))
            labels.append(str(getattr(seg, "SegmentLabel", "")))
            descriptions.append(str(getattr(seg, "SegmentDescription", "")))

            try:
                type_seq = seg.SegmentedPropertyTypeCodeSequence
                property_types.append(str(type_seq[0].CodeMeaning))
            except Exception:
                property_types.append("")

            try:
                cat_seq = seg.SegmentedPropertyCategoryCodeSequence
                property_categories.append(str(cat_seq[0].CodeMeaning))
            except Exception:
                property_categories.append("")

        out["dicom_seg_n_segments"] = len(seg_seq)
        out["dicom_seg_numbers"] = ";".join(numbers)
        out["dicom_seg_labels"] = ";".join(labels)
        out["dicom_seg_descriptions"] = ";".join(descriptions)
        out["dicom_seg_property_types"] = ";".join(property_types)
        out["dicom_seg_property_categories"] = ";".join(property_categories)

    except Exception as e:
        out["dicom_seg_labels"] = "READ_FAILED: " + str(e)

    return out

def get_image_info(path):
    img = sitk.ReadImage(path)

    info = {
        "size_x": img.GetSize()[0],
        "size_y": img.GetSize()[1],
        "size_z": img.GetSize()[2],
        "spacing_x": img.GetSpacing()[0],
        "spacing_y": img.GetSpacing()[1],
        "spacing_z": img.GetSpacing()[2],
        "origin": str(img.GetOrigin()),
        "direction": str(img.GetDirection())
    }

    return img, info

def mask_qc(image_path, mask_path):
    image, image_info = get_image_info(image_path)
    mask, mask_info = get_image_info(mask_path)

    arr = sitk.GetArrayFromImage(mask)

    labels = np.unique(arr)
    labels_nonzero = labels[labels != 0]

    nonzero = int(np.count_nonzero(arr))

    spacing = mask.GetSpacing()
    voxel_volume_mm3 = float(spacing[0] * spacing[1] * spacing[2])
    volume_cm3 = nonzero * voxel_volume_mm3 / 1000.0

    if nonzero > 0:
        coords = np.argwhere(arr != 0)

        zmin, ymin, xmin = coords.min(axis=0)
        zmax, ymax, xmax = coords.max(axis=0)

        bbox_vox_x = int(xmax - xmin + 1)
        bbox_vox_y = int(ymax - ymin + 1)
        bbox_vox_z = int(zmax - zmin + 1)

        bbox_mm_x = bbox_vox_x * spacing[0]
        bbox_mm_y = bbox_vox_y * spacing[1]
        bbox_mm_z = bbox_vox_z * spacing[2]
    else:
        bbox_vox_x = bbox_vox_y = bbox_vox_z = 0
        bbox_mm_x = bbox_mm_y = bbox_mm_z = 0

    same_size = image.GetSize() == mask.GetSize()
    same_spacing = all([
        abs(image.GetSpacing()[i] - mask.GetSpacing()[i]) < 1e-6
        for i in range(3)
    ])

    qc = {
        "image_size_x": image_info["size_x"],
        "image_size_y": image_info["size_y"],
        "image_size_z": image_info["size_z"],
        "mask_size_x": mask_info["size_x"],
        "mask_size_y": mask_info["size_y"],
        "mask_size_z": mask_info["size_z"],
        "image_spacing_x": image_info["spacing_x"],
        "image_spacing_y": image_info["spacing_y"],
        "image_spacing_z": image_info["spacing_z"],
        "mask_spacing_x": mask_info["spacing_x"],
        "mask_spacing_y": mask_info["spacing_y"],
        "mask_spacing_z": mask_info["spacing_z"],
        "same_size": same_size,
        "same_spacing": same_spacing,
        "mask_min": int(np.min(arr)),
        "mask_max": int(np.max(arr)),
        "mask_nonzero_voxels": nonzero,
        "mask_nonzero_labels": ";".join([str(int(x)) for x in labels_nonzero]),
        "n_nonzero_labels": int(len(labels_nonzero)),
        "voxel_volume_mm3": voxel_volume_mm3,
        "mask_volume_cm3": volume_cm3,
        "bbox_vox_x": bbox_vox_x,
        "bbox_vox_y": bbox_vox_y,
        "bbox_vox_z": bbox_vox_z,
        "bbox_mm_x": bbox_mm_x,
        "bbox_mm_y": bbox_mm_y,
        "bbox_mm_z": bbox_mm_z
    }

    return qc

############################################################
# 3. Load tables
############################################################

pairs = pd.read_csv(pair_csv)
pairs["patient_id"] = pairs["patient_id"].astype(str).str.upper().str.strip()

pairs = pairs[pairs["pair_status"] == "CANDIDATE_PAIR_CREATED"].copy()
pairs = pairs.sort_values("patient_id").reset_index(drop=True)

print("\n===== Loaded candidate pairs =====")
print("Candidate pairs:", len(pairs))

if os.path.exists(conv_log_csv):
    conv_log = pd.read_csv(conv_log_csv)
    print("\n===== DICOM SEG conversion log status =====")
    print(conv_log["status"].value_counts(dropna=False))
else:
    conv_log = None
    print("\nDICOM SEG conversion log not found.")

############################################################
# 4. Main QC loop
############################################################

rows = []

for i, row in pairs.iterrows():
    pid = row["patient_id"]

    if (i + 1) % 25 == 0 or i == 0:
        print(f"QC {i+1}/{len(pairs)}: {pid}")

    image_path = os.path.join(nrrd_root, pid, f"{pid}_image.nrrd")
    mask_path = os.path.join(nrrd_root, pid, f"{pid}_mask.nrrd")
    seg_file = row["seg_file_path"]

    one = {
        "patient_id": pid,
        "status": "FAILED",
        "message": "",
        "image_path": image_path,
        "mask_path": mask_path,
        "seg_file_path": seg_file,
        "image_exists": os.path.exists(image_path),
        "mask_exists": os.path.exists(mask_path),
        "image_file_size": os.path.getsize(image_path) if os.path.exists(image_path) else None,
        "mask_file_size": os.path.getsize(mask_path) if os.path.exists(mask_path) else None
    }

    try:
        if not os.path.exists(image_path):
            raise FileNotFoundError(f"Image NRRD not found: {image_path}")

        if not os.path.exists(mask_path):
            raise FileNotFoundError(f"Mask NRRD not found: {mask_path}")

        qc = mask_qc(image_path, mask_path)
        seg_meta = parse_seg_metadata(seg_file)

        one.update(qc)
        one.update(seg_meta)

        if one["mask_nonzero_voxels"] <= 0:
            one["status"] = "FAILED"
            one["message"] = "Empty mask"
        elif not one["same_size"]:
            one["status"] = "WARNING"
            one["message"] = "Image and mask size mismatch"
        elif not one["same_spacing"]:
            one["status"] = "WARNING"
            one["message"] = "Image and mask spacing mismatch"
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
# 5. Summary
############################################################

summary_items = []

def add_summary(item, value):
    summary_items.append({"item": item, "value": value})

add_summary("candidate_pairs", len(pairs))
add_summary("qc_PASS", int((qc_df["status"] == "PASS").sum()))
add_summary("qc_WARNING", int((qc_df["status"] == "WARNING").sum()))
add_summary("qc_FAILED", int((qc_df["status"] == "FAILED").sum()))
add_summary("image_exists_count", int(qc_df["image_exists"].sum()))
add_summary("mask_exists_count", int(qc_df["mask_exists"].sum()))
add_summary("empty_mask_count", int((qc_df.get("mask_nonzero_voxels", pd.Series(dtype=float)) <= 0).sum()))
add_summary("same_size_count", int(qc_df.get("same_size", pd.Series(dtype=bool)).sum()))
add_summary("same_spacing_count", int(qc_df.get("same_spacing", pd.Series(dtype=bool)).sum()))

if "mask_volume_cm3" in qc_df.columns:
    add_summary("mask_volume_cm3_min", float(qc_df["mask_volume_cm3"].min()))
    add_summary("mask_volume_cm3_q25", float(qc_df["mask_volume_cm3"].quantile(0.25)))
    add_summary("mask_volume_cm3_median", float(qc_df["mask_volume_cm3"].median()))
    add_summary("mask_volume_cm3_mean", float(qc_df["mask_volume_cm3"].mean()))
    add_summary("mask_volume_cm3_q75", float(qc_df["mask_volume_cm3"].quantile(0.75)))
    add_summary("mask_volume_cm3_max", float(qc_df["mask_volume_cm3"].max()))

if "n_nonzero_labels" in qc_df.columns:
    add_summary("n_cases_with_one_nonzero_label", int((qc_df["n_nonzero_labels"] == 1).sum()))
    add_summary("n_cases_with_multiple_nonzero_labels", int((qc_df["n_nonzero_labels"] > 1).sum()))

if "dicom_seg_n_segments" in qc_df.columns:
    add_summary("dicom_SEG_n_segments_min", float(qc_df["dicom_seg_n_segments"].min()))
    add_summary("dicom_SEG_n_segments_median", float(qc_df["dicom_seg_n_segments"].median()))
    add_summary("dicom_SEG_n_segments_max", float(qc_df["dicom_seg_n_segments"].max()))

summary_df = pd.DataFrame(summary_items)

print("\n===== QC status table =====")
print(qc_df["status"].value_counts(dropna=False))

print("\n===== Mask label count table =====")
if "n_nonzero_labels" in qc_df.columns:
    print(qc_df["n_nonzero_labels"].value_counts(dropna=False).sort_index())

print("\n===== Mask volume cm3 summary =====")
if "mask_volume_cm3" in qc_df.columns:
    print(qc_df["mask_volume_cm3"].describe())

print("\n===== DICOM SEG segment labels table =====")
if "dicom_seg_labels" in qc_df.columns:
    print(qc_df["dicom_seg_labels"].value_counts(dropna=False).head(20))

print("\n===== DICOM SEG n_segments table =====")
if "dicom_seg_n_segments" in qc_df.columns:
    print(qc_df["dicom_seg_n_segments"].value_counts(dropna=False).sort_index())

print("\n===== Potential suspicious masks =====")
suspicious = qc_df[
    (qc_df["status"] != "PASS") |
    (qc_df.get("mask_volume_cm3", 0) <= 0)
].copy()

print(suspicious[["patient_id", "status", "message"]].head(50))

############################################################
# 6. Save outputs
############################################################

qc_df.to_csv(out_csv, index=False, encoding="utf-8-sig")
summary_df.to_csv(summary_csv, index=False, encoding="utf-8-sig")

with pd.ExcelWriter(xlsx_out, engine="openpyxl") as writer:
    summary_df.to_excel(writer, sheet_name="summary", index=False)
    qc_df.to_excel(writer, sheet_name="case_QC", index=False)

    if "dicom_seg_labels" in qc_df.columns:
        label_tab = qc_df["dicom_seg_labels"].value_counts(dropna=False).reset_index()
        label_tab.columns = ["dicom_seg_labels", "n_cases"]
        label_tab.to_excel(writer, sheet_name="segment_label_table", index=False)

    if "n_nonzero_labels" in qc_df.columns:
        label_count_tab = qc_df["n_nonzero_labels"].value_counts(dropna=False).reset_index()
        label_count_tab.columns = ["n_nonzero_labels", "n_cases"]
        label_count_tab.to_excel(writer, sheet_name="mask_label_count", index=False)

print("\n===== Lung1 NRRD mask quality control completed =====")
print("Main output files:")
print(out_csv)
print(summary_csv)
print(xlsx_out)

print("\nKey output files:")
print("1. QC status table")
print("2. Mask label count table")
print("3. Mask volume cm3 summary")
print("4. DICOM SEG segment labels table")
print("5. DICOM SEG n_segments table")
