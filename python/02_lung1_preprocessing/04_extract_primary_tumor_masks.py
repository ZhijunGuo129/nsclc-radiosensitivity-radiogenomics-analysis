# -*- coding: utf-8 -*-

# 04_extract_primary_tumor_masks.py
#
# Extract primary-tumor masks from multi-label segmentations.
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

nrrd_full_root = os.path.join(
    lung1_root,
    "NRRD_full"
)

tumor_mask_root = os.path.join(
    lung1_root,
    "NRRD_tumor_only_fast"
)

metadata_dir = os.path.join(
    lung1_root,
    "metadata"
)

os.makedirs(tumor_mask_root, exist_ok=True)
os.makedirs(metadata_dir, exist_ok=True)

out_csv = os.path.join(
    metadata_dir,
    "lung1_primary_tumor_mask_qc.csv"
)

summary_csv = os.path.join(
    metadata_dir,
    "lung1_primary_tumor_mask_summary.csv"
)

xlsx_out = os.path.join(
    metadata_dir,
    "lung1_primary_tumor_mask_qc.xlsx"
)

############################################################
# 2. Helper functions
############################################################

def find_tumor_segment_from_dicom_seg(seg_file):
    """
    Return tumor segment number and label from DICOM SEG metadata.
    Priority:
      1. Exact SegmentLabel == Neoplasm, Primary
      2. Label contains both neoplasm and primary
      3. Label contains tumor or gtv
    """

    ds = pydicom.dcmread(seg_file, stop_before_pixels=True, force=True)

    if not hasattr(ds, "SegmentSequence"):
        raise RuntimeError("No SegmentSequence found in DICOM SEG.")

    segments = []

    for idx, seg in enumerate(ds.SegmentSequence):
        seg_number = getattr(seg, "SegmentNumber", idx + 1)
        seg_label = str(getattr(seg, "SegmentLabel", ""))

        segments.append({
            "index_1based": idx + 1,
            "segment_number": int(seg_number),
            "segment_label": seg_label
        })

    # exact match
    for x in segments:
        if x["segment_label"].lower().strip() == "neoplasm, primary":
            return x, segments

    # contains neoplasm and primary
    for x in segments:
        name = x["segment_label"].lower().strip()
        if ("neoplasm" in name) and ("primary" in name):
            return x, segments

    # fallback tumor / gtv
    for x in segments:
        name = x["segment_label"].lower().strip()
        if ("tumor" in name) or ("gtv" in name):
            return x, segments

    raise RuntimeError(
        "Tumor segment not found. Segment labels: " +
        ";".join([x["segment_label"] for x in segments])
    )

def write_sitk_image_uint8_like_reference(arr_uint8, reference_img, out_path):
    out_img = sitk.GetImageFromArray(arr_uint8)
    out_img.CopyInformation(reference_img)

    writer = sitk.ImageFileWriter()
    writer.SetFileName(out_path)
    writer.UseCompressionOn()
    writer.Execute(out_img)

    if not os.path.exists(out_path):
        raise RuntimeError("Failed to save tumor mask: " + out_path)

    return os.path.getsize(out_path)

############################################################
# 3. Load pair table
############################################################

pairs = pd.read_csv(pair_csv)
pairs["patient_id"] = pairs["patient_id"].astype(str).str.upper().str.strip()
pairs = pairs[pairs["pair_status"] == "CANDIDATE_PAIR_CREATED"].copy()
pairs = pairs.sort_values("patient_id").reset_index(drop=True)

print("\n===== Primary tumor-mask extraction =====")
print("Candidate pairs:", len(pairs))
print("Input NRRD full root:", nrrd_full_root)
print("Output tumor mask root:", tumor_mask_root)

############################################################
# 4. Main loop
############################################################

rows = []

for i, row in pairs.iterrows():

    pid = row["patient_id"]
    seg_file = row["seg_file_path"]

    image_path = os.path.join(
        nrrd_full_root,
        pid,
        f"{pid}_image.nrrd"
    )

    all_mask_path = os.path.join(
        nrrd_full_root,
        pid,
        f"{pid}_mask.nrrd"
    )

    out_dir = os.path.join(tumor_mask_root, pid)
    os.makedirs(out_dir, exist_ok=True)

    tumor_mask_path = os.path.join(
        out_dir,
        f"{pid}_tumor_mask.nrrd"
    )

    one = {
        "patient_id": pid,
        "status": "FAILED",
        "message": "",
        "image_path": image_path,
        "all_label_mask_path": all_mask_path,
        "tumor_mask_path": tumor_mask_path,
        "seg_file_path": seg_file,
        "tumor_segment_number": "",
        "tumor_segment_index_1based": "",
        "tumor_segment_label": "",
        "all_segment_labels": "",
        "mask_unique_labels": "",
        "tumor_label_used": "",
        "tumor_nonzero_voxels": "",
        "tumor_volume_cm3": "",
        "all_mask_nonzero_voxels": "",
        "all_mask_volume_cm3": "",
        "tumor_to_all_volume_ratio": "",
        "image_exists": os.path.exists(image_path),
        "all_mask_exists": os.path.exists(all_mask_path),
        "tumor_mask_file_size": ""
    }

    try:
        if (i + 1) % 25 == 0 or i == 0:
            print(f"Processing {i+1}/{len(pairs)}: {pid}")

        if not os.path.exists(image_path):
            raise FileNotFoundError("Image NRRD not found: " + image_path)

        if not os.path.exists(all_mask_path):
            raise FileNotFoundError("All-label mask NRRD not found: " + all_mask_path)

        if not os.path.exists(seg_file):
            raise FileNotFoundError("DICOM SEG file not found: " + seg_file)

        tumor_seg, all_segments = find_tumor_segment_from_dicom_seg(seg_file)

        one["tumor_segment_number"] = tumor_seg["segment_number"]
        one["tumor_segment_index_1based"] = tumor_seg["index_1based"]
        one["tumor_segment_label"] = tumor_seg["segment_label"]
        one["all_segment_labels"] = ";".join([x["segment_label"] for x in all_segments])

        mask_img = sitk.ReadImage(all_mask_path)
        arr = sitk.GetArrayFromImage(mask_img)

        labels = np.unique(arr)
        labels_int = [int(x) for x in labels.tolist()]
        labels_nonzero = [x for x in labels_int if x != 0]

        one["mask_unique_labels"] = ";".join([str(x) for x in labels_int])

        spacing = mask_img.GetSpacing()
        voxel_volume_mm3 = float(spacing[0] * spacing[1] * spacing[2])

        all_nonzero = int(np.count_nonzero(arr))
        all_volume_cm3 = all_nonzero * voxel_volume_mm3 / 1000.0

        one["all_mask_nonzero_voxels"] = all_nonzero
        one["all_mask_volume_cm3"] = all_volume_cm3

        # Primary assumption: label value equals DICOM SegmentNumber
        tumor_label = int(tumor_seg["segment_number"])

        # Fallback: if SegmentNumber not present in NRRD labels, use export order index
        if tumor_label not in labels_nonzero:
            fallback_label = int(tumor_seg["index_1based"])
            if fallback_label in labels_nonzero:
                tumor_label = fallback_label
            else:
                raise RuntimeError(
                    "Tumor label not found in all-label mask. "
                    f"SegmentNumber={tumor_seg['segment_number']}, "
                    f"index={tumor_seg['index_1based']}, "
                    f"mask labels={labels_nonzero}"
                )

        tumor_arr = (arr == tumor_label).astype(np.uint8)

        tumor_nonzero = int(np.count_nonzero(tumor_arr))
        tumor_volume_cm3 = tumor_nonzero * voxel_volume_mm3 / 1000.0

        one["tumor_label_used"] = tumor_label
        one["tumor_nonzero_voxels"] = tumor_nonzero
        one["tumor_volume_cm3"] = tumor_volume_cm3

        if all_volume_cm3 > 0:
            one["tumor_to_all_volume_ratio"] = tumor_volume_cm3 / all_volume_cm3

        if tumor_nonzero <= 0:
            raise RuntimeError("Extracted tumor-only mask is empty.")

        fsize = write_sitk_image_uint8_like_reference(
            tumor_arr,
            mask_img,
            tumor_mask_path
        )

        one["tumor_mask_file_size"] = fsize

        # QC flag
        # >1500 cm3 is unusually huge for primary tumor and should be checked.
        if tumor_volume_cm3 > 1500:
            one["status"] = "WARNING_LARGE_TUMOR"
            one["message"] = "Tumor mask extracted, but volume >1500 cm3. Needs visual check."
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

############################################################
# 5. Summary
############################################################

qc_df = pd.DataFrame(rows)

summary_rows = []

def add_summary(item, value):
    summary_rows.append({"item": item, "value": value})

add_summary("candidate_pairs", len(pairs))
add_summary("PASS", int((qc_df["status"] == "PASS").sum()))
add_summary("WARNING_LARGE_TUMOR", int((qc_df["status"] == "WARNING_LARGE_TUMOR").sum()))
add_summary("FAILED", int((qc_df["status"] == "FAILED").sum()))
add_summary("tumor_masks_created", int(qc_df["tumor_mask_path"].apply(os.path.exists).sum()))

vol = pd.to_numeric(qc_df["tumor_volume_cm3"], errors="coerce")
ratio = pd.to_numeric(qc_df["tumor_to_all_volume_ratio"], errors="coerce")

add_summary("tumor_volume_cm3_min", float(vol.min()))
add_summary("tumor_volume_cm3_q25", float(vol.quantile(0.25)))
add_summary("tumor_volume_cm3_median", float(vol.median()))
add_summary("tumor_volume_cm3_mean", float(vol.mean()))
add_summary("tumor_volume_cm3_q75", float(vol.quantile(0.75)))
add_summary("tumor_volume_cm3_max", float(vol.max()))

add_summary("tumor_to_all_volume_ratio_median", float(ratio.median()))
add_summary("tumor_to_all_volume_ratio_max", float(ratio.max()))

summary_df = pd.DataFrame(summary_rows)

print("\n===== Tumor-only FAST extraction status table =====")
print(qc_df["status"].value_counts(dropna=False))

print("\n===== Tumor segment label table =====")
print(qc_df["tumor_segment_label"].value_counts(dropna=False))

print("\n===== Tumor label used table =====")
print(qc_df["tumor_label_used"].value_counts(dropna=False).sort_index())

print("\n===== Tumor volume cm3 summary =====")
print(vol.describe())

print("\n===== Tumor-to-all-mask volume ratio summary =====")
print(ratio.describe())

print("\n===== Failed or warning cases =====")
print(qc_df.loc[
    qc_df["status"] != "PASS",
    ["patient_id", "status", "message", "tumor_segment_label", "tumor_label_used", "tumor_volume_cm3"]
].head(100))

############################################################
# 6. Save outputs
############################################################

qc_df.to_csv(out_csv, index=False, encoding="utf-8-sig")
summary_df.to_csv(summary_csv, index=False, encoding="utf-8-sig")

try:
    with pd.ExcelWriter(xlsx_out, engine="openpyxl") as writer:
        summary_df.to_excel(writer, sheet_name="summary", index=False)
        qc_df.to_excel(writer, sheet_name="case_QC", index=False)

        label_tab = qc_df["tumor_segment_label"].value_counts(dropna=False).reset_index()
        label_tab.columns = ["tumor_segment_label", "n_cases"]
        label_tab.to_excel(writer, sheet_name="tumor_label_table", index=False)

        label_value_tab = qc_df["tumor_label_used"].value_counts(dropna=False).reset_index()
        label_value_tab.columns = ["tumor_label_used", "n_cases"]
        label_value_tab.to_excel(writer, sheet_name="label_value_table", index=False)

except Exception as e:
    print("Excel writing failed, CSV outputs are still saved.")
    print(str(e))

print("\n===== Primary tumor-mask extraction completed =====")
print("Main output files:")
print(out_csv)
print(summary_csv)
print(xlsx_out)
print("Tumor mask root:")
print(tumor_mask_root)

print("\nKey output files:")
print("1. Tumor-only FAST extraction status table")
print("2. Tumor segment label table")
print("3. Tumor label used table")
print("4. Tumor volume cm3 summary")
print("5. Failed or warning cases")
