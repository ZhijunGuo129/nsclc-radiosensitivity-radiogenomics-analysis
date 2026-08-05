# -*- coding: utf-8 -*-

# 06_export_selected_masks_in_slicer.py
#
# Re-export selected primary-tumor masks in 3D Slicer.
#
# Paths are configured through environment variables.

import os
import csv
import traceback
import numpy as np
import vtk

import slicer
from DICOMLib import DICOMUtils


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
target_csv = os.path.join(
    lung1_root,
    "metadata",
    "lung1_selected_tumor_mask_targets.csv"
)

log_csv = os.path.join(
    lung1_root,
    "metadata",
    "lung1_targeted_mask_export_log.csv"
)

skip_existing_success = True

############################################################
# 2. Helper functions
############################################################

def read_targets(target_csv):
    rows = []
    with open(target_csv, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


def clear_scene():
    slicer.mrmlScene.Clear(0)


def get_first_scalar_volume():
    volume_nodes = slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode")

    if len(volume_nodes) == 0:
        return None

    best_node = None
    best_size = -1

    for node in volume_nodes:
        try:
            img = node.GetImageData()
            if img is None:
                continue

            dims = img.GetDimensions()
            size = int(dims[0]) * int(dims[1]) * int(dims[2])

            if size > best_size:
                best_size = size
                best_node = node

        except Exception:
            pass

    return best_node


def get_first_segmentation():
    seg_nodes = slicer.util.getNodesByClass("vtkMRMLSegmentationNode")

    if len(seg_nodes) == 0:
        return None

    best_node = None
    best_count = -1

    for node in seg_nodes:
        try:
            nseg = node.GetSegmentation().GetNumberOfSegments()

            if nseg > best_count:
                best_count = nseg
                best_node = node

        except Exception:
            pass

    return best_node


def get_segment_ids_and_names(seg_node):
    segmentation = seg_node.GetSegmentation()
    nseg = segmentation.GetNumberOfSegments()

    out = []

    for i in range(nseg):
        seg_id = segmentation.GetNthSegmentID(i)
        seg = segmentation.GetSegment(seg_id)
        seg_name = seg.GetName()

        out.append({
            "segment_index": i,
            "segment_id": seg_id,
            "segment_name": seg_name
        })

    return out


def find_tumor_segment_id(seg_node):
    seg_info = get_segment_ids_and_names(seg_node)

    # 1. Exact match to Neoplasm, Primary.
    for x in seg_info:
        name_low = x["segment_name"].lower().strip()

        if name_low == "neoplasm, primary":
            return x["segment_id"], x["segment_name"], seg_info

    # 2. Match labels containing both neoplasm and primary.
    for x in seg_info:
        name_low = x["segment_name"].lower().strip()

        if ("neoplasm" in name_low) and ("primary" in name_low):
            return x["segment_id"], x["segment_name"], seg_info

    # 3. Fallback match to tumor or GTV.
    for x in seg_info:
        name_low = x["segment_name"].lower().strip()

        if ("tumor" in name_low) or ("gtv" in name_low):
            return x["segment_id"], x["segment_name"], seg_info

    return None, None, seg_info


def export_one_segment_to_labelmap(seg_node, volume_node, tumor_segment_id, labelmap_node):
    try:
        seg_node.SetReferenceImageGeometryParameterFromVolumeNode(volume_node)
    except Exception:
        pass

    segment_ids = vtk.vtkStringArray()
    segment_ids.InsertNextValue(tumor_segment_id)

    logic = slicer.modules.segmentations.logic()

    try:
        logic.ExportSegmentsToLabelmapNode(
            seg_node,
            segment_ids,
            labelmap_node,
            volume_node,
            slicer.vtkSegmentation.EXTENT_REFERENCE_GEOMETRY
        )

        return "ExportSegmentsToLabelmapNode success"

    except Exception as e1:

        try:
            display_node = seg_node.GetDisplayNode()
            all_info = get_segment_ids_and_names(seg_node)

            for x in all_info:
                display_node.SetSegmentVisibility(x["segment_id"], False)

            display_node.SetSegmentVisibility(tumor_segment_id, True)

            logic.ExportVisibleSegmentsToLabelmapNode(
                seg_node,
                labelmap_node,
                volume_node,
                slicer.vtkSegmentation.EXTENT_REFERENCE_GEOMETRY
            )

            return "ExportVisibleSegmentsToLabelmapNode success"

        except Exception as e2:
            raise RuntimeError(
                "Tumor-only export failed. Method1: {} | Method2: {}".format(
                    str(e1),
                    str(e2)
                )
            )


def binarize_labelmap(labelmap_node):
    arr = slicer.util.arrayFromVolume(labelmap_node)
    bin_arr = (arr > 0).astype(np.uint8)

    slicer.util.updateVolumeFromArray(labelmap_node, bin_arr)

    try:
        slicer.util.arrayFromVolumeModified(labelmap_node)
    except Exception:
        pass

    return bin_arr


def save_node_checked(node, path):
    ok = slicer.util.saveNode(node, path)

    if not ok:
        raise RuntimeError("Failed to save node to: {}".format(path))

    if not os.path.exists(path):
        raise RuntimeError("Output file not found after save: {}".format(path))

    return os.path.getsize(path)


def write_log(log_rows):
    if len(log_rows) == 0:
        return

    with open(log_csv, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(log_rows[0].keys()))
        writer.writeheader()
        writer.writerows(log_rows)


def existing_mask_looks_valid(path):
    if not os.path.exists(path):
        return False

    try:
        return os.path.getsize(path) > 1000
    except Exception:
        return False


############################################################
# 3. Load target cases
############################################################

if not os.path.exists(target_csv):
    raise FileNotFoundError("Target CSV not found: {}".format(target_csv))

targets = read_targets(target_csv)
targets = sorted(targets, key=lambda x: x.get("patient_id", ""))

print("\n===== Targeted Slicer tumor-mask export =====")
print("Target cases:", len(targets))
print("Target CSV:", target_csv)
print("Log CSV:", log_csv)

############################################################
# 4. Main loop
############################################################

log_rows = []

for i, row in enumerate(targets, start=1):

    pid = row["patient_id"].upper().strip()
    ct_series_dir = row["ct_series_dir"]
    seg_series_dir = row["seg_series_dir"]
    final_mask_path = row["final_mask_path"]

    os.makedirs(os.path.dirname(final_mask_path), exist_ok=True)

    one_log = {
        "patient_id": pid,
        "status": "FAILED",
        "message": "",
        "case_index": i,
        "total_cases": len(targets),
        "ct_series_dir": ct_series_dir,
        "seg_series_dir": seg_series_dir,
        "final_mask_path": final_mask_path,
        "tumor_mask_size": "",
        "tumor_mask_shape": "",
        "tumor_mask_min": "",
        "tumor_mask_max": "",
        "tumor_mask_nonzero": "",
        "tumor_segment_id": "",
        "tumor_segment_name": "",
        "all_segment_names": "",
        "n_segments": "",
        "fast_status": row.get("fast_status", ""),
        "fast_message": row.get("fast_message", "")
    }

    print("\n==================================================")
    print("Processing {}/{}: {}".format(i, len(targets), pid))
    print("CT :", ct_series_dir)
    print("SEG:", seg_series_dir)
    print("Final mask:", final_mask_path)

    try:
        if skip_existing_success and existing_mask_looks_valid(final_mask_path):
            one_log["status"] = "SKIPPED_EXISTING"
            one_log["message"] = "Existing targeted tumor mask found; skipped"
            one_log["tumor_mask_size"] = os.path.getsize(final_mask_path)

            print("SKIPPED_EXISTING:", pid)

            log_rows.append(one_log)
            write_log(log_rows)
            continue

        clear_scene()

        if not os.path.isdir(ct_series_dir):
            raise FileNotFoundError("CT series dir not found: {}".format(ct_series_dir))

        if not os.path.isdir(seg_series_dir):
            raise FileNotFoundError("SEG series dir not found: {}".format(seg_series_dir))

        with DICOMUtils.TemporaryDICOMDatabase() as db:

            print("Importing CT DICOM...")
            DICOMUtils.importDicom(ct_series_dir, db)

            print("Importing SEG DICOM...")
            DICOMUtils.importDicom(seg_series_dir, db)

            patient_uids = db.patients()

            print("DICOM database patient UID count:", len(patient_uids))

            if len(patient_uids) == 0:
                raise RuntimeError("No DICOM patients imported.")

            loaded_any = False

            for patient_uid in patient_uids:
                print("Loading patient UID:", patient_uid)

                loaded_node_ids = DICOMUtils.loadPatientByUID(patient_uid)

                print("Loaded node IDs:", loaded_node_ids)

                if loaded_node_ids:
                    loaded_any = True

            if not loaded_any:
                raise RuntimeError("No nodes loaded.")

        volume_node = get_first_scalar_volume()
        seg_node = get_first_segmentation()

        if volume_node is None:
            raise RuntimeError("No scalar volume node loaded.")

        if seg_node is None:
            raise RuntimeError("No segmentation node loaded.")

        seg_info = get_segment_ids_and_names(seg_node)
        all_names = [x["segment_name"] for x in seg_info]

        tumor_segment_id, tumor_segment_name, seg_info = find_tumor_segment_id(seg_node)

        one_log["n_segments"] = len(seg_info)
        one_log["all_segment_names"] = ";".join(all_names)

        print("All segment names:")
        print(all_names)

        if tumor_segment_id is None:
            raise RuntimeError("Tumor segment not found. Segment names: {}".format(all_names))

        one_log["tumor_segment_id"] = tumor_segment_id
        one_log["tumor_segment_name"] = tumor_segment_name

        print("Selected tumor segment:", tumor_segment_name)

        labelmap_node = slicer.mrmlScene.AddNewNodeByClass(
            "vtkMRMLLabelMapVolumeNode",
            "{}_tumor_mask".format(pid)
        )

        export_msg = export_one_segment_to_labelmap(
            seg_node,
            volume_node,
            tumor_segment_id,
            labelmap_node
        )

        arr = binarize_labelmap(labelmap_node)

        one_log["tumor_mask_shape"] = str(arr.shape)
        one_log["tumor_mask_min"] = int(np.min(arr))
        one_log["tumor_mask_max"] = int(np.max(arr))
        one_log["tumor_mask_nonzero"] = int(np.count_nonzero(arr))

        print(export_msg)
        print("Tumor mask shape:", one_log["tumor_mask_shape"])
        print("Tumor mask min:", one_log["tumor_mask_min"])
        print("Tumor mask max:", one_log["tumor_mask_max"])
        print("Tumor mask nonzero:", one_log["tumor_mask_nonzero"])

        if one_log["tumor_mask_nonzero"] <= 0:
            raise RuntimeError("Tumor-only mask is empty.")

        mask_size = save_node_checked(labelmap_node, final_mask_path)

        one_log["tumor_mask_size"] = mask_size
        one_log["status"] = "SUCCESS"
        one_log["message"] = "OK"

        print("SUCCESS:", pid)
        print("Tumor mask:", final_mask_path, mask_size)

    except Exception as e:
        one_log["status"] = "FAILED"
        one_log["message"] = str(e)

        print("FAILED:", pid)
        print(str(e))
        traceback.print_exc()

    log_rows.append(one_log)
    write_log(log_rows)

############################################################
# 5. Summary
############################################################

print("\n===== Summary =====")

status_counts = {}

for r in log_rows:
    status = r["status"]
    status_counts[status] = status_counts.get(status, 0) + 1

print("Status counts:")

for k, v in status_counts.items():
    print(k, v)

print("\nFailed patients:")

for r in log_rows:
    if r["status"] == "FAILED":
        print(r["patient_id"], r["message"])

ok_n = sum([
    1 for r in log_rows
    if r["status"] in ["SUCCESS", "SKIPPED_EXISTING"]
])

print("\nSuccessful or skipped targeted count:")
print(ok_n)

print("\nLog saved to:")
print(log_csv)

print("\n===== Targeted Slicer tumor-mask export completed =====")
