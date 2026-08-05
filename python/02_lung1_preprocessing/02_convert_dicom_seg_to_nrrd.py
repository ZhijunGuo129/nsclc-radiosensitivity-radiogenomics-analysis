# -*- coding: utf-8 -*-

# 02_convert_dicom_seg_to_nrrd.py
#
# Convert Lung1 CT and DICOM SEG series to NRRD in 3D Slicer.
#
# Paths are configured through environment variables.

import os
import csv
import traceback
import numpy as np

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
pair_csv = os.path.join(
    lung1_root,
    "metadata",
    "lung1_ct_seg_candidate_pairs.csv"
)

out_root = os.path.join(
    lung1_root,
    "NRRD_full"
)

log_csv = os.path.join(
    lung1_root,
    "metadata",
    "lung1_dicomseg_conversion_log.csv"
)

os.makedirs(out_root, exist_ok=True)

############################################################
# 2. Settings
############################################################

# Skip cases with valid existing image and mask files to support restart.
skip_existing_success = True

# Number of processed cases between log writes.
save_log_every_n_cases = 1

############################################################
# 3. Helper functions
############################################################

def read_pairs(pair_csv):
    rows = []
    with open(pair_csv, "r", encoding="utf-8-sig", newline="") as f:
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

def export_segmentation_to_labelmap(seg_node, volume_node, labelmap_node):
    try:
        seg_node.SetReferenceImageGeometryParameterFromVolumeNode(volume_node)
    except Exception:
        pass

    logic = slicer.modules.segmentations.logic()

    try:
        logic.ExportAllSegmentsToLabelmapNode(
            seg_node,
            labelmap_node,
            slicer.vtkSegmentation.EXTENT_REFERENCE_GEOMETRY
        )
        return "Export method 1 success"
    except Exception as e1:
        try:
            logic.ExportAllSegmentsToLabelmapNode(
                seg_node,
                labelmap_node
            )
            return "Export method 2 success"
        except Exception as e2:
            raise RuntimeError(
                "Export segmentation failed. Method1: {} | Method2: {}".format(
                    str(e1),
                    str(e2)
                )
            )

def save_node_checked(node, path):
    ok = slicer.util.saveNode(node, path)
    if not ok:
        raise RuntimeError("Failed to save node to: {}".format(path))

    if not os.path.exists(path):
        raise RuntimeError("Output file not found after save: {}".format(path))

    return os.path.getsize(path)

def write_log(log_rows, log_csv):
    if len(log_rows) == 0:
        return

    with open(log_csv, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(log_rows[0].keys()))
        writer.writeheader()
        writer.writerows(log_rows)

def existing_files_look_valid(image_out, mask_out):
    if not os.path.exists(image_out):
        return False

    if not os.path.exists(mask_out):
        return False

    try:
        image_size = os.path.getsize(image_out)
        mask_size = os.path.getsize(mask_out)

        if image_size > 10000 and mask_size > 1000:
            return True
        else:
            return False

    except Exception:
        return False

############################################################
# 4. Load pair table
############################################################

if not os.path.exists(pair_csv):
    raise FileNotFoundError("Pair CSV not found: {}".format(pair_csv))

all_pairs = read_pairs(pair_csv)

pairs = [
    r for r in all_pairs
    if r.get("pair_status", "") == "CANDIDATE_PAIR_CREATED"
]

pairs = sorted(pairs, key=lambda x: x.get("patient_id", ""))

print("\n===== Lung1 DICOM SEG conversion =====")
print("Candidate pairs:", len(pairs))
print("Output root:", out_root)
print("Log CSV:", log_csv)

############################################################
# 5. Main loop
############################################################

log_rows = []

for i, row in enumerate(pairs, start=1):

    pid = row["patient_id"].upper().strip()

    ct_series_dir = row["ct_series_dir"]
    seg_series_dir = row["seg_series_dir"]

    out_dir = os.path.join(out_root, pid)
    os.makedirs(out_dir, exist_ok=True)

    image_out = os.path.join(out_dir, "{}_image.nrrd".format(pid))
    mask_out = os.path.join(out_dir, "{}_mask.nrrd".format(pid))

    one_log = {
        "patient_id": pid,
        "status": "FAILED",
        "message": "",
        "ct_series_dir": ct_series_dir,
        "seg_series_dir": seg_series_dir,
        "image_out": image_out,
        "mask_out": mask_out,
        "image_size": "",
        "mask_size": "",
        "mask_shape": "",
        "mask_min": "",
        "mask_max": "",
        "mask_nonzero": "",
        "n_segments": "",
        "case_index": i,
        "total_cases": len(pairs)
    }

    print("\n==================================================")
    print("Processing {}/{}: {}".format(i, len(pairs), pid))
    print("CT :", ct_series_dir)
    print("SEG:", seg_series_dir)

    try:
        if skip_existing_success and existing_files_look_valid(image_out, mask_out):
            one_log["status"] = "SKIPPED_EXISTING"
            one_log["message"] = "Existing image and mask files found; skipped"
            one_log["image_size"] = os.path.getsize(image_out)
            one_log["mask_size"] = os.path.getsize(mask_out)

            print("SKIPPED_EXISTING:", pid)
            print("Image:", image_out, one_log["image_size"])
            print("Mask :", mask_out, one_log["mask_size"])

            log_rows.append(one_log)

            if i % save_log_every_n_cases == 0:
                write_log(log_rows, log_csv)

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
                raise RuntimeError("DICOMUtils.loadPatientByUID loaded no nodes.")

        volume_node = get_first_scalar_volume()
        seg_node = get_first_segmentation()

        if volume_node is None:
            raise RuntimeError("No scalar volume node loaded.")

        if seg_node is None:
            raise RuntimeError("No segmentation node loaded.")

        n_segments = seg_node.GetSegmentation().GetNumberOfSegments()
        one_log["n_segments"] = n_segments

        print("Selected volume node:", volume_node.GetName())
        print("Selected segmentation node:", seg_node.GetName())
        print("Number of segments:", n_segments)

        if n_segments < 1:
            raise RuntimeError("Segmentation has zero segments.")

        labelmap_node = slicer.mrmlScene.AddNewNodeByClass(
            "vtkMRMLLabelMapVolumeNode",
            "{}_mask".format(pid)
        )

        export_msg = export_segmentation_to_labelmap(
            seg_node,
            volume_node,
            labelmap_node
        )

        arr = slicer.util.arrayFromVolume(labelmap_node)

        one_log["mask_shape"] = str(arr.shape)
        one_log["mask_min"] = int(np.min(arr))
        one_log["mask_max"] = int(np.max(arr))
        one_log["mask_nonzero"] = int(np.count_nonzero(arr))

        print(export_msg)
        print("Mask shape:", one_log["mask_shape"])
        print("Mask min:", one_log["mask_min"])
        print("Mask max:", one_log["mask_max"])
        print("Mask nonzero:", one_log["mask_nonzero"])

        if one_log["mask_nonzero"] <= 0:
            raise RuntimeError("Exported mask is empty.")

        image_size = save_node_checked(volume_node, image_out)
        mask_size = save_node_checked(labelmap_node, mask_out)

        one_log["image_size"] = image_size
        one_log["mask_size"] = mask_size
        one_log["status"] = "SUCCESS"
        one_log["message"] = "OK"

        print("SUCCESS:", pid)
        print("Image:", image_out, image_size)
        print("Mask :", mask_out, mask_size)

    except Exception as e:
        one_log["status"] = "FAILED"
        one_log["message"] = str(e)

        print("FAILED:", pid)
        print(str(e))
        traceback.print_exc()

    log_rows.append(one_log)

    if i % save_log_every_n_cases == 0:
        write_log(log_rows, log_csv)

############################################################
# 6. Final log
############################################################

write_log(log_rows, log_csv)

############################################################
# 7. Summary
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

print("\nSuccessful or skipped output count:")
ok_n = sum([1 for r in log_rows if r["status"] in ["SUCCESS", "SKIPPED_EXISTING"]])
print(ok_n)

print("\nLog saved to:")
print(log_csv)

print("\n===== Lung1 DICOM SEG conversion completed =====")
