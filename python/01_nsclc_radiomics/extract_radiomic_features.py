# -*- coding: utf-8 -*-

# extract_radiomic_features.py
#
# Extract original PyRadiomics features for NSCLC-Radiogenomics.
#
# Paths are configured through environment variables.

import os
import traceback
import pandas as pd
from radiomics import featureextractor

def required_env_directory(name):
    value = os.environ.get(name, "").strip()
    if not value:
        raise RuntimeError(f"Environment variable {name} is not set.")
    path = os.path.abspath(os.path.expanduser(value))
    if not os.path.isdir(path):
        raise FileNotFoundError(f"Directory specified by {name} does not exist: {path}")
    return path

project_dir = required_env_directory("RADIOGENOMICS_PROJECT_DIR")
nrrd_root = os.path.join(
    project_dir,
    "01_raw_data",
    "NSCLC_Radiogenomics",
    "NRRD_full"
)

pair_csv = os.path.join(
    project_dir,
    "02_metadata",
    "ct_seg_candidate_pairs.csv"
)

out_dir = os.path.join(
    project_dir,
    "07_results",
    "radiomics_features"
)

os.makedirs(out_dir, exist_ok=True)

out_csv = os.path.join(out_dir, "nsclc_radiomics_features.csv")
out_xlsx = os.path.join(out_dir, "nsclc_radiomics_features.xlsx")
log_csv = os.path.join(out_dir, "nsclc_radiomics_features_log.csv")

############################################################
# Load patient list
############################################################

pairs = pd.read_csv(pair_csv)

if "pair_status" in pairs.columns:
    pairs = pairs[pairs["pair_status"] == "CANDIDATE_PAIR_CREATED"].copy()

patient_ids = list(pairs["patient_id"].astype(str))

print("\n===== Patient list =====")
print("Number of patients from pair table:", len(patient_ids))
print(patient_ids[:10])

############################################################
# PyRadiomics settings
############################################################

settings = {
    "binWidth": 25,
    "label": 1,
    "interpolator": "sitkBSpline",
    "resampledPixelSpacing": [1, 1, 1],
    "correctMask": True,
    "geometryTolerance": 1e-4,
    "normalize": False,
    "removeOutliers": None
}

extractor = featureextractor.RadiomicsFeatureExtractor(**settings)

# Extract original-image features only.
extractor.disableAllImageTypes()
extractor.enableImageTypeByName("Original")

extractor.disableAllFeatures()
extractor.enableFeatureClassByName("firstorder")
extractor.enableFeatureClassByName("shape")
extractor.enableFeatureClassByName("glcm")
extractor.enableFeatureClassByName("glrlm")
extractor.enableFeatureClassByName("glszm")
extractor.enableFeatureClassByName("gldm")
extractor.enableFeatureClassByName("ngtdm")

print("\n===== PyRadiomics extractor settings =====")
print(extractor.settings)

print("\nEnabled image types:")
print(extractor.enabledImagetypes)

print("\nEnabled features:")
print(extractor.enabledFeatures)

############################################################
# Helper function
############################################################

def clean_value(v):
    """
    Convert PyRadiomics output values to CSV/Excel-friendly values.
    """
    try:
        if hasattr(v, "item"):
            return v.item()
    except Exception:
        pass

    try:
        if isinstance(v, (list, tuple)):
            return ";".join([str(x) for x in v])
    except Exception:
        pass

    return v

def save_current(feature_rows, log_rows):
    feature_df = pd.DataFrame(feature_rows)
    log_df = pd.DataFrame(log_rows)

    feature_df.to_csv(out_csv, index=False, encoding="utf-8-sig")
    log_df.to_csv(log_csv, index=False, encoding="utf-8-sig")

############################################################
# Main loop
############################################################

feature_rows = []
log_rows = []

print("\n===== NSCLC-Radiogenomics PyRadiomics extraction started =====")

for i, pid in enumerate(patient_ids, start=1):
    print("\n==================================================")
    print(f"Processing {i}/{len(patient_ids)}: {pid}")

    image_path = os.path.join(nrrd_root, pid, f"{pid}_image.nrrd")
    mask_path = os.path.join(nrrd_root, pid, f"{pid}_mask.nrrd")

    one_log = {
        "patient_id": pid,
        "status": "FAILED",
        "message": "",
        "image_path": image_path,
        "mask_path": mask_path,
        "n_features_total": None,
        "n_radiomics_features": None
    }

    try:
        if not os.path.exists(image_path):
            raise FileNotFoundError(f"Image not found: {image_path}")

        if not os.path.exists(mask_path):
            raise FileNotFoundError(f"Mask not found: {mask_path}")

        result = extractor.execute(image_path, mask_path)

        row = {"patient_id": pid}

        for k, v in result.items():
            row[k] = clean_value(v)

        radiomics_feature_names = [
            k for k in row.keys()
            if k.startswith("original_")
        ]

        feature_rows.append(row)

        one_log["status"] = "SUCCESS"
        one_log["message"] = "OK"
        one_log["n_features_total"] = len(row)
        one_log["n_radiomics_features"] = len(radiomics_feature_names)

        print("SUCCESS:", pid)
        print("Total output columns:", len(row))
        print("Radiomics feature columns:", len(radiomics_feature_names))

    except Exception as e:
        msg = str(e)
        one_log["status"] = "FAILED"
        one_log["message"] = msg
        print("FAILED:", pid)
        print(msg)
        traceback.print_exc()

    log_rows.append(one_log)

    # Save after each case to support safe restart.
    save_current(feature_rows, log_rows)

############################################################
# Save final outputs
############################################################

feature_df = pd.DataFrame(feature_rows)
log_df = pd.DataFrame(log_rows)

feature_df.to_csv(out_csv, index=False, encoding="utf-8-sig")
log_df.to_csv(log_csv, index=False, encoding="utf-8-sig")

with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
    feature_df.to_excel(writer, sheet_name="features", index=False)
    log_df.to_excel(writer, sheet_name="log", index=False)

print("\n===== NSCLC-Radiogenomics PyRadiomics extraction completed =====")
print("Feature table:", out_csv)
print("Excel table  :", out_xlsx)
print("Log table    :", log_csv)

print("\nStatus table:")
print(log_df["status"].value_counts(dropna=False))

print("\nFeature dataframe shape:")
print(feature_df.shape)

if len(feature_df) > 0:
    radiomics_cols = [c for c in feature_df.columns if c.startswith("original_")]
    diagnostics_cols = [c for c in feature_df.columns if c.startswith("diagnostics_")]

    print("\nNumber of diagnostics columns:")
    print(len(diagnostics_cols))

    print("\nNumber of original radiomics feature columns:")
    print(len(radiomics_cols))

    print("\nFirst 20 radiomics columns:")
    print(radiomics_cols[:20])

print("\nFailed patients:")
if "status" in log_df.columns:
    print(log_df[log_df["status"] == "FAILED"][["patient_id", "message"]])
