# -*- coding: utf-8 -*-

# extract_radiomic_features.py
#
# Extract original PyRadiomics features for Lung1.
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

############################################################
# 1. Paths
############################################################

lung1_root = required_env_directory("LUNG1_DATA_DIR")
pair_csv = os.path.join(
    lung1_root,
    "metadata",
    "lung1_ct_seg_candidate_pairs.csv"
)

image_root = os.path.join(
    lung1_root,
    "NRRD_full"
)

mask_root = os.path.join(
    lung1_root,
    "NRRD_tumor_only_final"
)

out_dir = os.path.join(
    lung1_root,
    "radiomics_features"
)

os.makedirs(out_dir, exist_ok=True)

out_csv = os.path.join(
    out_dir,
    "lung1_radiomics_features.csv"
)

out_xlsx = os.path.join(
    out_dir,
    "lung1_radiomics_features.xlsx"
)

log_csv = os.path.join(
    out_dir,
    "lung1_radiomics_features_log.csv"
)

############################################################
# 2. Load patient list
############################################################

pairs = pd.read_csv(pair_csv)
pairs["patient_id"] = pairs["patient_id"].astype(str).str.upper().str.strip()
pairs = pairs[pairs["pair_status"] == "CANDIDATE_PAIR_CREATED"].copy()
pairs = pairs.sort_values("patient_id").reset_index(drop=True)

patient_ids = pairs["patient_id"].tolist()

print("\n===== Lung1 patient list =====")
print("Number of patients:", len(patient_ids))
print("First 10:", patient_ids[:10])

############################################################
# 3. PyRadiomics settings
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
# 4. Helper functions
############################################################

def clean_value(v):
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
# 5. Main loop
############################################################

feature_rows = []
log_rows = []

print("\n===== Lung1 PyRadiomics extraction started =====")

for i, pid in enumerate(patient_ids, start=1):

    print("\n==================================================")
    print(f"Processing {i}/{len(patient_ids)}: {pid}")

    image_path = os.path.join(
        image_root,
        pid,
        f"{pid}_image.nrrd"
    )

    mask_path = os.path.join(
        mask_root,
        pid,
        f"{pid}_tumor_mask.nrrd"
    )

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
            raise FileNotFoundError(f"Tumor mask not found: {mask_path}")

        result = extractor.execute(image_path, mask_path)

        row = {
            "patient_id": pid
        }

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
# 6. Final save
############################################################

feature_df = pd.DataFrame(feature_rows)
log_df = pd.DataFrame(log_rows)

feature_df.to_csv(out_csv, index=False, encoding="utf-8-sig")
log_df.to_csv(log_csv, index=False, encoding="utf-8-sig")

with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
    feature_df.to_excel(writer, sheet_name="features", index=False)
    log_df.to_excel(writer, sheet_name="log", index=False)

############################################################
# 7. Summary
############################################################

print("\n===== Lung1 PyRadiomics extraction completed =====")
print("Feature table:", out_csv)
print("Excel table  :", out_xlsx)
print("Log table    :", log_csv)

print("\nStatus table:")
print(log_df["status"].value_counts(dropna=False))

print("\nFeature dataframe shape:")
print(feature_df.shape)

if len(feature_df) > 0:
    diagnostics_cols = [
        c for c in feature_df.columns
        if c.startswith("diagnostics_")
    ]

    radiomics_cols = [
        c for c in feature_df.columns
        if c.startswith("original_")
    ]

    print("\nNumber of diagnostics columns:")
    print(len(diagnostics_cols))

    print("\nNumber of original radiomics feature columns:")
    print(len(radiomics_cols))

    print("\nFirst 20 radiomics columns:")
    print(radiomics_cols[:20])

print("\nFailed patients:")
if "status" in log_df.columns:
    print(log_df[log_df["status"] == "FAILED"][["patient_id", "message"]])
