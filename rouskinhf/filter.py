import numpy as np
from .util import UKN

import pandas as pd
from sklearn.metrics import roc_auc_score
from .list_datapoints import ListofDatapoints
from .datapoint import Datapoint


def filter(listofdatapoints: ListofDatapoints, min_AUROC: int = 0.8):
    """Filters out duplicate sequences.
        Only keep the first occurence of a sequence if all the other structures are the same.

        Examples:
            >>> datapoints = ListofDatapoints([ Datapoint(reference='ref1', sequence='AACCGG', structure=[[1, 2], [3, 4]], dms=[1,0,0,0,0,1]),\
                                Datapoint(reference='ref2', sequence='AACCGG', structure=[[1, 2], [3, 4]], dms=[1,0,0,0,0,1]),\
                                Datapoint(reference='ref2', sequence='AACCGG', structure=[[1, 2], [3, 4]], dms=[1,0,0,0,0,1]),\
                                Datapoint(reference='ref4', sequence='AUGGC', structure=[[1, 2]], dms=[0,0,0,0,1]),\
                                Datapoint(reference='ref5', sequence='AUGGC', structure=[[0, 4], [2, 3]], dms=[0,0,0,0,0]),\
                                Datapoint(reference='ref6', sequence='not a regular sequence', structure=[[0, 4]], dms=[0,0,0,0,0]) ])
            >>> print(filter(datapoints, min_AUROC=0))
            Over a total of 6 datapoints, there are:
            ### OUTPUT
            - 3 valid datapoints
            ### MODIFIED
            - 1 multiple sequences with the same reference (renamed reference)
            - 1 duplicate sequences with different structure / dms / shape
            ### FILTERED OUT
            - 1 invalid datapoints (ex: sequence with non-regular characters)
            - 0 datapoints with bad structures
            - 2 duplicate sequences with the same structure / dms / shape
            - 0 datapoints removed because of low AUROC (<0)
        """

    # Remove None datapoints
    n_input_datapoints = len(listofdatapoints)
    datapoints, n_unvalid_datapoints = listofdatapoints.drop_none_dp()

    # Remove bad structures
    bad_structures_idx = []
    for idx, datapoint in enumerate(datapoints):
        if not datapoint._assert_structure():
            bad_structures_idx.append(idx)
    datapoints = [
        datapoint
        for idx, datapoint in enumerate(datapoints)
        if idx not in bad_structures_idx
    ]
    n_bad_structures_datapoints = len(bad_structures_idx)

    # Convert to pandas
    df = listofdatapoints.to_pandas(datapoints)

    # Remove duplicate or conflicting datapoints and keep track of the number of datapoints removed
    def drop_duplicates(df: pd.DataFrame, **kwargs):
        len_df_before = len(df)
        df.drop_duplicates(**kwargs)
        return len_df_before - len(df)

    # If multiple sequences with the same reference, rename the reference
    refs = dict()
    n_same_ref_datapoints = 0
    for idx, row in df.iterrows():
        if row["reference"] in refs:
            refs[row["reference"]] += 1
            df.at[idx, "reference"] = f"{row['reference']}_{refs[row['reference']]}"
            n_same_ref_datapoints += 1
        else:
            refs[row["reference"]] = 0

    # Keep only one datapoint per sequence and structure
    if "structure" in df.columns:
        n_duplicates_datapoints = drop_duplicates(
            df,
            subset=["sequence", "structure"],
            inplace=True,
            ignore_index=True,
            keep="first",
        )

    # Keep only one datapoint per sequence and signal
    for signal in ["dms", "shape"]:
        if signal in df.columns:
            if not "structure" in df.columns:
                n_duplicates_datapoints = 0
            n_duplicates_datapoints += drop_duplicates(
                df,
                subset=["sequence", signal],
                inplace=True,
                ignore_index=True,
                keep="first",
            )

    # Count how many multiple structures / dms with the same sequence
    n_same_seq_datapoints = 0
    for _, group in df.groupby("sequence"):
        if len(group) > 1:
            n_same_seq_datapoints += len(group) - 1

    ## Filter out references with low AUROC
    mask_high_AUROC = None
    if "structure" in df.columns and "dms" in df.columns or "shape" in df.columns:
        signal = "dms" if "dms" in df.columns else "shape"

        def calculate_auroc(row):
            # best signal is DMS, if unavailable use SHAPE
            signal = (
                "dms" if "dms" in row and not type(row["dms"]) == float else "shape"
            )
            if (
                (not "structure" in row)
                or (row["structure"] is None)
                or (type(row["structure"]) == float and np.isnan(row["structure"]))
                or len(row["structure"]) == 0
            ):
                return 1  #  if you can't calculate the AUROC, don't filter it out

            sig = np.array(row[signal])
            if signal == "dms":  # mask G/T/U bases for DMS
                mask_UKN = np.array([s in "GTU" for s in row["sequence"].upper()])
                sig[mask_UKN] = UKN

            isUnpaired = np.ones_like(sig)
            isUnpaired[np.array(row["structure"]).flatten()] = 0
            if set(isUnpaired[sig != UKN]) != set([0, 1]):
                return 0
            return roc_auc_score(isUnpaired[sig != UKN], sig[sig != UKN])

        # Create a boolean mask for rows with auroc score greater than or equal to a threshold
        mask_high_AUROC = df.apply(
            lambda row: calculate_auroc(row) >= min_AUROC, axis=1
        )
        df = df[mask_high_AUROC]

    # Convert back to list of datapoints
    listofdatapoints.datapoints = listofdatapoints.from_pandas(df)

    # Write report
    report = f"""Over a total of {n_input_datapoints} datapoints, there are:
### OUTPUT
- ALL: {len(listofdatapoints)} valid datapoints
- INCLUDED: {n_same_seq_datapoints} duplicate sequences with different structure / dms / shape
### MODIFIED
- {n_same_ref_datapoints} multiple sequences with the same reference (renamed reference)
### FILTERED OUT
- {n_unvalid_datapoints} invalid datapoints (ex: sequence with non-regular characters)
- {n_bad_structures_datapoints} datapoints with bad structures"""
    if "structure" in df.columns or "dms" in df.columns or "shape" in df.columns:
        report += f"""
- {n_duplicates_datapoints} duplicate sequences with the same structure / dms / shape"""
    if mask_high_AUROC is not None:
        report += f"""
- {np.sum(~mask_high_AUROC)} datapoints removed because of low AUROC (<{min_AUROC})"""

    return report
