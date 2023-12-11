from .list_datapoints import ListofDatapoints
from .path import Path
from .filter import filter as filter_datapoints

def convert(
    format: str,
    file_or_folder: str,
    name: str = None,
    path_out: str = "data",
    predict_structure: bool = False,
    filter: bool = True,
    min_AUROC=0.8,
    verbose: bool = True,
):
    """Converts a file or folder into a json file. Different formats are supported.
    
    Args:
        format (str): Format of the input file. Can be 'ct', 'seismic', 'json', 'bpseq' or 'fasta'.
        file_or_folder (str): Path to the file or folder to convert.
        name (str, optional): Name of the dataset. Defaults to None, in which case the name of the file or folder will be used.
        path_out (str, optional): Path to the output folder. Defaults to 'data'.
        predict_structure (bool, optional): Whether to predict the structure or not using RNAstructure. Defaults to False.
        filter (bool, optional): Whether to filter the datapoints or not. Defaults to True. Datapoints with no sequence or reference will be dropped anyways.
        min_AUROC (float, optional): Minimum AUROC to keep a datapoint. Defaults to 0.8.
        verbose (bool, optional): Whether to print the conversion report or not. Defaults to True.
    """
    assert format in ["ct", "seismic", "json", "bpseq", "fasta"], "Format not supported"

    if name is None:
        name = file_or_folder.split("/")[-1].split(".")[0]
    path = Path(name=name, root=path_out)
    path.make()

    if format == "ct":
        datapoints = ListofDatapoints.from_ct(
            file_or_folder, tqdm=True, verbose=verbose
        )

    elif format == "seismic":
        datapoints = ListofDatapoints.from_dreem_output(
            file_or_folder, predict_structure, tqdm=True, verbose=verbose
        )

    elif format == "json":
        datapoints = ListofDatapoints.from_json(
            file_or_folder, predict_structure, tqdm=True, verbose=verbose
        )

    elif format == "bpseq":
        datapoints = ListofDatapoints.from_bpseq(
            file_or_folder, tqdm=True, verbose=verbose
        )

    elif format == "fasta":
        datapoints = ListofDatapoints.from_fasta(
            file_or_folder, predict_structure, tqdm=True, verbose=verbose
        )

    if filter:
        report = filter_datapoints(datapoints, min_AUROC=min_AUROC)
    else:
        _, report = datapoints.drop_none_dp()
        report = (
            f"Drop {report} datapoints with None values (null sequence or reference)"
        )

    if verbose:
        print(report)

    with open(path.get_conversion_report(), "w") as f:
        f.write("# Conversion report \n\n")
        f.write(report)

    if path_out is not None:
        datapoints.to_json(path.get_data_json())

    return datapoints.to_dict()
