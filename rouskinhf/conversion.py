from .list_datapoints import ListofDatapoints
from .path import Path


def convert(
    format:str,
    file_or_folder:str,
    name:str=None,
    path_out:str='data',
    predict_structure:bool=False,
    filter:bool=True,
    min_AUROC=0.8,
    verbose:bool=True,
):
    assert format in ['ct', 'seismic', 'json', 'bpseq', 'fasta'], 'Format not supported'
    
    if name is None:
        name = file_or_folder.split('/')[-1].split('.')[0]
    path = Path(name=name, root=path_out)
    path.make()
    
    if format == 'ct':
        datapoints = ListofDatapoints.from_ct(file_or_folder, tqdm=True, verbose=verbose)
    
    elif format == 'seismic':
        datapoints = ListofDatapoints.from_dreem_output(file_or_folder, predict_structure, tqdm=True, verbose=verbose)
        
    elif format == 'json':
        datapoints = ListofDatapoints.from_json(file_or_folder, predict_structure, tqdm=True, verbose=verbose)
        
    elif format == 'bpseq':
        datapoints = ListofDatapoints.from_bpseq(file_or_folder, tqdm=True, verbose=verbose)
        
    elif format == 'fasta':
        datapoints = ListofDatapoints.from_fasta(file_or_folder, predict_structure, tqdm=True, verbose=verbose)
        
    if filter:
        report = datapoints.filter(min_AUROC=min_AUROC)
    else:
        _, report = datapoints.drop_none_dp()
        report = f"Drop {report} datapoints with None values (null sequence or reference)"
    
    if verbose:
        print(report)
        
    with open(path.get_conversion_report(), 'w') as f:
        f.write("# Conversion report \n\n")
        f.write(report)
    
    if path_out is not None:
        datapoints.to_json(path.get_data_json())
        
    return datapoints.to_dict()
