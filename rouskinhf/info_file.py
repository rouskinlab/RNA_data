import datetime, os, json
from typing import Any
from .path import PathDatafolder

STRUCTURE_FROM_RNASTRUCTURE = 'RNAstructure'
STRUCTURE_FROM_SOURCE = 'source'
DMS_FROM_RNASTRUCTURE = 'RNAstructure'
DMS_FROM_SOURCE = 'source'

class InfoFileWriterTemplate(PathDatafolder):

    def __init__(self, name, root, has_dms, has_structure) -> None:
        super().__init__(name, root)
        self.has_dms = has_dms
        self.has_structure = has_structure

        self.info = {
            'name': self.name,
            'upload_date': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'user': os.environ['USER'],
        }

    def write(self):
        """Write the info file and the readme file."""
        with open(self.get_info_file(), 'w') as f:
            json.dump(self.info, f, indent=4)

        with open(self.get_readme_file(), 'w') as f:
            f.write(f"""
---
license: mit
language:
  - en
tags:
  - chemistry
  - biology`
author: Silvi Rouskin
pretty_name: {self.name}
---

""")
            f.write(f"# {self.name}\n")
            for k, v in self.info.items():
                if k not in ['structure', 'about structure', 'DMS', 'about DMS','name', 'filtering_report']:
                    f.write(f"\t{k}: {v}\n")
            f.write(f"\n\tData types:")
            if self.has_structure:
                f.write(f"\n\t- structure ({self.info['about structure']})")
            if self.has_dms:
                f.write(f"\n\t- DMS ({self.info['about DMS']})")
            if 'filtering_report' in self.info:
                f.write(f"\n\n\tFiltering report: \n\t{self.info['filtering_report']}")
        return self
    
    def add_filtering_report(self, filtering_report):
        self.info['filtering_report'] = filtering_report
        return self


class InfoFileWriterFromDREEMoutput(InfoFileWriterTemplate):

    def __init__(self, name, root, path_in, predict_dms, predict_structure) -> None:
        super().__init__(name, root, has_dms=True, has_structure=predict_structure)

        file = json.load(open(path_in, 'r'))
        if 'user' in file:
            file['experimenter'] = file.pop('user')
        for k, v in file.items():
            if type(v) != dict:
                self.info[k] = v

        self.info['source'] = "`{}` (DREEM output format)".format(os.path.basename(path_in))
        if predict_structure:
            self.info['structure'] = STRUCTURE_FROM_RNASTRUCTURE
            self.info['about structure'] = 'Predicted using RNAstructure and the DMS reactivity data.'
        self.info['DMS'] = DMS_FROM_SOURCE
        self.info['about DMS'] = 'Measured experimentally using the attached experimental conditions.'


class InfoFileWriterFromCT(InfoFileWriterTemplate):

    def __init__(self, name, root, path_in, predict_dms, predict_structure) -> None:
        super().__init__(name, root, has_dms=predict_dms, has_structure=True)

        self.info['source'] = "`{}` (CT files format)".format(os.path.basename(path_in))
        self.info['structure'] = STRUCTURE_FROM_SOURCE
        self.info['about structure'] = 'Read from source.'
        if predict_dms:
            self.info['DMS'] = DMS_FROM_RNASTRUCTURE
            self.info['about DMS'] = 'Predicted using RNAstructure at 310K.'


class InfoFileWriterFromFasta(InfoFileWriterTemplate):

    def __init__(self, name, root, path_in, predict_dms, predict_structure) -> None:
        super().__init__(name, root, has_dms=predict_dms, has_structure=predict_structure)

        self.info['source'] = "`{}` (fasta file format)".format(os.path.basename(path_in))
        if predict_structure:
            self.info['structure'] = STRUCTURE_FROM_RNASTRUCTURE
            self.info['about structure'] = 'Predicted using RNAstructure at 310K.'
        if predict_dms:
            self.info['DMS'] = DMS_FROM_RNASTRUCTURE
            self.info['about DMS'] = 'Predicted using RNAstructure at 310K.'


def infoFileWriter(source, datafolder, predict_dms, predict_structure):

    if source == 'dreem_output':
        return InfoFileWriterFromDREEMoutput(datafolder.name, datafolder.path_out, datafolder.path_in, predict_dms=predict_dms, predict_structure=predict_structure)
    if source == 'ct':
        return InfoFileWriterFromCT(datafolder.name, datafolder.path_out, datafolder.path_in, predict_dms=predict_dms, predict_structure=predict_structure)
    if source == 'fasta':
        return InfoFileWriterFromFasta(datafolder.name, datafolder.path_out, datafolder.path_in, predict_dms=predict_dms, predict_structure=predict_structure)
    raise ValueError(f"Unknown source: {source}")

