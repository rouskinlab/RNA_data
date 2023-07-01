import datetime, os, json
from typing import Any
from .path import PathDatafolder

STRUCTURE_FROM_RNASTRUCTURE = 'RNAstructure'
STRUCTURE_FROM_SOURCE = 'source'
DMS_FROM_RNASTRUCTURE = 'RNAstructure'
DMS_FROM_SOURCE = 'source'

class InfoFileWriterTemplate(PathDatafolder):

    def __init__(self, name, root) -> None:
        super().__init__(name, root)

        self.info = {
            'name': self.name,
            'date': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'user': os.environ['USER'],
        }

    def write(self):
        """Write the info file."""
        with open(self.get_info_file(), 'w') as f:
            json.dump(self.info, f, indent=4)


class InfoFileWriterFromDREEMoutput(InfoFileWriterTemplate):

    def __init__(self, name, root, path_in) -> None:
        super().__init__(name, root)

        file = json.load(open(path_in, 'r'))
        if 'user' in file:
            file['experimenter'] = file.pop('user')
        for k, v in file.items():
            if type(v) != dict:
                self.info[k] = v

        self.info['structure'] = STRUCTURE_FROM_RNASTRUCTURE
        self.info['about structure'] = 'Predicted using RNAstructure and the DMS reactivity data.'
        self.info['DMS'] = DMS_FROM_SOURCE
        self.info['about DMS'] = 'Measured experimentally using the attached experimental conditions.'


class InfoFileWriterFromCT(InfoFileWriterTemplate):

    def __init__(self, name, root) -> None:
        super().__init__(name, root)

        self.info['structure'] = STRUCTURE_FROM_SOURCE
        self.info['about structure'] = 'Read from source.'
        self.info['DMS'] = DMS_FROM_RNASTRUCTURE
        self.info['about DMS'] = 'Predicted using RNAstructure at 310K.'


class InfoFileWriterFromFasta(InfoFileWriterTemplate):

    def __init__(self, name, root) -> None:
        super().__init__(name, root)

        self.info['structure'] = STRUCTURE_FROM_RNASTRUCTURE
        self.info['about structure'] = 'Predicted using RNAstructure at 310K.'
        self.info['DMS'] = DMS_FROM_RNASTRUCTURE
        self.info['about DMS'] = 'Predicted using RNAstructure at 310K.'


def infoFileWriter(source, datafolder):
    if source == 'dreem_output':
        return InfoFileWriterFromDREEMoutput(datafolder.name, datafolder.path_out, datafolder.path_in)
    if source == 'ct':
        return InfoFileWriterFromCT(datafolder.name, datafolder.path_out)
    if source == 'fasta':
        return InfoFileWriterFromFasta(datafolder.name, datafolder.path_out)
    raise ValueError(f"Unknown source: {source}")

