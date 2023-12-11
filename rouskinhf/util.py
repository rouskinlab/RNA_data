import os

# Define the one-hot encodings for the sequences and structures
seq2int = {
    "X": 0,
    "A": 1,
    "C": 2,
    "G": 3,
    "U": 4,
}

int2seq = {v: k for k, v in seq2int.items()}

dot2int = {".": 1, "(": 2, ")": 3, "X": 0}
int2dot = {v: k for k, v in dot2int.items()}

UKN = -1000


class DreemUtils:
    def flatten_json(data):
        out, row = [], {}

        for k, v in data.copy().items():
            if type(v) != dict:
                row[k] = v
            else:
                row["reference"] = k
                for k2, v2 in v.items():
                    if type(v2) != dict:
                        row[k2] = v2
                    else:
                        row["section"] = k2
                        for k3, v3 in v2.items():
                            row["cluster"] = k3
                            if type(v3) != dict:
                                row[k3] = v3
                            else:
                                for k4, v4 in v3.items():
                                    row[k4] = v4
                                out.append(row.copy())
        return out

    def sort_dict(mut_profiles):
        sorting_key = lambda item: 1000 * {int: 0, str: 0, float: 0, list: 1, dict: 2}[
            type(item[1])
        ] + ord(item[0][0])
        mut_profiles = {k: mut_profiles[k] for k in sorted(mut_profiles)}
        mut_profiles = dict(sorted(mut_profiles.items(), key=sorting_key))
        for k, v in mut_profiles.items():
            if type(v) is not dict:
                continue
            mut_profiles[k] = DreemUtils.sort_dict(v)
        return mut_profiles


def standardize_sequence(sequence):
    sequence = sequence.upper()
    sequence = sequence.replace("T", "U")
    # trim the sequence
    sequence = sequence.replace(" ", "").replace("\n", "").replace("\t", "")
    return sequence


def sequence_has_regular_characters(sequence):
    return not (set(sequence) - set("ACGU"))
