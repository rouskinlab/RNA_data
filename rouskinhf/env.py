import os


class Env:
    def get_hf_token() -> str:
        if "HUGGINGFACE_TOKEN" in os.environ:
            return os.environ["HUGGINGFACE_TOKEN"]
        raise Exception("HUGGINGFACE_TOKEN not found in environment variables")

    def get_data_folder() -> str:
        if "DATA_FOLDER" in os.environ:
            return os.environ["DATA_FOLDER"]
        return "data/datafolders"

    def get_rnastructure_path() -> str:
        if "RNASTRUCTURE_PATH" in os.environ:
            return os.environ["RNASTRUCTURE_PATH"]
        return ""

    def get_rnastructure_temp_path() -> str:
        if "RNASTRUCTURE_TEMP_PATH" in os.environ:
            return os.environ["RNASTRUCTURE_TEMP_PATH"]
        return os.path.join(os.path.dirname(__file__), "..", "temp")
