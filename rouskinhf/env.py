import os


class Env:
    """Class to get environment variables."""
    def get_hf_token() -> str:
        if "HUGGINGFACE_TOKEN" in os.environ:
            return os.environ["HUGGINGFACE_TOKEN"]
        raise Exception("HUGGINGFACE_TOKEN not found in environment variables")

    def get_rnastructure_path() -> str:
        if "RNASTRUCTURE_PATH" in os.environ:
            return os.environ["RNASTRUCTURE_PATH"]
        return ""

    def get_rnastructure_temp_path() -> str:
        if "RNASTRUCTURE_TEMP_PATH" in os.environ:
            return os.environ["RNASTRUCTURE_TEMP_PATH"]
        return os.path.join(os.path.dirname(__file__), "..", "temp")
