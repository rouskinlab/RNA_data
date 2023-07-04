Usage
=====

Create, push and pull data with datafolders
-------------------------------------------

The datafolder object is a wrapper around your local folder and HuggingFace API, to keep a consistent datastructure across your datasets. It contains multiple methods to create datasets from various input formats, store the data and metadata in a systematic way, and push / pull from HuggingFace.

On HuggingFace, the datafolder stores the data under the following structure:


.. code-block:: bash

    HUGGINGFACE DATAFOLDER
    - [datafolder name]
        - source
            - whichever file(s) you used to create the dataset (fasta, set of CTs, etc.).
        - data.json # the data under a human readable format.
        - info.json # the metadata of the dataset. This file indicates how we got the DMS signal and the structures (directly from the source or from a prediction).
        - README.md # the metadata of the dataset in a human readable format.

Locally, we have the same structure with the addition of `.npy` files which contain the data in a machine readable format. Each `.npy` file contains a numpy array of the data, and the name of the file is the name of the corresponding key in the `data.json` file. The source file won't be downloaded by default. Hence, the local structure is:

.. code-block:: bash

    LOCAL DATAFOLDER
    - [datafolder name]
        ...
        - README.md # the metadata of the dataset in a human readable format
        - references.npy
        - sequences.npy
        - base_pairs.npy
        - dms.npy


Note: some of the steps below will require getting through the HuggingFace tutorial first.

Tutorials
---------

Log in to HuggingFace
~~~~~~~~~~~~~~~~~~~~~

See the `tutorial <https://github.com/rouskinlab/RNA_data/blob/main/tutorials/huggingface.ipynb>`_.

Create a datafolder, push and pull data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See the `tutorial <https://github.com/rouskinlab/RNA_data/blob/main/tutorials/create_push_pull.ipynb>`_.

Use the data for your models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See the `tutorial <https://github.com/rouskinlab/RNA_data/blob/main/tutorials/use_for_models.ipynb>`_.
