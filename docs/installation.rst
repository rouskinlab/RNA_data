
Installation
============

1. Download the code by cloning the repo:
-----------------------------------------

.. code-block:: bash

    git clone https://github.com/rouskinlab/RNA_data

2. Create a conda environment with the required packages:
----------------------------------------------------------------------------------

.. code-block:: bash

    conda create -n rna_data python=3.10
    conda activate rna_data
    conda install -c bioconda rnastructure # you can also install this manually
    pip install -r requirements.txt


3. Create a ``env`` file
-------------------------

Rename ``env_template`` to ``env`` and edit the token to your own huggingface token.

.. code-block:: bash

    mv env_template env

You can get a token at `huggingface.co/settings/tokens <https://huggingface.co/settings/tokens>`_.
Check out the `HuggingFace tutorial <https://github.com/rouskinlab/RNA_data/blob/main/tutorials/huggingface.ipynb>`_.


- ``export HUGGINGFACE_TOKEN="your token here"``
- ``DATA_FOLDER``, which is the default path for storing the datafolders.
- ``RNASTRUCTURE_PATH``, you must change this to the path of your RNAstructure executable.
- ``RNASTRUCTURE_TEMP_FOLDER``, you can change this to the path of your choice.


4. Source the environment:
-----------------------------------------

.. code-block:: bash

    source env

5. Test the installation:
-----------------------------------------

.. code-block:: bash

    python -m pytest

6. Use the tutorial notebooks to get started.
----------------------------------------------------------------------------------

.. code-block:: bash

    jupyter notebook



