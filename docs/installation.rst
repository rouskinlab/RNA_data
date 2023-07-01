
Installation
============

1. Download the code by cloning the repo:
-----------------------------------------

.. code-block:: bash

    git clone https://github.com/rouskinlab/RNA_data


2. Rename ``env_template`` to ``env`` and edit the token to your own huggingface token.
---------------------------------------------------------------------------------------

You can get a token at `huggingface.co/rouskinlab <https://huggingface.co/rouskinlab>`_.

- ``export HUGGINGFACE_TOKEN="your token here"``
- ``DATA_FOLDER``, which is the default path for storing the datafolders.
- ``RNASTRUCTURE_PATH``, you must change this to the path of your RNAstructure executable.
- ``RNASTRUCTURE_TEMP_FOLDER``, you can change this to the path of your choice.

3. Create a conda environment with the required packages:
----------------------------------------------------------------------------------

.. code-block:: bash

    conda env create -n rna_data python=3.8
    conda activate rna_data
    conda install -c bioconda rnastructure
    pip install -r requirements.txt


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



