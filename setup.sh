# install python 3.10
sudo apt-get update
sudo apt-get install python3.10 -y
set alias python /opt/conda/bin/python3.10
python --version

# download and install dmsensei
git clone https://github.com/rouskinlab/DMSensei 
python -m pip install -r DMSensei/requirements.txt
export PYTHONPATH=DMSensei

# install rouskinhf
python -m pip install rouskinhf
export HUGGINGFACE_TOKEN="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" # get your own token!

# setup wandb
python -m pip install wandb
wandb login e880fa9907144e0820e2c725b19ad1c5238793ee # get your own token!