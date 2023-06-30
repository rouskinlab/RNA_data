import yaml

config = yaml.load(open('config.yaml', 'r'), Loader=yaml.FullLoader)
for k, v in config.items():
    exec(f"{k} = '{v}'")
