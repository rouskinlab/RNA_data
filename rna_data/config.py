import yaml

config = yaml.load(open('my_config.yaml', 'r'), Loader=yaml.FullLoader)
for k, v in config.items():
    exec(f"{k} = '{v}'")
