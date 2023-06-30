import yaml, os

my_config_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'my_config.yaml')

config = yaml.load(open(my_config_path, 'r'), Loader=yaml.FullLoader)
for k, v in config.items():
    exec(f"{k} = '{v}'")
