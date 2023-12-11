from setuptools import setup

requirements = open('./requirements.txt').read().splitlines()

setup(
    name="rouskinhf",
    packages=["rouskinhf"],
    license="MIT",
    install_requires=requirements,
)