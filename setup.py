from setuptools import setup

setup(
    name="rouskinhf",
    packages=["rouskinhf"],
    license="MIT",
    install_requires=[
        'pytest',
        'pytest_cov',
        'numpy',
        'pyarrow',
        'pandas',
        'huggingface',
        'huggingface_hub',
        'jupyter',
        'ipykernel',
        'ipython',
        'scipy',
        'pydantic',
        'pytest-dotenv',
        'python-dotenv',
        'tqdm',
        'pytest-env'
        ],
)