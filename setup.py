from setuptools import find_packages
from setuptools import setup

setup(
    name = "ToppCellPy",
    version = "0.1.0",
    description = "This package is used for ToppCell gene module generation and downstream analysis.",
    author = "Kang Jin",
    author_email = "Kang.Jin@cchmc.org",
    url = "https://github.com/KANG-BIOINFO/ToppCellPy",
    packages = find_packages(),
    install_requires = [
        "numpy",
        "pandas",
        "scanpy>=1.6.0",
        "matplotlib",
        "seaborn",
        "gseapy",
        "scipy",
        "pylab"
    ]
)