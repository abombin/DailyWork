
import importlib
import subprocess

def check_and_install_packages(packages):
    for package in packages:
        try:
            importlib.import_module(package)
            print(f"{package} is already installed.")
        except ImportError:
            print(f"{package} is not installed. Installing...")
            subprocess.check_call(['pip', 'install', package])
            print(f"{package} has been successfully installed.")

# Usage example
packages_to_check = ['pandas', 'argparse', 'numpy', 'itertools']

check_and_install_packages(packages_to_check)
