from setuptools import setup, find_packages

config = {
    'description': 'Predicts breakdown of protein sequences by proteases',
    'author': 'Soren Wacker',
    'url': 'https://github.com/sorenwacker',
    'download_url': 'https://github.com/sorenwacker/digest_simulator',
    'author_email': 'swacker@ucalgary.ca',
    'version': '0.0.1',
    'install_requires': [],
    'packages': find_packages(),
    'scripts': [],
    'name': 'digest_simulator'
}

setup(**config)
