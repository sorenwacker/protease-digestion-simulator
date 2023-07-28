from setuptools import setup, find_packages

config = {
    'description': '__DESC__',
    'author': '__AUTHOR__',
    'url': 'https://github.com/__GITHUB_USER__',
    'download_url': 'https://github.com/__GITHUB_USER__/__PROJECT__',
    'author_email': '__EMAIL__',
    'version': '0.0.1',
    'install_requires': [],
    'packages': find_packages(),
    'scripts': [],
    'name': '__PROJECT__'
}

setup(**config)
