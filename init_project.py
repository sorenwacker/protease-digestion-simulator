#!/usr/bin/env python

import os

if __name__ == '__main__':
    project_name = input('Project name: ').replace('-', '_')
    desc = input('Project description: ')
    author_name = input('Author name: ')
    email = input('Contact email address: ')
    git_user = input('Github username: ')

    assert project_name == project_name.replace(' ', ''), 'Project name cannot have spaces.'

    with open('setup.py', 'r') as file:
        setup_py = ''.join(file.readlines())
    
    setup_py = setup_py.replace('__PROJECT__', project_name)
    setup_py = setup_py.replace('__AUTHOR__', author_name)
    setup_py = setup_py.replace('__EMAIL__', email)
    setup_py = setup_py.replace('__DESC__', desc)
    setup_py = setup_py.replace('__GITHUB_USER__', git_user)

    with open('setup.py', 'w') as file:
        file.write(setup_py)

    with open('setup.cfg', 'r') as file:
        setup_cfg = ''.join(file.readlines())

    setup_cfg = setup_cfg.replace('__PROJECT__', project_name)

    with open('setup.cfg', 'w') as file:
        file.write(setup_cfg)
    
    os.system('rm -r template; rm -r tests/test__hello_world.py;')
    os.system(f'mkdir {project_name}; touch {project_name}/__init__.py')
    os.system('pip install versioneer; versioneer install')
    os.system('git add .gitattributes -f')
    os.system(f"echo '# {project_name}\n{desc}' > README.md")    

