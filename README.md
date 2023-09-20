# PyTools

This is the collection of Python scripts for pre- and post-processing related to [ICON-ART](https://www.icon-art.kit.edu/) model.

* [Plotscripts](./Plotscripts):
Scripts that can be executed or submitted as jobs on Mistral. It is recommended to use the latest Python module (python3/unstable). Scripts can also be used on other HPC systems if the required Python modules are installed. 

* [Plotscripts_with_Jupyter](./Plotscripts_with_Jupyter):
Scripts with Jupyter notebook.


Please submit an Issue if a problem occurs regarding one of the scripts.
If you want to make a modification to an existing script or upload your own script, please create a new branch and create a merge request with the main branch. 


# Working with the repository
To start working with this repository please `git clone` the repository
```bash
git clone  https://github.com/alihoshy/art_pytools.git
```

In case you have set up an ssh key, you can `git clone` via
```bash
git clone git@gitlab.dkrz.de:art/pytools.git
```

If you want to create your own branch, do the following
```bash
git checkout -b <branch name>
git push origin -u <branch name>
```

Do the following for committing work
```bash
git add <file name>
git commit -m 'TEXT about the changes done'
git push origin <branch name>


