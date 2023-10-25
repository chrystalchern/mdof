# `mdof`

<hr>
<br>

<img align="left" src="https://raw.githubusercontent.com/BRACE2/mdof/master/docs/_static/images/logos/mdof_readmefig.svg" width="250px" alt="mdof logo">

Fast and friendly system identification for structures. *[Learn More](https://brace2.github.io/mdof/)*

<div style="align:center">

[![Latest PyPI version](https://img.shields.io/pypi/v/mdof?logo=pypi&style=for-the-badge)](https://pypi.python.org/pypi/mdof)
[![Downloads per Month](https://img.shields.io/pypi/dm/mdof?style=for-the-badge)]((https://pypi.python.org/pypi/mdof))

</div>

<hr>
<br>

**Conventional dynamic eigenanalysis** is well-known to structural engineers -- from system properties and a given excitation, the goal is to determine the system's dynamic response.

```python
output_motion = eigen(M,C,K, input_motion)
```

**The `mdof` system id package** allows structural engineers to solve **inverse eigenanalysis** and related problems -- from structural vibrations, the goal is to identify the system properties.

Inverse eigenanalysis:
```python
eigvecs, eigvals = eigid(input_motion, output_motion)
```

State space system identification:
```python
A,B,C,D = sysid(input_motion, output_motion)
```


## Getting Started

Click [**JupyterLab on DataHub**](https://datahub.berkeley.edu/hub/user-redirect/git-pull?repo=https%3A%2F%2Fgithub.com%2FBRACE2%2Fmdof&urlpath=lab%2Ftree%2Fmdof%2Fnotebooks%2FREADME.ipynb&branch=master) (UC Berkeley users) or  [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/BRACE2/mdof/HEAD?labpath=notebooks%2FREADME.ipynb) (non-UC Berkeley users) to access and experiment with example Jupyter notebooks.

-------------------------------------------------

