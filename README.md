# `mdof`


<img align="left" src="https://raw.githubusercontent.com/BRACE2/OpenSeesRT/master/docs/figures/spectrum.svg" width="250px" alt="PEER Logo">

Fast and friendly system identification for structures.

<br>

**Conventional dynamic eigenanalysis** is well-known to structural engineers -- from system properties and a given excitation, the goal is to determine a system's dynamic responsess.

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


## Try it out!
Click [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/BRACE2/mdof/HEAD?labpath=notebooks%2FREADME.ipynb) to access and experiment with example Jupyter notebooks.

<div style="align:center">

[![Latest PyPI version](https://img.shields.io/pypi/v/mdof?logo=pypi&style=for-the-badge)](https://pypi.python.org/pypi/mdof)


</div>

-------------------------------------------------

