# `mdof`

<hr>
<br>

<img align="left" src="https://raw.githubusercontent.com/chrystalchern/mdof/master/docs/_static/images/logos/mdof_readmefig.svg" width="250px" alt="mdof logo">

Fast and friendly system identification for structures. *[Learn More](https://chrystalchern.github.io/mdof/)*

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

**Create beautiful visuals for historical modal property documentation** with integration of the [`osmg`](https://github.com/ioannis-vm/OpenSees_Model_Generator) and [`opensees`](https://pypi.org/project/opensees/) packages.

<table align="center">
<tr>
  <td>
  <img src="https://raw.githubusercontent.com/chrystalchern/mdof/master/docs/_static/images/gallery/LA_modes_core.png" width="225px" alt="historical mode shape documentation">
  </td>
  <td>
  <img src="https://raw.githubusercontent.com/chrystalchern/mdof/master/docs/_static/images/gallery/LA_FDD_02.png" width="650px" alt="historical spectral density documentation">
  </td>
</tr>
</table>




## Getting Started

Click [**JupyterLab on DataHub**](https://datahub.berkeley.edu/hub/user-redirect/git-pull?repo=https%3A%2F%2Fgithub.com%2FBRACE2%2Fmdof&urlpath=lab%2Ftree%2Fmdof%2Fnotebooks%2FREADME.ipynb&branch=master) (UC Berkeley users) or  [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/BRACE2/mdof/HEAD?labpath=notebooks%2FREADME.ipynb) (non-UC Berkeley users) to access and experiment with example Jupyter notebooks.

-------------------------------------------------

## Support

<table align="center">
<tr>
  <td>
    <a href="https://peer.berkeley.edu">
    Pacific Earthquake Engineering Research Center (PEER)
    </a>
  </td>

  <td>
    <a href="https://dot.ca.gov/">
    California Department of Transportation (Caltrans)
    </a>
  </td>

  <td>
    <a href="https://peer.berkeley.edu">
    Bridge Rapid Assessment Center for Extreme Events (BRACE2)
    </a>
  </td>

  <td>
    <a href="https://www.nsfgrfp.org/">
    National Science Foundation (NSF) Graduate Research Fellowship Program (GRFP)
    </a>
  </td>

</tr>

<tr>

  <td align="center">
    <a href="https://peer.berkeley.edu">
    <img src="https://raw.githubusercontent.com/chrystalchern/mdof/master/docs/assets/PEER_logo_old.svg"
         alt="PEER Logo" height="120px"/>
    </a>
  </td>

  <td align="center">
    <a href="https://dot.ca.gov/">
    <img src="https://raw.githubusercontent.com/claudioperez/sdof/master/docs/assets/Caltrans.svg.png"
         alt="Caltrans Logo" height="120px"/>
    </a>
  </td>

  <td align="center">
    <a href="https://peer.berkeley.edu">
    <img src="https://raw.githubusercontent.com/claudioperez/sdof/master/docs/assets/brace2_logo-new3_ungrouped.svg"
         alt="BRACE2 Logo" height="120px"/>
    </a>
  </td>

  <td align="center">
    <a href="https://www.nsfgrfp.org/">
    <img src="https://raw.githubusercontent.com/chrystalchern/mdof/master/docs/_static/images/logos/nsf_logo.jpg"
         alt="NSF Logo" height="120px"/>
    </a>
  </td>
 
</tr>
</table>
