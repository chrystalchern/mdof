# `mdof`

<hr>
<br>

<div>
 
  <img align="left" src="https://raw.githubusercontent.com/chrystalchern/mdof/master/docs/_static/images/logos/mdof_readmefig.svg" width="250px" alt="mdof logo">

  Fast and friendly system identification for structures.

  [**Documentation**](https://chrystalchern.github.io/mdof/)

  <div style="align:center">

  [![Latest PyPI version](https://img.shields.io/pypi/v/mdof?logo=pypi&style=for-the-badge)](https://pypi.python.org/pypi/mdof)
  [![Downloads per Month](https://img.shields.io/pypi/dm/mdof?style=for-the-badge)]((https://pypi.python.org/pypi/mdof))

  </div>

</div>

<br>
<hr>

**The `mdof` package** solves **inverse problems**. It is tailored for the identification of system properties from structural vibrations.

Modal identification:
```python
periods, modeshapes = modes(input_motion, output_motion, dt)
```

<!-- Output-only modal identification:
```python
periods, modeshapes = modes(output_motion, dt)
``` -->

State space system identification:
```python
A,B,C,D = sysid(input_motion, output_motion)
```

Response reconstruction:
```python
output = reconstruct(realization, dt, input_motion)
```
<!-- 
Inverse eigenanalysis:
```python
eigvecs, eigvals = eigid(input_motion, output_motion)
``` -->

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
- [**Documentation**](https://chrystalchern.github.io/mdof/)
- Examples:
  - [**JupyterLab on DataHub**](https://datahub.berkeley.edu/hub/user-redirect/git-pull?repo=https%3A%2F%2Fgithub.com%2Fchrystalchern%2Fmdof&urlpath=lab%2Ftree%2Fmdof%2Fnotebooks%2FREADME.ipynb&branch=master) (UC Berkeley users)
  - [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/chrystalchern/mdof/HEAD?labpath=notebooks%2FREADME.ipynb) (non-UC Berkeley users)

-------------------------------------------------

## Support

<table align="center">
<tr>
  <td>
    <a href="https://github.com/stairlab">
    STAIRlab (STructural Artificial Intelligence Research Lab)
    </a>
  </td>

  <td>
    <a href="https://peer.berkeley.edu">
    PEER (Pacific Earthquake Engineering Research Center)
    </a>
  </td>

  <td>
    <a href="https://dot.ca.gov/">
    Caltrans (California Department of Transportation)
    </a>
  </td>

  <td>
    <a href="https://peer.berkeley.edu">
    BRACE2 (Bridge Rapid Assessment Center for Extreme Events)
    </a>
  </td>

  <td>
    <a href="https://www.nsfgrfp.org/">
    NSF (National Science Foundation) GRFP (Graduate Research Fellowship Program)
    </a>
  </td>

</tr>

<tr>
  <td align="center">
    <a href="https://github.com/stairlab">
    <img src="https://raw.githubusercontent.com/chrystalchern/mdof/master/docs/_static/images/logos/stairlab.svg"
         alt="PEER Logo" height="120px"/>
    </a>
  </td>

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
