# Contact Angle Analyser for MD simulations 

This Python library is designed to calculate contact angles for a given trajectory file, written as part of my MPhil thesis. The algorithm is broadly based off https://pubs.acs.org/doi/abs/10.1021/la990171l.

Overview
---

1. Aggregate selected simulation frames
2. Evaluates the number density for each cell
3. For each row of cells, the following $\tanh$ function is fitted to the number density profile as a function of radial distance from the centre of mass:
$$
\rho(r)=A\left(1-\tanh \left(\frac{2\left(r-r_{\mathrm{i}}\right)}{d}\right)\right),
$$
of which $r_{\mathrm{i}}, $A$, $d$ are free parameters. In particular, $r_\mathrm{i}$ is the Gibbs equimolar dividing surface, denoting the boundary point for each row of cells.
4. Fits a circle to the all $r_\mathrm{i}$ for each row, from which the contact angle is evaluated 


Installation
---
Easiest method is to clone the repository and install the package via pip. Ensure the correct Python environment has been selected before installing the package. Since the library is still under active development, we recommend creating a new virtual environment for this package.
```
git clone https://github.com/andrew-liang-2001/contact_angle.git
pip install /path/to/repository
```
Additional Features
---
A few helper scripts are included. For example, one script automatically obtains the equilibrium frame after supplying a CP2K energy file. Other scripts exist for handling IO functions.

Licence
---
This project is licensed under the MIT License - see the LICENSE file for details.

