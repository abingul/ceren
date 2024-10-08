# ceren

Simulation of Photons Generated from Cherenkov Radiation

This program simulates generation of Cherenkov photons when a particle
is passing through an optical material (Radiator) having rectangular prism shape. 
These photons are generated inside the radiator when particle's speed is greater than 
that of speed of light in radiator. Hence, Cherenkov photons are created while the following 
inequality is satisfied:
```math
     \beta n > 1
```
where $\beta = v_{particle}/c$ and $n = c/v_{light}$, and c = speed of light in vacuum.

Production points of photons and their directions are written to a file called <b>cherenkov.dat</b> 
which can be used in NSC mode of Ansys Zemax OpticStudio.
Note that in order to use these Cherenkov Photons, you need to fill <b>Inside Of</b> column in 
Zemax NSC Mode. Then, Zemax can do Ray Tracing and model absorption, reflection and scattering effects 
for each photon if you specify properties of optical material in Zemax.

In the program, we assume that
  - all distance units are in mm
  - photon wavelength units are in um (micrometer)
  - particle energy/momentum units are in GeV
  - Energy lost and multiple scattering effects are ignored

EOF




