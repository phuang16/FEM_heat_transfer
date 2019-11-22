# FEM_heat_transfer

**Problem statement:**
Consider a vascularized tissue that involves the blood flow. The heat transfer analysis here will include not only conduction but also
convection (which is due to the fluid motion [2]). The aim is to solve for the temperature distribution over time under this condition.

**Method:**
- Finite element method (FEM) modified from [1] is implemented, which solves the overall heat transfer problem by assuming a system combined of a thin-vessel (1-D, transient convectio-diffusion) and the surrounding tissue without blood flow (2-D conduction). Therefore,  both the spatial and temporal characteristics of the tissue temperature will be simulated.
- Geometry & FEM discretization:
![Figure1]( https://github.com/phuang16/FEM_heat_transfer/blob/master/FEM_heat_transfer_geometry.png)

**Results:**
- Spatiotemporal distribution of temperature:
![Figure2]( https://github.com/phuang16/FEM_heat_transfer/blob/master/FEM_heat_transfer_results.png)

**References:**
- [1] V. De Santis, M. Feliziani, F. Maradei, and C. Buccella, "Finite-element analysis of temperature increase in vascularized
biological tissues exposed to RF sources," IEEE Transactions on Magnetics, vol. 45, pp. 1682-1685, 2009.
- [2] W. L. Roland, N. Perumal, and N. Kankanhalli, "Fundamentals of the finite element method for heat and fluid flow," ed: John Wiley
and Sons Ltd. England, 2004.
- Part of the code was adapted from: https://www.mathworks.com/matlabcentral/fileexchange/44296-fem-diffusion-convection-solution
