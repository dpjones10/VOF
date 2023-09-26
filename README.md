# VOF
Volume of Fluids Algorithm with a PLIC and a split advection scheme

This algorithm only tracks the interface of a two phase flow. It is not a Navier Stokes solver

Summary:
Square staggered grid structure
Forward Euler time integration
Mass conserving to machine precision (E_m < 10^(-12)) for incompressible flows
Geometric error converges roughly second order

Interface reconstruction options:
  1) Centered Columns Method
  2) Youngs' Method
  3) Mixed Youngs Centered
  4) ELVIRA

Boundary condition options:
  1) Periodic
  2) No penetration
  Note: No Pentration BC not implimented on ELVIRA
  Note: Boundary conditions controlled insdie interface normal computation functions

2D advection Options
  1) Eulerian Implcit - Eulerian Implicit
  2) Lagrangian Explicit - Lagrangian Explicit
  3) Lagrangian Explicit - Eulerian Implicit
  4) Eulerian Implicit - Lagrangian Explicit

