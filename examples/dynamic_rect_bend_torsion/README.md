# Dynamic rectangular beam: bending + torsion

This case is a cantilever beam with a rectangular section subjected to a constant transverse tip force and a constant tip torque.

Model:
- Beam axis: global `X`
- Length: `1.0 m`
- Elements: `2`
- Section size: `0.04 m x 0.02 m`
- Area: `8.0e-4 m^2`
- Torsion constant `J`: `7.323481481e-8 m^4`
- Principal inertias:
  - `Iy = 2.666666667e-8 m^4`
  - `Iz = 1.066666667e-7 m^4`
- Material:
  - `E = 2.10e11 Pa`
  - `G = 8.076923077e10 Pa`
  - `rho = 7850 kg/m^3`

Loads at node 3:
- `Fy = -100 N`
- `Mx = 50 N*m`

Time integration:
- `dt = 5.0e-4 s`
- total time `0.1 s`
- constant load during the whole step

Expected comparison quantities:
- node 2: `Uy`, `Rx`
- node 3: `Uy`, `Rx`

For this repository:
- `StructureSolver` writes node 2 and node 3 displacements to `disp_ele_2_StructureSolver.dat`
- `BeamStructure` writes the second element nodal displacement vector to `disp_ele_2_BeamStructure.dat`, which corresponds to node 2 and node 3

Abaqus file:
- `abaqus_rect_bend_torsion.inp`

Suggested Abaqus comparison:
1. Run the input file with `B31` elements.
2. Extract history output for node 2 and node 3.
3. Compare:
   - `U2` with your solver `YTRA`
   - `UR1` with your solver `XROT`
