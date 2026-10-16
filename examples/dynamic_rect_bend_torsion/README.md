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

Boundary and loads:
- Node `1` is fully fixed.
- Node `3` carries `Fy = -100 N` and `Mx = 50 N*m`.

Time integration:
- `dt = 5.0e-4 s`
- total time `0.1 s`
- constant load during the whole step

Expected comparison quantities:
- node 2: `Uy`, `Rx`
- node 3: `Uy`, `Rx`

Displacement history output:
- The current root-level solvers have their displacement history output statements commented out; they write `fieldstat` files for regression testing.
- To export histories from `StructureSolver`, enable the output block for `disp_ele_2_StructureSolver.dat`, which selects node `2` and node `3`.
- To export histories from `BeamStructure`, enable and adapt the output block for `disp_ele_2_BeamStructure.dat` to select node `2` and node `3`. The current commented example selects the tip of element `50`, which does not exist in this two-element case.
- Stored baseline histories are `old_disp_ele_2_StructureSolver.dat` and `old_disp_ele_2_BeamStructure.dat`.

Abaqus file:
- `abaqus_rect_bend_torsion.inp`

Suggested Abaqus comparison:
1. Run the input file with `B31` elements.
2. Extract history output for node 2 and node 3.
3. Compare:
   - `U2` with your solver `YTRA`
   - `UR1` with your solver `XROT`
