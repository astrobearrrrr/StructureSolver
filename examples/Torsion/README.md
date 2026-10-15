# Dynamic torsion

This case is a cantilever beam aligned with global `Z`, subjected to a constant tip torque about the beam axis.

Geometry:
- Beam length: `5.0 m`
- Nodes: `51`
- Elements: `50`
- Element length: `0.1 m`
- `Lspan = 0.5 m`
- `Rspan = 0.5 m`
- Total span extension length: `1.0 m`
- `dir` direction: `(0.965925826289067979, -0.258819045102521017, 0.0)`

Current `BeamStructure` triad convention:
- Local `x` is the beam axis (global `+Z`).
- Local `y` is the projected `dir` direction.
- Local `z = x x y`.
- `Segment_InitTriad_D` is the unique initialization entry for the element triad.
- `Segment_UpdateTriad_D` is the unique update entry for the element triad.
- `Segment_RotateMatrix` reads the current orientation directly from `triad_ee`.

Current material and section data from `Torsion.dat`:
- `E = 288.86836`
- `G = 115.47344`
- `A = 0.03464`
- `rho = 28.86836`
- `gamma = 0.0`
- `IP = 1.3855187E-05`
- `IA = 3.463796E-06`
- `IB = 2.886666E-03`

Interpretation:
- `IP` is the Saint-Venant torsion constant used in `G * IP / L`.
- `IA` is aligned with local `y`.
- `IB` is aligned with local `z`.

Boundary and load:
- Node `1` is fully fixed.
- Node `51` carries a constant torque `Mz = 2.8E-05` (the `ZROT` load component).
- All other force and moment components are zero.

Time integration from `inFlow.dat`:
- Input file: `Torsion.dat`
- Newmark: `gamma = 0.5`, `beta = 0.25`
- `dt = 0.06 s`
- `5000` steps
- Total time: `300 s`
- Maximum Newton iterations per step: `12`
- Convergence tolerance: `dtol = 1.0E-08`
- `dampM = 0`
- `dampK = 0`
- `gamma1 = 1.0`
- `alphaf1 = 0.0`

Regression:
- This directory contains baseline results for `BeamStructure`.
- Field statistics: `old_fieldstat_BeamStructure.dat`
- Displacement history: `old_disp_ele_2_BeamStructure.dat`

Abaqus comparison:
- Input file: `abaqus_Torsion.inp`
- Element type: `B31`
- Section type: `*Beam General Section`
- Beam axis is global `+Z`.
- Section axis `1` follows the current `dir` vector `(0.965925826289067979, -0.258819045102521017, 0.0)`.
- Property mapping is:
  - `I11 = IA`
  - `I22 = IB`
  - `J = IP`
- Node `1` is fully fixed; node `51` carries `Cload` component `6 = 2.8E-05`.
- Dynamic step uses:
  - `nlgeom=YES`
  - `alpha = 0.0`
  - `dt = 0.06 s`
  - Total time: `300 s`
- History output requests `U` and `UR` at node `26` (`NMID`) and node `51` (`NTIP`), every `50` increments (`3 s`).
- Compare Abaqus `UR3` with the solver's global `ZROT` response at the same node and time.
- Exported Abaqus results: `disp_ele_2_abaqus.dat`
- Origin comparison project: `Torsion-BeamStructure-Abaqus.opju`
