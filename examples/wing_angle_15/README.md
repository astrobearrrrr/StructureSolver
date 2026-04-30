# Wing Angle 15

This case is a cantilever wing-like beam aligned with global `Z`.

Geometry:
- Beam length: `5.0 m`
- Nodes: `51`
- Elements: `50`
- `Lspan = 0.5 m`
- `Rspan = 0.5 m`
- Total span extension length: `1.0 m`
- `dir` direction: `(0.0, 1.0, 0.0)`

Current `BeamStructure` triad convention:
- Local `x` is the beam axis.
- Local `y` is the projected `dir` direction.
- Local `z = x x y`.
- `Segment_InitTriad_D` is the unique initialization entry for the element triad.
- `Segment_UpdateTriad_D` is the unique update entry for the element triad.
- `Segment_RotateMatrix` reads the current orientation directly from `triad_ee`.

Current material and section data from `Wing.dat`:
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
- Node `51` carries:
  - `Fx = -0.008`
  - `Fy = 0.008`
  - `Fz = 0.008`
  - `Mx = 0.0`
  - `My = 0.0`
  - `Mz = 0.008`

Time integration from `inFlow.dat`:
- Newmark: `gamma = 0.5`, `beta = 0.25`
- `dt = 1.0e-4 s`
- `150` steps
- `12` Newton iterations per step
- `dampM = 0`
- `dampK = 0`
- `gamma1 = 1.0`
- `alphaf1 = 0.0`

Regression:
- This example is intended to be checked with `BeamStructure` only.
- Baseline file: `old_fieldstat_BeamStructure.dat`

Abaqus comparison:
- Input file: `abaqus_wing_angle_15.inp`
- Element type: `B31`
- Section type: `*Beam General Section`
- Beam axis is global `+Z`
- Section axis `1` is set to global `+Y`
- Property mapping is:
  - `I11 = IA`
  - `I22 = IB`
  - `J = IP`
- Dynamic step uses:
  - `nlgeom=YES`
  - `alpha = 0.0`
  - `dt = 1.0e-4 s`
  - total time `1.5e-2 s`
- The Abaqus file is kept consistent with the current `Wing.dat` geometry, material data, section data, and tip loading.
