Geometry v2 workflow verification: actual gprMax 4.0.1 outputs, 37 burial and 37 matching background traces, 500 MHz, 20 mm electromagnetic cells, 10 cm trace spacing. This is explicitly under-resolved and verifies the pipeline rather than numerical convergence or detectability. The shared detailed anatomy and portable soil field agree between Python and JavaScript. Three-dimensional geometry construction and selected transverse-plane export were also checked.


## General cemetery geometry (2.0)

`cemetery_solver_checks.json` records actual gprMax 4.0.1 single-trace checks for
a 27-degree rotated stepped stone cist with skeletal contents, in 2D and 3D.
Reproducible configs are `../examples/smoke_2d.json` and `smoke_3d.json`.
The 20/40 mm meshes and 12 ns windows are deliberately inadequate for skeletal
resolution and complete burial travel times. These check input compatibility
and finite receiver samples, not detectability or numerical convergence.

```bash
python burial_gpr.py generate examples/smoke_2d.json --out smoke2d --allow-underresolved
python burial_gpr.py run smoke2d --limit-traces 1 --threads 4
```

The old arrays and figures above retain their original single-coffin model;
they have not been relabelled as cemetery simulation results.
