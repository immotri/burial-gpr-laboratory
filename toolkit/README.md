# Burial GPR Laboratory v1.1

A reproducible scenario builder and gprMax 4.0.1 export workflow for a reference-inspired Viking Age coffin burial. The browser previews geometry; gprMax solves Maxwell's equations locally. This is a working research starting point, not a calibrated reconstruction or a validated prediction of detectability.

## Start with a small experiment

Use Python 3.11–3.13 (3.12 recommended).

```bash
python -m venv .venv
```

Activate the environment:

- Linux/macOS: `source .venv/bin/activate`
- Windows PowerShell: `.venv\Scripts\Activate.ps1`

Then:

```bash
python -m pip install -r requirements.txt
python burial_gpr.py estimate quickstart.json
python burial_gpr.py generate quickstart.json --out first_coffin
python burial_gpr.py run first_coffin --limit-traces 3 --threads 4
python burial_gpr.py collect first_coffin
```

The three-trace limit checks the workflow. It does not cover the burial: the full profile begins at the left domain margin. Omit `--limit-traces 3` for the full survey. The quickstart has a coffin and backfill, homogeneous layers, a 250 MHz source and materials frozen at the centre frequency. It excludes the smaller objects to keep the initial mesh manageable. It is still a full-wave 2D electromagnetic calculation, with the constitutive approximation stated explicitly.

For your own model, use the browser's **Export model**, save it as `my_model.json` beside the script, and substitute that name. The toolkit includes `model.json` with the complete layout and a 100/250/500 MHz selection. All values are SI except frequencies in MHz, time windows in ns, and ages in years.

`full_frequency_sweep.json` selects 100, 250, 500, 1000 and 2000 MHz with a 0.75 mm cell. This is a costly experiment even in 2D; estimate it before solving. Do not begin with a full 3D raster.

## What is represented

The photographs show a rectangular wooden coffin, supine skeleton, sword-like object and vessels near the feet. Dimensions are assumed, editable values; they were not measured from the photographs.

The default original coffin is 2.15 × 0.72 × 0.45 m, with 35 mm boards and its top 0.75 m below ground. The model includes topsoil, subsoil, an enlarged grave shaft and backfill, residual wood, collapsed lid panels, simplified bone elements, hollow sediment-filled vessels and a corroding metal object. There is no soft tissue at any age. Setting age zero gives an initial *skeleton-bearing* geometry, not a freshly deceased person.

Coordinates in the scenario are x/y horizontal and depth positive downward. gprMax inputs convert depth to elevation z, with the origin at the domain bottom. Coffin length is parallel to x; survey profiles can follow x or y. Arbitrary coffin azimuth and imported anatomical meshes are not implemented in this version.

Geometry v2 uses a shared procedural anatomy template: a hollow cranial vault and facial cavities, mandible and teeth, 24 vertebral bodies, 12 pairs of curved ribs, clavicles and scapulae, pelvic shells, joint ends, hollow long-bone shafts, hands and feet. Coffin boards are separate, the grave shaft narrows in stepped lifts, and vessels have hollow profiles with openings. The cortical long-bone wall is approximately 4 mm; the cranial shell is 8 mm. These are geometric approximations rather than an osteometric reconstruction or scan of a particular individual. The photographs contain more detailed objects than these primitives represent. The 3D view is transparent for inspection, whereas the solver receives the complete ordered geometry.

## Preservation is a hypothesis

For elapsed time t and a user-defined half-life H, the remaining index is:

`r(t) = 2^(-t/H)`

This convenient sensitivity law is not an archaeological decay calibration. Temperature, oxygen/redox state, groundwater history, pH, species, microbial conditions and soil chemistry are not inferred from age.

| Browser preset | Wood H | Bone H | Iron H | Void H | Subsoil water fraction | Subsoil σDC |
|---|---:|---:|---:|---:|---:|---:|
| Drained | 200 yr | 1800 yr | 1600 yr | 100 yr | 0.10 | 0.003 S/m |
| Damp | 600 yr | 1800 yr | 2000 yr | 180 yr | 0.18 | 0.006 S/m |
| Waterlogged | 5000 yr | 1000 yr | 4000 yr | 50 yr | 0.40 | 0.025 S/m |

These numbers deliberately produce contrasting scenarios; they are illustrative priors, not literature-derived half-lives. The presets also change target pore saturation, topsoil moisture and backfill properties. Once a preset is applied, parameters can be edited independently. Editing `environment` in JSON changes the environment label and the no-air-void rule; it does **not** automatically rewrite the explicitly stored numbers. Use the browser preset to update all its values together.

- Wood: the residual index controls effective mixing between the porous wood endmember and replacement infill. It is not a direct estimate of surviving wood mass. Nominal shell geometry subsides as wood is lost.
- Collapse: the maximum collapse index is 0.85 by default. A rigid-contents floor limits remaining coffin height to 0.28 m with bones and 0.24 m with vessels, preventing the simplified lid from cutting through the contents. This is not a geomechanical calculation.
- Bone: each whole simplified element is retained if its seeded random number is below the survival probability. Bone loss is element-wise, without gradual dissolution or fragmentation. Changing the seed changes the retained anatomy.
- Iron: a PEC core shrinks in width and thickness in proportion to sqrt(r), while a lossy corrosion halo remains. This approximates high-conductivity scattering, without skin-depth resolution, magnetic dispersion or chemical corrosion kinetics.
- Void: the surviving cavity height follows its own half-life. No coffin air void is assumed in the waterlogged preset. Vessels and the skull are sediment-filled. An independently surviving chest cavity is not included.

For 1200–1400-year burial scenarios, compare several preservation hypotheses and several seeds. Do not treat one generated geometry as the expected preservation state.

## Electromagnetic constitutive model

Every non-metal material has relative magnetic permeability 1 and no magnetic loss. The tool uses a complex refractive index mixing approximation (CRIM):

`sqrt(εmix*) = (1 - φ) sqrt(εsolid) + θ sqrt(εwater*) + (φ - θ)`

Here φ is porosity and θ is volumetric water content per bulk volume; 0 ≤ θ ≤ φ. For target materials, θ = φ × pore saturation. The pure-water dielectric endmember at approximately 20°C is:

`εwater*(ω) = 4.9 + 75.2 / (1 + iω × 9.231 ps)`

DC conductivity is an independent bulk parameter, added once:

`εtotal* = εmix* - i σDC / (ω ε0)`

This avoids counting DC conduction again in the Debye dielectric loss. The material table's effective conductivity includes the dielectric-loss contribution at the selected frequency. DC conductivity is not calculated from water fraction or salinity.

Initial isotropic solid permittivities are mineral 4.5, wood 3, bone 7, pottery 5.5 and corrosion material 8. These and target DC conductivities are editable in the material calibration section or JSON. They are plausible starting priors, not measured properties of archaeological samples. In particular, fresh medical bone data should not be substituted uncritically for diagenetic, porous archaeological bone.

In **dispersive mode**, a constrained non-negative fit represents each complex CRIM curve using ε∞ and up to three positive Debye poles. Fixed relaxation times are 9.231, 50 and 200 ps; fits use 50 MHz–5 GHz. The manifest records each fit's maximum relative complex error. The default material set fits within 0.4%; this is a constitutive fit check, not validation of the physical endmembers. The timestep is capped conservatively below the shortest relaxation time, in addition to the spatial CFL condition.

In **frozen mode**, ε′ and total effective conductivity are evaluated at the selected centre frequency and held constant through that run. This is a faster, narrowband approximation. It does not reproduce the full frequency-dependent loss across a broad pulse. Compare against dispersive runs before using it quantitatively across 100 MHz–2 GHz.

CRIM omits a calibrated clay/bound-water response, anisotropy and explicit salinity/temperature dependence. For clay-rich or saline soils, replace or extend the constitutive model with measured spectra or a suitable calibrated soil model. A field conductivity measured at low frequency is not automatically the correct total GHz loss.

## Heterogeneity and controls

In geometry v2, a seeded 12-mode band-limited sinusoidal field is normalised over a fixed physical lattice and quantised into 13 moisture bins. The browser and Python use the same coefficients and normalisation. The legacy geometry v1 implementation retains its older smoothed Gaussian field for reproducibility. The moisture standard deviation and nominal correlation scale are editable. The actual representation is a piecewise-constant tile model; its nominal scale is not an exactly calibrated covariance length. Conductivity changes modestly with the same moisture-bin index as an explicit modelling prior.

The field is identical across frequencies, survey lines and the corresponding background control. A background model contains the same natural layers and heterogeneity, with the grave and its contents removed. Subtracting its raw field histories isolates the complete burial perturbation, including the grave cut; it is an ideal synthetic control that generally cannot be measured at the real site. To isolate skeleton or metal contributions, create additional model variants with one target class disabled.

Use a new output folder for a changed project. The exporter rejects reuse of a folder with a different project hash, preventing old traces being mistaken for the new acquisition.

## Numerical and acquisition design

The cell size, along-profile trace interval, cross-profile line interval and solver time step are independent. Survey positions and offsets snap to the electromagnetic mesh. The manifest records actual trace spacing, line spacing and section coordinates. Use these values in sampling tests.

The planning screen assumes a Ricker pulse with upper working bandwidth fmax = 2.5 f0 (approximately a −30 dB amplitude tail). It recommends the smaller of:

- λmin / 15 in the slowest represented dielectric, including the wettest heterogeneity bin and the corrosion halo when enabled;
- one third of the smallest represented feature thickness.

This is a starting screen. Verify numerical dispersion, geometry rasterisation, mesh convergence, PML distance/thickness and time-window adequacy on the actual models. The shortest wavelength calculation uses a low-loss phase-velocity approximation; strongly conducting conditions require additional analysis.

Generation checks **every selected frequency**. `--allow-underresolved` explicitly permits a labelled screening experiment. Thin snapped boxes may disappear on a coarse mesh. In 3D, cylinders/spheres with subcell radii and coincident snapped cylinder endpoints are omitted rather than enlarged; the manifest records those primitive omissions. They must not be interpreted as simulated resolved features. The 2D voxelisation samples the true geometry at cell centres, so small objects may likewise be missed by the grid.

The spatial sampling screen uses v/(4 fmax). λ/4 at the centre frequency is shown as a heuristic reference, not a guaranteed vertical resolution. Intrinsic two-way soil loss is a plane-wave calculation, excluding geometrical spreading, reflections, antenna efficiency and noise.

A uniform 3D mesh becomes very expensive at the high-frequency end. Halving the cell size increases a fixed-volume 3D cell count by about 8 and usually doubles the timestep count. RAM figures are broad planning ranges, including a nominal fixed overhead; actual allocation depends on materials, solver options and platform. Runtime requires a measured pilot trace. `--threads` sets CPU OpenMP threads. No GPU is required; GPU/MPI orchestration is not part of this wrapper. Generated `.in` files can be used directly with an appropriately configured gprMax installation.

## 2D versus 3D

2D TM simulations sample an x/depth or y/depth section and imply invariance normal to that section. Their spreading, amplitudes and object shape differ from a finite 3D burial. The electric dipole must be normal to the 2D plane. They can test a sectional response economically but cannot establish the detectability of finite bones or cross-line sampling.

For 3D, set `survey.mode` to `3D`. Use one profile first:

```bash
python burial_gpr.py generate my_model.json --out three_d_profile
```

For a regular raster:

```bash
python burial_gpr.py generate my_model.json --out three_d_raster --raster
```

The runner processes every generated burial/background pair. For a full raster this can be a very large job. The antenna is a 1 A peak Ricker-driven ideal horizontal electric dipole, with an ideal receiver at the defined separation. Changing centre frequency holds peak source current fixed, not radiated energy or system sensitivity. Frequency comparisons do not model a real 100 MHz, 500 MHz or 2 GHz commercial antenna. Add a validated antenna model and measured system response for instrument-specific conclusions.

## Outputs and sampling tests

`collect` verifies matching trace timing, receiver positions, shapes and numbered acquisition order. It saves:

- raw burial, background and difference arrays as compressed NPZ;
- time in seconds and Tx profile coordinates in metres;
- three-panel radargrams with one common linear amplitude scale;
- complete 3D rasters as `f500_raster.npz`, with time × line × trace ordering.

NPZ files can be loaded with NumPy. No AGC, migration, noise, dewow or temporal resampling is applied. Source delay and the direct wave remain in the raw traces. Set the analysis time gate deliberately; it is not automatically converted from burial depth.

Thin a complete dense profile, using integer multiples of its actual spacing:

```bash
python burial_gpr.py sampling runs/f500_line000_burial_scan.npz \
  --spacings 0.01 0.02 0.05 0.10 \
  --gate-start-ns 15 --gate-end-ns 60 --out trace_comparison
```

Those example spacings match the default 1 cm interval on the 1 mm mesh. For a different cell size, change them to integer multiples of the actual interval in your manifest. Thinning uses no spatial antialias filter, intentionally reproducing undersampled acquisition. Each possible starting phase is assessed, with linear reconstruction only over the sampled support. Relative L2 error measures information loss under that reconstruction, not archaeological detection probability.

For a complete 3D raster:

```bash
python burial_gpr.py grid-sampling three_d_raster/f500_raster.npz \
  --trace-spacings 0.01 0.02 0.05 \
  --line-spacings 0.02 0.04 0.10 \
  --gate-start-ns 15 --gate-end-ns 60 --out grid_comparison
```

Again use integer multiples of the actual dense spacings. This computes a gated RMS map, thins both spatial axes, checks every possible grid phase, and measures bilinear map reconstruction error. The RMS map is not a migrated depth slice. Only complete, matched raster lines are assembled. Large phase sweeps and cubes need additional RAM and runtime.

## Verification included

Run:

```bash
python -m unittest -v
```

The tests cover CRIM endmembers, material passivity and Debye fits, decay limits, deterministic geometry, invalid moisture, background invariance, geometry intersections, sweep-wide mesh checks, mesh-aligned acquisition, and phase-aware 1D/2D spatial sampling metrics.

The development verification also ran actual gprMax 4.0.1 burial/background 2D cases, a dispersive single trace, and a 3D geometry construction check. Included `verification/` outputs come from the explicitly **under-resolved** `smoke_test.json`: their purpose is demonstrating file generation, solver execution, collection and sampling analysis. They are not numerical-convergence results or evidence that a skeleton can be identified. The smaller `quickstart.json` also ran paired three-trace calculations with a mesh that passes the initial planning screen; convergence remains to be established.

## Suggested research sequence

1. Calibrate local layer depths, seasonal water contents, bulk conductivities and target spectra where feasible.
2. Compare drained, damp and waterlogged hypotheses at 1200, 1300 and 1400 years; vary preservation parameters independently of age.
3. Begin with background + grave shaft; add residual wood, bones, vessels, metal and optional voids one class at a time.
4. Establish mesh and boundary convergence on a limited set of profiles.
5. Run orthogonal 3D profiles and both horizontal polarisations for the most informative scenarios.
6. Acquire one sufficiently dense 3D dataset per physical scenario and source band; thin it to compare trace and line spacing with all grid phases.
7. Add measured antenna response and realistic noise before interpreting detection performance. Keep equal display scales and document source normalisation.

## Sources

The numerical method and material command syntax follow the official gprMax documentation. The archaeological literature motivates preservation-dependent contrasts, rather than the illustrative decay constants.

- Damiata, B.N., Steinberg, J.M., Bolender, D.J. & Zoëga, G. (2013). Imaging skeletal remains with ground-penetrating radar: comparative results over two graves from Viking Age and Medieval churchyards on the Stóra-Seyla farm, northern Iceland. *Journal of Archaeological Science* 40, 268–278. https://www.fiskecenter.umb.edu/Staff/Steinberg/Damiata_etal_2013.pdf
- *Physical Parameters and Contrasts of Wooden Objects in Lacustrine Environment: Ground Penetrating Radar and Geoelectrics* (2020), *Geosciences* 10(4), 146. https://www.mdpi.com/2076-3263/10/4/146
- Birchak, J.R., Gardner, C.G., Hipp, J.E. & Victor, J.M. (1974). High dielectric constant microwave probes for sensing soil moisture. *Proceedings of the IEEE* 62, 93–98. https://doi.org/10.1109/PROC.1974.9388
- gprMax 4.0.1 installation and version: https://github.com/gprMax/gprMax/releases
- Material/input reference, including the illustrative single-Debye water constants: https://docs.gprmax.com/en/latest/input_hash_cmds.html
- Modelling guidance: https://docs.gprmax.com/en/latest/gprmodelling.html
- Output reference: https://docs.gprmax.com/en/latest/output.html

## Optional local browser interface

The `web/` folder contains the same geometry/material editor. From the toolkit root run `python -m http.server 8000` and open `http://localhost:8000/web/`. Keep that server local to your computer. The local interface exports JSON; use the Python commands above to generate and solve models. Its guide links to the included README rather than downloading a second toolkit.

## Linked section explorer (v1.1)

The Model space page contains an orbitable 3D view and simultaneous XY, XZ and YZ material sections. Drag the 3D view to orbit; use a wheel or touch gesture to zoom. The view uses WebGL where available and a projected software renderer otherwise.

- **Horizontal XY:** `sections.depth_m`, measured downward from the ground surface.
- **Longitudinal XZ:** `sections.y_m`, measured from the domain y origin.
- **Transverse YZ:** `sections.x_m`, measured from the domain x origin.

Sliders, numerical entries and draggable crosshairs update the linked coordinates. Coloured 3D handles can also move the planes. Changes are stored in the exported project JSON. Use the per-section +/− controls to magnify around the crosshair, and disable **Focus on burial** to inspect the entire domain. PNG exports contain the plotted section; the filename records its fixed coordinate. The sections are material maps, not simulated GPR images.

Use **Use for GPR** on either vertical section to select the profile direction and section offset. The 2D polarisation is updated to remain normal to its plane. Moving the currently selected vertical section also updates the survey offset. Sections outside the usable lateral survey margin remain available for inspection; they are not silently moved to a different simulation position. Antenna and cell-grid snapping still apply during export and are recorded in the manifest. The horizontal section is an inspection plane rather than a horizontal 2D GPR propagation model.

**Show lid** and **Section planes** control display only. Hidden infill in the 3D cutaway remains present in the material sections and simulator. Anatomical volume retention still follows the explicit preservation scenario; anatomical detail is not evidence of preservation at the real site.

The v1.1 verification includes exact browser/Python comparisons of all volume primitives, the portable soil lattice and all three material-section rasters. The shared template is in `anatomy.json`; `scripts/build_anatomy.py` in the Site source regenerates it. `geometry_version: 2` selects the detailed model. The packaged verification example is regenerated with v2. The original v1 geometry remains in the Python module for explicitly versioned legacy projects.
