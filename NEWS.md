# iSensors 1.3.0 (in development)

## Breaking changes

* The `"mean_normed"` and `"median_normed"` signals and the `normBy` argument of
  `CalcSensors()` were removed. They belonged to an earlier scoring approach and,
  with `normBy = "cols"`, divided by the wrong cell's mean or median. Asking for
  them gives an error that names the replacement.
* `CalcSensors()` now computes `signals = "mean"` by default (was
  `"mean_normed"`).
* Panels must contain at least 3 genes. `iSensorsTransPanelCreate()` and
  `iSensorsCisTransPanelCreate()` refuse to create a smaller panel,
  `LoadSensors()` refuses to load one, and random panels must be at least 3
  genes. `CalcSensors()` skips panels with fewer than 3 genes detected in the
  data, with a warning that names them.
* Five shipped panels with 2 genes were removed: `AAG-aux-trans-IAA`,
  `CBR-aux-trans-IAA`, `MPO-aux-trans-IAA`, `SHI-aux-trans-IAA` and
  `SHI-aux-trans-PAT`. 676 panels remain; all 100 species keep other panels.

## Panel fixes

* `ATH-aux-trans-ARF` (Arabidopsis *trans* ARF panel) now contains all 23
  members of the Arabidopsis ARF family. ARF1 (AT1G59750) was missing from the
  1.2.3 panel (22 genes); `gene_metadata` gains the corresponding row. The
  activator-only `ATH-aux-trans-A-ARF` panel and all non-Arabidopsis
  `*-aux-trans-ARF` panels are unchanged. This fix was first built locally as
  1.2.4 and is released here.

## Bug fixes

* A panel of which only one gene was found in the data (after genes with zero
  variance are left out) got that gene's mean over all cells, the same value in
  every cell. Such panels are now skipped (see the 3-gene minimum above).
* `iSensorsTransPanelCreate()` and `iSensorsCisTransPanelCreate()` work when
  called from inside a function (a script's own function, a loop in a function,
  a test). They previously failed with "object not found" and saved no file.
  They no longer overwrite and then delete an object of the same name as the
  panel in the calling environment, and they return the panel invisibly.
* `iSensorsCisTransPanelCreate()` works without attaching `universalmotif`
  first; it gives a clear error if the package is not installed.

## Documentation

* `main` is the default branch again and holds the current release, so
  `install_github("MironovaLab/iSensors")` installs it (it previously installed
  the old 1.1 release).
* README rewritten, with the iSensors logo and a link to the online tutorial:
  working installation, quick start, how the scores are calculated, a catalogue
  of the 676 panels and 100 species codes.
* `CalcSensors()`: the `seurLayer` default is documented as `"data"` (it was
  described as `"RNA"`) and the examples run.
* `LoadSensors()`: filter examples use the codes in the panel file names
  (`species = "ATH"`, `hormone = "aux"`); they previously used `"AT"` and
  `"auxin"`, which match no panels.

## Clean-up

* Comments translated to English; commented-out code removed. No change in
  behaviour.
* The test suite now runs with `R CMD check` (the runner was in the wrong folder
  and disabled), with new tests for the scores and the panel fixes.
* Removed developer scratch scripts (`scripts/`) and Jupyter checkpoints; the
  tutorial script is now `tutorial/tutorial.R`.

## Known issues

* Random control panels are re-drawn for each requested signal type.
