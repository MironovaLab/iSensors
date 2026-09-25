# iSensors 1.3.0 (in development)

## Panel fixes

* `ATH-aux-trans-ARF` (Arabidopsis *trans* ARF panel) now contains all 23
  members of the Arabidopsis ARF family. ARF1 (AT1G59750) was missing from the
  1.2.3 panel (22 genes); `gene_metadata` gains the corresponding row. The
  activator-only `ATH-aux-trans-A-ARF` panel and all non-Arabidopsis
  `*-aux-trans-ARF` panels are unchanged. This fix was first built locally as
  1.2.4 and is released here.

## Documentation

* `main` is the default branch again and holds the current release, so
  `install_github("MironovaLab/iSensors")` installs it (it previously installed
  the old 1.1 release).
* README rewritten: working installation, quick start, how the
  scores are calculated, a catalogue of the 681 panels and 100 species codes.
* `CalcSensors()`: the `seurLayer` default is documented as `"data"` (it was
  described as `"RNA"`), `normBy` is documented, and the examples run.
* `LoadSensors()`: filter examples use the codes in the panel file names
  (`species = "ATH"`, `hormone = "aux"`); they previously used `"AT"` and
  `"auxin"`, which match no panels.

## Clean-up

* Comments translated to English; commented-out code removed. No change in
  behaviour.
* Removed developer scratch scripts (`scripts/`) and Jupyter checkpoints; the
  tutorial script is now `tutorial/tutorial.R`.

## Known issues (to be addressed in 1.3.0)

* `"mean_normed"` and `"median_normed"` with `normBy = "cols"` divide by the
  wrong cell's mean or median (the corrected version is on `iSensors-dev`).
* Random control panels are re-drawn for each requested signal type.
