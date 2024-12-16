# RomanHLWASTools

This is the set of tools used for tiling the Roman HLWAS in the survey definition committee report.

## Building the tiling

The `surveyscript/` directory contains an iPython notebook, `tiling.ipynb`, that builds the list of pointings. It first builds maps of the sky properties, then builds a survey footprint, then tiles it, and finally outputs the survey information. The internal map projection used for `tiling.ipynb` is Equatorial, Cylindrical Equal Area.

Inputs used:
* The SFD dust map (`ebv_1e4.dat.gz`: this is a flat gzipped text file in the required projection, in units of 1e-4 magnitudes).
* The table of slew and settle times provided by the Roman project, `slew-settle-table.txt` (first column is the slew angle in degrees, the second is the slew and settle time in seconds).
* `bscat.dat` is a table of stars to avoid (colums are RA in degrees; Dec in degrees; J band AB mag from 2MASS) — the tiling code moves pointings that contain these objects.
* `letters.dat` is a bitmap for labeling plots.

The outputs include various values and plots in the iPython notebook itself. It also includes three output files:
* `skychart.png`: a map of the sky showing the Medium and Wide regions
* `footprint.fits`: a FITS image of the sky with WCS (suitable for being imported into another application)
* `hlwas_tiling.txt`: a table of the pointings that were generated

**The output tiling file used for the report can be found in** `tilingdata/hlwas_tiling_241206.txt.gz`**.**

## Making coverage maps

The simple module `coveragemap.py` is provided for making detailed FITS coverage maps from the list of tilings and a WCS for the coverage map. This is best done following the example in `example_coverage.py`, which calls the tiling file `tilingdata/hlwas_tiling_241206.txt.gz` (as used in the HLWAS survey definition committee report). See also the instructions in `coverage_multiexp` (in `coveragemap.py`).

These maps take into account field distortions and are good for, e.g., survey window functions and effects of chip gaps, but the distortion map is a cubic approximation and has errors of a few arcsec.

Note that the example (covering 8x8 degrees at 1 sample every 7.2 arcsec) is designed to run in ~ 1 hour on a laptop; mosaiced maps of the whole survey will take a lot longer.

## References

The dust map used is a reprojection of the SFD map:  
Schlegel, Finkbeiner, Davis, Astrophys. J. 500:525 (1998)

The Two Micron All Sky Survey (2MASS)  
M.F. Skrutskie, R.M. Cutri, R. Stiening, M.D. Weinberg, S. Schneider, J.M. Carpenter, C. Beichman, R. Capps, T. Chester, J. Elias, J. Huchra, J. Liebert, C. Lonsdale, D.G. Monet, S. Price, P. Seitzer, T. Jarrett, J.D. Kirkpatrick, J. Gizis, E. Howard, T. Evans, J. Fowler, L. Fullmer, R. Hurt, R. Light, E.L. Kopan, K.A. Marsh, H.L. McCallon, R. Tam, S. Van Dyk, and S. Wheelock, 2006, AJ, 131, 1163.
