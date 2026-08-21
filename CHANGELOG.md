# Changelog

See Git commit history for more detailed list of changes.

## v3.0.0/upcoming release

* Real-axis implementation of normal-state equations
* Logarithmic sampling of real axis
* Coulomb/static contribution to energy shift
* Separation of self-consistency steps to enable non-self-consistent solutions
* Improved rescaling of Coulomb pseudopotential using matrix inversion
* Calculation of quasiparticle density of states using analytic continuation
* OpenMP parallelization
* Model input data for one-dimensional chain
* Improved documentation including website with examples
* Case-insensitive input parameters
* New inputs: `cutoffP`, `diag`, `divdos`, `restoren`, `toln`, `chiC`, `align0`,
  `logscale`, `apade`, `realgw`, `krakro` (aliases: `steps`, `points`, `gap`)
* Changed inputs: `upper` replaces `clip`, `cdosfix` replaces `imitate`
* New outputs: `domega`, `Re[Sigma]`, `Im[Sigma]`, `phi`, `Re[phi]`, `Im[phi]`,
  `chiC`, `states`, `inspect`

## v2.0.0/2024-11-29

This is the version of the code used for the paper *Understanding the origin of
superconducting dome in electron-doped MoS&#x2082; monolayer* by N. Girotto
Erhardt, J. Berges, S. Ponc&eacute;, and D. Novko, npj 2D Mater. Appl. **9**, 44
(2024).

* Support for arbitrary Eliashberg spectral functions (beyond Einstein model)
* Improved "unscaling" of Coulomb pseudopotential
* Improved project structure with license, logo, and single source directory
* New inputs: `a2f`/`a2F`, `bands` (new alias: `DOS`)
* New outputs: `omegaE`, `omegaLog`, `omega2nd`

## v1.1.0/2020-08-10

This is the version of the code used for my dissertation *Many-body
instabilities in two-dimensional materials*, https://doi.org/10.26092/elib/250
(2020).

* Possibility to provide non-rescaled Coulomb parameter
* Calculation of quasiparticle density of states
* Improved Python wrapper (Python-3 support, box-shaped and step-like DOS)
* New inputs: `muC`, `unscale`, `lower`, `upper`, `eta`/`0+`

## v1.0.0/2016-10-09

This is the original version of the code used for my Master's thesis *On the
scope of McMillan's formula*, https://scipost.org/theses/132 (2016).
