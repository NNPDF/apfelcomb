# apfelcomb
Generates FK tables for NNPDF fits.

## Project summary and aim

The aim of `apfelcomb` is to provide a set of tools for the generation
of FK tables / APFELgrid tables for DIS, DYP and Hadronic processes.

These tools are linked to the `libnnpdf` API for the
data/kinematics/process elaboration, `APFEL` for the evolution
kernels, and `nnpdfcpp` for the theory definitions in
`data/theory.db`. Each dataset is configurable via databases located
in the `db` folder.

### Release and Tag policy

The library is tagged and released when a major and stable status is achieved. 
Tags and releases do not necessarily follows the NNPDF releases.

### Code development policy/rules

Developers must never commit code structure modifications to master. The development pattern should follow these rules:
- Open an issue explaining your bug or feature request. If you report a bug, post information to reproduce it.
- The resolution of issues must be performed in a new branch through a pull request.
- If you have already a local version of the code that you would like to merge in the master, open a pull request.
- The pull request must be reviewed by at least 2 core developers.

### Code style

Uses C++11 features.

### Continuous integration (CI)

CI is actually not implemented in the current repository.

### Testing

Testing is actually not implemented in the current repository.

## Installation

`apfelcomb` depends on the following libraries:

- APFEL
- libnnpdf

### Prerequisites

* Download and compile latest version of APFEL trunk
```Shell
git clone https://github.com/scarrazza/apfel.git
./configure --prefix=/your/location/
make
make clean
```
* Add paths to your .bashrc

### Install libnnpdf

```Shell
cd /nnpdfcpp/trunk/libnnpdf/
# Configure and install
./configure --prefix=/your/location --enable-safemode
# (Note the latest option is necessary)
make
make install
```
* Add paths to your .bashrc

### APPLgrid database
* For each new observable add info on APPLgrid to the applgrid.db in the /nnpdfcpp/trunk/apfelcomb directory. All entries are trivial but three
- maxprec = maximum relative precision of data in the new dataset
- xmin = this is obtained by running
```Shell
  ./appl_gridinfo <APPLGRID_location>
```
- nxpt = this is obtained by running
```Shell
  ./appl_optgrid <entry number in the database> <entry number in the theory database>
```
(The theory database.nb is in /data/theory.db. For example 11 is the
entry for NLO as=0.118 FONLLB...)  Nxpt is the minimum number of
points in the FK grid to obtain a max precision better than the
precision of the data.  NB: if two tables need to be merged nxpt must
be the same

### GENFK scripts
* For each new set, you need to add a corresponding genFK script into the /genFK_scripts/ folder
  In this script you should perform any post-processing required (merging FK tables, normalising, etc).
  Minimally this script should perform a check that the FK table file is present.
  Lots of examples are provided in /genFK_scripts/

### Running APFELcomb
* Once the applgrid.db is all set, run
```Shell
./appl_comb <entry number in the database> <entry number in the theory database>
```
* In a short time you should get your FK tables

## Documentation

### Code documentation

The code is documented with Doxygen, if you find methods or classes
not fully documented open a issue request.

### Layout documentation

For specifications about data please check the `nnpdfcpp` repository in `data/doc`.
