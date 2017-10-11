# APFELcomb 
The APFELcomb project consists of a set of tools for the generation
of FK tables, which provide the mechanism for computing predictions
from theory in the NNPDF framework. Broadly speaking this is acheived
by taking DGLAP evolution kernels from `APFEL` and combining them with
interpolated parton-level observable kernels of various forms.

The mechanism behind APFELcomb is documented in [1605.02070].
The various data formats used in APFELcomb are described in `nnpdfcpp/data/doc/'.

## Prerequisites
APFELcomb depends on the following libraries

* **APFEL** *github.com/scarrazza/apfel.git*
* **libnnpdf** *github.com/NNPDF/libnnpdf*
* **APPLgrid 1.4.70-nnpdf** *github.com/NNPDF/external/applgrid-1.4.70-nnpdf*

And datafiles from
* **nnpdfcpp** *github.com/NNPDF/nnpdfcpp*
* **nnpdf-applgrids** *github.com/NNPDF/applgrids*

## Compilation and setup 

Compilation flags and various paths are defined in `Makefile.inc`.
These are mostly inferred from package-config files with the exception of

- **RESULTS_PATH** *(default ./results)*
  Defines the path results are written to
- **DATA_PATH** *(default ../nnpdfcpp/data)*
  Defines the path to the COMMONDATA repository
- **APPL_PATH** *(default ../applgrids)*
  Defines the path to the nnpdf-applgrid repository
- **DB_PATH** *(default ./db)*
  Defines the path to the APFELcomb database

The defaults are configured assuming that both the nnpdfcpp and applgrid repositories are
located at `../`.

With these paths set, compiling the APFELgrid code should be as simple as
```Shell
make
```

## Structure and generation process

Each FK table is generated piecewise in one or more `subgrids`. The subgrids
implemented in APFELcomb can be displayed by running the script

```Shell
./scripts/disp_grids.py
```

Typically DIS and FKGenerator Drell-Yan tables are made of only one subgrid, whereas
FK tables generated from APPLgrids have one subgrid per APPLgrid file. How subgrids
are merged into grids, and the generation parameters of each subgrid, is specified in
the `db/apfelcomb.db` database.

Generating an individual subgrid is performed by running

```Shell
./apfel_comb <source=app/dis/dyp> <subgrid id> <theory id>
```

where <app/dis/dyp> specifies whether the subgrid is in the applgrid, DIS or DYP subgrid categories in the database,
<subgrid id> is the corresponding ID in that database (visible in the disp\_grids script) and <theory id> specifies
the NNPDF theory index desired (the entry in nnpdfcpp/data/theory.db). As an example:

```Shell
./apfel_comb app 500 65 
```
Will generate the subgrid for CDFZRAP and theory 65 (NNPDF3.1 NNLO perturbative charm). The resulting FK subgrid
will be written out to 

```Shell
$RESULTS\_PATH/theory\_<theoryID>/subgrids/FK\_<setname>\_<subgridID>.dat.
```

Once all the relevant subgrids for the desired dataset(s) are generated, you should run

```Shell
./merge_allgrids.py <theory id>
```

which will loop over all datasets and attempt to merge their subgrids into a complete FK table. The resulting final
FK table should be stored at

```Shell
$RESULTS\_PATH/theory\_<theoryID>/fastkernel/FK\_<setname>.dat.
```

### Helper scripts

TODO

### C-factor scaling

TODO

### Implementing a new dataset

TODO
