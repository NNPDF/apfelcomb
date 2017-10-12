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
  Defines the path results are written to.
- **DATA_PATH** *(default ../nnpdfcpp/data)*
  Defines the path to the COMMONDATA repository.
- **APPL_PATH** *(default ../applgrids)*
  Defines the path to the nnpdf-applgrid repository.
- **DB_PATH** *(default ./db)*
  Defines the path to the APFELcomb database.

The defaults are configured assuming that both the nnpdfcpp and applgrid repositories are
located at `../`

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

where `<app/dis/dyp>` specifies whether the subgrid is in the applgrid, DIS or DYP subgrid categories in the database,
`<subgrid id>` is the corresponding ID in that database (visible in the disp\_grids script) and `<theory id>` specifies
the desired NNPDF theory index (the entry in nnpdfcpp/data/theory.db). As an example:

```Shell
./apfel_comb app 500 65 
```
Will generate the subgrid for CDFZRAP and theory 65 (NNPDF3.1 NNLO perturbative charm). The resulting FK subgrid
will be written out to 

```Shell
$RESULTS_PATH/theory_<theoryID>/subgrids/FK_<setname>_<subgridID>.dat.
```

Once all the relevant subgrids for the desired dataset(s) are generated, you should run

```Shell
./merge_allgrids.py <theory id>
```

which will loop over all datasets and attempt to merge their subgrids into a complete FK table. The resulting final
FK table should be stored at

```Shell
$RESULTS_PATH/theory_<theoryID>/fastkernel/FK_<setname>.dat.
```
### Implementing a new FK table

To implement a new FK table you must first add a corresponding entry into the apfelcomb database under the `grids` table.
These entries are comprised of the following fields.

- **id**		- The primary key identifier of the FK table
- **setname**		- The COMMONDATA SetName of the corresponding dataset
- **name**		- The name of the FK table
- **description**	- A one-line description of the FK table
- **nx**		- The number of x-grid interpolation points
- **positivity**	- A flag specifying if the FK table is a positivity set
- **source**		- Specifies if the corresponding subgrids are [APP/DIS/DYP]

Here it is important to note that setname and name may be different in the case of compound observables such
as ratios, where multiple FK tables are required to compute predictions for a single dataset. The `nx` parameter
specified the interpolation accuracy of the dataset (this must currently be tuned by hand). The `positiviy` parameter
restricts the observable to NLO matrix elements and disables target-mass corrections.

Once this entry is complete, you must move on to adding entries in the corresponding subgrid table.

### Implementing a new APPLgrid subgrid 

To add a new APPLgrid-based subgrid, you must add a corresponding entry into the `app\_subgrids` table of the
apfelcomb database. One entry should be added for each APPLgrid making up the final target FK table.
The entries have the following fields:

- **id** 	- The primary key identifier of the subgrid. 
- **fktarget**	- The name of the FK table this subgrid belongs to. 
- **applgrid**	- The filename of the corresponding APPLgrid.
- **fnlobin**   - The fastNLO index if the table is a fastNLO grid, or -1 if not.
- **ptmin**	- The minimum perturbative order (1 when the LO is zero, 0 if not).
- **pdfwgt**	- A boolean flag, 1 if the APPLgrid has PDF weighting, 0 if not.
- **ppbar**	- A boolean flag, 1 if the APPLgrid should be transformed to ppbar beams, 0 if not.
- **mask**	- A boolean mask, specifiying which APPLgrid entries should be considered datapoints.
- **operators** - A list of operators to handle certain special cases (see below).

MOREDESCRIPTIONHERE

### Operators

### Helper scripts

TODO

### C-factor scaling

TODO


TODO
