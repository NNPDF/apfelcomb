# APFELcomb 
The APFELcomb project consists of a set of tools for the generation
of FK tables, which provide the mechanism for computing predictions
from theory in the NNPDF framework. Broadly speaking this is achieved
by taking DGLAP evolution kernels from `APFEL` and combining them with
interpolated parton-level observable kernels of various forms.

The mechanism behind APFELcomb is documented in `[1605.02070]`.

The various data formats used in APFELcomb are described in `nnpdfcpp/data/doc/`.

**Table of Contents**
- [Prerequisites](#prerequisites)
- [Compilation and setup](#compilation-and-setup)
- [Structure and generation process](#structure-and-generation-process)
- [Implementing a new FK table](#implementing-a-new-fk-table)
  - [Implementing a new APPLgrid subgrid](#implementing-a-new-applgrid-subgrid)
  - [Implementing a new DIS or DYP subgrid](#implementing-a-new-dis-or-dyp-subgrid)
  - [Subgrid operators](#subgrid-operators)
  - [Compound files and C-factors](#compound-files-and-c-factors)
  - [Important note on subgrid ordering](#important-note-on-subgrid-ordering)
  - [Important note on committing changes](#important-note-on-committing-changes)
- [Helper scripts](#helper-scripts)
- [Generating a complete theory](#generating-a-complete-theory)

## Prerequisites
APFELcomb depends on the following libraries

* **APFEL** *github.com/scarrazza/apfel.git*
* **nnpdf** *github.com/NNPDF/nnpdf*
* **APPLgrid 1.4.70-nnpdf** *github.com/NNPDF/external/applgrid-1.4.70-nnpdf*

And data files from
* **nnpdf-applgrids** *github.com/NNPDF/applgrids*

## Compilation and setup 
Compilation flags and various paths are defined in `Makefile.inc`.
These are mostly inferred from package-config files with the exception of

- **RESULTS_PATH** *(default ./results)*
  Defines the path results are written to.
- **APPL_PATH** *(default ../applgrids)*
  Defines the path to the nnpdf-applgrid repository.
- **DB_PATH** *(default ./db)*
  Defines the path to the APFELcomb database.

The defaults are configured assuming that both the nnpdfcpp and applgrid repositories are
located at `../`

With these paths set, compiling the APFELcomb code should be as simple as
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
the `db/apfelcomb.db` database. The database itself is not stored in the repository,
but it is built from the sqlite dump at `db/apfelcomb.dat`. This is done automatically
by the APFELcomb makefile.


Generating an individual subgrid is performed by running

```Shell
./apfel_comb <source=app/dis/dyp> <subgrid id> <theory id>
```

where `<app/dis/dyp>` specifies whether the subgrid is in the applgrid, DIS or DYP subgrid categories in the database,
`<subgrid id>` is the corresponding ID in that database (visible in the `disp\_grids` script) and `<theory id>` specifies
the desired NNPDF theory index (the entry in nnpdfcpp/data/theory.db). As an example:

```Shell
./apfel_comb app 500 65 
```
Will generate the subgrid for CDFZRAP and theory 65 (NNPDF3.1 NNLO perturbative charm). The resulting FK subgrid
will be written out to 

```Shell
$RESULTS_PATH/theory_<theory id>/subgrids/FK_<setname>_<subgrid id>.dat.
```

Once all the relevant subgrids for the desired dataset(s) are generated, you should run

```Shell
./merge_allgrids.py <theory id>
```

which will loop over all datasets and attempt to merge their subgrids into a complete FK table. The resulting final
FK table should be stored at

```Shell
$RESULTS_PATH/theory_<theory id>/fastkernel/FK_<setname>.dat.
```
## Implementing a new FK table
To implement a new FK table you must first add a corresponding entry into the apfelcomb database under the `grids` table.
These entries are comprised of the following fields.

- **id**		- The primary key identifier of the FK table.
- **setname**		- The COMMONDATA set name of the corresponding dataset.
- **name**		- The name of the FK table.
- **description**	- A one-line description of the FK table.
- **nx**		- The number of x-grid interpolation points.
- **positivity**	- A flag specifying if the FK table is a positivity set.
- **source**		- Specifies if the corresponding subgrids are [APP/DIS/DYP].

Here it is important to note that **setname** and **name** may be different in the case of compound observables such
as ratios, where multiple FK tables are required to compute predictions for a single dataset. The `nx` parameter
specifies the interpolation accuracy of the dataset (this must currently be tuned by hand). The `positivity` parameter
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
- **ppbar**	- A boolean flag, 1 if the APPLgrid should be transformed to *ppbar* beams, 0 if not.
- **mask**	- A boolean mask, specifying which APPLgrid entries should be considered data points.
- **operators** - A list of operators to handle certain special cases (see below).

Here the mask should have as many entries as APPLgrid bins and each boolean value should be separated by
a space. For example, for an applgrid with five bins where we want to exclude the penultimate bin, the mask
would be:

```
1 1 1 0 1
```
The applgrid filename assumes that the grid can be found at

```Shell
$APPL_PATH/<setname>/<applgrid>
```
where `APPL_PATH` is defined in Makefile.am, `<setname>` is the corresponding COMMONDATA set name specified in the grids table,
and `<applgrid>` is specified in the field described above.

### Implementing a new DIS or DYP subgrid 
New DIS or DYP subgrids should be entered respectively into the `dis_subgrids` or `dyp_subgrids` tables of the apfelcomb database.
Typically only one subgrid is needed per DIS or DYP FK table. Each subgrid entry has the following fields:

- **id**	- The primary key identifier of the subgrid
- **fktarget**	- The name of the FK table this subgrid belongs to
- **operators** - A list of operators to handle certain special cases (see below).

For DIS there is one additional field:
- **process**	- The process string of the observable (e.g DIS\_F2P, see APFEL)

### Subgrid operators
Subgrid operators are used to provide certain subgrid-wide transformations that can be useful in certain circumstances.
They are formed by a key-value pair with syntax:

```Shell
<KEY>:<VALUE>
```
If using multiple operators, they should be comma-separated. Currently these operators are implemented:

- \*:*V* - Duplicate the subgrid data point (there must be only one for this operator) *V* times.
- +:*V*  - Increment the starting data point index of this subgrid by *V*.
- N:*V*  - Normalise all data points in this subgrid by *V*.

The \* operator is typically used for normalised cross-sections, where the total cross-section computation (a single data point)
must be duplicated *N\_dat* times to correspond to the size of the COMMONDATA file. The + operator is typically used to compensate
for missing subgrids, for example when a COMMONDATA file begins with several data points that cannot yet be computed from theory,
the + operator can be used to 'skip' those points. The N operator is used to perform unit conversions or the like.

### Compound files and C-factors
If your new dataset is a compound observable (that is, theory predictions are a function of more than one FK-product) then you
should write a corresponding COMPOUND file as described in the data specifications at `nnpdfcpp/data/doc/`. This compound
file should be stored in the APFELcomb repository under the `compound` directory.

C-factors should be in the format once again specified in `nnpdfcpp/data/doc/`, and stored in the nnpdfcpp repo under
the `nnpdfcpp/data/N*LOCFAC/` directory.

### Important note on subgrid ordering
If your FK table consists of more than one subgrid to be merged into a single table, then the ordering
of the subgrids in their subgrid **id** is vital. The `merge_allgrids.py` script will merge the subgrids
in order of their **id**. So if you are constructing an FK table for a merged W+/W-/Z dataset, it is crucial
that the ordering of the corresponding W+/W-/Z subgrids in id matches the ordering in COMMONDATA.

### Important note on committing changes
If you have made a modification to the apfelcomb.db database, once you are happy with it you *must* export it to the 
plain-text dump file at `db/apfelcomb.dat`. This file must then be committed. It is important to note that the binary
sqlite database is not stored in the repository.

A helper script is provided to do this. If you want to convert your binary database to the text dump, run
`db/generate_dump.sh` and then commit the resulting `apfelcomb.dat` file.

## Helper scripts

Several helper scripts are provided to make using APFELcomb easier (particularly when generating a full set of FK tables for a particular theory).
- `scripts/disp_grids.py` displays a full list of APPLgrid, DIS or DYP subgrids implemented in APFELcomb.
- `run_allgrids.py [theoryID] [job script]` scans the results directory and submits jobs for all missing subgrids for the specified theory.
- `test_submit.py` is an example [job script] to be used for `run\_allgrids.py`. These scripts specify how jobs are launched on a given cluster.
- `hydra_submit.py` is the [job script] for the HYDRA cluster in Oxford.
- `merge_allgrids.py [theoryID]` merges all subgrids in the results directory for a specified theory into final FK tables. This does not delete subgrids.
- `finalise.sh [theoryID]` runs C-factor scaling, copies COMPOUND files, deletes the subgrids, and finally compresses the result into a theory.tgz file ready for upload.


## Generating a complete theory

The general workflow for generating a complete version of theory 53 on the hydra cluster is then.
```Shell
./run_allgrids.py 53 ./hydra_submit.sh # Submit all APFELcomb subgrid-jobs
# Once all subgrid jobs have successfully finished
./merge_allgrids.py 53 # Merge subgrids into FK tables
# If merging is successful
./finalise.sh 53
# Results in a final theory at ./results/theory_53.tgz
```
