# NeuroC

A collections of programs morphology cloning
The following packages can be found
- axon shrinker:

For each morphology of the FILES_FOLDER, remove the axon splice described by the corresponding annotation (ie. located between the end of the dendritic annotation and the start of the axonal annotation) and replace it by an intermediate vertical segment. For each input morphology,
    NSAMPLES output morphologies are generated, each with a different length of the replaced segment. Lengths spans from 0 to the length of initially spliced segment

- jitter:
Create clones of a given morphology with some jitter to make them all different. There are two kinds of jitter: rotation and scaling.
Rotational jitter rotates each section around its parent axis or around the PCA ([Principal Component Analysis](https://en.wikipedia.org/wiki/Principal_component_analysis)) of all descendant points.




# Installation

In a fresh virtualenv:
```bash
pip install  --index-url  https://bbpteam.epfl.ch/repository/devpi/bbprelman/dev/+simple/ neuroc[plotly]
```

# Usage
In a shell, do:

```bash
neuroc --help
```
to list all functionalities.


## Axon shrinker


```bash
neuroc axon_shrinker files_dir annotations_dir output_dir
```
to shrink axons.


## Rat to human scaling
```bash
neuroc scale rat-to-human HUMAN_DIR RAT_DIR MTYPE_MAPPING OUTPUT_DIR
```

Will scale the rat cells in RAT_DIR to human cells dimensions.
HUMAN\_DIR should be a dir with the following structure:
- Must be **only** composed of sub-folders whose filename is a layer name
- Each sub folder should be composed of morphology files whose first part of the filename before the '_' is considered as the **mtype**

RAT\_DIR should be a directory containing rat morphology files **and a neuronDB.xml file.

MTYPE\_MAPPING\_FILE is a YAML file containing a dictionary where:
- a key is a human mtype or **all**
- the value is a list of rat mtypes to associate with the key. Or a list of one 'all' element
```
ls
$ RAT_DIR

$ RAT_DIR\L1
$ RAT_DIR\L2
$ ...
$ RAT_DIR\L6

$ RAT_DIR\L1\AC_cell_name.swc
$ RAT_DIR\L1\BTC_cell_name.swc
$ ...
```
