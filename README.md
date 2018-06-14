# NeuroC

A collections of programs morphology cloning
  The following packages can be found
    - axon shrinker:

      For each morphology of the FILES_FOLDER, remove the axon splice described by the corresponding annotation (ie. located between the end of the dendritic annotation and the start of the axonal annotation) and replace it by an intermediate vertical segment. For each input morphology,
    NSAMPLES output morphologies are generated, each with a different length of the replaced segment. Lengths spans from 0 to the length of initially spliced segment


# Installation

In a fresh virtualenv:
```bash
git clone  --index-url  https://bbpteam.epfl.ch/repository/devpi/bbprelman/dev/+simple/ ssh://bbpcode.epfl.ch/nse/NeuroC
cd NeuroC
pip install .
```

# Usage
In a shell, do:

```bash
neuroc --help
```
to list all functionalities, or:


```bash
neuroc axon_shrinker files_dir annotations_dir output_dir
```
to shrink axons.
