## Using MDAnalysis on MARCC

### Why?
`MDAnalysis` is available on MARCC. However, the only version available is 0.15.0. This version does not support the file format versions which Gromacs 2016 uses, so we must use a more up-to-date version of `MDAnalysis`, specifically version 0.18.0. 

### How?
This is not entirely easy to do, as Python will prioritize globally installed packages over locally installed packages, so installing `MDAnalysis` ver. 0.18.0 using `pip install --user` doesn't seem to work (I also tried changing the order of Python's `site-packages` directory in the `PATH`, but it didn't seem to work consistently). Knowing this, as far as I can tell, the easiest/most convenient way to do this is use Python's `virtualenv` package, which is installed on MARCC. Also, since this sets up an entirely new environment which is self-contained, it prevents any of the packages from interfering with MARCC's global packages.

### `Virtualenv`
`virtualenv` sets up a virtual python environment, consisting of its own Python and pip installations and downloaded modules. The set up is fairly straightforward.

  1) Create the virtual environment directory using `virtualenv $YOURPATH/md_env`, where `$YOURPATH` is wherever you want the directory to be located (it really doesn't matter).
  2) Activate the virtual environment using `source $YOURPATH/md_env/bin/activate`. This activates the environment and changes the shell's `$PATH` accordingly. The terminal should change to reflect this environment being sourced by putting `(md_env)`, or something similar, at the beginning of the terminal line.
  3) Ensure everything is set up properly by using `which python` and `which pip` show `$YOURPATH/md_env/bin/python` and `$YOURPATH/md_env/bin/pip`. If everything worked properly, proceed.
  4) Install `MDAnalysis` and its dependencies:
     ```bash
     pip install backports.functools-lru-cache biopython certifi chardet citeproc-py cycler decorator duecredit funcsigs GridDataFormats gsd idna joblib kiwisolver lxml matplotlib MDAnalysis mmtf-python mock msgpack networkx numpy pbr pyparsing python-dateutil pytz requests scipy six subprocess32 urllib3
     ```
  5) Use `pip freeze | less |  grep "MDAnalysis"` to make sure the installed version of `MDAnalysis` is 0.18.0.
