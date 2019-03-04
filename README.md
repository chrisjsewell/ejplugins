# ejplugins

[![Build_Status](https://travis-ci.org/chrisjsewell/ejplugins.svg?branch=master)](https://travis-ci.org/chrisjsewell/ejplugins)
[![PyPI](https://img.shields.io/pypi/v/ejplugins.svg)](https://pypi.python.org/pypi/ejplugins/)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/ejplugins/badges/version.svg)](https://anaconda.org/conda-forge/ejplugins)

Parser plugins for the [jsonextended package](https://jsonextended.readthedocs.io) and validation schema, to convert
output files from materials simulation packages to a JSON format.

## Usage

### Conda (recommended)

    >> conda install -c conda-forge ejplugins

### Pip

    >> pip install ejplugins

or, to use Pymatgen and ASE parsers:

    >> pip install ejplugins[science]

Either use independently:

```python
from ejplugins.qespresso import QEmainPlugin
with open("path/to/my.qe.out") as f:
    output = QEmainPlugin.read_file(f)
```

Or with `jsonextended`:

```python
from ejplugins.qespresso import QEmainPlugin
from jsonextended import plugins, ejson
plugins.load_plugin_classes([QEmainPlugin])

ejson.to_dict("path/to/qespresso/outputs")
```

See ejplugins/test_files for example input/outputs.