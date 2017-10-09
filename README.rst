ejplugins
=========

|Build_Status| |PyPI|

Parser plugins for the `jsonextended
package <https://jsonextended.readthedocs.io>`__ and validation schema,
to convert output files from materials simulation packages to a JSON
format.

Usage
-----

::

    pip install ejplugins
    pip install jsonextended

Either use independently:

.. code:: python

    from ejplugins.qespresso import QEmainPlugin
    with open("path/to/my.qe.out") as f:
        output = QEmainPlugin.read_file(f)

Or with ``jsonextended``:

.. code:: python

    from ejplugins.qespresso import QEmainPlugin
    from jsonextended import plugins, ejson
    plugins.load_plugin_classes([QEmainPlugin])

    ejson.to_dict("path/to/qespresso/outputs")

See ejplugins/test\_files for example input/outputs.

.. |Build_Status| image:: https://travis-ci.org/chrisjsewell/ejplugins.svg?branch=master
   :target: https://travis-ci.org/chrisjsewell/ejplugins
.. |PyPI| image:: https://img.shields.io/pypi/v/ejplugins.svg
   :target: https://pypi.python.org/pypi/ejplugins/
