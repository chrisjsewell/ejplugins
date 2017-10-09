
from jsonextended import plugins

from ejplugins.crystal import (CrystalOutputPlugin, DOSSPlugin, BANDPlugin, ECH3CubePlugin, ECH3OutPlugin)
from ejplugins.qespresso import QEmainPlugin, QEbandPlugin
from ejplugins.cif import CIFPlugin
from ejplugins.gulp import GULPOutPlugin
from ejplugins.lammps import LAMMPSAtomDumpPlugin, LAMMPSSysDumpPlugin

__version__ = "0.2.0"


def load_all_parsers():
    """ load all parsers

    Returns
    -------
    errors: list

    """
    parsers = [CrystalOutputPlugin, DOSSPlugin, BANDPlugin, ECH3CubePlugin, ECH3OutPlugin,
               QEmainPlugin, QEbandPlugin, CIFPlugin, GULPOutPlugin, LAMMPSAtomDumpPlugin, LAMMPSSysDumpPlugin]

    return plugins.load_plugin_classes(parsers, "parsers")
