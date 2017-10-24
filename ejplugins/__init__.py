from jsonextended import plugins

from ejplugins.crystal import (CrystalOutputPlugin, CrystalSCFLogPlugin, DOSSPlugin, BANDPlugin,
                               ECH3CubePlugin, ECH3OutPlugin)
from ejplugins.qespresso import QEmainPlugin, QEbandPlugin, QEChargeDensityPlugin
from ejplugins.cif import CIFPlugin
from ejplugins.gulp import GULPOutPlugin
from ejplugins.lammps import LAMMPSAtomDumpPlugin, LAMMPSSysDumpPlugin

from ejplugins.pymatgen_decode import Encode_Pymatgen
from ejplugins.ase_decode import Encode_ASE

from ejplugins.utils import load_test_file, validate_against_schema


def load_all_parsers():
    """ load all parsers

    Returns
    -------
    errors: list

    """
    parsers = [CrystalOutputPlugin, CrystalSCFLogPlugin, DOSSPlugin, BANDPlugin, ECH3CubePlugin, ECH3OutPlugin,
               QEmainPlugin, QEbandPlugin, QEChargeDensityPlugin, CIFPlugin, GULPOutPlugin, LAMMPSAtomDumpPlugin, LAMMPSSysDumpPlugin]

    return plugins.load_plugin_classes(parsers, "parsers")


def load_all_encoders():
    return plugins.load_plugin_classes([Encode_Pymatgen, Encode_ASE])


__version__ = "0.4.4"
