from jsonextended import plugins

from ejplugins.crystal import (CrystalOutputPlugin, CrystalSCFLogPlugin, CrystalDOSPlugin, BANDPlugin,
                               ECH3CubePlugin, ECH3OutPlugin)
from ejplugins.qespresso import (QEmainPlugin, QEbandPlugin, QEChargeDensityPlugin,
                                 QELowdinPlugin, QEnscfPlugin, QEdosPlugin)
from ejplugins.gulp import GULPOutPlugin
from ejplugins.lammps import LAMMPSAtomDumpPlugin, LAMMPSSysDumpPlugin

# optional dependencies
try:
    from ejplugins.pymatgen_decode import Encode_Pymatgen
    from ejplugins.ase_decode import Encode_ASE
    from ejplugins.cif import CIFPlugin
except ImportError:
    pass

from ejplugins.utils import load_test_file, validate_against_schema


def load_all_parsers():
    """ load all parsers

    Returns
    -------
    errors: list

    """
    parsers = [CrystalOutputPlugin, CrystalSCFLogPlugin, CrystalDOSPlugin,
               BANDPlugin, ECH3CubePlugin, ECH3OutPlugin,
               QEmainPlugin, QEbandPlugin, QEChargeDensityPlugin, QELowdinPlugin,
               GULPOutPlugin,
               LAMMPSAtomDumpPlugin, LAMMPSSysDumpPlugin]
    try:
        from ejplugins.cif import CIFPlugin
        parsers.append(CIFPlugin)
    except ImportError:
        pass

    return plugins.load_plugin_classes(parsers, "parsers")


def load_all_encoders():
    try:
        from ejplugins.pymatgen_decode import Encode_Pymatgen
        from ejplugins.ase_decode import Encode_ASE
    except ImportError:
        raise ImportError("pymatgen and/or ase not installed")
    return plugins.load_plugin_classes([Encode_Pymatgen, Encode_ASE])


__version__ = "0.10.0"
