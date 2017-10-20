import os
import pytest
import json
from jsonextended import plugins, edict, ejson
from jsonextended.encoders.ndarray import Encode_NDArray

from ejplugins.crystal import (CrystalOutputPlugin, CrystalSCFLogPlugin,
                               DOSSPlugin, BANDPlugin, ECH3CubePlugin, ECH3OutPlugin)
from ejplugins.qespresso import QEmainPlugin, QEChargeDensityPlugin
from ejplugins.cif import CIFPlugin
from ejplugins.gulp import GULPOutPlugin
from ejplugins.lammps import LAMMPSAtomDumpPlugin, LAMMPSSysDumpPlugin
from ejplugins.utils import _get_all_schema_name

import jsonschema

file_folder = os.path.join(os.path.dirname(__file__), "test_files")
schema_folder = os.path.join(os.path.dirname(__file__), "schema")


@pytest.mark.parametrize("testplugin,filename", [
    (CrystalOutputPlugin, "scf_only.crystal.out"),
    (CrystalOutputPlugin, "scf_and_opt.crystal.out"),
    (CrystalOutputPlugin, "scf_and_opt_slab.crystal.out"),
    (CrystalSCFLogPlugin, "scf_and_opt.crystal.scflog"),
    (DOSSPlugin, "crystal.doss.f25"),
    (BANDPlugin, "crystal.band.f25"),
    (ECH3CubePlugin, "crystal.ech3_dat.prop3d"),
    (ECH3OutPlugin, "crystal.ech3.out"),
    (QEmainPlugin, "scf.qe.out"),
    (QEmainPlugin, "scf_with_fermi.qe.out"),
    (QEmainPlugin, "vcrelax.qe.out"),
    (QEmainPlugin, "band.qe.out"),
    (QEChargeDensityPlugin, "scf.qe.charge"),
    (CIFPlugin, "FeS_troilite.cif"),
    (GULPOutPlugin, "reaxf_noopt.gulp.out"),
    (GULPOutPlugin, "reaxf_opt.gulp.out")
])
def test_plugins(testplugin, filename):
    plugins.unload_all_plugins()
    plugins.load_plugin_classes([Encode_NDArray])
    plugins.load_plugin_classes([testplugin], "parsers")

    inpath = os.path.join(file_folder, filename)
    output = plugins.parse(inpath)
    # print(json.dumps(output, indent=2, default=plugins.encode))
    outpath = os.path.join(file_folder, filename + ".json")

    # if "scf.qe.charge" in filename:
    #     with open(outpath, "w") as f:
    #         json.dump(output, f, indent=2, default=plugins.encode)

    expected = ejson.to_dict(outpath)

    assert edict.diff(output, expected, np_allclose=True) == {}


def test_opt_scflog_merge():
    plugins.load_plugin_classes([CrystalOutputPlugin, CrystalSCFLogPlugin], "parsers")
    opt = plugins.parse(os.path.join(file_folder, "scf_and_opt.crystal.out"))
    scflog = plugins.parse(os.path.join(file_folder, "scf_and_opt.crystal.scflog"))
    opt_all = edict.merge([opt, scflog], list_of_dicts=True)
    expected = ejson.to_dict(os.path.join(file_folder, "opt_merge_scflog.crystal.out.json"))
    assert edict.diff(opt_all, expected, np_allclose=True) == {}

@pytest.mark.parametrize(
    "fname,sname", [
        ("scf_only.crystal.out.json", "crystal_out.json"),
        ("scf_and_opt.crystal.out.json", "crystal_out.json"),
        ("scf_and_opt_slab.crystal.out.json", "crystal_out.json"),
        ("opt_merge_scflog.crystal.out.json", "crystal_out.json"),
        ("crystal.band.f25.json", "crystal_band.json"),
        ("scf.qe.out.json", "qe_out.json"),
        ("scf_with_fermi.qe.out.json", "qe_out.json"),
        ("vcrelax.qe.out.json", "qe_out.json"),
        ("band.qe.out.json", "qe_out.json"),
        ("FeS_troilite.cif.json", "cif.json")
    ]
)
def test_against_schema(fname, sname):

    with open(os.path.join(file_folder, fname)) as f:
        injson = json.load(f)
    with open(os.path.join(schema_folder, sname)) as f:
        schema = json.load(f)

    jsonschema.validate(injson, schema)


def test_get_schema():
    assert _get_all_schema_name() == ['cif', 'crystal_band', 'crystal_out', 'qe_out']


