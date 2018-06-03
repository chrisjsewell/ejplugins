import json
import os

import jsonschema
import pytest
from ejplugins.cif import CIFPlugin
from ejplugins.crystal import (CrystalOutputPlugin, CrystalSCFLogPlugin,
                               BANDPlugin, ECH3CubePlugin, ECH3OutPlugin, CrystalDOSPlugin)
from ejplugins.gcube import GaussianCube, ejdict_to_gcube
from ejplugins.gulp import GULPOutPlugin
from ejplugins.qespresso import QEmainPlugin, QEChargeDensityPlugin, QELowdinPlugin, QEdosPlugin
from ejplugins.utils import _get_all_schema_name
from jsonextended import plugins, edict, ejson
from jsonextended.encoders.ndarray import Encode_NDArray

file_folder = os.path.join(os.path.dirname(__file__), "test_files")
schema_folder = os.path.join(os.path.dirname(__file__), "schema")


@pytest.mark.parametrize("testplugin,filename", [
    (CrystalOutputPlugin, "scf_only.crystal.out"),
    (CrystalOutputPlugin, "scf_and_opt.crystal.out"),
    (CrystalOutputPlugin, "scf_and_opt_slab.crystal.out"),
    (CrystalOutputPlugin, "crystal17_spin_opt.crystal.out"),
    (CrystalSCFLogPlugin, "scf_and_opt.crystal.scflog"),
    # (DOSSPlugin, "crystal.doss.f25"),
    (CrystalDOSPlugin, "crystal.doss.f25"),
    (CrystalDOSPlugin, "crystal_pdos_spin.doss.f25"),
    (BANDPlugin, "crystal.band.f25"),
    (ECH3CubePlugin, "crystal.ech3_dat.prop3d"),
    (ECH3OutPlugin, "multi_atom_bulk.ech3.out"),
    (ECH3CubePlugin, "multi_atom_bulk.ech3_dat.prop3d"),
    (ECH3OutPlugin, "crystal.ech3.out"),
    (QEmainPlugin, "scf.qe.out"),
    (QEmainPlugin, "scf_with_fermi.qe.out"),
    (QEmainPlugin, "vcrelax.qe.out"),
    (QEmainPlugin, "band.qe.out"),
    (QELowdinPlugin, "relax.qe.pdos.out"),
    (QEChargeDensityPlugin, "scf.qe.charge"),
    (QEChargeDensityPlugin, "surface.qe.charge"),
    (QEdosPlugin, "test.qe.dos"),
    (CIFPlugin, "FeS_troilite.cif"),
    (GULPOutPlugin, "reaxf_noopt.gulp.out"),
    (GULPOutPlugin, "reaxf_opt.gulp.out"),
    (GaussianCube, "CHG_H2O_small.cube"),
])
def test_plugins(testplugin, filename):
    plugins.unload_all_plugins()
    plugins.load_plugin_classes([Encode_NDArray])
    plugins.load_plugin_classes([testplugin], "parsers")

    inpath = os.path.join(file_folder, filename)
    output = plugins.parse(inpath)
    # print(json.dumps(output, indent=2, default=plugins.encode))
    outpath = os.path.join(file_folder, filename + ".json")

    # if ".crystal.out" in filename:
    #     with open(outpath, "w") as f:
    #         json.dump(output, f, indent=1, default=plugins.encode)
    print("reading expected")
    expected = ejson.to_dict(outpath)

    assert edict.diff(output, expected, np_allclose=True) == {}


def test_opt_scflog_merge():
    plugins.load_plugin_classes([CrystalOutputPlugin, CrystalSCFLogPlugin], "parsers")
    opt = plugins.parse(os.path.join(file_folder, "scf_and_opt.crystal.out"))
    scflog = plugins.parse(os.path.join(file_folder, "scf_and_opt.crystal.scflog"))
    opt_all = edict.merge([opt, scflog], list_of_dicts=True)
    # with open(os.path.join(file_folder, "opt_merge_scflog.crystal.out.json"), "w") as f:
    #     json.dump(opt_all, f, indent=1, default=plugins.encode)

    expected = ejson.to_dict(os.path.join(file_folder, "opt_merge_scflog.crystal.out.json"))
    assert edict.diff(opt_all, expected, np_allclose=True) == {}


@pytest.mark.parametrize(
    "fname,sname", [
        ("scf_only.crystal.out.json", "crystal_out.json"),
        ("scf_and_opt.crystal.out.json", "crystal_out.json"),
        ("scf_and_opt_slab.crystal.out.json", "crystal_out.json"),
        ("opt_merge_scflog.crystal.out.json", "crystal_out.json"),
        ("crystal.band.f25.json", "crystal_band.json"),
        ("crystal.doss.f25.json", "crystal_doss.json"),
        ("crystal_pdos_spin.doss.f25.json", "crystal_doss.json"),
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
    assert _get_all_schema_name() == ['cif', 'crystal_band', 'crystal_doss', 'crystal_out', 'edensity', 'qe_out']


def test_ejdict_to_gcube():
    plugins.unload_all_plugins()
    plugins.load_plugin_classes([Encode_NDArray])
    plugins.load_plugin_classes([GaussianCube], "parsers")

    inpath = os.path.join(file_folder, "CHG_H2O_small.cube")
    data = plugins.parse(inpath)

    gcube_path = ejdict_to_gcube(data)

    assert gcube_path.file_content[0:10] == ['H2O Density Calculation.', 'MP2 Total Density',
                                             '     3  -7.596700  -7.596700  -7.596700',
                                             '    50   0.303868   0.000000   0.000000',
                                             '    50   0.000000   0.303868   0.000000',
                                             '    50   0.000000   0.000000   0.303868',
                                             '     1   1.000000  -0.000006   1.430423   0.984133',
                                             '     8   8.000000  -0.000006  -0.000006  -0.123023',
                                             '     1   1.000000  -0.000006  -1.430435   0.984133',
                                        ' 1.19286E-11  1.53229E-11  1.96312E-11  2.51042E-11  3.20622E-11  4.09056E-11']

    assert gcube_path.file_content[-1] == " 6.77960E-11  5.61461E-11"


def test_ejdict_to_gcube_crystal():
    plugins.unload_all_plugins()
    plugins.load_plugin_classes([Encode_NDArray])
    plugins.load_plugin_classes([ECH3CubePlugin, ECH3OutPlugin], "parsers")

    inpath1 = os.path.join(file_folder, "multi_atom_bulk.ech3_dat.prop3d")
    data = plugins.parse(inpath1)
    inpath2 = os.path.join(file_folder, "multi_atom_bulk.ech3.out")
    adata = plugins.parse(inpath2)

    # expected_path = os.path.join(file_folder, "multi_atom_bulk.ech3.cube")
    # expected_data = plugins.parse(expected_path)

    gcube_path = ejdict_to_gcube(data, density=0, adata=adata)
    #new_data = plugins.parse(gcube_path)

    assert gcube_path.file_content[0:10] == ['CRYSTAL', 'charge',
                                             '    30   5.468942  -6.513792 -10.875708',
                                             '    65   0.174544   0.100212  -0.000685',
                                             '    65  -0.174544   0.100212   0.000685',
                                             '   120  -0.091149   0.000000   0.181262',
                                             '    16  16.021000  -1.131871   4.830045  -2.704532',
                                             '    16  16.021000   6.598999   4.830045  -8.171509',
                                             '    16  16.021000   1.131832  -4.830046   2.704526',
                                             '    16  16.021000  -6.599038  -4.830046   8.171504']

    assert gcube_path.file_content[-1] == " 2.60932E-02  2.52211E-02  2.43080E-02  2.35170E-02  2.29831E-02  2.27950E-02"




