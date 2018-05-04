import numpy as np
from jsonextended import units, mockpath, edict
from ejplugins.utils import split_numbers, codata, validate_against_schema

try:
    import pathlib
except ImportError:
    import pathlib2 as pathlib


class GaussianCube(object):
    """ parser plugin for jsonextended

    as specified at http://h5cube-spec.readthedocs.io/en/latest/cubeformat.html

    """
    plugin_name = 'gaussian_cube'
    plugin_descript = 'read gaussian cube charge/spin density file'
    file_regex = '*.cube'

    def read_file(self, f, **kwargs):
        comments1 = f.readline().strip()
        comments2 = f.readline().strip()
        line = f.readline().strip()
        inputs = split_numbers(line)
        if len(inputs) > 4 and inputs[4] != 1:
            # TODO implement NVAL != 1
            raise NotImplementedError("not yet implemented NVAL != 1: {0}".format(line))
        natoms = inputs[0]
        centre = -1 * np.array(inputs[1:4]) * codata[("Bohr", "Angstrom")]
        if natoms < 0:
            # TODO implement DSET_IDS
            raise NotImplementedError("not yet implemented DSET_IDS")
        an, ax, ay, az = split_numbers(f.readline().strip())
        bn, bx, by, bz = split_numbers(f.readline().strip())
        cn, cx, cy, cz = split_numbers(f.readline().strip())

        if an <= 0 or bn <= 0 or cn <= 0:
            raise ValueError("an, bn and cn must be positive integers")

        avec = [a * an * codata[("Bohr", "Angstrom")] for a in [ax, ay, az]]
        bvec = [b * bn * codata[("Bohr", "Angstrom")] for b in [bx, by, bz]]
        cvec = [c * cn * codata[("Bohr", "Angstrom")] for c in [cx, cy, cz]]
        #centre = 0.5 * (np.array(avec) + np.array(bvec) + np.array(cvec))

        atomic_numbers = []
        nuclear_charges = []
        ccoords = []
        for _ in range(int(natoms)):
            anum, ncharge, x, y, z = split_numbers(f.readline().strip())
            atomic_numbers.append(int(anum))
            nuclear_charges.append(ncharge)
            ccoord = (np.asarray([x, y, z]) * codata[("Bohr", "Angstrom")]) + centre
            ccoords.append(ccoord.tolist())

        values = []
        line = f.readline().strip()
        while line:
            values += line.split()
            line = f.readline().strip()

        return {
            "title": comments1,
            #"na": int(an), "nb": int(bn), "nc": int(cn),
            "cell_vectors": {
                "a": {"units": "angstrom", "magnitude": avec},
                "b": {"units": "angstrom", "magnitude": bvec},
                "c": {"units": "angstrom", "magnitude": cvec}
            },
            #"centre": [0, 0, 0],
            "densities": [{
                "type": comments2,
                "magnitude": np.array(values, dtype=float).reshape((int(an), int(bn), int(cn)))
            }],
            "atoms": {"ccoords": {"units": "angstrom",
                                  "magnitude": ccoords},
                      "nuclear_charge": nuclear_charges,
                      "atomic_number": atomic_numbers}
        }


def ejdict_to_gcube(data, fpath=None, density=0,
                    include_atoms=True, adata=None, cell_tol=1E-3):
    """

    Parameters
    ----------
    data: dict
    fpath: str or None
        output file path or, if None, write to MockPath
    density: int
        take density from data["densities"][density]
    include_atoms: bool
        include atoms in gaussian cube
    adata: dict or None
        separate atom data (for instance for Crystal output)
    cell_tol: float or None
        if not None, raise and error if the data and adata cell vectors are not within this tolerance

    Returns
    -------
    fpath: pathlib.Path or jsonextended.mockpath.MockPath

    """
    if include_atoms and adata is not None:
        if cell_tol:
            cdiff = edict.diff(data["cell_vectors"], adata["cell_vectors"],
                               np_allclose=True, rtol=cell_tol, atol=cell_tol)
            if cdiff:
                raise ValueError("data and adata have different cell vectors: {}".format(cdiff))
        data["atoms"] = adata["atoms"]

    validate_against_schema(data, "edensity")

    if fpath is None:
        fpath = mockpath.MockPath("test.cube", is_file=True)
    else:
        fpath = pathlib.Path(fpath)

    with fpath.open("w") as f:
        data = units.combine_quantities(data)
        data = units.apply_unitschema(data, {"a": "angstrom", "b": "angstrom", "c": "angstrom", "ccoords": "angstrom"},
                                      as_quantity=False)
        natoms = 0 if "atoms" not in data or not include_atoms else len(data["atoms"]["ccoords"])

        avec = np.asarray(data["cell_vectors"]["a"]) / codata[("Bohr", "Angstrom")]
        bvec = np.asarray(data["cell_vectors"]["b"]) / codata[("Bohr", "Angstrom")]
        cvec = np.asarray(data["cell_vectors"]["c"]) / codata[("Bohr", "Angstrom")]
        centre = 0.5 * (avec + bvec + cvec)
        centre_offset = -1*centre

        f.write(data["title"] + "\n")
        f.write(data["densities"][density]["type"] + "\n")
        dense = np.asarray(data["densities"][density]["magnitude"])
        na, nb, nc = dense.shape
        f.write("{0:6d} {1:10.6f} {2:10.6f} {3:10.6f}\n".format(natoms, *centre_offset.tolist()))
        f.write("{0:6d} {1:10.6f} {2:10.6f} {3:10.6f}\n".format(na, *(avec/na).tolist()))
        f.write("{0:6d} {1:10.6f} {2:10.6f} {3:10.6f}\n".format(nb, *(bvec/nb).tolist()))
        f.write("{0:6d} {1:10.6f} {2:10.6f} {3:10.6f}\n".format(nc, *(cvec/nc).tolist()))

        if data.get("atoms", False):
            for i, c in enumerate(data["atoms"]["ccoords"]):
                atomic_number = data["atoms"]["atomic_number"][i]
                nuclear_charge = data["atoms"]["nuclear_charge"][i]
                ccoord = (np.array(c) / codata[("Bohr", "Angstrom")]) + centre_offset

                f.write("{0:6d} {1:10.6f} {2:10.6f} {3:10.6f} {4:10.6f}\n".format(atomic_number, nuclear_charge,
                                                                                  *ccoord.tolist()))
        dense = dense.flatten().tolist()
        dlength = len(dense)
        output = []
        for i in range(int(dlength/6.)+1):
            if dlength > i*6 + 6:
                output.append("{0:12.5E} {1:12.5E} {2:12.5E} {3:12.5E} {4:12.5E} {5:12.5E}".format(*dense[i*6: i*6+6]))
            else:
                output.append(" ".join(["{0:12.5E}".format(v) for v in dense[i*6: dlength]]))

        f.write("\n".join(output))

    return fpath
