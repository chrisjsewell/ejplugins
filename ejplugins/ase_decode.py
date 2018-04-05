import warnings
from pymatgen.io.ase import AseAtomsAdaptor
# from pymatgen.core.structure import Structure
with warnings.catch_warnings(record=True):
    warnings.filterwarnings("ignore", category=ImportWarning)
    import pymatgen as pym
from ase import Atoms


class Encode_ASE(object):
    """

    Examples
    --------
    >>> from ase.build import bulk
    >>> struct_ase = bulk("Fe", a=1, orthorhombic=True)
    >>> struct_dict = Encode_ASE().to_json(struct_ase)
    >>> print(sorted([str(s) for s in struct_dict.keys()]))
    ['@ase', '@class', '@module', 'charge', 'lattice', 'sites']
    >>> isinstance(Encode_ASE().from_json(struct_dict), Atoms)
    True

    """
    plugin_name = 'ase.Atoms'
    plugin_descript = 'encode/decode ase.Atoms'
    objclass = Atoms
    dict_signature = ['@ase', '@class', '@module', 'lattice', 'sites']

    def to_json(self, obj):
        struct = AseAtomsAdaptor.get_structure(obj)
        dct = struct.as_dict()
        dct["@ase"] = True
        return dct

    def from_json(self, obj):
        struct = pym.Structure.from_dict(obj)
        return AseAtomsAdaptor.get_atoms(struct)


