import warnings
# from pymatgen.core.structure import Structure
# from pymatgen.core.lattice import Lattice
with warnings.catch_warnings(record=True):
    warnings.filterwarnings("ignore", category=ImportWarning)
    import pymatgen as pym


class Encode_Pymatgen(object):
    """

    Examples
    --------
    >>> struct = pym.Structure.from_spacegroup(229, pym.Lattice.from_parameters(1,1,1,90,90,90), ['Fe'], [[0, 0, 0]])
    >>> struct_dict = Encode_Pymatgen().to_json(struct)
    >>> print(sorted([str(s) for s in struct_dict.keys()]))
    ['@class', '@module', 'charge', 'lattice', 'sites']
    >>> isinstance(Encode_Pymatgen().from_json(struct_dict), pym.Structure)
    True

    """
    plugin_name = 'pymatgen.Structure'
    plugin_descript = 'encode/decode pymatgen.Structure'
    objclass = pym.Structure
    dict_signature = ['@class', '@module', 'lattice', 'sites']

    def to_json(self, obj):
        return obj.as_dict()

    def from_json(self, obj):
        return pym.Structure.from_dict(obj)

