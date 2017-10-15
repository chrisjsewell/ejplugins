from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice


class Encode_Pymatgen(object):
    """

    Examples
    --------
    >>> struct = Structure.from_spacegroup(229, Lattice.from_parameters(1,1,1,90,90,90), ['Fe'], [[0, 0, 0]])
    >>> struct_dict = Encode_Pymatgen().to_json(struct)
    >>> print(sorted(struct_dict.keys()))
    ['@class', '@module', 'lattice', 'sites']
    >>> isinstance(Encode_Pymatgen().from_json(struct_dict), Structure)
    True

    """
    plugin_name = 'pymatgen.Structure'
    plugin_descript = 'encode/decode pymatgen.Structure'
    objclass = Structure
    dict_signature = ['@class', '@module', 'lattice', 'sites']

    def to_json(self, obj):
        return obj.as_dict()

    def from_json(self, obj):
        return Structure.from_dict(obj)

