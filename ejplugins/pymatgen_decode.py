import warnings

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
    >>> Encode_Pymatgen().to_str(struct)
    'Structure(Fe2, Comp=Fe, SpGrp=229)'

    """
    plugin_name = 'pymatgen.Structure'
    plugin_descript = 'encode/decode pymatgen.Structure'
    objclass = pym.Structure
    dict_signature = ['@class', '@module', 'lattice', 'sites']
    allow_other_keys = True

    def to_json(self, obj):
        return obj.as_dict()

    def from_json(self, obj):

        modname = obj["@module"]
        classname = obj["@class"]
        mod = __import__(modname, globals(), locals(), [classname], 0)
        if hasattr(mod, classname):
            cls_ = getattr(mod, classname)
            data = {k: v for k, v in obj.items()
                    if k not in ["@module", "@class"]}
            if hasattr(cls_, "from_dict"):
                return cls_.from_dict(data)
            else:
                raise ValueError("the class {0}.{1} does not have a from_dict method".format(modname, classname))
        else:
            raise ValueError("the module {0} does not have the required class {1}".format(modname, classname))

    def to_str(self, obj):
        name = obj.__class__.__name__
        return "{0}({1}, Comp={2}, SpGrp={3})".format(name,
            obj.formula.replace(" ", ""), obj.composition.reduced_formula, obj.get_space_group_info()[1])


class Encode_SymmOp(object):
    """

    Examples
    --------
    >>> symop = pym.SymmOp.from_rotation_and_translation()
    >>> symop_dict = Encode_SymmOp().to_json(symop)
    >>> print(sorted([str(s) for s in symop_dict.keys()]))
    ['@class', '@module', 'matrix', 'tolerance']
    >>> isinstance(Encode_SymmOp().from_json(symop_dict), pym.SymmOp)
    True
    >>> Encode_SymmOp().to_str(symop)
    'SymOp(R=[[1. 0. 0.],[0. 1. 0.],[0. 0. 1.]], T=[0. 0. 0.])'

    """
    plugin_name = 'pymatgen.SymmOp'
    plugin_descript = 'encode/decode pymatgen.SymmOp'
    objclass = pym.SymmOp
    dict_signature = ['@class', '@module', 'matrix']
    allow_other_keys = True

    def to_json(self, obj):
        return obj.as_dict()

    def from_json(self, obj):

        modname = obj["@module"]
        classname = obj["@class"]
        mod = __import__(modname, globals(), locals(), [classname], 0)
        if hasattr(mod, classname):
            cls_ = getattr(mod, classname)
            data = {k: v for k, v in obj.items()
                    if k not in ["@module", "@class"]}
            if hasattr(cls_, "from_dict"):
                return cls_.from_dict(data)
            else:
                raise ValueError("the class {0}.{1} does not have a from_dict method".format(modname, classname))
        else:
            raise ValueError("the module {0} does not have the required class {1}".format(modname, classname))

    def to_str(self, obj):
        outstr = str(obj).replace("\n", " ").replace("Rot: ", "R=").replace(" tau ", ", T=").replace("]  [", "],[")
        outstr = outstr.replace("  ", " ").replace("[ ", "[")
        return "SymOp({})".format(outstr)

