import glob
import json
import os
import re
from decimal import Decimal

import numpy as np
import jsonschema
from jsonextended import ejson, plugins
from jsonextended.encoders.ndarray import Encode_NDArray

try:
    basestring
except NameError:
    basestring = str


def split_numbers(string, as_decimal=False):
    """ get a list of numbers from a string (even with no spacing)

    Parameters
    ----------
    string: str
    as_decimal: bool
        if True return floats as decimal.Decimal objects

    Returns
    --------
    float_list: List

    Examples
    --------
    >>> split_numbers("1")
    [1.0]

    >>> split_numbers("1 2")
    [1.0, 2.0]

    >>> split_numbers("1.1 2.3")
    [1.1, 2.3]

    >>> split_numbers("1e-3")
    [0.001]

    >>> split_numbers("-1-2")
    [-1.0, -2.0]

    >>> split_numbers("1e-3-2")
    [0.001, -2.0]

    """
    _match_number = re.compile('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[+-]?\ *[0-9]+)?')
    string = string.replace(" .", " 0.")
    string = string.replace("-.", "-0.")
    return [Decimal(s) if as_decimal else float(s) for s in re.findall(_match_number, string)]


# CODATA 2014 taken from
# http://arxiv.org/pdf/1507.07956.pdf
codata = {
    ("Rydberg", "eV"): 13.605693012183622,
    ("Hartree", "eV"): 27.21138602,
    ("Bohr", "Angstrom"): 0.5291772105638411,
    ("Debye", "C.m"): 3.33564E-30
}

_test_file_folder = os.path.join(os.path.dirname(__file__), "test_files")


def load_test_file(fname):
    """ load a test json file

    Parameters
    ----------
    fname: str

    Returns
    -------

    """
    if not fname.endswith(".json"):
        fname += ".json"
    with plugins.plugins_context([Encode_NDArray]):
        data = ejson.to_dict(os.path.join(_test_file_folder, fname))
    return data


_schema_folder = os.path.join(os.path.dirname(__file__), "schema")


def _get_all_schema_name():
    """return a list off all the schame names"""
    schemas = []
    for path in glob.glob(os.path.join(_schema_folder, "*.json")):
        schemas.append(os.path.splitext(os.path.basename(path))[0])
    return schemas


def nddim_validator(validator, value, instance, schema):
    """ validators for n-dimensional data shape

    Parameters
    ----------
    validator
    value
    instance
    schema

    Returns
    -------

    """
    dim = len(np.asarray(instance).shape)
    if value != dim:
        yield jsonschema.ValidationError(
            "object is of dimension {} not {}".format(dim, value))


def ndtype_validator(validator, value, instance, schema):
    """ validators for n-dimensional data dtype

    Parameters
    ----------
    validator
    value
    instance
    schema

    Returns
    -------

    """
    try:
        np.asarray(instance, dtype=value)
    except (ValueError, TypeError):
        yield jsonschema.ValidationError(
            "object cannot be coerced to type %s" % value)


def validate_against_schema(data, schema):
    """ validate json-type data against a schema

    Parameters
    ----------
    data: dict
    schema: str or dict
        either a reference to a schema in ejplugins.schema (not including extension) or a full schema
    """
    if isinstance(schema, basestring):
        schemas = _get_all_schema_name()
        if schema not in schemas:
            raise ValueError("{0} not in available schema: {1}\nschema folder: {2}".format(schema, schemas, _schema_folder))

        with open(os.path.join(_schema_folder, "{}.json".format(schema))) as f:
            schema = json.load(f)

    # add some validators for n-dimensional data
    # an example of custom validators: https://lat.sk/2017/03/custom-json-schema-type-validator-format-python/
    validator = jsonschema.validators.extend(
        jsonschema.Draft4Validator,
        validators={
            'nddim': nddim_validator,
            'ndtype': ndtype_validator})

    # by default, only validates lists
    validator(schema, types={"array": (list, tuple)}).validate(data)


_SYM2ANUM = dict(Ac=89, Ag=47, Al=13, Am=95, Ar=18, As=33, At=85, Au=79, B=5, Ba=56, Be=4, Bh=107, Bi=83, Bk=97, Br=35,
                 C=6, Ca=20, Cd=48, Ce=58, Cf=98, Cl=17, Cm=96, Cn=112, Co=27, Cr=24, Cs=55, Cu=29, Db=105, Ds=110,
                 Dy=66, Er=68, Es=99, Eu=63, F=9, Fe=26, Fl=114, Fm=100, Fr=87, Ga=31, Gd=64, Ge=32, H=1, He=2, Hf=72,
                 Hg=80, Ho=67, Hs=108, I=53, In=49, Ir=77, K=19, Kr=36, La=57, Li=3, Lr=103, Lu=71, Lv=116, Md=101,
                 Mg=12, Mn=25, Mo=42, Mt=109, N=7, Na=11, Nb=41, Nd=60, Ne=10, Ni=28, No=102, Np=93, O=8, Os=76, P=15,
                 Pa=91, Pb=82, Pd=46, Pm=61, Po=84, Pr=59, Pt=78, Pu=94, Ra=88, Rb=37, Re=75, Rf=104, Rg=111, Rh=45,
                 Rn=86, Ru=44, S=16, Sb=51, Sc=21, Se=34, Sg=106, Si=14, Sm=62, Sn=50, Sr=38, Ta=73, Tb=65, Tc=43,
                 Te=52, Th=90, Ti=22, Tl=81, Tm=69, U=92, Uuh=117, Uup=115, Uut=113, V=23, W=74, Xe=54, Y=39, Yb=70,
                 Zn=30, Zr=40)


def symbol2anum(symbol):
    """
    Parameters
    ----------
    symbol: str
    """
    return _SYM2ANUM.get(symbol, None)
