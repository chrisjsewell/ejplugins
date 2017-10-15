import json
import re
import os
import glob
from decimal import Decimal

import jsonschema
from jsonextended import ejson


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
    return [Decimal(s) if as_decimal else float(s) for s in re.findall(_match_number, string)]

# CODATA 2014 taken from
# http://arxiv.org/pdf/1507.07956.pdf
codata = {
    ("Rydberg", "eV"): 13.605693012183622,
    ("Hartree", "eV"): 27.21138602,
    ("Bohr", "Angstrom"): 0.5291772105638411,
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
    return ejson.to_dict(os.path.join(_test_file_folder, fname))

_schema_folder = os.path.join(os.path.dirname(__file__), "schema")


def _get_all_schema_name():
    """return a list off all the schame names"""
    schemas = []
    for path in glob.glob(os.path.join(_schema_folder, "*.json")):
        schemas.append(os.path.splitext(os.path.basename(path))[0])
    return schemas


def validate_against_schema(data, sname):
    """get a validation schema by name (not including extension) and validate the data against it"""
    schemas = _get_all_schema_name()
    if sname not in schemas:
        raise ValueError("{0} not in available schema: {1}\nschema folder: {2}".format(sname, schemas, _schema_folder))

    with open(os.path.join(_schema_folder, "{}.json".format(sname))) as f:
        schema = json.load(f)

    # by default, only validates lists
    validator = jsonschema.Draft4Validator(schema, types={"array": (list, tuple)})
    validator.validate(data)

