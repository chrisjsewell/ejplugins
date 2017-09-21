from collections import OrderedDict
import warnings
with warnings.catch_warnings(record=True):
    warnings.filterwarnings("ignore", category=ImportWarning)
    from pymatgen.io.cif import CifParser as cifparser

class CIFPlugin(object):
    """ cif parser plugin for jsonextended
    """
    plugin_name = 'crytallographic_info_file'
    plugin_descript = 'read crytallographic information files, assumes single structure per file'
    file_regex = '*.cif'

    def read_file(self, file_obj, **kwargs):
        cif = cifparser.from_string(file_obj.read())
        dic = cif.as_dict()
        dic = dic[list(dic.keys())[0]]
        # remove underscores for key names
        new_dic = OrderedDict()
        for key in dic.keys():
            newkey = key[1:] if key.startswith('_') else key
            new_dic[newkey] = dic[key]
        new_dic['structures'] = [s.as_dict() for s in cif.get_structures()]

        return new_dic
