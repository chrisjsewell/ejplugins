# -*- coding: utf-8 -*-
import copy


# TODO split into sections " Input for Configuration" and "Output for configuration"
class GULPOutPlugin(object):
    """ parser plugin for jsonextended
    """
    plugin_name = 'gulp_out'
    plugin_descript = 'read GULP output'
    file_regex = '*.gulp.out'

    @staticmethod
    def reaxff_map(type):
        gulp_ename_map = {
            'bond': 'Bond',
            'bpen': 'Double-Bond Valence Angle Penalty',
            'lonepair': 'Lone-Pair',
            'over': 'Coordination (over)',
            'under': 'Coordination (under)',
            'val': 'Valence Angle',
            'pen': 'Double-Bond Valence Angle Penalty',
            'coa': 'Valence Angle Conjugation',
            'tors': 'Torsion',
            'conj': 'Conjugation',
            'hb': 'Hydrogen Bond',
            'vdw': 'van der Waals',
            'coulomb': 'Coulomb',
            'self': 'Charge Equilibration'}
        return gulp_ename_map.get(type, type)

    def read_file(self, f, **kwargs):
        data = {}
        line = f.readline()
        while line:
            fields = line.strip().split()

            if ' '.join(fields[:4]) == 'Total number atoms/shells =':
                data['natoms'] = int(fields[4])
            elif ' '.join(fields[:2]) == 'Formula =':
                data['formula'] = fields[2]

            elif ' '.join(fields[:4]) == '**** Optimisation achieved ****':
                data['optimised'] = True
            elif "No variables to optimise - single point performed" in line:
                data['optimised'] = True
            elif ' '.join(fields[:4]) == '**** Too many failed' and len(fields) > 5:
                if fields[6] == 'optimise':
                    data['optimised'] = False
            elif ' '.join(fields[:2]) == '**** Maximum' and len(fields) > 7:
                if ' '.join(fields[4:5] + [fields[8]]) == 'function calls reached':
                    data['optimised'] = False

            elif ' '.join(fields[:4]) == 'Total lattice energy =':
                units = ' '.join(fields[5:])
                if units == 'eV':
                    etype = 'initial'
                    if 'initial' in data:
                        if 'lattice_energy' in data['initial']:
                            if not 'final' in etype:
                                data['final'] = {}
                            etype = 'final'
                    else:
                        data['initial'] = {}
                    data[etype]['lattice_energy'] = {}
                    data[etype]['lattice_energy']['primitive'] = {'magnitude': float(fields[4]), 'units': units}

            elif ' '.join(fields[:4]) == 'Total lattice energy :':
                etype = 'initial'
                if 'initial' in data:
                    if 'lattice_energy' in data['initial']:
                        if not 'final' in etype:
                            data['final'] = {}
                        etype = 'final'
                else:
                    data['initial'] = {}
                data[etype]['lattice_energy'] = {}
                fields = f.readline().strip().split()
                assert fields[0] == 'Primitive' and fields[5] == 'eV'
                data[etype]['lattice_energy']['primitive'] = {'magnitude': float(fields[4]), 'units': fields[5]}
                fields = f.readline().strip().split()
                assert fields[0] == 'Non-primitive' and fields[5] == 'eV'
                data[etype]['lattice_energy']['conventional'] = {'magnitude': float(fields[4]), 'units': fields[5]}

            elif ' '.join(fields[:4]) == 'ReaxFF : Energy contributions:':
                assert 'reaxff_ebreakdown' not in data
                data['reaxff_ebreakdown'] = {}
                f.readline()
                line = f.readline()
                while line.strip():
                    fields = line.strip().split()
                    data['reaxff_ebreakdown'][self.reaxff_map(fields[0][2:-1])] = {
                        'magnitude': float(fields[2]), 'units': fields[3]}
                    line = f.readline()

            line = f.readline()

        if 'optimised' not in data:
            data['optimised'] = False
        if 'initial' in data and not 'final' in data and data['optimised']:
            data['final'] = copy.deepcopy(data['initial'])

        data["creator"] = {"program": "GULP"}

        return data
