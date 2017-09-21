# -*- coding: utf-8 -*-
import ast


# TODO add lammps unit tests
class LAMMPSSysDumpPlugin(object):
    """ parser plugin for jsonextended
    """
    plugin_name = 'lammps_sys_out'
    plugin_descript = 'read LAMMPS final system output (with real units)'
    file_regex = '*.lammps.sys.real.dump'

    # to convert to per atom
    _MOLES = 6.022140857e23
    _lammps_ename_map = {
        'eb': 'Bond',
        'ea': 'Coordination (over)',  ## or technically over/under, or what is Atom???
        'elp': 'Lone-Pair',
        'ev': 'Valence Angle',
        'epen': 'Double-Bond Valence Angle Penalty',
        'ecoa': 'Valence Angle Conjugation',
        'ehb': 'Hydrogen Bond',
        'et': 'Torsion',
        'eco': 'Conjugation',
        'ew': 'van der Waals',
        'ep': 'Coulomb',
        'eqeq': 'Charge Equilibration'}

    def read_file(self, f, **kwargs):
        line = f.readline().strip()
        while line.startswith('#'):
            line = f.readline().strip()
        headers = line.split()
        line = f.readline().strip()
        while line.startswith('#'):
            line = f.readline().strip()
        values = line.split()
        data = {}
        for name, value in zip(headers, values):
            if name == 'natoms':
                data['natoms'] = int(value)
            elif name == 'temp':
                data['temperature'] = {'magnitude': float(value), 'units': 'kelvin'}
            elif name == 'vol':
                data['volume'] = {'magnitude': float(value), 'units': 'angstrom^3'}
            elif name == 'press':
                data['pressure'] = {'magnitude': float(value), 'units': 'atmospheres'}
            elif name == 'teng':
                data['lattice_energy'] = {}
                data['lattice_energy']['primitive'] = {'magnitude': float(value) / self._MOLES, 'units': 'kcal'}
            elif name == 'peng':
                data['potential_energy'] = {'magnitude': float(value) / self._MOLES, 'units': 'kcal'}
            elif name == 'keng':
                data['kinetic_energy'] = {'magnitude': float(value) / self._MOLES, 'units': 'kcal'}
            elif name in list(self._lammps_ename_map.keys()):
                if not 'reaxff_ebreakdown' in data:
                    data['reaxff_ebreakdown'] = {}
                sdata = {'magnitude': float(value) / self._MOLES, 'units': 'kcal'}
                data['reaxff_ebreakdown'][self._lammps_ename_map[name]] = sdata
        return data


class LAMMPSAtomDumpPlugin(object):
    """ parser plugin for jsonextended
    """
    plugin_name = 'lammps_atom_out'
    plugin_descript = 'read LAMMPS final system output (with real units)'
    file_regex = '*.lammps.atom.real.dump'

    @staticmethod
    def tryeval(val):
        try:
            val = ast.literal_eval(val)
        except ValueError:
            pass
        return val

    def read_file(self, f, **kwargs):
        steps = []
        step = None

        line = f.readline()
        while line:
            fields = line.strip().split()
            if not len(fields) >= 2:
                line = f.readline()
                continue
            if not fields[0] == 'ITEM:':
                line = f.readline()
                continue
            if fields[1] == 'TIMESTEP':
                # new step
                if step is not None:
                    steps.append(step)
                step = {}
                step['timestep'] = int(f.readline().strip().split()[0])
            elif fields[1] == 'BOX':
                # add 0. to ensure tilt factors are set
                xfields = [float(s) for s in f.readline().strip().split()] + [0.]
                xlo_bound, xhi_bound, xy = xfields[0:3]
                yfields = [float(s) for s in f.readline().strip().split()] + [0.]
                ylo_bound, yhi_bound, xz = yfields[0:3]
                zfields = [float(s) for s in f.readline().strip().split()] + [0.]
                zlo_bound, zhi_bound, yz = zfields[0:3]

                xlo, xhi = xlo_bound - min(0.0, xy, xz, xy + xz), xhi_bound - max(0.0, xy, xz, xy + xz)
                ylo, yhi = ylo_bound - min(0.0, yz), yhi_bound - max(0.0, yz)
                zlo, zhi = zlo_bound, zhi_bound
                step['cell_vectors'] = {'origin': (xlo, ylo, zlo),
                                        'a': (xhi - xlo, 0., 0.),
                                        'b': (xy, yhi - ylo, 0.),
                                        'c': (xz, yz, zhi - zlo)}
            elif fields[1] == 'ATOMS':
                vhead = fields[2:]
                vval = []
                while line:
                    line = f.readline()
                    fields = line.strip().split()
                    if not fields:
                        break
                    if fields[0] == 'ITEM':
                        break
                    vval.append([self.tryeval(field) for field in fields])

                step['atoms'] = {v[0]: v[1:] for v in zip(vhead, *vval)}

                continue

            line = f.readline()
        if step and step is not None:
            steps.append(step)
        return {'configs': steps}

