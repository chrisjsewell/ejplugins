# -*- coding: utf-8 -*-
"""
adapted from: https://github.com/lmmentel/ase-espresso

Original work Copyright (C) 2013-2015 SUNCAT
Modified work Copyright 2015-2017 Lukasz Mentel

                    GNU GENERAL PUBLIC LICENSE
                       Version 3, 29 June 2007

"""
from fnmatch import fnmatch
import warnings

import numpy as np
# from jsonextended import edict

from ejplugins.utils import codata, split_numbers


def raise_error(msg, line, i, start_line):
    """ raise error when reading line

    Parameters
    ----------
    msg: str
    line: str
    i: int
    start_line: int

    Returns
    -------

    """
    raise IOError("{0}, for line #{1}: {2}".format(msg, i+start_line, line.strip()))


# TODO separate contributions
def read_forces(lines, start_line):
    """
    Read the forces from the PWSCF output file

    Parameters
    ----------
    lines: List

    Returns
    -------
    forces: dict
    units: str

    """
    units_conv = codata[("Rydberg", "eV")] / codata[("Bohr", "Angstrom")]

    forces = {}
    in_forces = False

    for i, line in enumerate(lines):
        line = line.strip()
        if line.startswith("Forces acting on atoms"):
            if in_forces:
                raise_error("found multiple forces data", line, i, start_line)
            if "Ry/au" not in line:
                raise_error("expecting forces to be in Ry/au", line, i, start_line)
            in_forces = True
            if lines[i+1].strip():
                raise_error("Expecting blank line after initial forces line", lines[i+1], i+1, start_line)
            forces = {"peratom": []}
            j = 2
            while fnmatch(lines[i+j].strip(), 'atom *'):
                aid, atype, fx, fy, fz = split_numbers(lines[i+j])
                forces["peratom"].append({"units": 'eV/angstrom',
                                          "magnitude": (fx * units_conv, fy * units_conv, fz * units_conv)})
                j += 1

    return forces if forces else None


def read_energies(lines, start_line):
    """
    Read the energies from the PWSCF output file

    The method read total energies (zero point) and smearing contributions
    (-TS) and returns either a tuple ``(energy, free_energy)`` or a list
    of such tuples

    Parameters
    ----------
    lines: List
    start_line: int

    Returns
    -------
    energy: dict
    units: str

    """
    total_energy = None
    contributions = {}
    fermi = None
    for i, line in enumerate(lines):
        line = line.strip()
        if fnmatch(line, "!*total energy*=*"):
            if "Ry" not in line:
                raise_error("Expecting energy to end Ry", line, i, start_line)
            if total_energy:
                raise_error("found multiple energies in section", line, i, start_line)
            total_energy = split_numbers(line)[0] * codata[("Rydberg", "eV")]

        if "The total energy is the sum of the following terms:" in line:
            if contributions:
                raise_error("already found energy contributions for this section", line, i, start_line)
            if lines[i+1].strip():
                raise_error("Expecting blank line after initial energy contributions line", lines[i+1], i+1, start_line)
            j = 2
            while lines[i+j].strip():
                if lines[i+j].strip().startswith("->"):
                    # These are sub component contributions, so ignore
                    j += 1
                    continue
                if not fnmatch(lines[i+j].strip(), "*=*Ry"):
                    raise_error("Expecting energy contribution to be in the form '*=*Ry'",
                                lines[i+j], i+j, start_line)
                name = lines[i+j].strip().split("=")[0]  # type: str
                name = name.replace("contribution", "").replace(" contrib.", "").replace("  ", " ").strip()
                value = split_numbers(lines[i+j])[0] * codata[("Rydberg", "eV")]
                if name in contributions:
                    raise_error("Multiple instances of energy contribution {} found in section".format(name),
                                lines[i+j], i+j, start_line)
                contributions[name] = {"units": "eV", "magnitude": value}
                j += 1

        if "the Fermi energy is" in line:
            if fermi:
                raise_error("already found a fermi energy for this section", line, i, start_line)
            if "ev" not in line:
                raise_error("expected fermi in eV", line, i, start_line)
            fermi = {"magnitude": split_numbers(line)[0], "units": "eV"}

    out = {}
    if total_energy:
        out["total"] = {"magnitude": total_energy, "units": "eV"}
    if contributions:
        out["components"] = contributions
    if fermi is not None:
        out["fermi"] = fermi

    return None if not out else out


def read_scf(lines, start_line):
    """ read scf data

    Parameters
    ----------
    lines: List

    Returns
    -------

    """
    scf = None
    scf_cyc = None
    in_scf = False
    last_cyc_num = None
    atomic_charges_peratom = None
    spin_density_peratom = None

    for i, line in enumerate(lines):
        line = line.strip()
        if "Self-consistent Calculation" in line:
            scf = []
            in_scf = True
            continue
        if "End of self-consistent calculation" in line:
            break
        if not in_scf:
            continue
        if line.startswith("iteration #"):
            # start new cycle
            if scf_cyc is not None:
                scf.append(scf_cyc)
            scf_cyc = {}
            # check we are adding them in sequential order
            cur_cyc_num = split_numbers(line)[0]
            if last_cyc_num is not None:
                if cur_cyc_num != last_cyc_num + 1:
                    raise_error("was expecting the SCF iteration number to be "
                                "{0}".format(int(last_cyc_num + 1)), line, i, start_line)
            last_cyc_num = cur_cyc_num

        if line.startswith("total energy"):
            if not line.endswith("Ry"):
                raise_error("was expecting total energy to be in Rydberg (end Ry)", line, i, start_line)
            if "energy" not in scf_cyc:
                scf_cyc["energy"] = {}
            scf_cyc["energy"]["total"] = {"units": "eV",
                                          "magnitude": split_numbers(line)[0] * codata[("Rydberg", "eV")]}

        # The total magnetization is the integral of the magnetization in the cell:
        #     MT=∫ (nup-ndown) d3 r
        #
        # The absolute magnetization is the integral of the absolute value of the magnetization in the cell:
        #     MA=∫ |nup-ndown| d3 r
        #
        # In a simple ferromagnetic material they should be equal (except possibly for an overall sign).
        # In simple antiferromagnets (like FeO) MT is zero and MA is twice the magnetization of each of the two atoms.
        if line.startswith("total magnetization"):
            scf_cyc["spin_density_total"] = split_numbers(line)[0]
        if line.startswith("absolute magnetization"):
            scf_cyc["spin_density_absolute"] = split_numbers(line)[0]

        # TODO make note that the sum of these spins/charges will not exactly equal final total_spin_density
        # see; https://www.mail-archive.com/pw_forum@pwscf.org/msg24862.html
        if "Magnetic moment per site:" in line:
            if atomic_charges_peratom or spin_density_peratom:
                raise_error("found multiple magnetic moments per site in same scf run", line, i, start_line)
            atomic_charges_peratom = []
            spin_density_peratom = []
            nxtline = i+1
            while lines[nxtline].strip():
                if not fnmatch(lines[nxtline].strip(), "atom:*charge:*magn:*constr:*"):
                    raise_error("was expecting magnetic moment fields; atom:, charge:, magn:, constr:",
                                lines[nxtline], nxtline, start_line)
                atom, charge, magn, constr = split_numbers(lines[nxtline])
                atomic_charges_peratom.append(charge)
                spin_density_peratom.append(magn)
                nxtline += 1

    # add last scf cycle
    if scf_cyc is not None:
        scf.append(scf_cyc)

    if atomic_charges_peratom is None:
        magmom = None
    else:
        magmom = {"atomic_charges_peratom": atomic_charges_peratom, "spin_density_peratom": spin_density_peratom}

    return scf, magmom


def read_cell(lines, start_line):
    """
    Read unit cell parameters from the PWSCF output file

    Parameters
    ----------
    lines: list

    Returns
    -------
    cell: dict

    """

    alat = None
    cell = {}
    # TODO you can get a section after the optimisation ("Begin final coordinates") but then a final scf calculation
    # followed by a recalculated final unit cell
    final_coords = False

    # additional unit cell information
    #bli_lines = [line for line in lines if 'bravais-lattice index' in line]
    #brav_latt_indices = [int(line.split('=')[1].strip()) for line in bli_lines]

    for i, line in enumerate(lines):
        line = line.strip()
        if 'lattice parameter (alat)' in line:
            if "a.u." not in line:
                raise_error("expecting alat in a.u.", line, i, start_line)
            if alat is not None and not final_coords:
                raise_error("found multiple alat values", line, i, start_line)
            alat = split_numbers(line)[0] * codata[("Bohr", "Angstrom")]

        # initial and final
        if 'crystal axes: (cart. coord. in units of alat)' in line:
            if cell and not final_coords:
                raise_error("found multiple cell coordinates", line, i, start_line)
            if not alat:
                raise_error("expecting to have found alat before", line, i, start_line)
            for key, newline in zip(["a", "b", "c"], lines[i + 1: i + 4]):
                x, y, z = split_numbers(newline.split('=')[1])
                cell[key] = {"units": "angstrom", "magnitude": (x*alat, y*alat, z*alat)}

        if "Begin final coordinates" in line:
            final_coords = True

        # optimisation step
        if fnmatch(line, "CELL_PARAMETERS*(*alat*=*)"):
            if cell:
                raise_error("found multiple cell coordinates", line, i, start_line)
            alat = split_numbers(line)[0] * codata[("Bohr", "Angstrom")]
            for key, newline in zip(["a", "b", "c"], lines[i + 1: i + 4]):
                x, y, z = split_numbers(newline)
                cell[key] = {"units": "angstrom", "magnitude": (x*alat, y*alat, z*alat)}

    return cell if cell else None


def read_atoms(lines, start_line):
    """
    Read ion positions from the PWSCF output file

    Parameters
    ----------
    lines: List

    Returns
    -------
    ids: List
    symbols: List
    positions: List

    """
    fcoords = []
    ccoords = []
    final_symbols = None
    final_ids = None
    alat = None

    # TODO you can get a section after the optimisation ("Begin final coordinates") but then a final scf calculation
    # followed by a recalculated final unit cell
    final_coords = False

    for i, line in enumerate(lines):
        line = line.strip()

        if 'lattice parameter (alat)' in line:
            if "a.u." not in line:
                raise_error("expecting alat in a.u.", line, i, start_line)
            if alat is not None and not final_coords:
                raise_error("found multiple alat values", line, i, start_line)
            alat = split_numbers(line)[0] * codata[("Bohr", "Angstrom")]

        # initial and final
        if 'Crystallographic axes' in line:
            if fcoords and not final_coords:
                raise_error("already found fractional atomic positions in section", line, i, start_line)
            if lines[i+1].strip():
                raise_error("Expecting blank line after initial Crystallographic axes line",
                            lines[i + 1], i+1, start_line)
            if not fnmatch(lines[i+2].strip(), "site*atom*positions*"):
                raise_error("Expecting 2nd line after Crystallographic axes to have headers: site, atom and positions",
                            lines[i + 2], i+2, start_line)
            symbols = []
            ids = []
            fcoords = []
            j = 3
            while lines[i+j].strip():
                symbols.append(lines[i+j].strip().split()[1].strip('0123456789'))
                ids.append(int(split_numbers(lines[i+j])[0]))
                fcoords.append(split_numbers(lines[i+j])[-3:])
                j += 1

            if not final_symbols:
                final_symbols = symbols
            if not final_ids:
                final_ids = ids
            if ids != final_ids:
                raise_error("difference in ids: {0}, from previously read {1}".format(ids, final_ids),
                            line, i, start_line)
            if symbols != final_symbols:
                raise_error("difference in symbols: {0}, from previously read {1}".format(symbols, final_symbols),
                            line, i, start_line)

        if line == 'Cartesian axes':
            if ccoords and not final_coords:
                raise_error("already found cartesian atomic positions in section", line, i, start_line)
            if lines[i+1].strip():
                raise_error("Expecting blank line after initial Cartesian axes line",
                            lines[i + 1], i+1, start_line)
            if not fnmatch(lines[i+2].strip(), "site*atom*positions*alat*"):
                raise_error("Expecting 2nd line after Crystallographic axes to have "
                            "headers: site, atom and positions (in alat units)", lines[i+2], i+2, start_line)
            symbols = []
            ids = []
            ccoords = []
            j = 3
            while split_numbers(lines[i+j]):
                symbols.append(lines[i+j].strip().split()[1].strip('0123456789'))
                ids.append(int(split_numbers(lines[i+j])[0]))
                x, y, z = split_numbers(lines[i+j])[-3:]
                ccoords.append((x*alat, y*alat, z*alat))
                j += 1

            if not final_symbols:
                final_symbols = symbols
            if not final_ids:
                final_ids = ids
            if ids != final_ids:
                raise_error("difference in ids: {0}, from previously read {1}".format(ids, final_ids),
                            line, i, start_line)
            if symbols != final_symbols:
                raise_error("difference in symbols: {0}, from previously read {1}".format(symbols, final_symbols),
                            line, i, start_line)

        if "Begin final coordinates" in line:
            final_coords = True

        # optimization step
        if "ATOMIC_POSITIONS" in line:
            coords = []
            symbols = []
            ids = []
            j = 1
            while split_numbers(lines[i + j]):
                symbols.append(lines[i + j].strip().split()[0].strip('0123456789'))
                ids.append(j)
                coords.append(split_numbers(lines[i + j])[-3:])
                j += 1

            if not final_symbols:
                final_symbols = symbols
            if not final_ids:
                final_ids = ids
            if ids != final_ids:
                raise_error("difference in ids: {0}, from previously read {1}".format(ids, final_ids),
                            line, i, start_line)
            if symbols != final_symbols:
                raise_error("difference in symbols: {0}, from previously read {1}".format(symbols, final_symbols),
                            line, i, start_line)
            if "crystal_sg" in line:
                warnings.warn("haven't yet implemented reading crystal_sg")
                continue
            if "crystal" in line:
                if fcoords:
                    raise_error("fractional coordinates already set for this section", line, i, start_line)
                fcoords = coords
            if "alat" in line:
                if ccoords:
                    raise_error("cartesian coordinates already set for this section", line, i, start_line)
                ccoords = [(x*alat, y*alat, z*alat) for x, y, z in coords]
            if "bohr" in line:
                if ccoords:
                    raise_error("cartesian coordinates already set for this section", line, i, start_line)
                ccoords = [(x*codata[("Bohr", "Angstrom")], y*codata[("Bohr", "Angstrom")],
                            z*codata[("Bohr", "Angstrom")]) for x, y, z in coords]
            if "angstrom" in line:
                if ccoords:
                    raise_error("cartesian coordinates already set for this section", line, i, start_line)
                ccoords = coords

    return (final_ids, final_symbols,
            fcoords if fcoords else None, {"magnitude": ccoords, "units": "angstrom"} if ccoords else None)


def read_stress(lines, start_line):
    """
    Read the stress from the PWSCF output file

    Parameters
    ----------
    lines: List

    Returns
    -------
    stress: dict


    """
    stress = []
    # ASE convention for the stress tensor appears to differ
    # from the PWscf one by a factor of -1
    conv_factor = -1.0 * codata[("Rydberg", "eV")] / codata[("Bohr", "Angstrom")] ** 3

    for i, line in enumerate(lines):
        line = line.strip()
        if fnmatch(line, "total*stress*"):
            if "Ry/bohr**3" not in line:
                raise_error("expecting stress in Ry/bohr**3", line, i, start_line)

            for nrow, newline in enumerate(lines[i + 1: i + 4]):
                stress.append([s*conv_factor for s in split_numbers(newline)[:3]])

    if not stress:
        return None
    else:
        return {"magnitude": stress, "units": 'eV/angstrom^3'}


def get_band_mapping(lines):
    """

    Parameters
    ----------
    lines

    Returns
    -------
    crystal_coord_map: list

    """
    # by default they are given in cartesian*2pi/alat, which is not helpful
    crystal_coord_map = []
    in_kpoints = False
    in_crystal_coord = False
    for i, line in enumerate(lines):
        line = line.strip()
        # if fnmatch(line, "End of*calculation"):
        #     break
        if "number of k points=" in line:
            in_kpoints = True
        if not in_kpoints:
            continue
        if not line and not fnmatch(lines[i+2], "*k*(*)*=*(*)*wk*=*"):
            break
        if "cryst. coord." in line:
            in_crystal_coord = True
        if not in_crystal_coord:
            continue
        if fnmatch(line, "*k*(*)*=*(*)*wk*=*"):
            i, x, y, z, wk = split_numbers(line)
            crystal_coord_map.append((x, y, z))

    return crystal_coord_map


# TODO occupation numbers (space below energies, if present)
# TODO maybe in cartesian and crystal
def read_bands(lines, start_line, crystal_coord_map=None):
    """ read band structure

    Parameters
    ----------
    lines: List
    start_line: int
    crystal_coord_map: List
        mapping (by index) of k-point coordinate to crystal coordinates

    Returns
    -------

    """
    data = {"coord_type": "crystal" if crystal_coord_map else "cartesian*2pi/alat",
            "all": {"kpoints": [], "pws": [], "energies_per_band": []},
            "spinup": {"kpoints": [], "pws": [], "energies_per_band": []},
            "spindown": {"kpoints": [], "pws": [], "energies_per_band": []}}
    key = "all"
    in_band_data = False
    for i, line in enumerate(lines):
        line = line.strip()
        if fnmatch(line, "End of*calculation"):
            in_band_data = True
        if not in_band_data:
            continue

        if fnmatch(line, "*SPIN UP*"):
            key = "spinup"
        elif fnmatch(line, "*SPIN DOWN*"):
            key = "spindown"
        elif fnmatch(line, "*k*=*PWs*bands*"):
            if "ev" not in line:
                raise_error("expecting units in ev", line, i, start_line)
            if not lines[i+1].strip() == "":
                raise_error("expecting empty line after k-point", lines[i+1], i+1, start_line)

            if not crystal_coord_map:
                coords = tuple(split_numbers(line)[0:3])
            else:
                coords = crystal_coord_map[len(data[key]["kpoints"])]

            data[key]["kpoints"].append(coords)
            data[key]["pws"].append(int(split_numbers(line)[3]))

            j = 2
            energies = []
            while lines[i+j].strip():
                energies += split_numbers(lines[i+j])
                j += 1
            data[key]["energies_per_band"].append(energies)

    some_data = False
    for key in ["all", "spinup", "spindown"]:
        if not data[key]["kpoints"]:
            data.pop(key)
            continue
        some_data = True
        data[key]["energies_per_band"] = [{"units": "eV", "magnitude": e}
                                          for e in np.array(data[key]["energies_per_band"]).T.tolist()]
    return data if some_data else None


# TODO get metadata, like version numbers, cpu/wall time(s)
def first_parse(lines):
    """ first parse of file to fin start/end of sections and get metadata

    Parameters
    ----------
    lines: List

    Returns
    -------

    """
    scf_start_first = None
    scf_end_last = None
    opt_start = None
    opt_end = None
    steps_num = None
    steps = []
    warnings = []
    for i, line in enumerate(lines):
        line = line.strip()
        if "DEPRECATED" in line:
            warnings.append(line)
        if "Geometry Optimization" in line:
            if "End of" in line:
                if not opt_start:
                    raise IOError("found end of geometry optimization before start in line #{0}: {1}".format(i, line))
                if opt_end:
                    raise IOError("found multiple geometry optimization ends in line #{0} and #{1}".format(opt_end, i))
                opt_end = i
            else:
                if opt_start:
                    raise IOError("found multiple geometry optimization starts in "
                                  "line #{0} and #{1}".format(opt_start, i))
                opt_start = i
        if fnmatch(line, "number of * steps *=*") and opt_start and not opt_end:
            new_step = split_numbers(line)[0]
            if not steps_num is None:
                if new_step != steps_num + 1:
                    pass
                    # NOTE: this isn't strictly true since the step history can be reset
                    # TODO maybe use the line above "number of scf cycles = "
                    # raise IOError("expecting geometry optimization step {0} at "
                    #               "line #{1}: {2}".format(steps_num+1, i, line))
                steps[-1] = (steps[-1][0], i)
                steps_num += 1
            else:
                steps_num = 0
            # steps_num = new_step
            steps.append((i, None))

        if "Self-consistent Calculation" in line and not scf_start_first:
            scf_start_first = i
        if "End of self-consistent calculation" in line:
            scf_end_last = i

    if opt_end and steps:
        steps[-1] = (steps[-1][0], opt_end)

    if opt_start and not opt_end:
        raise IOError("found start of geometry optimization but no end")

    if scf_start_first and not scf_end_last:
        raise IOError("found start of SCF but no end")

    return opt_start, opt_end, steps, scf_start_first, scf_end_last, warnings


def get_data_section(lines, start_line, crystal_coord_map=None):
    """ get data for a particular section

    Parameters
    ----------
    lines: list
    start_line: int
    crystal_coord_map: list or None

    Returns
    -------

    """
    out = {}

    ids, symbols, fcoords, ccoords = read_atoms(lines, start_line)
    out["ids"] = ids
    out['symbols'] = symbols
    out['fcoords'] = fcoords
    out['ccoords'] = ccoords
    out['cell_vectors'] = read_cell(lines, start_line)
    out['energy'] = read_energies(lines, start_line)
    out["scf"], out["sphere_integration"] = read_scf(lines, start_line)
    out['forces'] = read_forces(lines, start_line)
    out['stress'] = read_stress(lines, start_line)
    out["bands"] = read_bands(lines, start_line, crystal_coord_map)

    return out


class QEmainPlugin(object):
    """ quantum espresso scf output parser plugin for jsonextended
    """
    plugin_name = 'quantum_espresso_scf_output'
    plugin_descript = 'read quantum espresso scf output'
    file_regex = '*.qe.out'

    def read_file(self, file_obj, **kwargs):
        """

        Parameters
        ----------
        file_obj: file_object
        kwargs

        Returns
        -------

        """

        lines = file_obj.readlines()  # type: List

        all = {}

        opt_start, opt_end, opt_steps, scf_start_first, scf_end_last, warnings = first_parse(lines)

        all["warnings"] = warnings if warnings else None
        all["errors"] = None  # TODO check for non terminating errors (e.g. reaching max scf steps)
        if opt_start:
            crystal_coord_map = get_band_mapping(lines[0:opt_start])
            all["initial"] = get_data_section(lines[0:opt_start], 0, crystal_coord_map)
            all["optimisation"] = [get_data_section(lines[i:j], i, crystal_coord_map) for i, j in opt_steps]
            all["final"] = get_data_section(lines[opt_end:], opt_end, crystal_coord_map)
        elif scf_start_first:
            crystal_coord_map = get_band_mapping(lines[0:scf_start_first])
            all["initial"] = get_data_section(lines[0:scf_start_first], 0, crystal_coord_map)
            all["optimisation"] = [get_data_section(lines[scf_start_first:scf_end_last+1],
                                                    scf_start_first, crystal_coord_map)]
            all["final"] = get_data_section(lines[scf_end_last:], scf_end_last, crystal_coord_map)
        else:
            crystal_coord_map = get_band_mapping(lines)
            all["initial"] = None
            all["optimisation"] = None
            all["final"] = get_data_section(lines, 0, crystal_coord_map)

        # TODO copy data from initial or final optimisation to final section (if not already present)?
        # (for instance if the geometry does not change) or leave for user to deal with?

        all["creator"] = {"program": "Quantum Espresso"}

        return all


class QEbandPlugin(QEmainPlugin):
    """ quantum espresso output parser plugin for jsonextended
    """
    plugin_name = 'quantum_espresso_band_output'
    plugin_descript = 'read quantum espresso output (identical to main)'
    file_regex = '*.qe.band.out'
