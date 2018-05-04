# -*- coding: utf-8 -*-
"""
adapted from: https://github.com/lmmentel/ase-espresso

Original work Copyright (C) 2013-2015 SUNCAT
Modified work Copyright 2015-2017 Lukasz Mentel

                    GNU GENERAL PUBLIC LICENSE
                       Version 3, 29 June 2007

"""
import warnings
from datetime import datetime
from fnmatch import fnmatch

import numpy as np
from ase import Atoms
from ejplugins.utils import codata, split_numbers, symbol2anum
from jsonextended import edict


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
    raise IOError("{0}, for line #{1}: {2}".format(msg, i + 1 + start_line, line.strip()))


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
            if lines[i + 1].strip():
                raise_error("Expecting blank line after initial forces line", lines[i + 1], i + 1, start_line)
            forces = {"peratom": []}
            j = 2
            while fnmatch(lines[i + j].strip(), 'atom *'):
                aid, atype, fx, fy, fz = split_numbers(lines[i + j])
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
            if lines[i + 1].strip():
                raise_error("Expecting blank line after initial energy contributions line", lines[i + 1], i + 1,
                            start_line)
            j = 2
            while lines[i + j].strip():
                if lines[i + j].strip().startswith("->"):
                    # These are sub component contributions, so ignore
                    j += 1
                    continue
                if not fnmatch(lines[i + j].strip(), "*=*Ry"):
                    raise_error("Expecting energy contribution to be in the form '*=*Ry'",
                                lines[i + j], i + j, start_line)
                name = lines[i + j].strip().split("=")[0]  # type: str
                name = name.replace("contribution", "").replace(" contrib.", "").replace("  ", " ").strip()
                value = split_numbers(lines[i + j])[0] * codata[("Rydberg", "eV")]
                if name in contributions:
                    raise_error("Multiple instances of energy contribution {} found in section".format(name),
                                lines[i + j], i + j, start_line)
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
                if last_magmom_num % 100 == 0.:
                    pass  # it appears that the magnetic moment is printed every 100 iterations
                else:
                    raise_error("found multiple magnetic moments per site in same scf run", line, i, start_line)
            last_magmom_num = cur_cyc_num
            atomic_charges_peratom = []
            spin_density_peratom = []

            nxtline = i + 1
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
    # followed by a recalculated final unit cell (see ejplugins/test_files/qe_vc_relax_with_multiple_final_coords.out)
    final_coords = False

    # additional unit cell information
    # bli_lines = [line for line in lines if 'bravais-lattice index' in line]
    # brav_latt_indices = [int(line.split('=')[1].strip()) for line in bli_lines]

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
                cell[key] = {"units": "angstrom", "magnitude": (x * alat, y * alat, z * alat)}

        if "Begin final coordinates" in line:
            final_coords = True

        # optimisation step
        if fnmatch(line, "CELL_PARAMETERS*(*alat*=*)"):
            alat = split_numbers(line)[0] * codata[("Bohr", "Angstrom")]
            if cell and final_coords:
                # there can be a 'Begin/End final coordinates' in the final opt step
                new_cell = {}
                for key, newline in zip(["a", "b", "c"], lines[i + 1: i + 4]):
                    x, y, z = split_numbers(newline)
                    new_cell[key] = {"units": "angstrom", "magnitude": (x * alat, y * alat, z * alat)}
                if edict.diff(cell, new_cell, np_allclose=True):
                    raise_error("normal and final cell coordinates are different", line, i, start_line)
            elif cell:
                raise_error("found multiple cell coordinates", line, i, start_line)
            else:

                for key, newline in zip(["a", "b", "c"], lines[i + 1: i + 4]):
                    x, y, z = split_numbers(newline)
                    cell[key] = {"units": "angstrom", "magnitude": (x * alat, y * alat, z * alat)}

    return cell if cell else None


def read_atoms(lines, start_line, section="n/a"):
    """
    Read ion positions from the PWSCF output file

    Parameters
    ----------
    lines: List
    start_line: int
    section: str

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
                raise_error("already found fractional atomic positions in section ({})".format(section),
                            line, i, start_line)
            if lines[i + 1].strip():
                raise_error("Expecting blank line after initial Crystallographic axes line",
                            lines[i + 1], i + 1, start_line)
            if not fnmatch(lines[i + 2].strip(), "site*atom*positions*"):
                raise_error("Expecting 2nd line after Crystallographic axes to have headers: site, atom and positions",
                            lines[i + 2], i + 2, start_line)
            symbols = []
            ids = []
            fcoords = []
            j = 3
            while lines[i + j].strip():
                symbols.append(lines[i + j].strip().split()[1].strip('0123456789'))
                ids.append(int(split_numbers(lines[i + j])[0]))
                fcoords.append(split_numbers(lines[i + j])[-3:])
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
                raise_error("already found cartesian atomic positions in section ({})".format(section),
                            line, i, start_line)
            if lines[i + 1].strip():
                raise_error("Expecting blank line after initial Cartesian axes line",
                            lines[i + 1], i + 1, start_line)
            if not fnmatch(lines[i + 2].strip(), "site*atom*positions*alat*"):
                raise_error("Expecting 2nd line after Crystallographic axes to have "
                            "headers: site, atom and positions (in alat units)", lines[i + 2], i + 2, start_line)
            symbols = []
            ids = []
            ccoords = []
            j = 3
            while split_numbers(lines[i + j]):
                symbols.append(lines[i + j].strip().split()[1].strip('0123456789'))
                ids.append(int(split_numbers(lines[i + j])[0]))
                x, y, z = split_numbers(lines[i + j])[-3:]
                ccoords.append((x * alat, y * alat, z * alat))
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
            # TODO add fixed positions to output
            fixed_pos = []
            symbols = []
            ids = []
            j = 1
            while split_numbers(lines[i + j]):
                symbols.append(lines[i + j].strip().split()[0].strip('0123456789'))
                ids.append(j)
                nvalues = [float(f) for f in lines[i + j].split()[1:]]
                coords.append(nvalues[0:3])
                if len(nvalues) == 3:
                    fixed_pos.append((False, False, False))
                elif len(nvalues) == 6:
                    fixed_pos.append(tuple([not bool(i) for i in nvalues[3:6]]))
                else:
                    raise_error("ATOMIC_POSITIONS: expecting either 'symbol x y z' or 'symbol x y z fx fy fz'",
                                line, i, start_line)

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
                if fcoords and fcoords != coords:
                    raise_error("fractional coordinates already set for this section ({})".format(section),
                                line, i, start_line)
                fcoords = coords
            if "alat" in line:
                coords = [(x * alat, y * alat, z * alat) for x, y, z in coords]
                if ccoords and ccoords != coords:
                    raise_error("cartesian coordinates already set for this section ({})".format(section),
                                line, i, start_line)
                ccoords = coords
            if "bohr" in line:
                coords = [(x * codata[("Bohr", "Angstrom")], y * codata[("Bohr", "Angstrom")],
                           z * codata[("Bohr", "Angstrom")]) for x, y, z in coords]
                if ccoords and ccoords != coords:
                    raise_error("cartesian coordinates already set for this section ({})".format(section),
                                line, i, start_line)
                ccoords = coords
            if "angstrom" in line:
                if ccoords and ccoords != coords:
                    raise_error("cartesian coordinates already set for this section ({})".format(section),
                                line, i, start_line)
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
                stress.append([s * conv_factor for s in split_numbers(newline)[:3]])

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
        if not line and not fnmatch(lines[i + 2], "*k*(*)*=*(*)*wk*=*"):
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
            if not lines[i + 1].strip() == "":
                raise_error("expecting empty line after k-point", lines[i + 1], i + 1, start_line)

            if not crystal_coord_map:
                coords = tuple(split_numbers(line)[0:3])
            else:
                coords = crystal_coord_map[len(data[key]["kpoints"])]

            data[key]["kpoints"].append(coords)
            data[key]["pws"].append(int(split_numbers(line)[3]))

            j = 2
            energies = []
            while lines[i + j].strip():
                energies += split_numbers(lines[i + j])
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
    non_terminating_errors = []
    found_job_done = False
    start_time, end_time, elapsed_time, nprocs = None, None, None, None
    for i, line in enumerate(lines):
        line = line.strip()
        if line.startswith("JOB DONE."):
            found_job_done = True
        # convergence NOT achieved after *** iterations: stopping
        if "convergence NOT achieved" in line:
            non_terminating_errors.append(line)
        if "Error in routine cdiaghg" in line:
            non_terminating_errors.append(line)
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

            # Can get:
            # lsda relaxation :  a final configuration with zero
            #                    absolute magnetization has been found
            #
            # the program is checking if it is really the minimum energy structure
            # by performing a new scf iteration without any "electronic" history
            #
        # TODO should this be a separate step?
        if fnmatch(line, '*performing a new scf iteration without any "electronic" history'):
            steps[-1] = (steps[-1][0], i)
            steps_num += 1
            steps.append((i, None))

        if "Self-consistent Calculation" in line and not scf_start_first:
            scf_start_first = i
        if "End of self-consistent calculation" in line:
            scf_end_last = i

        if fnmatch(line, "Parallel version *, running on * processors"):
            if nprocs is not None:
                raise IOError("found second nprocs on line {0}: {1}".format(i, line))
            nprocs = int(line.split()[-2])
        elif fnmatch(line, "Parallel version *, running on * processor cores"):
            if nprocs is not None:
                raise IOError("found second nprocs on line {0}: {1}".format(i, line))
            nprocs = int(line.split()[-3])

        # NB: time uses spaces instead of 0, e.g. This run was terminated on:  12:43: 3   6Sep2017
        if fnmatch(line, "Program*starts on * at *"):
            if start_time is not None:
                raise IOError("found second start time on line {0}: {1}".format(i, line))
            line_segs = line.split("at")
            atime = line_segs[-1].replace(" ", "")
            adate = line_segs[-2].split()[-1]
            start_time = datetime.strptime(atime + ' ' + adate, '%H:%M:%S %d%b%Y')
        if line.startswith("This run was terminated on:"):
            if end_time is not None:
                raise IOError("found second end time on line {0}: {1}".format(i, line))
            line_segs = line.split()
            adate = line_segs[-1]
            atime = ":".join(" ".join(line_segs[:-1]).split(":")[-3:]).replace(" ", "")
            end_time = datetime.strptime(atime + ' ' + adate, '%H:%M:%S %d%b%Y')

    if not found_job_done:
        non_terminating_errors.append("Did not find indicator: JOB DONE.")

    if opt_end and steps:
        steps[-1] = (steps[-1][0], opt_end)

    if opt_start and not opt_end:
        raise IOError("found start of geometry optimization but no end")

    if scf_start_first and not scf_end_last:
        raise IOError("found start of SCF but no end")

    if start_time and end_time:
        delta_time = end_time - start_time
        m, s = divmod(delta_time.total_seconds(), 60)
        h, m = divmod(m, 60)
        elapsed_time = "%d:%02d:%02d" % (h, m, s)

    meta = {"elapsed_time": elapsed_time, "nprocs": nprocs}

    return opt_start, opt_end, steps, scf_start_first, scf_end_last, warnings, non_terminating_errors, meta


def get_data_section(lines, start_line, crystal_coord_map=None, section="n/a"):
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

    ids, symbols, fcoords, ccoords = read_atoms(lines, start_line, section)
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
    plugin_descript = 'read quantum espresso main output (from pw.x)'
    file_regex = '*.qe.out'

    def read_file(self, file_obj, skip_opt=False, **kwargs):
        """

        Parameters
        ----------
        file_obj: file_object
        skip_opt: bool
            if True, skip retrieval of optimisation step data

        Returns
        -------

        """

        lines = file_obj.readlines()  # type: List

        all = {}

        opt_start, opt_end, opt_steps, scf_start_first, scf_end_last, warnings, non_terminating_errors, meta = first_parse(
            lines)

        all["warnings"] = warnings if warnings else None
        all[
            "errors"] = non_terminating_errors  # TODO check for more non terminating errors (e.g. reaching max opt steps)
        if opt_start:
            crystal_coord_map = get_band_mapping(lines[0:opt_start])
            all["initial"] = get_data_section(lines[0:opt_start], 0, crystal_coord_map, "initial")
            if skip_opt:
                all["optimisation"] = None
            else:
                all["optimisation"] = [get_data_section(lines[i:j], i, crystal_coord_map, "opt{}".format(step + 1))
                                       for step, (i, j) in enumerate(opt_steps)]
            all["final"] = get_data_section(lines[opt_end:], opt_end, crystal_coord_map, "final")
        elif scf_start_first:
            crystal_coord_map = get_band_mapping(lines[0:scf_start_first])
            all["initial"] = get_data_section(lines[0:scf_start_first], 0, crystal_coord_map, "initial")
            if skip_opt:
                all["optimisation"] = None
            else:
                all["optimisation"] = [get_data_section(lines[scf_start_first:scf_end_last + 1],
                                                        scf_start_first, crystal_coord_map, "opt1")]
            all["final"] = get_data_section(lines[scf_end_last:], scf_end_last, crystal_coord_map, "final")
        else:
            crystal_coord_map = get_band_mapping(lines)
            all["initial"] = None
            all["optimisation"] = None
            all["final"] = get_data_section(lines, 0, crystal_coord_map, "final")

        # TODO copy data from initial or final optimisation to final section (if not already present)?
        # (for instance if the geometry does not change) or leave for user to deal with?

        all["creator"] = {"program": "Quantum Espresso"}
        all["meta"] = meta

        return all


class QEnscfPlugin(QEmainPlugin):
    """ quantum espresso output parser plugin for jsonextended
    """
    plugin_name = 'quantum_espresso_nscf_output'
    plugin_descript = 'read quantum espresso output (identical to main)'
    file_regex = '*.qe.nscf.out'


class QEbandPlugin(QEmainPlugin):
    """ quantum espresso output parser plugin for jsonextended
    """
    plugin_name = 'quantum_espresso_band_output'
    plugin_descript = 'read quantum espresso output (identical to main)'
    file_regex = '*.qe.band.out'


_PNUM2TYPE = {0: "charge (pseudo)",
              1: "potential (total)",
              2: "potential (ionic)",
              3: "states (at specific energy)",
              4: "electronic entropy",
              5: "STM",
              6: "spin",
              7: "wavefuction contribution",
              8: "electron localisation function",
              9: "charge (minus superposition of atomic densities)",
              10: "states (integrated)",
              11: "potential (V_bare + V_h)",
              12: "sawtooth electric field potential",
              13: "noncollinear magnetization",
              17: "charge (all-electron valence)",
              18: "exchange and correlation magnetic field",
              19: "Reduced density gradient",
              20: "charge * 2nd eig Hessian",
              21: "charge (all-electron)"}


class QEChargeDensityPlugin(object):
    """quantum espresso output parser plugin for jsonextended

    """
    plugin_name = 'quantum_espresso_charge_density_output'
    plugin_descript = 'read quantum espresso charge density output (from pp.x), ' \
                      'in Gaussian Cube Format (output_format=6)'
    file_regex = '*.qe.charge'

    def read_file(self, f, **kwargs):
        comment = f.readline().strip()  # first line blank
        line = f.readline().strip()
        try:
            na, nb, nc, _, _, _, natoms, ntyp = [int(i) for i in line.split()]
        except:
            raise IOError(
                "file format incorrect, expected 8 fields; na, nb, nc, _, _, _, natoms, ntyp: {0}".format(line))

        line = f.readline().strip()
        try:
            bravais_lattice_index, alat, blat, clat, alpha, beta, gamma = [float(i) for i in line.split()]
            alat *= codata[("Bohr", "Angstrom")]
            blat *= codata[("Bohr", "Angstrom")]
            clat *= codata[("Bohr", "Angstrom")]
            bravais_lattice_index = int(bravais_lattice_index)
        except:
            raise IOError("file format incorrect, expected ibrav, a, b, c, alpha, beta, gamma: {0}".format(line))

        if bravais_lattice_index == 0:  # free
            try:
                line = f.readline().strip()
                avec = [float(i) * alat for i in line.split()]
                line = f.readline().strip()
                bvec = [float(i) * alat for i in line.split()]
                line = f.readline().strip()
                cvec = [float(i) * alat for i in line.split()]
            except:
                raise IOError("file format incorrect, expected fields; x, y, z: {0}".format(line))
        elif bravais_lattice_index == 1:  # cubic P (sc)
            avec = (alat, 0., 0.)
            bvec = (0., alat, 0.)
            cvec = (0., 0., alat)
        elif bravais_lattice_index == 2:  # cubic F (fcc)
            avec = (-alat / 2., 0., alat / 2)
            bvec = (0., alat / 2., alat / 2.)
            cvec = (-alat / 2., alat / 2., 0.)
        elif bravais_lattice_index == 3:  # cubic I (bcc)
            avec = (alat / 2., alat / 2., alat / 2)
            bvec = (-alat / 2., alat / 2., alat / 2)
            cvec = (-alat / 2., -alat / 2., alat / 2)
        elif bravais_lattice_index == -3:  # cubic I (bcc), more symmetric axis
            avec = (-alat / 2., alat / 2., alat / 2)
            bvec = (alat / 2., -alat / 2., alat / 2)
            cvec = (alat / 2., alat / 2., -alat / 2)
        else:
            # TODO implemented ibrav > 3
            raise NotImplementedError("haven't yet implemented ibrav > 3")

        line = f.readline()  # TODO find out what numbers in this line are (third is Ecutoff)
        plotnum = int(line.strip().split()[3])
        dtype = _PNUM2TYPE[plotnum]

        typ_lookup = {}
        for _ in range(ntyp):
            line = f.readline().strip()
            try:
                i, symbol, valence_electrons = line.split()
            except:
                raise IOError("file format incorrect, expected fields; i, symbol, valence_electrons: {0}".format(line))
            symbol = "".join([s for s in symbol if s.isalpha()])
            atomic_number = symbol2anum(symbol)
            typ_lookup[i] = (symbol, atomic_number, float(valence_electrons))

        ccoords = []
        symbols = []
        atomic_numbers = []
        valence_charges = []
        for _ in range(natoms):
            line = f.readline().strip()
            try:
                i, a, b, c, atyp = line.split()
                ccoord = np.array([a, b, c], dtype=float) * alat
                ccoords.append(ccoord.tolist())
                symbols.append(typ_lookup[atyp][0])
                atomic_numbers.append(typ_lookup[atyp][1])
                valence_charges.append(typ_lookup[atyp][2])
            except:
                raise IOError("file format incorrect, expected fields; i, a, b, c, atyp: {0}".format(line))

        atoms = Atoms(symbols=symbols, positions=ccoords, cell=[avec, bvec, cvec], pbc=[True, True, True])
        atoms.wrap()
        ccoords = atoms.positions.tolist()

        charge_density = []
        line = f.readline().strip().split()
        while line:
            charge_density += line
            line = f.readline().strip().split()
        dense = np.array(charge_density, dtype=float).reshape((nc, nb, na)).transpose()

        return {
            "title": "Quantum Espresso " + comment,
            # "na": na, "nb": nb, "nc": nc,
            "cell_vectors": {
                "a": {"units": "angstrom", "magnitude": avec},
                "b": {"units": "angstrom", "magnitude": bvec},
                "c": {"units": "angstrom", "magnitude": cvec}
            },
            # "centre": [0, 0, 0],
            "densities": [{
                "type": dtype,
                "magnitude": dense
            }],
            "atoms": {"ccoords": {"units": "angstrom",
                                  "magnitude": ccoords},
                      # "nuclear_charge": nuclear_charges,
                      "nuclear_charge": valence_charges,
                      "symbols": symbols, "atomic_number": atomic_numbers},
            "creator": {"program": "Quantum Espresso"}
        }


class QELowdinPlugin(object):
    """quantum espresso output parser plugin for jsonextended

    https://www.researchgate.net/post/What_is_the_origin_of_difference_between_Lowdin_charge_when_calculated_with_plane_wave_code_and_gaussian

    """
    plugin_name = 'quantum_espresso_pdos_output'
    plugin_descript = 'read quantum espresso charge density output (from pw.x), mainly for Lodwin charges partitioning'
    file_regex = '*.qe.pdos.out'

    def read_file(self, f, **kwargs):
        dic = {"final": {"lowdin": {"charge": [], "spin": []}}, "creator": {"program": "Quantum Espresso"}}
        in_lowdin = False
        natoms = 0
        for line in f:
            if "Lowdin Charges" in line:
                in_lowdin = True
                continue
            if not in_lowdin:
                continue
            if "Spilling Parameter" in line:
                in_lowdin = False
                dic["final"]["lowdin"]["spilling"] = split_numbers(line)[0]
                continue
            if "Atom #" in line:
                natoms = split_numbers(line)[0]
            if "total charge" in line:
                charge = split_numbers(line.split("total charge")[1])[0]
                dic["final"]["lowdin"]["charge"].append(charge)
                continue
            if "polarization" in line:
                spin = split_numbers(line.split("polarization")[1])[0]
                dic["final"]["lowdin"]["spin"].append(spin)
                continue

        if len(dic["final"]["lowdin"]["charge"]) != natoms:
            raise IOError("number of atoms ({0}) is different to the number of total charge values found ({1})".format(
                len(dic["final"]["lowdin"]["charge"]), natoms))
        if len(dic["final"]["lowdin"]["spin"]) != natoms:
            raise IOError("number of atoms ({0}) is different to the number of total spin values found ({1})".format(
                len(dic["final"]["lowdin"]["spin"]), natoms))

        return dic


class QEdosPlugin(object):
    """quantum espresso output parser plugin for jsonextended

    https://www.researchgate.net/post/What_is_the_origin_of_difference_between_Lowdin_charge_when_calculated_with_plane_wave_code_and_gaussian

    """
    plugin_name = 'quantum_espresso_dos'
    plugin_descript = 'read quantum espresso total density of states output'
    file_regex = '*.qe.dos'

    def read_file(self, f, **kwargs):
        first_line = f.readline()
        # TODO extend to use with no spin polarised calculation output
        if not fnmatch(first_line, "*E (eV)*dosup(E)*dosdw(E)*Int dos(E)*"):
            raise IOError("*.qe.dos header does not match; *E (eV)*dosup(E)*dosdw(E)*Int dos(E)*")
        line = f.readline().strip()
        es, ups, downs = [], [], []
        line = f.readline()  # first line has no energy
        while line:
            e, up, down, intd = line.split()
            es.append(float(e))
            ups.append(float(up))
            downs.append(float(down))
            line = f.readline().strip()

        return {"tdos": {"energy": {"units": "eV", "magnitude": es}, "up": ups, "down": downs}}
