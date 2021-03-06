# -*- coding: utf-8 -*-
"""
originally written for CRYSTAL14, but also tested against and extended for CRYSTAL17
"""
import copy
import logging
import re
import warnings
from fnmatch import fnmatch

import numpy as np
from jsonextended import edict
from ejplugins.utils import codata, split_numbers

try:
    from distutils.util import strtobool
except ImportError:
    from distutils import strtobool

logger = logging.getLogger(__name__)

_ihferm_map = {0: 'closed shell, insulating system',
               1: 'open shell, insulating system',
               2: 'closed shell, conducting system',
               3: 'open shell, conducting system'}


class BANDPlugin(object):
    """ band parser plugin for jsonextended
    """
    plugin_name = 'crystal_band'
    plugin_descript = 'read CRYSTAL band output'
    file_regex = '*.band.f25'

    def read_file(self, f, **kwargs):
        band_segments = []

        line = f.readline().strip()
        while line:
            if line.startswith('-%-'):
                segment = {}
                segment['system'] = _ihferm_map[int(line[3])]
                segment['calc_type'] = line[4:8]
                segment["coord_type"] = "crystal*ISS"

                nbands, nkpts, dummy, kdist, fermi = split_numbers(line[8:])
                nbands, nkpts, kdist = int(nbands), int(nkpts), float(kdist)
                segment['nkpts_per_band'] = nkpts
                segment['fermi'] = {"units": "eV", "magnitude": float(fermi) * codata[("Hartree", "eV")]}
                line = f.readline().strip()
                EMIN, EMAX = split_numbers(line)
                segment['erange'] = {"units": "eV", "magnitude": (EMIN * codata[("Hartree", "eV")],
                                                                  EMAX * codata[("Hartree", "eV")])}
                line = f.readline().strip()
                i1, i2, i3, j1, j2, j3 = [float(i) for i in line.split()]
                segment['krange'] = [(i1, i2, i3), (j1, j2, j3)]
                line = f.readline().strip()
                energy = []
                while not line.startswith('-%-') and line:
                    energy += split_numbers(line)
                    line = f.readline().strip()
                allenergies = np.array([energy[i * nbands:i * nbands + nbands] for i in range(nkpts)]).T
                allenergies *= codata[("Hartree", "eV")]
                segment["kdist"] = np.linspace(0, kdist * nkpts, nkpts).tolist()

                kstart = np.array((i1, i2, i3))
                kend = np.array((j1, j2, j3))
                segment["kpoints"] = []
                for i in np.linspace(0, nkpts, nkpts).tolist():
                    kpoint = kstart + (i / float(nkpts)) * (kend - kstart)
                    segment["kpoints"].append(kpoint.tolist())

                segment["energies_per_band"] = []
                for energies in allenergies.tolist():
                    segment["energies_per_band"].append({"units": "eV", "magnitude": energies})

                # segment['kspan'] = kdist * nkpts

                band_segments.append(segment)
            else:
                line = f.readline().strip()
        return {"segments": band_segments,
                "creator": {"program": "CRYSTAL"}}


class DOSSPlugin(object):
    """ doss parser plugin for jsonextended

    Notes
    -----
    File Format:

    1ST RECORD : -%-,IHFERM,TYPE,NROW,NCOL,DX,DY,COSXY (format : A3,I1,A4,2I5,1P,(3E12.5))
    2ND RECORD : X0,Y0 (format : 1P,6E12.5)
    3RD RECORD : I1,I2,I3,I4,I5,I6 (format : 6I3)
    4TH RECORD
    AND FOLLOWING : ((RDAT(I,J),I=1,NROW),J=1,NCOL) (format : 1P,6E12.5)

    Meaning of the variables:
    1   NROW            1 (DOSS are written one projection at a time)
        NCOL            number of energy points in which the DOS is calculated
        DX              energy increment (hartree)
        DY              not used
        COSXY           Fermi energy (hartree)
    2   X0              energy corresponding to the first point
        Y0              not used
    3   I1              number of the projection;
        I2              number of atomic orbitals of the projection;
        I3,I4,I5,I6     not used
    4   RO(J),J=1,NCOL  DOS: density of states ro(eps(j)) (atomic units).

    """
    plugin_name = 'crystal_doss'
    plugin_descript = 'read CRYSTAL density of states output'
    file_regex = '*.doss.f25'

    def __init__(self):
        warnings.warn("DOSSPlugin is deprecated in favour of CrystalDOSPlugin", FutureWarning)
        super(DOSSPlugin, self).__init__()

    def read_file(self, f, **kwargs):
        doss_segments = []

        line = f.readline().strip()
        while line:
            if line.startswith('-%-'):
                projection = {}
                projection['ihferm'] = _ihferm_map[int(line[3])]
                projection['type'] = line[4:8]

                nrows, ncols, _, denergy, fermi = split_numbers(line[8:])
                nrows, ncols, energy_delta = (int(nrows), int(ncols),
                                              float(denergy) * codata[("Hartree", "eV")])
                projection['fermi_energy'] = {"units": "eV", "magnitude": float(fermi) * codata[("Hartree", "eV")]}

                line = f.readline().strip()
                initial_energy = split_numbers(line)[1] * codata[("Hartree", "eV")]

                line = f.readline().strip()
                projection['projid'], projection['norbitals'], _, _, _, _ = [int(i) for i in line.split()]
                line = f.readline().strip()
                dos = []
                while not line.startswith('-%-') and line:
                    dos += split_numbers(line)
                    line = f.readline().strip()
                projection['dos'] = dos
                projection['energy'] = {"units": "eV",
                                        "magnitude": np.linspace(initial_energy,
                                                                 initial_energy + len(dos) * energy_delta,
                                                                 len(dos)).tolist()}
                doss_segments.append(projection)
            else:
                line = f.readline().strip()
        return {"projections": doss_segments,
                "creator": {"program": "CRYSTAL"}}


class CrystalDOSPlugin(object):
    """ doss parser plugin for jsonextended

    Notes
    -----
    File Format:

    1ST RECORD : -%-,IHFERM,TYPE,NROW,NCOL,DX,DY,COSXY (format : A3,I1,A4,2I5,1P,(3E12.5))
    2ND RECORD : X0,Y0 (format : 1P,6E12.5)
    3RD RECORD : I1,I2,I3,I4,I5,I6 (format : 6I3)
    4TH RECORD
    AND FOLLOWING : ((RDAT(I,J),I=1,NROW),J=1,NCOL) (format : 1P,6E12.5)

    Meaning of the variables:
    1   NROW            1 (DOSS are written one projection at a time)
        NCOL            number of energy points in which the DOS is calculated
        DX              energy increment (hartree)
        DY              not used
        COSXY           Fermi energy (hartree)
    2   X0              energy corresponding to the first point
        Y0              not used
    3   I1              number of the projection;
        I2              number of atomic orbitals of the projection;
        I3,I4,I5,I6     not used
    4   RO(J),J=1,NCOL  DOS: density of states ro(eps(j)) (atomic units).

    """
    plugin_name = 'crystal_doss'
    plugin_descript = 'read CRYSTAL density of states output'
    file_regex = '*.doss.f25'

    def read_file(self, f, **kwargs):

        system_type = None
        fermi_energy = None
        energy_delta = None
        initial_energy = None
        len_dos = None
        alpha_projections = {}
        beta_projections = {}
        proj_number = 0

        line = f.readline().strip()
        while line:
            if line.startswith('-%-'):
                proj_number += 1

                if system_type is None:
                    system_type = line[3]
                elif not system_type == line[3]:
                    raise IOError("projection {0} has different system type ({1}) to previous ({2})".format(
                        proj_number, line[3], system_type
                    ))

                if not line[4:8] == "DOSS":
                    raise IOError("projection {0} is not of type DOSS".format(proj_number))

                nrows, ncols, _, denergy, fermi = split_numbers(line[8:])
                # nrows, ncols = (int(nrows), int(ncols))

                if energy_delta is None:
                    energy_delta = denergy
                elif not energy_delta == denergy:
                    raise IOError("projection {0} has different delta energy ({1}) to previous ({2})".format(
                        proj_number, denergy, energy_delta
                    ))
                if fermi_energy is None:
                    fermi_energy = fermi
                elif not fermi_energy == fermi:
                    raise IOError("projection {0} has different fermi energy ({1}) to previous ({2})".format(
                        proj_number, fermi, fermi_energy
                    ))

                line = f.readline().strip()
                ienergy = split_numbers(line)[1]

                if initial_energy is None:
                    initial_energy = ienergy
                elif not initial_energy == ienergy:
                    raise IOError("projection {0} has different initial energy ({1}) to previous ({2})".format(
                        proj_number, ienergy, initial_energy
                    ))

                line = f.readline().strip()
                projid, norbitals, _, _, _, _ = [int(i) for i in line.split()]
                line = f.readline().strip()
                dos = []
                while not line.startswith('-%-') and line:
                    dos += split_numbers(line)
                    line = f.readline().strip()

                if len_dos is None:
                    len_dos = len(dos)
                elif not len_dos == len(dos):
                    raise IOError("projection {0} has different dos value lengths ({1}) to previous ({2})".format(
                        proj_number, len(dos), len_dos
                    ))

                if projid not in alpha_projections:
                    alpha_projections[projid] = {"id": projid, "norbitals": norbitals, "dos": dos}
                elif projid in beta_projections:
                    raise IOError("three data sets with same projid ({0}) were found".format(projid))
                else:
                    beta_projections[projid] = {"id": projid, "norbitals": norbitals, "dos": dos}
            else:
                line = f.readline().strip()

        system_type = _ihferm_map[int(system_type)]
        fermi_energy = {"units": "eV", "magnitude": float(fermi_energy) * codata[("Hartree", "eV")]}

        energy_delta = float(energy_delta) * codata[("Hartree", "eV")]
        initial_energy = float(initial_energy) * codata[("Hartree", "eV")]
        len_dos = int(len_dos)
        energies = {"units": "eV",
                    "magnitude": np.linspace(initial_energy, initial_energy + len_dos * energy_delta, len_dos).tolist()}

        total_alpha = None
        total_beta = None
        if alpha_projections:
            total_alpha = alpha_projections.pop(max(alpha_projections.keys()))
        if beta_projections:
            total_beta = beta_projections.pop(max(beta_projections.keys()))

        return {"energy": energies,
                "system_type": system_type,
                "fermi_energy": fermi_energy,
                "total_alpha": total_alpha,
                "total_beta": total_beta,
                "projections_alpha": list(alpha_projections.values()) if alpha_projections else None,
                "projections_beta": list(beta_projections.values()) if beta_projections else None,
                "creator": {"program": "CRYSTAL"}}


class ECH3OutPlugin(object):
    """ parser plugin for jsonextended
    """
    plugin_name = 'crystal_ech3_out'
    plugin_descript = 'read CRYSTAL14 charge density out data'
    file_regex = '*ech3.out'

    def read_file(self, file_obj, **kwargs):
        """
 DIRECT LATTICE VECTOR COMPONENTS (ANGSTROM)
         0.00000    2.71000    2.71000
         2.71000    0.00000    2.71000
         2.71000    2.71000    0.00000
...
 *******************************************************************************
   ATOM N.AT.  SHELL    X(A)      Y(A)      Z(A)      EXAD       N.ELECT.
 *******************************************************************************
    1   14 SI    4     0.678     0.678     0.678   1.233E-01      14.000
    2   14 SI    4    -0.677    -0.678    -0.677   1.233E-01      14.000
 *******************************************************************************
        """
        import pymatgen  # optional dependency

        sites, coords, nelect = [], [], []
        line = file_obj.readline()
        while line:
            if line.strip().startswith('DIRECT LATTICE VECTOR COMPONENTS'):
                avec = np.array(file_obj.readline().strip().split(), dtype=float)
                bvec = np.array(file_obj.readline().strip().split(), dtype=float)
                cvec = np.array(file_obj.readline().strip().split(), dtype=float)
            if line.strip().startswith('ATOM N.AT.'):
                line = file_obj.readline()
                line = file_obj.readline()
                while not line.strip().startswith('***'):
                    atom = line.strip().split()
                    sites.append(int(atom[1]))
                    nelect.append(float(atom[8]))
                    coords.append(np.array(atom[4:7], dtype=float))
                    line = file_obj.readline()
            line = file_obj.readline()

        lattice = pymatgen.Lattice([avec, bvec, cvec])
        struct = pymatgen.Structure(lattice, sites, coords,
                               to_unit_cell=True,
                               coords_are_cartesian=True,
                               site_properties={"nelect": nelect})

        return {
            "cell_vectors": {
                "a": {"units": "angstrom", "magnitude": avec.tolist()},
                "b": {"units": "angstrom", "magnitude": bvec.tolist()},
                "c": {"units": "angstrom", "magnitude": cvec.tolist()}
            },
            "atoms": {"ccoords": {"units": "angstrom",
                                  "magnitude": struct.cart_coords.tolist()},
                      "nuclear_charge": struct.site_properties["nelect"],
                      "atomic_number": struct.atomic_numbers
                      },
            "creator": {"program": "CRYSTAL"}
        }


class ECH3CubePlugin(object):
    """ parser plugin for jsonextended
    """
    plugin_name = 'crystal_ech3_cube'
    plugin_descript = 'read CRYSTAL charge density cube data'
    file_regex = '*ech3_dat.prop3d'

    def read_file(self, f, **kwargs):

        na, nb, nc = [int(i) for i in f.readline().strip().split()]
        origin = [float(i) for i in f.readline().strip().split()]
        da_vec = np.array([float(i) for i in f.readline().strip().split()])
        db_vec = np.array([float(i) for i in f.readline().strip().split()])
        dc_vec = np.array([float(i) for i in f.readline().strip().split()])
        # TODO na - 1 seems to work, but is it right?
        avec = (da_vec * (na - 1) * codata[("Bohr", "Angstrom")]).tolist()
        bvec = (db_vec * (nb - 1) * codata[("Bohr", "Angstrom")]).tolist()
        cvec = (dc_vec * (nc - 1) * codata[("Bohr", "Angstrom")]).tolist()

        name = f.readline().strip()
        densities = []
        assert name == 'Charge density'
        charge_density = []
        line = f.readline().strip().split()
        while line:
            if line[0] == 'Spin':
                break
            charge_density += [float(s) for s in line]
            line = f.readline().strip().split()
        cdense = np.array(charge_density).reshape((nc, nb, na)).transpose()
        densities.append({
            "type": "charge",
            "magnitude": cdense
        })
        spin_density = []
        line = f.readline().strip().split()
        while line:
            spin_density += [float(s) for s in line]
            line = f.readline().strip().split()
        if spin_density:
            sdense = np.array(spin_density).reshape((nc, nb, na)).transpose()
            densities.append({
                "type": "spin",
                "magnitude": sdense
            })

        return {
            "title": "CRYSTAL",
            # "na": na, "nb": nb, "nc": nc,
            "cell_vectors": {
                "a": {"units": "angstrom", "magnitude": avec},
                "b": {"units": "angstrom", "magnitude": bvec},
                "c": {"units": "angstrom", "magnitude": cvec}
            },
            # "centre": [0, 0, 0],
            "densities": densities,
            # "atoms": {"ccoords": {"units": "angstrom",
            #                       "magnitude": ccoords},
            #           #"nuclear_charge": nuclear_charges,
            #           "nuclear_charge": valence_charges,
            #           "symbols": symbols, "atomic_number": atomic_numbers},
            "creator": {"program": "CRYSTAL"}
        }


def split_output(lines):
    """ split up the crystal output into sections

    Parameters
    ----------
    lines: list of str

    Returns
    -------

    """
    start_line_no = None
    geom_input_end = None
    scf_init_start_no = None
    scf_init_end_no = None
    opt_start_no = None
    opt_end_no = None
    mulliken_starts = None
    final_opt = None
    run_warnings = []
    non_terminating_errors = []
    errors = []
    mpi_abort = False
    second_opt_line = False
    telapse_line = False
    band_gaps = None
    meta = {}
    for i, line in enumerate(lines):
        if "************************" in line and start_line_no is None:
            start_line_no = i
        elif "WARNING" in line.upper():
            run_warnings.append(line.strip())
        elif "ERROR" in line:
            errors.append(line.strip())
        elif "SCF abnormal end" in line:
            errors.append(line.strip())
        elif "MPI_Abort" in line:
            # only record one mpi_abort event (to not clutter output)
            if not mpi_abort:
                errors.append(line.strip())
                mpi_abort = True
        elif line.strip().startswith("* GEOMETRY EDITING"):
            if geom_input_end is None:
                geom_input_end = i
            else:
                raise IOError("found two lines starting '* GEOMETRY EDITING' in initial data:"
                              " {0} and {1}".format(geom_input_end, i))

        elif "CRYSTAL - SCF - TYPE OF CALCULATION :" in line:
            if opt_start_no is not None:
                continue
            if scf_init_start_no is not None:
                raise IOError("found two lines starting scf ('CRYSTAL - SCF - ') in initial data:"
                              " {0} and {1}".format(scf_init_start_no, i))
            scf_init_start_no = i
        elif "SCF ENDED" in line:
            if opt_start_no is not None:
                continue
            if "CONVERGE" not in line:
                non_terminating_errors.append(line.strip())
            if scf_init_end_no is not None:
                raise IOError("found two lines ending scf ('SCF ENDED') in initial data:"
                              " {0} and {1}".format(scf_init_end_no, i))
            scf_init_end_no = i
        # elif "STARTING GEOMETRY OPTIMIZATION" in line: #not the same in CRYSTAL17
        elif "OPTOPTOPTOPT" in line:
            if opt_start_no is not None:
                if second_opt_line:
                    raise IOError("found two lines starting opt ('STARTING GEOMETRY OPTIMIZATION'):"
                                  " {0} and {1}".format(opt_start_no, i))
                else:
                    second_opt_line = True
            opt_start_no = i
        elif "CONVERGENCE ON GRADIENTS SATISFIED AFTER THE FIRST OPTIMIZATION CYCLE" in line:
            if opt_start_no is not None:
                if second_opt_line:
                    raise IOError("found two lines starting opt ('STARTING GEOMETRY OPTIMIZATION'):"
                                  " {0} and {1}".format(opt_start_no, i))
                else:
                    second_opt_line = True
            opt_start_no = i
        elif "OPT END -" in line:
            if opt_end_no is not None:
                raise IOError("found two lines ending opt ('OPT END -'):"
                              " {0} and {1}".format(opt_end_no, i))
            opt_end_no = i
        elif opt_end_no and "BAND GAP" in line:
            # NB: this is new for CRYSTAL17
            band_gaps = {} if band_gaps is None else band_gaps
            if fnmatch(line.strip(), "ALPHA BAND GAP:*eV"):
                bgvalue = split_numbers(line)[0]
                bgtype = "alpha"
            elif fnmatch(line.strip(), "BETA BAND GAP:*eV"):
                bgvalue = split_numbers(line)[0]
                bgtype = "beta"
            elif fnmatch(line.strip(), "BAND GAP:*eV"):
                bgvalue = split_numbers(line)[0]
                bgtype = "all"
            else:
                raise IOError("found a band gap of unknown format at line {0}: {1}".format(i, line))
            if bgtype in band_gaps:
                raise IOError("band gap data already contains {0} value before line {1}: {2}".format(bgtype, i, line))
            band_gaps[bgtype] = {"magnitude": bgvalue, "units": "eV"}
        elif "CONVERGENCE TESTS UNSATISFIED" in line.upper():
            non_terminating_errors.append(line.strip())
        elif line.strip().startswith("MULLIKEN POPULATION ANALYSIS"):
            # can have ALPHA+BETA ELECTRONS and ALPHA-BETA ELECTRONS (denoted in line above mulliken_starts)
            if mulliken_starts is None:
                mulliken_starts = [i]
            else:
                mulliken_starts.append(i)
        elif "FINAL OPTIMIZED GEOMETRY" in line:
            if final_opt is not None:
                raise IOError("found two lines starting final opt geometry ('FINAL OPTIMIZED GEOMETRY'):"
                              " {0} and {1}".format(final_opt, i))
            final_opt = i
        elif "TELAPSE" in line:
            telapse_line = i

    if telapse_line:
        total_seconds = int(split_numbers(lines[telapse_line].split("TELAPSE")[1])[0])
        m, s = divmod(total_seconds, 60)
        h, m = divmod(m, 60)
        meta["elapsed_time"] = "%d:%02d:%02d" % (h, m, s)

    if errors:
        return (start_line_no, geom_input_end, scf_init_start_no, scf_init_end_no, opt_start_no, opt_end_no,
                mulliken_starts, final_opt, run_warnings, non_terminating_errors, errors, meta, band_gaps)

    if start_line_no is None:
        raise IOError("couldn't find start of program run (denoted *****)")
    if geom_input_end is None:
        raise IOError("couldn't find end of geometry input (denoted * GEOMETRY EDITING)")
    if scf_init_start_no is None:
        raise IOError("didn't find an SCF (as expected)")
    if scf_init_start_no is not None and scf_init_end_no is None:
        raise IOError("found start of scf but not end")
    if scf_init_end_no is not None and scf_init_start_no is None:
        raise IOError("found end of scf but not start")
    if opt_start_no is not None and opt_end_no is None:
        raise IOError("found start of optimisation but not end")
    if opt_end_no is not None and opt_start_no is None:
        raise IOError("found end of optimisation but not start")

    return (start_line_no, geom_input_end, scf_init_start_no, scf_init_end_no, opt_start_no, opt_end_no,
            mulliken_starts, final_opt, run_warnings, non_terminating_errors, errors, meta, band_gaps)


def read_before_run(lines):
    """ read metadata output in the file before the actual program start

    Parameters
    ----------
    lines: list of str

    Returns
    -------

    """
    meta = {}
    for line in lines:
        line = line.strip()
        if fnmatch(line, "date:*"):
            meta["date"] = line.replace("date:", "").strip()
        # TODO fo this better; perhaps "16 PROCESSORS WORKING" just above start of run
        if fnmatch(line, "resources_used.ncpus =*"):
            meta["nprocs"] = int(line.replace("resources_used.ncpus =", ""))
    return meta


def get_geometry(dct, i, line, lines, startline=0):
    """ update dict get geometry related variables

    Parameters
    ----------
    dct
    i
    line
    lines
    startline: int

    Returns
    -------

    """
    if fnmatch(line, "LATTICE PARAMETERS*(*)"):
        if not ("ANGSTROM" in line and "DEGREES" in line):
            raise IOError("was expecting lattice parameters in angstroms and degrees on line:"
                          " {0}, got: {1}".format(startline + i, line))

    for pattern, field, pattern2 in [('PRIMITIVE*CELL*', "primitive_cell", "ATOMS IN THE ASYMMETRIC UNIT*"),
                                     ('CRYSTALLOGRAPHIC*CELL*', "crystallographic_cell",
                                      "COORDINATES IN THE CRYSTALLOGRAPHIC CELL")]:
        if fnmatch(line, pattern):
            if not fnmatch(lines[i + 1].strip(), "A*B*C*ALPHA*BETA*GAMMA"):
                raise IOError("was expecting A B C ALPHA BETA GAMMA on line:"
                              " {0}, got: {1}".format(startline + i + 1, lines[i + 1]))
            dct[field] = edict.merge([dct.get(field, {}),
                                      {"cell_parameters":
                                           dict(zip(['a', 'b', 'c', 'alpha', 'beta', 'gamma'],
                                                    [{"units": "angstrom" if i < 3 else "degrees",
                                                      "magnitude": p} for i, p in enumerate(split_numbers(lines[i + 2]))
                                                     ]))}])
        if fnmatch(line, pattern2):
            periodic = [True, True, True]
            if not fnmatch(lines[i + 1].strip(), "ATOM*X/A*Y/B*Z/C"):
                # for 2d (slab) can get z in angstrom (and similar for 1d)
                if fnmatch(lines[i + 1].strip(), "ATOM*X/A*Y/B*Z(ANGSTROM)*"):
                    periodic = [True, True, False]
                elif fnmatch(lines[i + 1].strip(), "ATOM*X/A*Y(ANGSTROM)*Z(ANGSTROM)*"):
                    periodic = [True, False, False]
                elif fnmatch(lines[i + 1].strip(), "ATOM*X(ANGSTROM)*Y(ANGSTROM)*Z(ANGSTROM)*"):
                    periodic = [False, False, False]
                    cell_params = dict(zip(['a', 'b', 'c', 'alpha', 'beta', 'gamma'],
                                           [{"units": "angstrom" if i < 3 else "degrees",
                                             "magnitude": p} for i, p in enumerate([500., 500., 500., 90., 90., 90.])]))
                    dct[field] = edict.merge([dct.get(field, {}), {"cell_parameters": cell_params}])
                else:
                    raise IOError("was expecting ATOM X Y Z (in units of ANGSTROM or fractional) on line:"
                                  " {0}, got: {1}".format(startline + i + 1, lines[i + 1]))
            if not all(periodic) and "cell_parameters" not in dct.get(field, {}):
                raise IOError("require cell parameters to have been set for non-periodic directions in line"
                              " #{0} : {1}".format(startline + i + 1, lines[i + 1]))
            a, b, c, alpha, beta, gamma = [None] * 6
            if not all(periodic):
                cell = dct[field]["cell_parameters"]
                a, b, c, alpha, beta, gamma = [cell[p]["magnitude"] for p in ["a", "b", "c", "alpha", "beta", "gamma"]]

            nextindx = i + 3
            atom_data = {'ids': [], 'assymetric': [], 'atomic_numbers': [], 'symbols': [], "fcoords": []}
            atom_data["pbc"] = periodic
            while lines[nextindx].strip() and not lines[nextindx].strip()[0].isalpha():
                fields = lines[nextindx].strip().split()
                atom_data['ids'].append(fields[0])
                atom_data['assymetric'].append(bool(strtobool(fields[1])))
                atom_data['atomic_numbers'].append(int(fields[2]))
                atom_data['symbols'].append(fields[3].lower().capitalize())
                if all(periodic):
                    atom_data['fcoords'].append([float(fields[4]), float(fields[5]), float(fields[6])])
                elif periodic == [True, True, False] and alpha == 90 and beta == 90:
                    atom_data['fcoords'].append([float(fields[4]), float(fields[5]), float(fields[6]) / c])
                # TODO other periodic types (1D, 0D)
                nextindx += 1

            if not atom_data["fcoords"]:
                atom_data.pop("fcoords")
            dct[field] = edict.merge([dct.get(field, {}), atom_data])

    # TODO These ccoords DON'T work with lattice parameters (at least for final run)
    if fnmatch(line, "CARTESIAN COORDINATES - PRIMITIVE CELL*"):
        if not fnmatch(lines[i + 2].strip(), "*ATOM*X(ANGSTROM)*Y(ANGSTROM)*Z(ANGSTROM)"):
            raise IOError("was expecting ATOM X(ANGSTROM) Y(ANGSTROM) Z(ANGSTROM) on line:"
                          " {0}, got: {1}".format(startline + i + 2, lines[i + 2]))

        nextindx = i + 4
        atom_data = {'ids': [],
                     'atomic_numbers': [],
                     'symbols': [],
                     "ccoords": {"units": "angstrom", "magnitude": []}
                     }
        while lines[nextindx].strip() and not lines[nextindx].strip()[0].isalpha():
            fields = lines[nextindx].strip().split()
            atom_data['ids'].append(fields[0])
            atom_data['atomic_numbers'].append(int(fields[1]))
            atom_data['symbols'].append(fields[2].lower().capitalize())
            atom_data['ccoords']["magnitude"].append([float(fields[3]), float(fields[4]), float(fields[5])])
            nextindx += 1
        dct["primitive_cell"] = edict.merge([dct.get("primitive_cell", {}), atom_data])

    if fnmatch(line, "DIRECT LATTICE VECTORS CARTESIAN COMPONENTS*"):
        if "ANGSTROM" not in line:
            raise IOError("was expecting lattice vectors in angstroms on line:"
                          " {0}, got: {1}".format(startline + i, line))
        if not fnmatch(lines[i + 1].strip(), "X*Y*Z"):
            raise IOError("was expecting X Y Z on line:"
                          " {0}, got: {1}".format(startline + i + 1, lines[i + 1]))
        if "crystallographic_cell" not in dct:
            dct["crystallographic_cell"] = {}
        if "cell_vectors" in dct["crystallographic_cell"]:
            raise IOError("found multiple cell vectors on line:"
                          " {0}, got: {1}".format(startline + i + 1, lines[i + 1]))
        vectors = {"a": {"units": "angstrom", "magnitude": split_numbers(lines[i + 2])},
                   "b": {"units": "angstrom", "magnitude": split_numbers(lines[i + 3])},
                   "c": {"units": "angstrom", "magnitude": split_numbers(lines[i + 4])}}

        dct["primitive_cell"]["cell_vectors"] = vectors


def get_symmetry(dct, i, line, lines, startline=0):
    """ update dict with symmetry related variables

    Parameters
    ----------
    dct
    i
    line
    lines
    startline: int

    Returns
    -------

    """
    if fnmatch(line, "*SYMMOPS - TRANSLATORS IN FRACTIONAL UNITS*"):
        nums = split_numbers(line)
        if not len(nums) == 1:
            raise IOError("was expecting a single number, representing the number of symmops, on this line:"
                          " {0}, got: {1}".format(startline + i, line))
        nsymmops = int(nums[0])
        if not fnmatch(lines[i + 1], "*MATRICES AND TRANSLATORS IN THE CRYSTALLOGRAPHIC REFERENCE FRAME*"):
            raise IOError("was expecting CRYSTALLOGRAPHIC REFERENCE FRAME on this line"
                          " {0}, got: {1}".format(startline + i + 1, lines[i + 1].strip()))
        if not fnmatch(lines[i + 2], "*V*INV*ROTATION MATRICES*TRANSLATORS*"):
            raise IOError("was expecting symmetry headers on this line"
                          " {0}, got: {1}".format(startline + i + 2, lines[i + 2].strip()))
        symmops = []
        for j in range(nsymmops):
            values = split_numbers(lines[i + 3 + j])
            if not len(values) == 14:
                raise IOError("was expecting 14 values for symmetry data on this line"
                              " {0}, got: {1}".format(startline + i + 3 + j, lines[i + 3 + j].strip()))
            symmops.append(values[2:14])
        dct["primitive_symmops"] = symmops


def read_init(lines, startline):
    """ read initial setup data (starting after intital geometry input)

    Parameters
    ----------
    lines: List of str
    startline: int

    Returns
    -------

    """
    init = {"calculation": {"spin": False}}
    for i, line in enumerate(lines):
        line = line.strip()
        if line.startswith("TYPE OF CALCULATION :"):
            init["calculation"]["type"] = line.replace("TYPE OF CALCULATION :", "").strip().lower()
            if "HAMILTONIAN" in lines[i + 1]:
                regex = r"\(EXCHANGE\)\[CORRELATION\] FUNCTIONAL:\((.*)\)\[(.*)\]"
                string = lines[i + 3].strip()
                if re.match(regex, string):
                    init["calculation"]["functional"] = {"exchange": re.search(regex, string).group(1),
                                                         "correlation": re.search(regex, string).group(2)}
        if "SPIN POLARIZ" in line:
            init["calculation"]["spin"] = True

        get_geometry(init, i, line, lines, startline)
        get_symmetry(init, i, line, lines, startline)

    return init


def read_scf(lines, startline):
    """ read scf data

    Parameters
    ----------
    lines: List of str
    startline: int

    Returns
    -------

    """
    scf = []
    scf_cyc = None
    last_cyc_num = None
    for i, line in enumerate(lines):
        line = line.strip()

        if fnmatch(line, "CYC*"):

            # start new cycle
            if scf_cyc is not None:
                scf.append(scf_cyc)
            scf_cyc = {}

            # check we are adding them in sequential order
            cur_cyc_num = split_numbers(line)[0]
            if last_cyc_num is not None:
                if cur_cyc_num != last_cyc_num + 1:
                    raise IOError("was expecting the SCF cyle number to be "
                                  "{0} in line: {1}".format(int(last_cyc_num + 1), line))
            last_cyc_num = cur_cyc_num

            if fnmatch(line, "*ETOT*"):
                if not fnmatch(line, "*ETOT(AU)*"):
                    raise IOError("was expecting units in a.u. on line:"
                                  " {0}, got: {1}".format(startline + i, line))
                # this is the initial energy of the configuration and so actually the energy of the previous run
                if scf:
                    scf[-1]["energy"] = scf[-1].get("energy", {})
                    scf[-1]["energy"]["total"] = {"magnitude": split_numbers(line)[1] * codata[("Hartree", "eV")],
                                                  "units": "eV"}

        elif scf_cyc is None:
            continue

        # The total magnetization is the integral of the magnetization in the cell:
        #     MT=∫ (nup-ndown) d3 r
        #
        # The absolute magnetization is the integral of the absolute value of the magnetization in the cell:
        #     MA=∫ |nup-ndown| d3 r
        #
        # In a simple ferromagnetic material they should be equal (except possibly for an overall sign).
        # In simple antiferromagnets (like FeO) MT is zero and MA is twice the magnetization of each of the two atoms.

        if line.startswith("CHARGE NORMALIZATION FACTOR"):
            scf_cyc["CHARGE NORMALIZATION FACTOR".lower().replace(" ", "_")] = split_numbers(line)[0]
        if line.startswith("SUMMED SPIN DENSITY"):
            scf_cyc["spin_density_total"] = split_numbers(line)[0]

        if line.startswith("TOTAL ATOMIC CHARGES"):
            scf_cyc["atomic_charges_peratom"] = []
            j = i + 1
            while len(lines[j].strip().split()) == len(split_numbers(lines[j])):
                scf_cyc["atomic_charges_peratom"] += split_numbers(lines[j])
                j += 1
        if line.startswith("TOTAL ATOMIC SPINS"):
            scf_cyc["spin_density_peratom"] = []
            j = i + 1
            while len(lines[j].strip().split()) == len(split_numbers(lines[j])):
                scf_cyc["spin_density_peratom"] += split_numbers(lines[j])
                j += 1
            scf_cyc["spin_density_absolute"] = sum([abs(s) for s in split_numbers(lines[i + 1])])

    # add last scf cycle
    if scf_cyc:
        scf.append(scf_cyc)

    return scf


def read_post_scf(lines, startline):
    """ read post initial scf data

    Parameters
    ----------
    lines: list of str
    startline: int

    Returns
    -------

    """
    post_scf = {}
    for i, line in enumerate(lines):
        if fnmatch(line.strip(), "TOTAL ENERGY*DE*"):
            if not fnmatch(line.strip(), "TOTAL ENERGY*AU*DE*"):
                raise IOError("was expecting units in a.u. on line:"
                              " {0}, got: {1}".format(startline + i, line))
            post_scf["energy"] = post_scf.get("energy", {})
            if "total_corrected" in post_scf["energy"]:
                raise IOError("total corrected energy found twice, on line:"
                              " {0}, got: {1}".format(startline + i, line))
            post_scf["energy"]["total_corrected"] = {"magnitude": split_numbers(line)[1] * codata[("Hartree", "eV")],
                                                     "units": "eV"}

    return post_scf


def read_opt(lines, startline):
    """ read geometric optimisation

    Parameters
    ----------
    lines: list of str
    startline: int

    Returns
    -------

    """
    if "CONVERGENCE ON GRADIENTS SATISFIED AFTER THE FIRST OPTIMIZATION CYCLE" in lines[0]:
        if "OPT END -" not in lines[-1]:
            raise IOError("expecting OPT END in line {0}: {1}".format(startline + len(lines), lines[-1]))
        if not fnmatch(lines[-1], "*E(AU)*"):
            raise IOError("was expecting units in a.u. on line:"
                          " {0}, got: {1}".format(startline + len(lines), lines[-1]))
        return [{"energy": {"total_corrected": {"magnitude": split_numbers(lines[-1])[0] * codata[("Hartree", "eV")],
                                                "units": "eV"}}}]

    opt = []
    opt_cyc = None
    scf_start_no = None
    failed_opt_step = False

    for i, line in enumerate(lines):
        if i == 0:
            continue
        line = line.strip()
        if fnmatch(line, "*OPTIMIZATION*POINT*"):
            if opt_cyc is not None and not failed_opt_step:
                opt.append(opt_cyc)
            opt_cyc = {}
            scf_start_no = None
            failed_opt_step = False
        elif opt_cyc is None:
            continue

        # when using ONELOG optimisation key word
        if "CRYSTAL - SCF - TYPE OF CALCULATION :" in line:
            if scf_start_no is not None:
                raise IOError("found two lines starting scf ('CRYSTAL - SCF - ') in opt step {0}:".format(len(opt))
                              + " {0} and {1}".format(scf_start_no, i))
            scf_start_no = i
        elif "SCF ENDED" in line:
            if "CONVERGE" not in line:
                pass  # errors.append(line.strip())
            opt_cyc["scf"] = read_scf(lines[scf_start_no + 1:i + 1], startline + i + 1)

        get_geometry(opt_cyc, i, line, lines, startline)

        # TODO move to read_post_scf?
        if fnmatch(line, "TOTAL ENERGY*DE*"):
            if not fnmatch(line, "TOTAL ENERGY*AU*DE*AU*"):
                raise IOError("was expecting units in a.u. on line:"
                              " {0}, got: {1}".format(startline + i, line))
            opt_cyc["energy"] = opt_cyc.get("energy", {})
            opt_cyc["energy"]["total_corrected"] = {"magnitude": split_numbers(line)[1] * codata[("Hartree", "eV")],
                                                    "units": "eV"}

        for param in ["MAX GRADIENT", "RMS GRADIENT", "MAX DISPLAC", "RMS DISPLAC"]:
            if fnmatch(line, "{}*CONVERGED*".format(param)):
                if "convergence" not in opt_cyc:
                    opt_cyc["convergence"] = {}
                opt_cyc["convergence"][param.lower().replace(" ", "_")] = bool(strtobool(line.split()[-1]))

        if fnmatch(line, "*SCF DID NOT CONVERGE. RETRYING WITH A SMALLER OPT STEP*"):
            # TODO add failed optimisation steps with dummy energy and extra parameter?
            # for now discard this optimisation step
            failed_opt_step = True

    if opt_cyc and not failed_opt_step:
        opt.append(opt_cyc)

    return opt


def read_final(lines, startline):
    """ read final setup data

    Parameters
    ----------
    lines: list of str
    startline: int

    Returns
    -------

    """
    final = {}
    for i, line in enumerate(lines):
        line = line.strip()
        get_geometry(final, i, line, lines, startline)
        get_symmetry(final, i, line, lines, startline)

    return final


def read_mulliken(lines, mulliken_indices):
    """

    Parameters
    ----------
    lines: List of str
    mulliken_indices: List of int

    Returns
    -------

    """
    mulliken = {}

    for i, indx in enumerate(mulliken_indices):
        name = lines[indx - 1].strip().lower()
        if not (name == "ALPHA+BETA ELECTRONS".lower() or name == "ALPHA-BETA ELECTRONS".lower()):
            raise IOError("was expecting mulliken to be alpha+beta or alpha-beta on line:"
                          " {0}, got: {1}".format(indx - 1, lines[indx - 1]))

        mulliken[name.replace(" ", "_")] = {"ids": [], "symbols": [], "atomic_numbers": [], "charges": []}

        if len(mulliken_indices) > i + 1:
            searchlines = lines[indx + 1:mulliken_indices[i + 1]]
        else:
            searchlines = lines[indx + 1:]
        charge_line = None
        for j, line in enumerate(searchlines):
            if fnmatch(line.strip(), "*ATOM*Z*CHARGE*SHELL*POPULATION*"):
                charge_line = j + 2
                break
        if charge_line is None:
            continue

        while searchlines[charge_line].strip() and not searchlines[charge_line].strip()[0].isalpha():
            fields = searchlines[charge_line].strip().split()
            # shell population can wrap multiple lines, the one we want has the label in it
            if len(fields) != len(split_numbers(searchlines[charge_line])):
                mulliken[name.replace(" ", "_")]["ids"].append(int(fields[0]))
                mulliken[name.replace(" ", "_")]["symbols"].append(fields[1].lower().capitalize())
                mulliken[name.replace(" ", "_")]["atomic_numbers"].append(int(fields[2]))
                mulliken[name.replace(" ", "_")]["charges"].append(float(fields[3]))

            charge_line += 1

    return mulliken


class CrystalOutputPlugin(object):
    """ a class for parsing Crystal DFT Simulation Data

    """

    plugin_name = 'crystal.out'
    plugin_descript = 'read main crystal output file'
    file_regex = '*crystal.out'

    def read_file(self, file, log_warnings=True, **kwargs):

        lines = file.read().splitlines()

        (start_line_no, geom_input_end, scf_init_start_no, scf_init_end_no, opt_start_no, opt_end_no, mulliken_starts,
         final_opt, run_warnings, non_terminating_errors, errors, meta, band_gaps) = split_output(lines)

        if run_warnings and log_warnings:
            logger.warning("the following warnings were noted:\n  {}".format("\n  ".join(run_warnings)))

        errors_all = errors + non_terminating_errors
        if errors or non_terminating_errors:
            ## MPI aborts should be the result of another error, so remove them if this is the case
            errors_noabort = [e for e in errors_all if "MPI_Abort" not in e]
            errors_all = errors_noabort if errors_noabort else errors_all
            logger.warning("the following errors were found:\n  {}".format("\n  ".join(errors_all)))

        if start_line_no is None or start_line_no is None:
            initial = None
        else:
            initial = read_init(lines[geom_input_end:scf_init_start_no], geom_input_end)
            if scf_init_start_no is None or scf_init_end_no is None:
                initial["scf"] = None

            else:
                #  initial["scf_type"] = lines[scf_start_no].replace("CRYSTAL - SCF - TYPE OF CALCULATION :", "").strip(),
                initial["scf"] = read_scf(lines[scf_init_start_no + 1:scf_init_end_no + 1], scf_init_start_no + 1)
                if opt_start_no is not None:
                    initial = edict.merge([initial, read_post_scf(lines[scf_init_end_no + 1:opt_start_no],
                                                                  scf_init_start_no + 1)])
                elif final_opt is not None:
                    initial = edict.merge([initial, read_post_scf(lines[scf_init_end_no + 1:final_opt],
                                                                  scf_init_start_no + 1)])
                elif mulliken_starts is not None:
                    initial = edict.merge([initial, read_post_scf(lines[scf_init_end_no + 1:mulliken_starts[0]],
                                                                  scf_init_start_no + 1)])
                else:
                    initial = edict.merge([initial, read_post_scf(lines[scf_init_end_no + 1:], scf_init_start_no + 1)])

        # if errors:
        #     return {"warnings": run_warnings,
        #             "errors": errors_all,
        #             "meta": None if start_line_no is None else read_before_run(lines[:start_line_no]),
        #             "initial": initial,
        #             "creator": {"program": "Crystal14"}
        #             }

        output = {"warnings": run_warnings,
                  "errors": errors_all,
                  "meta": meta if start_line_no is None else edict.merge(
                      [meta, read_before_run(lines[:start_line_no])]),
                  "initial": initial,
                  "creator": {"program": "CRYSTAL"}
                  }

        if opt_start_no is not None and opt_end_no is not None:
            if errors:
                try:
                    output["optimisation"] = read_opt(lines[opt_start_no:opt_end_no + 1], opt_start_no)
                except:
                    pass
            else:
                output["optimisation"] = read_opt(lines[opt_start_no:opt_end_no + 1], opt_start_no)
        else:
            output["optimisation"] = None

        if final_opt is not None:
            if errors:
                try:
                    output["final"] = read_final(lines[final_opt:], final_opt)
                except:
                    pass
            else:
                output["final"] = read_final(lines[final_opt:], final_opt)
        else:
            output["final"] = {}
        if band_gaps is not None:
            output["final"]["band_gaps"] = band_gaps
        # output["final"] = output["final"] if output["final"] else None

        # we make sure that the final section holds the final energy and primitive geometry
        if "primitive_cell" not in output["final"]:
            if output["optimisation"] is not None:
                output["final"]["primitive_cell"] = copy.deepcopy(
                    output["optimisation"][-1].get("primitive_cell", None))
            else:
                output["final"]["primitive_cell"] = copy.deepcopy(output["initial"].get("primitive_cell", None))
        if "energy" not in output["final"]:
            if output["optimisation"] is not None:
                output["final"]["energy"] = copy.deepcopy(output["optimisation"][-1].get("energy", None))
            else:
                output["final"]["energy"] = copy.deepcopy(output["initial"].get("energy", None))

        if "primitive_symmops" not in output["final"]:
            if "primitive_symmops" in output["initial"] and output["optimisation"] is None:
                output["final"]["primitive_symmops"] = copy.deepcopy(output["initial"]["primitive_symmops"])

        if mulliken_starts is not None:
            if errors:
                try:
                    output["mulliken"] = read_mulliken(lines, mulliken_starts)
                except:
                    pass
            else:
                output["mulliken"] = read_mulliken(lines, mulliken_starts)

        return output


class CrystalSCFLogPlugin(object):
    """ a class for parsing Crystal DFT Simulation Data

    """

    plugin_name = 'crystal.scflog'
    plugin_descript = 'read crystal .SCFLOG file, created when running an optimisation'
    file_regex = '*crystal.scflog'

    def read_file(self, file, *args, **kwargs):

        lines = file.read().splitlines()
        opt = []
        scf_start_no = None
        for i, line in enumerate(lines):
            if "CRYSTAL - SCF - TYPE OF CALCULATION :" in line:
                if scf_start_no is not None:
                    raise IOError("found two lines starting scf ('CRYSTAL - SCF - ') in opt step {0}:".format(len(opt))
                                  + " {0} and {1}".format(scf_start_no, i))
                scf_start_no = i
            elif "SCF ENDED" in line:
                if "CONVERGE" not in line:
                    pass  # errors.append(line.strip())
                opt.append({"scf": read_scf(lines[scf_start_no + 1:i + 1], i + 1)})
                scf_start_no = None

        return {"optimisation": opt}
