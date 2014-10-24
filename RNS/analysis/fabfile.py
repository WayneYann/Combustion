"""Fabric (fabfile.org) tasks for RNS.  See README."""

import collections
import glob
import itertools
import os.path
import pickle
import jinja2

import numpy as np

from fabric.api import *
from fabric.colors import *
from fabric.utils import *
from fabric.contrib.files import exists

DiffSameDomainRefined = {
    '1d': os.path.abspath('../../../BoxLib/Tools/C_util/Convergence/DiffSameDomainRefined1d.Linux.g++.gfortran.ex'),
    '2d': os.path.abspath('../../../BoxLib/Tools/C_util/Convergence/DiffSameDomainRefined2d.Linux.g++.gfortran.ex'),
}

ErrorTuple = collections.namedtuple('ErrorTuple', [ 'method', 'max_level', 'nx', 'cfl', 'level', 'norm', 'error' ])

#### dme tests #####################################################################################

@task
def dme_flameball_convergence(force_run=False):

    force_run = False

    base   = os.path.abspath('../bin/FlameBall')
    exe    = os.path.join(base, 'RNS2d.SDC.Linux.gcc.gfortran.DEBUG.OMP.ex')
    inputs = os.path.join(base, 'inputs.dme')
    probin = os.path.join(base, 'probin.dme')
    diff   = DiffSameDomainRefined['2d']

    stop_time = 20.0e-9

    dts = [ 0.5e-9, 1.0e-9, 2.0e-9, 4.0e-9 ]

    with cd(base):
        run("make -j 8")

    for dt in dts:
        dname = os.path.join(base, 'dt_%.2g' % dt)
        args  = 'stop_time=%g rns.fixed_dt=%g mlsdc.max_iters=8' % (stop_time, dt)
        if not exists(dname) or force_run:
            RNS(dname, exe, inputs, args, probin)

    diffs = []
    for crse, fine in zip(dts[:-1], dts[1:]):
        puts(green("comparing dt %s and %s" % (crse, fine)))
        dname1 = 'dt_%.2g' % crse
        dname2 = 'dt_%.2g' % fine
        plot1  = sorted(glob.glob(os.path.join(base, dname1, "plt*")))[-1]
        plot2  = sorted(glob.glob(os.path.join(base, dname2, "plt*")))[-1]

        out = run(" ".join([diff, "infile1="+plot1, "reffile="+plot2, "norm=0"]))
        diffs.append(parse_diff(out))

    header = parse_first_header(os.path.join(base, dname))

    puts(green("convergence rates"))
    for fine, crse in zip(diffs[:-1], diffs[1:]):
        for lev in crse:
            e1 = np.asarray(crse[lev])
            e2 = np.asarray(fine[lev])

            rate = list(np.log10(e1/e2)/np.log10(2))
            puts(yellow("|-"), False)
            puts(yellow("| species | rate | err |"), False)
            puts(yellow("|-"), False)
            for species, rate, err in zip(header, rate, list(e1)):
                puts(yellow("| %30s | %6.2f | %12.2e |" % (species, rate, err)), False)
            puts(yellow("|-"), False)

            # puts(yellow(e2/e1), False)
            # puts(yellow(np.log10(e2/e1)/np.log10(2)), False)


#### mach2bump test ################################################################################

@task
def mach2bump(trial, force_run=False, run_reference=True):

    diff = os.path.abspath('../../../BoxLib/Tools/C_util/Convergence/DiffSameDomain1d.Linux.g++.gfortran.DEBUG.ex')
    forse_run = as_boolean(force_run)
    run_reference = as_boolean(run_reference)

    base   = os.path.abspath('../bin/Mach2Bump')
    inputs = os.path.join(base, 'inputs-bndry-test')
    exe = {
        'rk':  os.path.join(base, 'RNS1d.Linux.gcc.gfortran.ex'),
        'sdc': os.path.join(base, 'RNS1d.SDC.Linux.gcc.gfortran.ex'),
    }

    if not exists(exe['rk']):
        build(base, ["DIM=1", "USE_SDCLIB=FALSE", "DEBUG=FALSE", "USE_OMP=FALSE", "USE_MPI=FALSE"])

    if not exists(exe['sdc']):
        build(base, ["DIM=1", "USE_SDCLIB=TRUE", "DEBUG=FALSE", "USE_OMP=FALSE", "USE_MPI=FALSE"])

    name = lambda method, max_level, nx, cfl: os.path.join(
        base, trial + '.d', '%s_lev%d_nx%04d_cfl%0.2f' % (method, max_level, nx, cfl))

    nxs  = [ 64, 128, 256, 512 ]
    cfls = [ 0.2, 0.35, 0.5 ]

    # high-res reference run
    cfl = 0.25
    nx  = 4096
    reference = name('ref', 0, nx, cfl)
    if run_reference:
        if not exists(reference) or force_run:
            RNS(reference, exe['sdc'], inputs, "rns.cfl=%f amr.max_level=0 amr.n_cell=%d" % (cfl, nx))

    # sdc/rk runs
    for method, nx, cfl, max_level in itertools.product([ 'sdc', 'rk' ], nxs, cfls, [ 0, 1 ] ):
        dname = name(method, max_level, nx, cfl)
        if not exists(dname) or force_run:
            blocking_factor = nx/8
            max_grid_size = nx/8
            args = [ "rns.cfl=%.2f" % cfl,
                     "amr.n_cell=%d" % nx,
                     "amr.max_level=%d" % max_level ]
            RNS(dname, exe[method], inputs, " ".join(args))

    errors = []
    for method, nx, cfl, max_level, norm in itertools.product(
            [ 'rk', 'sdc' ], nxs, cfls, [ 0, 1 ], [ 0, 1, 2 ]):
        error = diff_reference([ name(method, max_level, nx, cfl) ], reference, diff, norm=norm)
        for level in range(max_level+1):
            errors.append(ErrorTuple(method, max_level, nx, cfl, level, norm, error[0][level][0]))

    with open('errors_%s.pkl' % trial, 'w') as f:
        pickle.dump(errors, f)


class NoRepeat:
    def __init__(self):
        self.last = None
    def __call__(self, new):
        if self.last == new:
            return ""
        else:
            self.last = new
        return new

@task
def mach2bump_results(trial):

    tex_table = jinja2.Template(r"""
\begin{table}
    \tbl{ {{caption}} }{
\begin{tabular}{lcccc} \toprule
  Scheme & CFL & Cells & Error & Rate \\ \midrule
{%- for r in rows %}
  {{r.method}} & {{r.cfl}} & {{r.nx}} & {{r.error}} & {{r.rate}} \\
{%- endfor %}
    \bottomrule
    \end{tabular}%
}
\end{table}
""")

    org_table = jinja2.Template(r"""
| scheme | cfl | cells | error | rate |
|-
{%- for r in rows %}
| {{r.method}} | {{r.cfl}} | {{r.nx}} | {{r.error}} | {{r.rate}} |
{%- endfor %}
""")

    from math import log10

    with open('errors_%s.pkl' % trial, 'r') as f:
        errors = pickle.load(f)

    unzip = lambda x: zip(*xy)

    cfls       = set([ x.cfl for x in errors ])
    nxs        = set([ x.nx for x in errors ])
    methods    = set([ x.method for x in errors ])
    max_levels = set([ x.max_level for x in errors ])
    norms      = set([ x.norm for x in errors ])

    norm = 0

    rows = []

    mfilt = NoRepeat()
    cfilt = NoRepeat()

    for method in methods:
        for cfl in cfls:

            l0 = sorted([ (x.nx, x.error) for x in errors
                          if x.method == method
                          and x.cfl == cfl
                          and x.level == 0
                          and x.max_level == 1
                          and x.norm == norm ])

            l1 = sorted([ (x.nx, x.error) for x in errors
                          if x.method == method
                          and x.cfl == cfl
                          and x.level == 1
                          and x.max_level == 1
                          and x.norm == norm ])

            rate0 = [''] + [ log10(l0[i-1][1] / l0[i][1]) / log10(2.0) for i in range(1, len(l0)) ]
            rate1 = [''] + [ log10(l1[i-1][1] / l1[i][1]) / log10(2.0) for i in range(1, len(l1)) ]
            l01   = zip(nxs, [ x[1] for x in l0 ], [ x[1] for x in l1 ], rate0, rate1 )

            for r in l01:
                nx = "%d/%d" % (r[0], 2*r[0])
                er = "%.02e/%.02e" % (r[1], r[2])
                ra = "%.02f/%.02f" % (r[3], r[4]) if r[3] else ''
                rows.append({'method': mfilt(method), 'cfl': cfilt(cfl), 'nx': nx, 'error': er, 'rate': ra })

    print tex_table.render(caption="ASDF", rows=rows)
    print org_table.render(caption="ASDF", rows=rows)




@task
def mach2bump_fixed_dt(trial, force_run=False):

    diff = os.path.abspath('../../../BoxLib/Tools/C_util/Convergence/DiffSameDomain1d.Linux.g++.gfortran.DEBUG.ex')
    forse_run = as_boolean(force_run)

    base   = os.path.abspath('../bin/Mach2Bump')
    inputs = os.path.join(base, 'inputs-bndry-test')
    exe = {
        'rk':  os.path.join(base, 'RNS1d.Linux.gcc.gfortran.ex'),
        'rk2': os.path.join(base, 'RNS1d.Linux.gcc.gfortran.ex'),
        'rk3': os.path.join(base, 'RNS1d.Linux.gcc.gfortran.ex'),
        'rk4': os.path.join(base, 'RNS1d.Linux.gcc.gfortran.ex'),
        'sdc': os.path.join(base, 'RNS1d.SDC.Linux.gcc.gfortran.ex'),
    }

    if not exists(exe['rk4']):
        build(base, ["DIM=1", "USE_SDCLIB=FALSE", "DEBUG=FALSE", "USE_OMP=FALSE", "USE_MPI=FALSE"])

    if not exists(exe['sdc']):
        build(base, ["DIM=1", "USE_SDCLIB=TRUE", "DEBUG=FALSE", "USE_OMP=FALSE", "USE_MPI=FALSE"])

    name = lambda method, max_level, nx, dt: os.path.join(
        base, trial + '.d', '%s_lev%d_nx%04d_dt%g' % (method, max_level, nx, dt))

    nxs  = [ 64, 128, 256, 512 ]
    dts  = [ 0.005, 0.0025, 0.00125 ] # relative to nx=64

    # high-res reference run
    nx = 4096
    dt = dts[0]*64/nx/2
    reference = name('ref', 0, nx, dt)
    if not exists(reference) or force_run:
        RNS(reference, exe['sdc'], inputs, "rns.fixed_dt=%g amr.max_level=0 amr.n_cell=%d" % (dt, nx))

    # reference = name('sdc', 0, 512, 0.5)
    header = parse_first_header(reference)

    # sdc/rk runs
    for method, nx, dt, max_level in itertools.product([ 'sdc', 'rk4' ], nxs, dts, [ 0, 1 ] ):
        dt = dt*64/nx
        dname = name(method, max_level, nx, dt)
        if not exists(dname) or force_run:
            args = [ "rns.fixed_dt=%g" % dt,
                     "amr.n_cell=%d" % nx,
                     "amr.max_level=%d" % max_level ]
            RNS(dname, exe[method], inputs, " ".join(args))

    errors = []
    for method, nx, dt, max_level, norm in itertools.product(
            [ 'rk4', 'sdc' ], nxs, dts, [ 0, 1 ], [ 0, 1, 2 ]):
        dt = dt*64/nx

    # for method, nx, cfl, max_level, norm in itertools.product(
    #         [ 'rk' ], nxs, [0.1, 0.2, 0.3], [ 0 ], [ 0, 1, 2 ]):
        error = diff_reference([ name(method, max_level, nx, dt) ], reference, diff, norm=norm)
        # error = diff_initial([ name(method, max_level, nx, cfl) ], diff, norm=norm)
        for level in range(max_level+1):
            errors.append(ErrorTuple(method, max_level, nx, dt, level, norm, error[0][level][0]))

    with open('errors_%s.pkl' % trial, 'w') as f:
        pickle.dump(errors, f)

    # sdc_conv = diff_reference([ name('sdc', cfl) for cfl in cfls[::-1] ], reference, diff, norm=0)
    # rk_conv = diff_reference([ name('rk', cfl) for cfl in cfls[::-1] ], reference, diff, norm=0)

    # echo_conergence_rates(sdc_conv, header)
    # echo_conergence_rates(rk_conv, header)


#### acoustic pulse tests ##########################################################################

@task
def acoustic_pulse_convergence(trial, sdc=True, force_run=False):

    diff = DiffSameDomainRefined['1d']

    sdc = as_boolean(sdc)
    forse_run = as_boolean(force_run)

    if sdc:
        base   = os.path.abspath('../bin/AcousticPulse')
        exe    = os.path.join(base, 'RNS1d.SDC.Linux.gcc.gfortran.DEBUG.OMP.ex')
        inputs = os.path.join(base, 'inputs-test-sdc')
        if not exists(exe):
            build(base, ["DIM=1", "USE_SDCLIB=TRUE", "DEBUG=TRUE", "USE_OMP=TRUE", "USE_MPI=FALSE"])
    else:
        base   = os.path.abspath('../bin/AcousticPulse')
        exe    = os.path.join(base, 'RNS1d.Linux.gcc.gfortran.DEBUG.OMP.ex')
        inputs = os.path.join(base, 'inputs-test')
        if not exists(exe):
            build(base, ["DIM=1", "USE_SDCLIB=FALSE", "DEBUG=TRUE", "USE_OMP=TRUE", "USE_MPI=FALSE"])

    # run convergence tests
    runs = []
    for nx in [ 128, 256, 512 ]:
        dname = os.path.join(base, trial, 'nx' + str(nx))
        if not sdc:
            dname += 'rk2'
        puts(red(dname))
        if not exists(dname) or force_run:
            args  = "amr.n_cell=%d amr.blocking_factor=%d" % (nx, nx/2)
            RNS(dname, exe, inputs, args)
        runs.append(dname)

    # echo convergence rates
    dname = os.path.join(base, trial, 'nx128')
    header = parse_first_header(os.path.join(base, dname))

    if sdc:
        puts(green("convergence rates, 0-norm"))
        diffs = diff_consecutive(runs, exe=diff, norm=0)
        echo_conergence_rates(diffs, header)

        puts(green("convergence rates, 2-norm"))
        diffs = diff_consecutive(runs, exe=diff, norm=2)
        echo_conergence_rates(diffs, header)
    else:
        reference = os.path.join(base, trial, 'nx512')

        puts(green("convergence rates, 0-norm"))
        diffs = diff_reference(runs[1:], reference, exe=diff, norm=0)
        echo_conergence_rates(diffs, header)

        puts(green("convergence rates, 2-norm"))
        diffs = diff_reference(runs[1:], reference, exe=diff, norm=2)
        echo_conergence_rates(diffs, header)


#### utils #########################################################################################

def as_boolean(x):
    return x in [ True, 'True', 1, '1', 'true', 't', 'yes', 'y' ]


def build(base, options):
    with cd(base):
        run("make " + " ".join(options))


def RNS(dname='', exe=None, inputs='inputs', args='', probin=None, dry_run=False):

    puts(green("running trial " + dname))

    run("mkdir -p " + dname)
    with cd(dname):
        run("cp " + inputs + " .")
        run("cp " + exe + " .")
        if probin:
            run("cp " + probin + " .")

        exe    = os.path.join('.', os.path.basename(exe))
        inputs = os.path.basename(inputs)
        tee    = "| tee stdout 2> stderr"
        cmd    = " ".join([exe, inputs, args, tee])

        run("echo '%s' > cmd.sh" % cmd)

        if dry_run:
            puts(red(cmd))
        else:
            run(cmd)


def parse_header(plt):

    import linecache
    import os
    hname = os.path.join(plt, 'Header')

    nspecies = int(linecache.getline(hname, 2))
    species = []
    for l in range(3, 3+nspecies):
        species.append(linecache.getline(hname, l).strip())

    nlevels = int(linecache.getline(hname, 3+nspecies).strip())
    time    = float(linecache.getline(hname, 4+nspecies).strip())

    return species, nlevels, time


def parse_first_header(dname):
    plt = sorted(glob.glob(os.path.join(dname, "plt*")))[0]
    species, _, _ = parse_header(plt)
    return species


def parse_diff(diff):
    errors = {}
    for l in diff.split("\n"):
        l = l.strip()
        if l.startswith(tuple('0 1 2 3 4 5 6 7 8 9'.split())):
            l = l.split()
            lev = int(l[0])
            err = map(float, l[1:])
            errors[lev] = err
    return errors


def diff_consecutive(runs, exe=None, norm=0):

    # test convergence
    diffs = []
    for crse, fine in zip(runs[:-1], runs[1:]):
        puts(green("comparing %s and %s" % (crse, fine)))
        plot1 = sorted(glob.glob(os.path.join(crse, "plt*")))[-1]
        plot2 = sorted(glob.glob(os.path.join(fine, "plt*")))[-1]

        out = run(" ".join([exe, "infile1="+plot1, "reffile="+plot2, "norm=%d" % norm]))
        diffs.append(parse_diff(out))

    return diffs


def diff_reference(runs, fine, exe=None, norm=0):

    plot2 = sorted(glob.glob(os.path.join(fine, "plt*")))[-1]

    # test convergence
    diffs = []
    for crse in runs:
        puts(green("comparing %s and %s" % (crse, fine)))
        try:
            plot1 = sorted(glob.glob(os.path.join(crse, "plt*")))[-1]
        except IndexError:
            print "ERROR: No plot filesls  found in: %s" % crse
            raise

        out = run(" ".join([exe, "infile1="+plot1, "reffile="+plot2, "norm=%d" % norm]))
        diffs.append(parse_diff(out))

    return diffs


def diff_initial(runs, exe=None, norm=0):

    # test convergence
    diffs = []
    for crse in runs:
        puts(green("comparing %s to itself" % crse))
        plot1 = sorted(glob.glob(os.path.join(crse, "plt*")))[-1]
        plot2 = sorted(glob.glob(os.path.join(crse, "plt*")))[0]

        out = run(" ".join([exe, "infile1="+plot1, "reffile="+plot2, "norm=%d" % norm]))
        diffs.append(parse_diff(out))

    return diffs


def echo_convergence_rates(diffs, variables):

    echo = lambda s: puts(yellow(s), False)

    for crse, fine in zip(diffs[:-1], diffs[1:]):
        for lev in crse:
            e1 = np.asarray(crse[lev])
            e2 = np.asarray(fine[lev])

            rate = list(np.log10(e1/e2)/np.log10(2))
            echo("")
            echo("| variable | rate | err |")
            echo("|-")
            for variable, rate, err in zip(variables, rate, list(e1)):
                echo("| %30s | %6.2f | %12.2e |" % (variable, rate, err))
            echo("")
