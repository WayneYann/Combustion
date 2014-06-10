"""Fabric (fabfile.org) tasks for RNS.  See README."""

import glob
import os.path

import numpy as np

from fabric.api import *
from fabric.colors import *
from fabric.utils import *
from fabric.contrib.files import exists


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

        if dry_run:
            puts(red(cmd))
        else:
            run(cmd)


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


def parse_header(plt):
    import linecache
    hname = os.path.join(plt, 'Header')

    nspecies = int(linecache.getline(hname, 2))
    species = []
    for l in range(3, 3+nspecies):
        species.append(linecache.getline(hname, l).strip())

    return species


@task
def dme_flameball_convergence(force_run=False):

    force_run = False

    base   = os.path.abspath('../bin/FlameBall')
    exe    = os.path.join(base, 'RNS2d.SDC.Linux.gcc.gfortran.DEBUG.OMP.ex')
    inputs = os.path.join(base, 'inputs.dme')
    probin = os.path.join(base, 'probin.dme')
    diff   = os.path.abspath('../../../BoxLib/Tools/C_util/Convergence/DiffSameDomainRefined2d.Linux.g++.gfortran.ex')

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

    header = parse_header(sorted(glob.glob(os.path.join(base, dname, "plt*")))[0])

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


@task
def acoustic_pulse_convergence(trial, force_run=False):

    base   = os.path.abspath('../bin/AcousticPulse')
    exe    = os.path.join(base, 'RNS1d.SDC.Linux.gcc.gfortran.DEBUG.OMP.ex')
    inputs = os.path.join(base, 'inputs-test-sdc')
    diff   = os.path.abspath('../../../BoxLib/Tools/C_util/Convergence/DiffSameDomainRefined1d.Linux.g++.gfortran.ex')

    # run convergence
    runs = []
    for nx in [ 128, 256, 512 ]:
        dname = os.path.join(base, trial, 'nx' + str(nx))
        if not exists(dname) or force_run:
            args  = "amr.n_cell=%d amr.blocking_factor=%d" % (nx, nx/2)
            RNS(dname, exe, inputs, args)
        runs.append(dname)

    # test convergence
    diffs = []
    for crse, fine in zip(runs[:-1], runs[1:]):
        puts(green("comparing %s and %s" % (crse, fine)))
        plot1 = sorted(glob.glob(os.path.join(crse, "plt*")))[-1]
        plot2 = sorted(glob.glob(os.path.join(fine, "plt*")))[-1]

        out = run(" ".join([diff, "infile1="+plot1, "reffile="+plot2, "norm=0"]))
        diffs.append(parse_diff(out))


    puts(green("convergence rates"))
    for crse, fine in zip(diffs[:-1], diffs[1:]):
        for lev in crse:
            e1 = np.asarray(crse[lev])
            e2 = np.asarray(fine[lev])

            puts(yellow(), False)


@task
def test_final_integrate():

    base   = os.path.abspath('../bin/FlameBall')
    exe    = os.path.join(base, 'RNS2d.SDC.Linux.gcc.gfortran.DEBUG.OMP.ex')
    inputs = os.path.join(base, 'inputs.dme')
    probin = os.path.join(base, 'probin.dme')

    with cd(base):
        run("make -j 8")
        RNS(os.path.join(base, 'test'), exe, inputs, probin=probin)

