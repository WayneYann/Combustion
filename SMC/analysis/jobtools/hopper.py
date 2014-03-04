"""Hopper scheduler (Cray PBS)."""

from StringIO import StringIO as strio

import base

from fabric.api import *


class HopperPBS(base.Container):

    cores_per_node = 24

    def submit(self, job, rwd=None, exe=None, inputs=None, queue=None,
               width=None, depth=None, pernode=None,
               walltime=None, stdout=None, stderr=None,
               pbs_opts=None, aprun_opts=None, verbose=False, dry_run=None,
               cmds=None, **kwargs):

        #
        # import defaults from fabric env
        #

        if width is None:
            width = getattr(env, 'width', None)

        if depth is None:
            depth = getattr(env, 'depth', None)

        if pernode is None:
            pernode = getattr(env, 'pernode', None)

        if queue is None:
            queue = getattr(env, 'queue', None)

        if pbs_opts is None:
            pbs_opts = getattr(env, 'pbs_opts', None)

        if aprun_opts is None:
            aprun_opts = getattr(env, 'aprun_opts', None)

        if walltime is None:
            walltime = getattr(env, 'walltime', None)

        if stdout is None:
            stdout = getattr(env, 'stdout', 'stdout')

        if stderr is None:
            stderr = getattr(env, 'stderr', 'stderr')

        verbose = verbose or getattr(env, 'verbose', False)

        #
        # sanity checks
        #

        if not rwd:
            raise ValueError("Invalid remote working directory (rwd): not set.")

        if not exe:
            raise ValueError("Invalid executable (exe): not set.")

        if pernode and width and pernode > width:
            pernode = width

        #
        # construct pbs script
        #
        opts = []               # aprun options

        pbs = [ '#!/bin/sh' ]
        pbs.append("#PBS -N " + job.name)
        if queue:
            pbs.append("#PBS -q " + queue)
        if width and depth:
            mppwidth = width*depth
            if pernode:
                mppwidth = self.cores_per_node * (width*depth/pernode)
            pbs.append("#PBS -l mppwidth=" + str(mppwidth))
            opts.append("-n " + str(width))
            opts.append("-d " + str(depth))
        elif width:
            pbs.append("#PBS -l mppwidth=" + str(width))
            opts.append("-n " + str(width))
        elif depth:
            pbs.append("#PBS -l mppdepth=" + str(depth))
            opts.append("-d " + str(depth))
        if pernode:
            # pbs.append("#PBS -l mppnppn=" + str(pernode))
            opts.append("-N " + str(pernode))
        if walltime:
            pbs.append("#PBS -l walltime=" + walltime)
        if stdout:
            pbs.append("#PBS -o {rwd}/" + stdout)
        if stderr:
            pbs.append("#PBS -e {rwd}/" + stderr)
        if pbs_opts:
            pbs.extend(pbs_opts)
        pbs.append("#PBS -V")

        pbs.append("")
        pbs.append("cd {rwd}")
        if cmds:
            pbs.extend(cmds)
        if depth:
            pbs.append("export OMP_NUM_THREADS=" + str(depth))
        pbs.append("aprun {opts} {exe} {inputs}")

        if aprun_opts:
            opts.extend(aprun_opts)

        opts = ' '.join(opts)
        pbs  = '\n'.join(pbs).format(opts=opts, rwd=rwd, exe=exe, inputs=inputs)

        # push to remote host
        run_script = rwd + '/pbs.sh'
        put(strio(pbs), run_script)

        if not dry_run:
            # submit to queue
            print 'SUBMITTING ', job.name, ' USING ', exe
            run('qsub ' + run_script)

        else:
            print pbs


