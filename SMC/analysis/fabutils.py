"""Fabric (fabfile.org) utilities for launching jobs."""

import glob
import os

from fabric.api import *
from fabric.colors import *
from fabric.utils import *
from fabric.contrib.files import *


env.mpirun   = 'mpirun'
env.nprocs   = 1
env.nthreads = 1


@task
def rsync():
  """Push (rsync) directories in env.rsync to env.host."""

  if env.host == 'localhost':
    return

  exclude = [ '*~', 't', '.git', '*.exe', 'plt*', '*.pyc', '*.so', 'staging' ]
  exclude = map(lambda x: "'" + x + "'", exclude)
  exclude = ' --exclude '.join(exclude)

  for d in env.rsync_dirs:
      command = "rsync -aPz --exclude {exclude} {src}/ {host}:{dst}".format(
          exclude=exclude, host=env.host, src=d, dst=d)
      local(command)


@task
def make(target=''):
  """Run *make* in the env.bin directory."""

  with cd(env.bin):
    run('make %s' % target)


class stage(object):
    """Context hander: create a staging directory upon entry and push
    to env.host:env.bin upon exit.
    """

    def __init__(self, stage='staging'):
        self.stage = stage

    def __enter__(self):
        self.mkdir()
        return self

    def __exit__(self, type, value, tb):
        import shutil
        puts(green('syncing stage'))
        local('rsync -auz {src}/* {host}:{dst}'.format(src=self.stage, host=env.host, dst=env.bin))
        shutil.rmtree(self.stage)

    def mkdir(self, *dirnames):
        path = os.path.join(self.stage, *dirnames)
        os.makedirs(path)
        return path + os.sep
