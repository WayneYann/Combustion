
from base import Container, NoInputError


class Job(Container):
    """Job container."""

    def add_to(self, jobs):
        jobs.add(self)

    def update_params(self, **kwargs):
        if hasattr(self, 'param_file'):
            with open(self.param_file, 'r') as f:
                template = f.read()
            self.inputs = template.format(**kwargs)
        else:
            raise NoInputError('update_params called but param_file not set.')
        return self

    def write_params(self, fname):
        with open(fname, 'w') as f:
            f.write(self.inputs)
        return self


