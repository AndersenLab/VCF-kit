import pytest
import vcfkit
from vcfkit import *
from subprocess import check_output, Popen, PIPE
import hashlib

from cStringIO import StringIO
import sys

class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        sys.stdout = self._stdout


class Capturing_err(list):
    def __enter__(self):
        self._stderr = sys.stderr
        sys.stderr = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        sys.stderr = self._stderr

def terminal(command):
    return Popen(command, stdout=PIPE, stderr=PIPE ).communicate()
