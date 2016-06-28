from clint.textui import colored, puts, indent, progress
import os
import re

def message(message, n_indent = 2):
    with indent(n_indent):
        puts(colored.blue('\nSearching...\n'))


def boolify(s):
    if s == 'True':
        return True
    if s == 'False':
        return False
    raise ValueError("huh?")

def autoconvert(s):
    for fn in (boolify, int, float):
        try:
            return fn(s)
        except ValueError:
            pass
    return s

def parse_region(region):
    return re.split("[:-]+", region)