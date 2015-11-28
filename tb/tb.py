#! /usr/bin/env python
"""
usage:
  tb <command> [<args>...] 
  tb setup
  tb -h | --help
  tb --version

commands:
  tajima
  filter
  call
  primer
  genome
  rename
  phylo
  freq
  geno
  vcf2tsv

"""
from docopt import docopt
from subprocess import call, check_output, CalledProcessError
from utils.vcf import *
from clint.textui import colored, puts, indent
import sys
import tb
import os


debug = None
if len(sys.argv) == 1:
    debug = [""]


def getScriptPath():
    return os.path.dirname(tb.__file__)


def main():
    args = docopt(__doc__,
                  argv=debug,
                  options_first=True)
    argv = [args['<command>']] + args['<args>']
    program_list = {"bwa": "bwa",
                    "samtools": "samtools",
                    "bcftools": "bcftools",
                    "blast": "blastn",
                    "muscle": "muscle"}
    if args["<command>"] == "setup":
        """
            Use Homebrew to install programs!
        """
        program_installed = program_list.keys()
        for install_name, program in program_list.items():
            check_output(["brew", "tap", "homebrew/science"])
            try:
                with indent(4):
                    puts(colored.blue("Installing " + install_name))
                check_output(["brew", "install", install_name])
                program_installed.remove(install_name)
            except CalledProcessError:
                try:
                    check_output(["which", program])
                    with indent(4):
                        puts(colored.blue(program + " previously installed"))
                    program_installed.remove(install_name)
                except CalledProcessError:
                    with indent(4):
                        puts(colored.red("Error installing " + install_name))
        if len(program_installed) == 0:
            with indent(4):
                puts(colored.blue("Programs successfully installed!"))
        else:
            with indent(4):
                puts(colored.red("Error: Not all programs successfully installed: " + ", ".join(program_installed)))
    elif args["<command>"] == "":
        print(__doc__)
        for prog in program_list.values():
            try:
                check_output(["which", prog])
            except CalledProcessError:
                with indent(4):
                    puts(
                        colored.red(prog + " not installed. Use a package manager to install or try using 'tb.py setup'\n"))
    elif args['<command>'] in ['tajima', 'filter', 'call', 'primer', 'genome', 'rename', 'phylo', 'freq', "geno", "vcf2tsv"]:
        comm = ['python', getScriptPath() + '/' + args["<command>"] + ".py"] + argv
        exit(call(comm))

if __name__ == '__main__':
    main()
