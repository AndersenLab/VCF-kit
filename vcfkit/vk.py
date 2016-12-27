#! /usr/bin/env python
"""
usage:
  vk <command> [<args>...] 
  vk setup
  vk -h | --help
  vk --version

commands:
  calc
  call
  filter
  geno
  genome
  hmm
  phylo
  primer
  rename
  tajima
  vcf2tsv

"""
from vcfkit import __version__
from utils import lev, message
from docopt import docopt
from subprocess import call, check_output, CalledProcessError
from utils.vcf import *
from clint.textui import colored, puts, indent
import sys
import vk
import os
import signal
signal.signal(signal.SIGINT, lambda x,y: sys.exit(0))
command_list = [x.strip() for x in filter(len, __doc__.splitlines()[8:])]



debug = None
if len(sys.argv) == 1:
    debug = [""]


def getScriptPath():
    return os.path.dirname(vk.__file__)


def main():
    args = docopt(__doc__,
                  argv=debug,
                  options_first=True,
                  version=__version__)
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
                        colored.red(prog + " not installed. Use a package manager to install or try using 'vk setup'\n"))
    elif args['<command>'] in command_list:
        comm = ['python', getScriptPath() + '/' + args["<command>"] + ".py"] + argv
        exit(call(comm))
    else:
        levs = [(x, lev(args['<command>'], x)) for x in command_list]
        closest =  min(levs, key = lambda x: x[1])[0]
        command = args['<command>']
        message("There is no command '{command}'. Did you mean 'vk {closest}'?".format(**locals()))

if __name__ == '__main__':
    main()
