from clint.textui import colored, puts, puts_err, indent
import os
from . import message


def get_genome_directory_file():
    """
        Return genome directory.
    """
    module_path = os.path.split(os.path.realpath(__file__))[0]
    return module_path + "/.genome_directory"


def get_genome_directory():
    genome_directory_file = get_genome_directory_file()
    if os.path.exists(genome_directory_file):
        with open(genome_directory_file, "r") as f:
            genome_directory = os.path.abspath(os.path.expanduser(f.read()))      
    else:
        with open(genome_directory_file, "w") as f:
            genome_directory = os.path.expanduser("~/.genome/")
            f.write(genome_directory)
    # Create directory if not exists
    if not os.path.isdir(genome_directory):
        os.mkdir(genome_directory)
    return genome_directory


def get_genome_list():
    """
        Return list of downloaded genomes.
    """
    return [os.path.split(x[0])[1] for x in os.walk(get_genome_directory())][1:]



def output_genome_list():
    """
        Outputs list of available genomes
    """
    message('\n'.join(get_genome_list()))


def resolve_reference_genome(loc):
    """
        Resolve location of reference genome file.
    """

    if loc is None:
        message("You must specify a genome:")
        output_genome_list()
        exit()

    if os.path.exists(loc):
        return loc
    else:
        if loc in get_genome_list():
            reference_location = "{gd}/{loc}/{loc}.fa.gz".format(gd = get_genome_directory(), loc = loc)
            with indent(4):
                puts_err(colored.green("\nUsing reference located at %s\n" % reference_location))
            return reference_location
        else:
            with indent(4):
                exit(puts_err(colored.red("\nGenome '%s' does not exist\n" % loc)))
