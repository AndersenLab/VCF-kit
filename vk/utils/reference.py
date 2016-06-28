from clint.textui import colored, puts, indent
import os


def get_genome_directory_file():
    """
        
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
    return genome_directory



def get_genome_list():
    """
        Return list of downloaded genomes.
    """
    return [os.path.split(x[0])[1] for x in os.walk(get_genome_directory())][1:]


def resolve_reference_genome(loc):
    """
        Resolve location of reference genome file.
    """
    if os.path.exists(loc):
        return loc
    else:
        if loc in get_genome_list():
            return "{gd}/{loc}/{loc}.fa.gz".format(gd = get_genome_directory(), loc = loc)
        else:
            with indent(4):
                puts_err(colored.red("\nGenome does not exist\n"))