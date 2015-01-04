from subprocess import PIPE, Popen
import os, sys

class bcolors:
    BOLD = "\033[1m"
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

def command(command, shell = False):
    print command
    comm, err = Popen(command, stdout=PIPE, stderr=PIPE, shell = shell).communicate()
    if err != "":
        if type(command) == list:
            command = " ".join(command)
        raise Exception(bcolors.BOLD + "\n" + command + bcolors.ENDC +  "\n" + bcolors.WARNING + "BCFtools Error " + bcolors.ENDC + err)
    else:
        return comm.strip()


def getScriptPath():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

def get_plot(plot_name):
    path = getScriptPath()
    base = open(path + "/R/base.R",'r').read()
    plot = open(path + "/R/" + plot_name + ".R",'r').read()
    return base.replace("{plot}", plot)

def boolify(s):
    """ http://stackoverflow.com/questions/7019283/automatically-type-cast-parameters-in-python """
    if s == 'True':
        return True
    if s == 'False':
        return False
    raise ValueError("huh?")

def set_type(s):
    for fn in (boolify, int, float):
        try:
            return fn(s)
        except ValueError:
            pass
    return s

def bcftools_version():
	"""
		Return the bcftools version
	"""
	version, err = Popen(["bcftools","--version"], stdout=PIPE, stderr=PIPE).communicate()
	if err is not "":
		return version.split("\n")[0]

def make_dir(dirname):
  if not os.path.exists(dirname):
    os.mkdir(dirname)


def remove_file(file):
    try:
        os.remove(file)
    except:
        pass

def error(text):
    """ Reports an error to the user and exits """
    print("")
    print(bcolors.FAIL + "VCF-Toolbox Error " + bcolors.ENDC + text)
    print("")
    sys.exit(0)

def replace_all(text, find, replace):
    for i in find:
        text = text.replace(i, replace)
    return text

def replace_all_at_end(text, find, replace):
    for i in find:
        if text.endswith(i):
            text = text.replace(i, replace)
    return text


def bc(text, color):
    return getattr(bcolors,color) + text + bcolors.ENDC

def within_range(x,lower,upper):
    return all([x >= lower, x <= upper])
