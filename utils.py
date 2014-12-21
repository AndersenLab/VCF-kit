from subprocess import PIPE, Popen
import os, sys


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

def command(command):
  comm, err = Popen(command, stdout=PIPE, stderr=PIPE).communicate()
  if err != "":
    raise Exception(bcolors.WARNING + "BCFtools Error " + bcolors.ENDC + err)
  else:
    return comm.strip()


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

class bcolors:
    BOLD = "\033[1m"
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

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


def bc(text, color):
    return getattr(bcolors,color) + text + bcolors.ENDC
