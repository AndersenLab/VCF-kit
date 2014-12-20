from subprocess import PIPE, Popen


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
    raise Exception(bcolors.WARNING + "BCFtools Error " + bcolors.ENDC + self.error)
  else:
    return comm.strip()


def bcftools_version():
	"""
		Return the bcftools version
	"""
	version, err = Popen(["bcftools","--version"], stdout=PIPE, stderr=PIPE).communicate()
	if err is not "":
		return version.split("\n")[0]



