from subprocess import PIPE, Popen


def bcftools_version():
	"""
		Return the bcftools version
	"""
	version, err = Popen(["bcftools","--version"], stdout=PIPE, stderr=PIPE).communicate()
	if err is not "":
		return version.split("\n")[0]



