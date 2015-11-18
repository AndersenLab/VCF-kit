from setuptools import setup
import glob

setup(name='vcf-toolbox',
      version='0.0.1',
      packages=['tb','tb.utils'],
      description='Tools for working with VCF files',
      url='https://github.com/AndersenLab/vcf-toolbox',
      author='Daniel Cook',
      author_email='danielecook@gmail.com',
      license='MIT',
      entry_points="""
      [console_scripts]
      tb = tb.tb:main
      """,
      install_requires=["cython","docopt", "cyvcf2", "biopython", "clint", "requests"],
      zip_safe=False)