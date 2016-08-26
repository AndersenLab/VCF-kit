from setuptools import setup
import glob

setup(name='vcfkit',
      version='0.0.1',
      packages=['vcfkit','vcfkit.utils'],
      description='Tools for working with VCF files',
      url='https://github.com/AndersenLab/vcf-toolbox',
      author='Daniel Cook',
      author_email='danielecook@gmail.com',
      license='MIT',
      entry_points="""
      [console_scripts]
      vk = vcfkit.vk:main
      """,
      install_requires=["cython","docopt", "cyvcf2", "biopython", "clint", "requests"],
      zip_safe=False,
      setup_requires=['pytest-runner'],
      tests_require=['pytest'])
