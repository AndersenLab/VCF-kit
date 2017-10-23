from setuptools import setup
import glob
import os
with open('requirements.txt') as f:
    required = f.read().splitlines()

from vcfkit import __version__

def gen_data_files(*dirs):
    results = []

    for src_dir in dirs:
        for root,dirs,files in os.walk(src_dir):
            results.append((root, map(lambda f:root + "/" + f, files)))
    return results


setup(name='VCF-kit',
      version=__version__,
      packages=['vcfkit','vcfkit.utils'],
      description='Assorted utilities for the variant call format',
      url='https://github.com/AndersenLab/VCF-kit',
      author='Daniel E. Cook',
      author_email='danielecook@gmail.com',
      license='MIT',
      entry_points="""
      [console_scripts]
      vk = vcfkit.vk:main
      """,
      install_requires=required,
      zip_safe=False,
      package_data={
        '': ['static/*', 'static/**/*'],
      },
      include_package_data=True,
      keywords=["VCF", "variant", "caller", "format", "snps", "genetic", "variation", "genetics"],
      data_files=gen_data_files("static"),
      setup_requires=['pytest-runner'],
      tests_require=['pytest', 'coveralls'],
      classifiers=["Development Status :: 4 - Beta","Operating System :: MacOS",
                   "Operating System :: Unix",
                   "Operating System :: POSIX :: Linux",
                   "License :: OSI Approved :: MIT License"])
