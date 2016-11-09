from setuptools import setup
import glob
with open('requirements.txt') as f:
    required = f.read().splitlines()

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
      install_requires=required,
      zip_safe=False,
      include_package_data=True,
      setup_requires=['pytest-runner'],
      tests_require=['pytest'])
