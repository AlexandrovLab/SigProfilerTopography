from setuptools import setup,find_packages
import os
import shutil

#remove the dist folder first if exists
if os.path.exists("dist"):
	shutil.rmtree("dist")

def readme():
	with open('README.rst') as f:
		return(f.read())

VERSION = '1.0.87'

def write_version_py(filename='SigProfilerTopography/version.py'):
	# Copied from numpy setup.py
	cnt = """
# THIS FILE IS GENERATED FROM SIGPROFILERTOPOGRAPHY SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
	
	"""
	fh = open(filename, 'w')
	fh.write(cnt % {'version': VERSION,})
	fh.close()

write_version_py()

setup(name="SigProfilerTopography",
    version=VERSION,
    author="Burcak Otlu",
    author_email="burcakotlu@eng.ucsd.edu",
    description="SigProfilerTopography provides topography analyses for substitutions, dinucleotides and indels for all given samples.",
    url="https://github.com/AlexandrovLab/SigProfilerTopography",
    license='UCSD',
    packages=find_packages(),
    install_requires=[
	"SigProfilerMatrixGenerator>=1.1.27",
	"SigProfilerSimulator>=1.1.2",
	"SigProfilerAssignment>=0.1.0",
	"XlsxWriter>=1.3.7",
        "pandas>=1.1.5",
        "numpy>=1.20.1",
        "matplotlib>=2.2.2",
        "scipy>=1.1.0",
        "statsmodels>=0.9.0",
        "fastrand>=1.2",
        "psutil>=5.6.3",
        "intervaltree>=3.1.0"],
    include_package_data=True,
	zip_safe=False)
