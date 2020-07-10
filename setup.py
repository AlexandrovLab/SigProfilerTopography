from setuptools import setup,find_packages
import os
import shutil

#remove the dist folder first if exists
if os.path.exists("dist"):
	shutil.rmtree("dist")

def readme():
	with open('README.rst') as f:
		return(f.read())

VERSION = '0.3.69'

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
    description="SigProfilerTopography provides topography analyses for substitutions, dinucleotides and indels for each sample and all samples pooled.",
    url="https://github.com/AlexandrovLab/SigProfilerTopography",
    license='UCSD',
    packages=find_packages(),
    install_requires=[
        "sigprofilermatrixgenerator>=1.0.21",
	    "sigprofilersimulator>=0.2.15",
        "matplotlib>=2.2.2",
        "scipy>=1.1.0",
        "pandas>=0.23.4",
        "numpy>=1.14.3",
        "statsmodels>=0.9.0",
        "fastrand>=1.2",
        "psutil>=5.6.3"],
    include_package_data=True,
	zip_safe=False)
