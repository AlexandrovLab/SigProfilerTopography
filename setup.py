from setuptools import setup,find_packages

setup(name="SigProfilerTopography",
    version="0.0.23",
    author="Burcak Otlu",
    author_email="burcakotlu@eng.ucsd.edu",
    description="SigProfilerTopography provides topography analyses for substitutions, dinucleotides and indels for each sample and all samples pooled.",
    url="https://github.com/AlexandrovLab/SigProfilerTopography",
	license='UCSD',
    packages=find_packages(),
    install_requires=[
        "matplotlib==2.2.2",
        "scipy==1.1.0",
        "pandas==0.23.4",
        "numpy==1.14.3",
        "statsmodels==0.9.0",
        "twobitreader"],
    include_package_data=True,
	zip_safe=False)
