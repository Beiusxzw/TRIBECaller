import os
import shutil
import sys
from cx_Freeze import setup, Executable
from TRIBECaller._version import __version__

os.environ['TCL_LIBRARY'] = r'C:\bin\Python37-32\tcl\tcl8.6'
os.environ['TK_LIBRARY'] = r'C:\bin\Python37-32\tcl\tk8.6'

__version__ = '1.0.0'
base = None
if sys.platform == 'win32':
    base = 'Win32GUI'

include_files = []
includes = ['numpy','scipy','pysam','tqdm']
excludes = []
packages = ['numpy']

setup(
    name='TRIBECaller',
    description='TRIBECaller',
    version=__version__,
    executables=[Executable('main.py', base=base)],
    options = {'build_exe': {
        'packages': packages,
        'includes': includes,
        'include_files': include_files,
        'include_msvcr': True,
        'excludes': excludes,
    }},
)