# import mantid algorithms, numpy and matplotlib

import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

# mantid crystallography
from mantid import config
config.setLogLevel(0, quiet=True)
from mantid.simpleapi import *
from mantid.geometry import CrystalStructure, ReflectionGenerator, ReflectionConditionFilter



#Usage: run script by typing:
# 1) ensure you are in an environment where mantid* is installed
# 2) at terminal prompy, type:
# > python cif2lau /myDir/cifname 
# 
# (n.b. the extension .cif is not needed)
# by default, the output .lau file will be written to the same directory and have the same
# filename as the input cif file. to override this behaviour, specify the full output path
# (again the extension .lau is not needed)
#interpret input instruction
parser = argparse.ArgumentParser(
    description='''cif2lau a tool using mantid to calculate Bragg reflection
    properties from a cif format file and output a peak list to lau format
    file

    The user can specify both input and output files. Input files must have the
    extension .cif and output files have extension .lau. When specifying file names
    or locations it is not necessary to include these extensions. If only filenames
    are provided, the path is assumed to be the current directory
    ''',
    epilog="",
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-i',
    type=str,
    help='input file, .cif extension assumed',
    required = True)
parser.add_argument('-o',
    type = str,
    help ='output file, .lau extension assumed. If not specified will write to same directory as input cif',
    required=False)
parser.add_argument('-dmin',
    type=float,
    help='if specified, default minimum d-spacing (0.75 Ang) will be overriden',
    required=False)
parser.add_argument('-dmax',
    type=float,
    help='if specified, calculated reflections with longer d-spacings will be excluded',
    required=False)
parser.add_argument('-v',
    action='store_true')

args = parser.parse_args()
# if extension not provided add it
if args.i[-3:].lower()!= 'cif':
    inputCifPath=args.i + '.cif'
else:
    inputCifPath=args.i

if args.o:
    outputLauPath = args.o
    if outputLauPath[-3:].lower()!='lau':
        outputLauPath = outputLauPath + '.lau'
else:
    outputLauPath = inputCifPath[:-4] +'.lau'

print(f'\ninput file: {inputCifPath}')
print(f'output lau path: {outputLauPath}')

#specify full path for input file
cifPath = inputCifPath
#specify full path for output file
lauPath = outputLauPath
#apply min and max d-spacing limits to process
if args.dmin:
    dMin=args.dmin
else:
    dMin = 0.75

if args.dmax:
    dMax=args.dmax
else:
    dMax = 100.0

##########################################################

#Check file exists before proceeding
if not os.path.exists(inputCifPath):
    print(f'ERROR: input cif file {inputCifPath} doesn\'t exist!')
    sys.exit()


CreateSampleWorkspace(OutputWorkspace='Phase1')
LoadCIF(Workspace='Phase1', InputFile=cifPath)
ws = mtd['Phase1']
crystal1 = ws.sample().getCrystalStructure()
#Generate reflections
generator = ReflectionGenerator(crystal1)
# Create list of unique reflections between 0.7 and 3.0 Angstrom
hkls = generator.getUniqueHKLsUsingFilter(dMin, dMax, ReflectionConditionFilter.StructureFactor)
# Calculate d and F^2
dValues = generator.getDValues(hkls)
fSquared = generator.getFsSquared(hkls)
pg = crystal1.getSpaceGroup().getPointGroup()

# Make list of tuples and sort by d-values, descending, include point group for multiplicity.
reflections = sorted([(hkl, d, fsq, len(pg.getEquivalents(hkl))) for hkl, d, fsq in zip(hkls, dValues, fSquared)],
                                key=lambda x: x[1] - x[0][0]*1e-6, reverse=True)


lines = ['H K L Multiplicity d_spacing |F|^2'] #header line    
for ref in reflections:
    stol = 1/(2*ref[1])
    h = int(ref[0][0])
    k = int(ref[0][1])
    l = int(ref[0][2])
    lines.append(f'{h:4d} {k:4d} {l:4d} {ref[3]:4d} {ref[1]:.4f} {ref[2]:.4f}')

allLines = '\n'.join(lines)
if args.v:
    print(allLines)

try:
    with open(lauPath,'w+') as f:
        f.writelines(allLines)
except:
    print('ERROR: can\'t open output lau file: {lauPath}')
    sys.exit()

print(f'Successfully created file: {lauPath}')
    