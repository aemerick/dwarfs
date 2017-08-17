from dwarfs.analysis.halogen_analysis import hgen_tools
import glob

output = glob.glob('*.tipsy.std')[0]
hgen_tools.convert_to_ascii(output, units='cgs')
