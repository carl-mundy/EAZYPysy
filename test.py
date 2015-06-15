import eazyPysy as ez
import numpy as np

eazy_path = "/Users/CJM2013/Desktop/eazy_test/src/eazy"

outpath = "/Users/CJM2013/Desktop/eazy_test/output/uds_dr8_15_may_2015_specz/"
outprefix = "uds_dr8_prior"

inpath = "/Users/CJM2013/Desktop/eazy_test/inputs/uds_specz/"
inpref = "uds.dr8"
initialzpfile = 'uds.dr8.zeropoint'


eazy = ez.eazy(outpath, outprefix)
eazy.calcZeropoints(eazy_path, inpath, inpref, initialzpfile, tol = 0.0075,
                    verbose=True, plot=False, exclude=[0])