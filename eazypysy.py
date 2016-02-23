import numpy as np
import os, shutil, sys
from astropy.table import Table
from subprocess import call
import matplotlib.pyplot as plt
from scipy.integrate import simps

class eazy(object):

    def __init__(self, OUTPATH, OUTPREFIX):
        """Initialise the eazy object.

        Arguments:
        OUTPATH -- path to directory containing EAZY output files
        OUTPREFIX -- prefix string of the .param, .translate and .zeropoint files
        """

        self.path = OUTPATH
        self.prefix = OUTPREFIX
        self.apply_prior = False

        if not self.path.endswith("/"):
            self.path = self.path + "/"

        self.params, self.apply_prior = self.getParams(self.path, self.prefix)
        self.zout = self.getOut(self.path, self.prefix)
        tfd = self.getTempfilt(self.path, self.prefix)
        self.NOBJ = tfd['NOBJ']
        self.NTEMP = tfd['NTEMP']
        self.NZ = tfd['NZ']
        self.NFILT = tfd['NFILT']
        self.zgrid = tfd['zgrid']

        self.z_best = self.getBin(self.path, self.prefix)

    def getParams(self, path, prefix, verbose = True):
        """Get the contents of the .param EAZY output file and return as a dictionary object.

        Arguments:
        path -- path to directory containing EAZY output files
        prefix -- prefix string of the .param, .translate and .zeropoint files
        verbose -- defines whether output is written to the terminal
        """

        fpath = path + prefix + ".param"

        if not os.path.isfile(fpath):
            if verbose: print '.params file not found - please check!'
            return False

        params = dict(np.loadtxt(fpath, dtype=str, comments="#"))

        if bool(params['APPLY_PRIOR']) or params['APPLY_PRIOR'].lower() == 'y':
            apply_prior = True
        else:
            apply_prior = False

        if verbose: print ".params file found and read in correctly!"
        return params, apply_prior

    def getCache(self, path, prefix, verbose = True):
        """Return the contents of the EAZY cache file in a dictionary

        Arguments:
        path -- path to directory containing EAZY output files
        prefix -- prefix string of the .param, .translate and .zeropoint files
        verbose -- defines whether output is written to the terminal
        """

        fpath = path + prefix + "_cache"

        if not os.path.isfile(fpath):
            if verbose: print "Cache file not found - please check!"
            return False

        with open(fpath, "rb") as f:
            s = np.fromfile(f, count=4, dtype=np.int32)
            NFILT, NTEMP, NZ, NOBJ = s[0], s[1], s[2], s[3]

            tempfilt = np.fromfile(f, count=NFILT*NTEMP*NZ, dtype=np.double).reshape((NZ,NTEMP,NFILT))
            lc = np.fromfile(f, count=NFILT, dtype=np.double)
            zgrid = np.fromfile(f, count=NZ, dtype=np.double)
            fnu = np.fromfile(f, count=NFILT*NOBJ, dtype=np.double).reshape((NOBJ,NFILT))
            efnu = np.fromfile(f, count=NFILT*NOBJ, dtype=np.double).reshape((NOBJ,NFILT))

        keys = ['NFILT','NTEMP','NZ','NOBJ','tempfilt','lc','zgrid','fnu','efnu']
        values = [NFILT, NTEMP, NZ, NOBJ, tempfilt, lc, zgrid, fnu, efnu]

        if verbose: print ".cache file found and read in correctly!"
        return dict(zip(keys, values))

    def getOut(self, path, prefix, verbose = True):
        """Open the ascii .zout file and return as an astropy.table Table object.

        Arguments:
        path -- path to directory containing EAZY output files
        prefix -- prefix string of the .param, .translate and .zeropoint files
        verbose -- defines whether output is written to the terminal
        """
        fpath = path + prefix + ".zout"

        if not os.path.isfile(fpath):
            if verbose: print ".zout file not found - please check!"
            return False

        try:
            zout = Table.read(fpath, format="ascii")
            if verbose: print ".zout file found and read in correctly!"
            return zout
        except:
            if verbose: print "Could not read .zout file - please check!"
            return False

    def getCoeff(self, path, prefix, verbose = True):
        """Open the .coeff EAZY output file and return contents as a dictionary object

        Arguments:
        path -- path to directory containing EAZY output files
        prefix -- prefix string of the .param, .translate and .zeropoint files
        verbose -- defines whether output is written to the terminal
        """

        fpath = path + prefix + ".coeff"
        # print 'get coeff path', fpath

        if not os.path.isfile(fpath):
            if verbose: print ".coeff file not found - please check!"
            return False

        with open(fpath, "rb") as f:
            s = np.fromfile(f, count=4, dtype=np.int32)
            NFILT, NTEMP, NZ, NOBJ = s[0], s[1], s[2], s[3]
            coeffs = np.fromfile(f, count = NOBJ*NTEMP, dtype = np.double).reshape((NOBJ, NTEMP))
            izbest = np.fromfile(f, count = NOBJ, dtype = np.int32)
            tnorm = np.fromfile(f, count = NTEMP, dtype = np.double)

        keys = ['NFILT','NTEMP','NZ','NOBJ','coeffs','izbest','tnorm']
        values = [NFILT, NTEMP, NZ, NOBJ, coeffs, izbest, tnorm]

        if verbose: print ".coeff file found and read in correctly!"
        return dict(zip(keys, values))

    def getTempfilt(self, path, prefix, verbose = True):
        """Open the .tempfilt EAZY output file and return contents as a dictionary object

        Arguments:
        path -- path to directory containing EAZY output files
        prefix -- prefix string of the .param, .translate and .zeropoint files
        verbose -- defines whether output is written to the terminal
        """

        fpath = path + prefix + '.tempfilt'

        if not os.path.isfile(fpath):
            if verbose: print ".tempfilt file not found - please check!"
            return False

        with open(fpath, "rb") as f:
            s = np.fromfile(f, count = 4, dtype = np.int32)
            NFILT, NTEMP, NZ, NOBJ = s[0], s[1], s[2], s[3]

            tempfilt = np.fromfile(f, count = NZ*NTEMP*NFILT, dtype = np.double).reshape((NZ, NTEMP, NFILT))
            lc = np.fromfile(f, count = NFILT, dtype = np.double)
            zgrid = np.fromfile(f, count = NZ, dtype = np.double)
            fnu = np.fromfile(f, count = NOBJ*NFILT, dtype = np.double).reshape((NOBJ, NFILT))
            efnu = np.fromfile(f, count = NOBJ*NFILT, dtype = np.double).reshape((NOBJ, NFILT))

        keys = ['NFILT','NTEMP','NZ','NOBJ','tempfilt','lc','zgrid','fnu','efnu']
        values = [NFILT, NTEMP, NZ, NOBJ, tempfilt, lc, zgrid, fnu, efnu]

        if verbose: print ".tempfilt file found and read in correctly!"
        return dict(zip(keys, values))

    def getPz(self, path, prefix, verbose = True):
        """Open the .pz EAZY output file and return contents as a dictionary object

        Arguments:
        path -- path to directory containing EAZY output files
        prefix -- prefix string of the .param, .translate and .zeropoint files
        verbose -- defines whether output is written to the terminal
        """

        fpath = path + prefix + ".pz"

        if not os.path.isfile(fpath):
            if verbose: print ".pz file not found - please check!"
            return False

        with open(fpath, "rb") as f:
            s = np.fromfile(f, count=2, dtype=np.int32)
            NZ, NOBJ = s[0], s[1]

            chi2 = np.fromfile(f, count=NZ*NOBJ, dtype=np.double).reshape((NOBJ,NZ))

            params, apply_prior = self.getParams(path, prefix)

            keys = ['NZ','NOBJ','chi2','apply_prior']
            values = [NZ, NOBJ, chi2, apply_prior]

            if apply_prior is True:
                s = np.fromfile(f, count=1, dtype=np.int32)
                NK = s[0]
                kbins = np.fromfile(f, count=NK, dtype=np.double)
                priorzk = np.fromfile(f, count=NZ*NK, dtype=np.double).reshape((NK,NZ))
                kidx = np.fromfile(f, count=NOBJ, dtype=np.int32)
                kidx[ np.logical_and(kidx < 0, kidx > NK-1) ] = -99
                self.kidx = kidx

                keys.extend(['NK','kbins','priorzk','kidx'])
                values.extend([NK, kbins, priorzk, kidx])

        if verbose: print ".pz file found and read in correctly!"
        return dict(zip(keys, values))

    def getBin(self, path, prefix, verbose = True):
        """Open the .zbin EAZY output file and return contents as a dictionary object

        Arguments:
        path -- path to directory containing EAZY output files
        prefix -- prefix string of the .param, .translate and .zeropoint files
        verbose -- defines whether output is written to the terminal
        """

        fpath = path + prefix + ".zbin"

        if not os.path.isfile(fpath):
            if verbose: print ".zbin file not found - please check!"
            return False

        with open(fpath, "rb") as f:
            s = np.fromfile(f, count=1, dtype=np.int32)
            NOBJ = s[0]

            z_a = np.fromfile(f, count=NOBJ, dtype=np.double)

            params, apply_prior = self.getParams(path, prefix)

            if apply_prior is True:
                z_p = np.fromfile(f, count=NOBJ, dtype=np.double)

            z_m1 = np.fromfile(f, count=NOBJ, dtype=np.double)
            z_m1 = np.fromfile(f, count=NOBJ, dtype=np.double)

            if apply_prior is True:
                z_peak = np.fromfile(f, count=NOBJ, dtype=np.double)
                z_best = z_peak
            else:
                z_best = z_a

        if verbose: print ".zbin file found and read in correctly!"
        return z_best

    def getPDF(self, idx, path = False, prefix = False, verbose = True):
        """Compute the normalised PDF for an object in the output

        Arguments:
        idx -- zero indexed position of object in the input catalogue
        path -- path to directory containing EAZY output files
        prefix -- prefix string of the .param, .translate and .zeropoint files
        verbose -- defines whether output is written to the terminal
        """

        if path is False:   path = self.path
        if prefix is False: prefix = self.prefix

        if not isinstance(idx, (list, np.ndarray)):
            idx = [idx]

        params, apply_prior = self.getParams(path, prefix)

        pz_d = self.getPz(path, prefix)
        chi2, priorzk, kidx = pz_d['chi2'], pz_d['priorzk'], pz_d['kidx']

        print priorzk.shape

        out = []

        for i in idx:
            g_pdf = np.exp(-0.5*chi2[i,:])
            if apply_prior is True and kidx[i] < priorzk.shape[0] and kidx[i] >= 0:
                g_pdf *= priorzk[kidx[i],:]
            # Normalise them!
            g_pdf /= simps(g_pdf, self.zgrid)
            # Save them for output
            out.append(g_pdf)

        if verbose: print 'PDFs returned successfully!'
        return np.array(out)

    def writeParams(self, outpath = False, newvals = False, clobber = False, verbose = True):
        """Write a new .params file, overwritten with new key values

        Arguments:
        outpath -- file you wish to create and overwrite
        newvals -- dictionary of values to change
        verbose -- defines whether output is written to the terminal
        """

        if newvals is not False and type(newvals) is not dict:
            if verbose: print "New parameters not in dictionary format - please check!"
            return False

        outparams = self.params

        if newvals is not False:
            for key in newvals:
                outparams[key] = newvals[key]

        if outpath is not False:
            if os.path.isfile(outpath) and clobber is False:
                if verbose: print "Out path already exists and clobber is not allowed - please check!"
                return False

            with open(outpath, "w") as f:
                for key in outparams:
                    f.write("{key}{tab}{tab}{value}{newline}".format(key=key, value=outparams[key],
                                                                     tab="\t", newline="\n"))

        if verbose: print ".params file written successfully!"
        return True

    def getZeropoints (self, indir, inpref, verbose = True):
        """Grab the current zeropoints from the .zeropoints file

        Arguments:
        indir -- path to directory containing the EAZY input files
        inpref -- prefix of the input files in indir
        """

        if os.path.isfile(indir+inpref+'.zeropoint'):
            zeropoints = np.genfromtxt(indir+inpref+'.zeropoint', usecols=1, dtype=np.float32)
            if verbose: print "zeropoints captured successfully!"
            return zeropoints
        else:
            if verbose: print ".zeropoint file not found - please check!"
            return False

    def calcZeropoints(self, eazypath, indir, inpref, tol = 1e-2, clobber = True,
                       maxiters = 100, exclude = False, verbose = True, plot = False):
        """Calculate the photometric zeropoints for a spectroscopic sample

        Arguments:
        eazypath -- absolute path to the eazy program file, e.g. "~/eazyfiles/src/eazy"
        indir -- path to directory containing the EAZY input files
        inpref -- prefix of the input files in indir
        tol -- tolerance for photometric zeropoint calculation, e.g tol = 0.01 means to within 1%
        clobber -- whether you wish to overwrite files
        maxiters -- if tol is set too low, it will continue indefinitely. Set limit here.
        exclude -- a zero indexed array of filters to keep constant/'anchor'
        verbose -- True if you want terminal output, False otherwise
        plot -- True if you want to show a plot of the process, False otherwise
        """

        zp_table = Table.read(indir+inpref+".zeropoint", format="ascii.no_header")

        newpvals = {"GET_ZP_OFFSETS": 1,
                    "FIX_ZSPEC": 1,
                    "BINARY_OUTPUT":1,
                    }
        self.writeParams(indir+inpref+".param", newvals=newpvals, clobber=clobber, verbose=verbose)

        zp_arr, diff_arr = [], []

        # zeropoints = np.ones_like(zp_table['col2'])
        zeropoints = self.getZeropoints(indir, inpref, verbose)
        print zeropoints
        if (zeropoints is False) or (len(zeropoints) == 0):
            zeropoints = np.ones_like(zp_table['col2'])
        any_above_tol = True
        iterations = 1
        filt_nums = zp_table['col1'].data
        newzppath = indir + inpref + ".zeropoint"

        if exclude is not False and isinstance(exclude, (list, np.ndarray)):
            if verbose: print "Excluding filters: ", exclude

        print 'exclude', exclude

        while any_above_tol:

            zp_arr.append( zeropoints )

            with open(newzppath, "w") as f:
                for i, med in enumerate(zeropoints):
                    f.write('{0}    {1:.4f} {2}'.format(filt_nums[i], med, '\n'))

            output = self.runEAZY(eazypath, indir, inpref)

            if output is not 0:
                if verbose: print '-----> EAZY output error. Code = {0}'.format(output)
                return False

            coeffs_d = self.getCoeff(self.path, self.prefix, verbose=verbose)
            tempfilt_d = self.getTempfilt(self.path, self.prefix, verbose=verbose)

            coeffs, tempfilt, izbest = coeffs_d['coeffs'], tempfilt_d['tempfilt'], coeffs_d['izbest']
            fnu, efnu = tempfilt_d['fnu'], tempfilt_d['efnu']

            fit_seds = (coeffs[:,:,None] * tempfilt[izbest]).sum(1)
            medians = np.copy(zeropoints)
            corrected_ratios = np.ones_like(medians)

            for i in range(len(medians)):
                cut = ((fnu > 2*efnu) * (efnu > 0.))[:,i]
                ratio = (fit_seds[cut,i]/fnu[cut,i])
                medians[i] = np.median(ratio)
                if exclude is not False and isinstance(exclude, (list, np.ndarray)):
                    if i in exclude:
                        medians[i] = 1.
                corrected_ratios[i] = np.median((fit_seds[cut,i]/(zeropoints[i]*fnu[cut,i])))
                if exclude is not False and isinstance(exclude, (list, np.ndarray)):
                    if i in exclude:
                        corrected_ratios[i] = 1.

            diff_arr.append(corrected_ratios)

            zeropoints = medians

            if verbose:
                print 'Iteration {0:.0f}: '.format(iterations),
                for i, zp in enumerate(zeropoints):
                    print '{0:.5f}'.format(corrected_ratios[i]),
                print '\n'

            iterations += 1
            if iterations > maxiters:
                if verbose: print "iterations exceeded max iterations ({0}) - please check!".format(maxiters)
                return False

            any_above_tol = (np.abs(corrected_ratios - 1.) > tol).any()

        zp_arr = np.array(zp_arr)
        diff_arr = np.array(diff_arr)

        print 'Final zeropoints:'
        for zp in zeropoints:
            print '{0:.4f}'.format(zp),
        print '\n'
        with open(newzppath,'w') as file:
            for i, med in enumerate(zeropoints):
                file.write('{0}    {1:.4f} {2}'.format(filt_nums[i], med, '\n'))
        print 'Written to file: {0}'.format(newzppath)

        if plot is True:
            fig, ax = plt.subplots(1,2, figsize=(10,4))

            for i in range(len(zeropoints)):
                ax[0].plot(zp_arr[:,i])
                ax[1].plot(diff_arr[:,i])

            ax[0].set_title("Zeropoints")
            ax[1].set_title("Observed Ratios")

            ax[1].plot([0, iterations], [1.+tol, 1+tol], '--k', alpha=0.7)
            ax[1].plot([0, iterations], [1.-tol, 1-tol], '--k', alpha=0.7)
            ax[1].plot([0, iterations], [1., 1.], '-k', lw=2, alpha=0.7)
            ax[1].set_xlim(0, iterations-2)

            plt.show()

        # Run eazy one final time now that zeropoints are written
        newpvals = {"FIX_ZSPEC": 0,
                    }
        self.writeParams(indir+inpref+".param", newvals=newpvals, clobber=clobber, verbose=verbose)
        output = self.runEAZY(eazypath, indir, inpref)

        self.__init__(self.path, self.prefix)

        return True

    def runEAZY(self, eazypath, inpath, inpref):
        """Run EAZY and collect the output.

        Arguments:
        eazypath -- Absolute or relative path to the eazy source file
        inpath -- path to directory containing input files
        inpref -- input file prefix
        """

        command = ("cd {directory}; {0} -p {1} -t {2} -z {3}".format(eazypath, inpref+".param", inpref+".translate",
                                                                inpref+".zeropoint",
                                                                directory=inpath))
        output = call(command, stdout = open('ezpysy.stdout','w'), stderr = open('ezpysy.stderr','w'), shell = True)

        return output

    def plotComparison(self, photoz = "z_peak", show = True):
        """Plot a comparison between photometric and spectroscopic redshifts

        Arguments:
        photoz -- which redshift estimator to use in the EAZY output. Default is z_peak.
        show -- True if you want a plot shown, False if not.
        """

        fig, ax = plt.subplots(1, 1, figsize=(5.5, 5))

        spec_z = self.zout["z_spec"]
        photo_z = self.zout[photoz]

        xy_max = np.max([spec_z.max(), photo_z.max()]) + 0.1
        prng = np.arange(0, xy_max+0.05, 0.05)

        ax.plot(spec_z, photo_z, 'ok', ms=4, alpha=0.6)
        ax.plot(prng, prng, '--r', lw=3, alpha=0.7)
        ax.plot(prng, prng + (1.+prng)*0.15, '-.r', lw=2, alpha=0.7)
        ax.plot(prng, prng - (1.+prng)*0.15, '-.r', lw=2, alpha=0.7)
        ax.set_ylim(0, xy_max), ax.set_xlim(0, xy_max)
        ax.set_xlabel("spec-z"), ax.set_ylabel("photo-z")

        if show:
            plt.show()

        return fig, ax

    def getCompStats(self, photoz = "z_peak", verbose = True):
        """Calculate some useful statistics for the photo-z/spec-z comparison.

        Arguments:
        photoz -- which redshift estimator to use in the EAZY output. Default is z_peak.
        verbose -- terminal output if True, False otherwise.
        """

        specz = self.zout["z_spec"]
        photoz = self.zout[photoz]

        dz = (photoz - specz)
        diff = (dz / (1.+specz))

        nmad = 1.4826 * np.median( np.abs( dz - np.median(dz) ) )
        mean_offset = np.mean(diff)
        median_offset = np.median(diff)
        dz1s = np.mean(np.abs(diff))

        outlier1 = ((np.abs(diff) > 0.15).sum(dtype = float) / self.NOBJ)
        outlier2 = ((np.abs(diff) > 3.*nmad).sum(dtype = float) / self.NOBJ)

        # print np.mean(np.abs(diff))

        # print nmad, outlier1, outlier2, mean_offset, median_offset

        if verbose:
            print "#"*35
            print "NMAD: \t\t\t{0:1.3f}".format(nmad)
            print "dz/1+z:\t\t\t{0:1.3f}".format(dz1s)
            print "nu 1: \t\t\t{0:1.1f}%".format(outlier1*100.)
            print "nu 2: \t\t\t{0:1.1f}%".format(outlier2*100.)
            print "mean offset: \t\t{0:1.3f}".format(mean_offset)
            print "median offset: \t\t{0:1.3f}".format(median_offset)
            print "#"*35

        keys = ["nmad", "nu1", "nu2", "mean_offset", "median_offset"]
        values = [nmad, outlier1, outlier2, mean_offset, median_offset]

        return dict(zip(keys, values))
