import numpy as np
import os, shutil, sys
from astropy.table import Table
from subprocess import call
import matplotlib.pyplot as plt

class eazy(object):

    def __init__(self, OUTPATH, OUTPREFIX):

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

    def getParams(self, path, prefix, verbose = True):

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
        """
        Open the EAZY cache file that contains various items such as the redshift grid
        and the catalogue fluxes.
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
        """
        Open the .zout file from EAZY and read in the contents as an astropy Table object.
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

        fpath = path + prefix + ".pz"

        if not os.path.isfile(fpath):
            if verbose: print ".pz file not found - please check!"
            return False

        with open(fpath, "rb") as f:
            s = np.fromfile(f, count=2, dtype=np.int32)
            NZ, NOBJ = s[0], s[1]

            chi2 = np.fromfile(f, count=NZ*NOBJ, dtype=np.double).reshape((NOBJ,NZ))

            params, apply_prior = self.getParams(fpath, prefix)

            keys = ['NZ','NOBJ','chi2','apply_prior']
            values = [NZ, NOBJ, chi2, apply_prior]

            if apply_prior is True:
                s = np.fromfile(f, count=1, dtype=np.int32)
                NK = s[0]
                kbins = np.fromfile(f, count=NK, dtype=np.double)
                priorzk = np.fromfile(f, count=NZ*NK, dtype=np.double).reshape((NZ,NK))
                kidx = np.fromfile(f, count=NOBJ, dtype=np.int32)
                kidx[ np.logical_and(self.kidx < 0, self.kidx > NK-1) ] = -99

                keys.append(['NK','kbins','priorzk','kidx'])
                values.append([NK, kbins, priorzk, kidx])

        if verbose: print ".pz file found and read in correctly!"
        return dict(zip(keys, values))

    def getBin(self, path, prefix, verbose = True):

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


        if path is False:   path = self.path
        if prefix is False: prefix = self.prefix

        if not isinstance(idx, (list, np.ndarray)):
            idx = [idx]

        params, apply_prior = self.getParams(path, prefix)
        pz_d = self.getPz(path, prefix)
        chi2, priorzk, kidx = pz_d['chi2'], pz_d['priorzk'], pz_d['kidx']

        out = []

        for i in idx:
            g_pdf = np.exp(-0.5*chi2[i,:])
            if apply_prior is True:
                g_pdf *= priorzk[kidx[i],:]
            # Normalise them!
            out.append(g_pdf)

        if verbose: print 'PDFs returned successfully!'
        return np.array(out)

    def writeParams(self, outpath = False, newvals = False, clobber = False, verbose = True):

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

    def calcZeropoints(self, eazypath, indir, inpref, tol = 1e-2, clobber = True,
                       maxiters = 100, exclude = False, verbose = True, plot = False):

        zp_table = Table.read(indir+inpref+".zeropoint", format="ascii.no_header")

        newpvals = {"GET_ZP_OFFSETS": 1,
                    "FIX_ZSPEC": 1,
                    "BINARY_OUTPUT":1,
                    }
        self.writeParams(indir+inpref+".param", newvals=newpvals, clobber=clobber, verbose=verbose)

        zp_arr, diff_arr = [], []

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
                        medians[i] = 1.000
                corrected_ratios[i] = np.median((fit_seds[cut,i]/(zeropoints[i]*fnu[cut,i])))

            diff_arr.append(corrected_ratios)

            # print '1', corrected_ratios
            # if exclude is not False and isinstance(exclude, (list, np.ndarray)):
            #     corrected_ratios[exclude] = 1.000
            # print '2', corrected_ratios

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

        command = ("cd {directory}; {0} -p {1} -t {2} -z {3}".format(eazypath, inpref+".param", inpref+".translate",
                                                                inpref+".zeropoint",
                                                                directory=inpath))
        output = call(command, stdout = open('ezpysy.stdout','w'), stderr = open('ezpysy.stderr','w'), shell = True)

        return output

    def plotComparison(self, photoz = "z_peak", show = True):

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

        plt.show() 
        return fig, ax

    def getCompStats(self, photoz = "z_peak", verbose = True):

        specz = self.zout["z_spec"]
        photoz = self.zout[photoz]

        dz = (photoz - specz)
        diff = (dz / (1.+specz))

        nmad = 1.4826 * np.median( np.abs( dz - np.median(dz) ) )
        mean_offset = np.mean(diff)
        median_offset = np.median(diff)

        outlier1 = ((np.abs(diff) > 0.15).sum(dtype = float) / self.NOBJ)
        outlier2 = ((np.abs(diff) > 3.*nmad).sum(dtype = float) / self.NOBJ)

        print nmad, outlier1, outlier2, mean_offset, median_offset

        if verbose:
            print "#"*35
            print "NMAD: \t\t\t{0:1.3f}".format(nmad)
            print "nu 1: \t\t\t{0:1.1f}%".format(outlier1*100.)
            print "nu 2: \t\t\t{0:1.1f}%".format(outlier2*100.)
            print "mean offset: \t\t{0:1.3f}".format(mean_offset)
            print "median offset: \t\t{0:1.3f}".format(median_offset)
            print "#"*35

        keys = ["nmad", "nu1", "nu2", "mean_offset", "median_offset"]
        values = [nmad, outlier1, outlier2, mean_offset, median_offset]

        return dict(zip(keys, values))
