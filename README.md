# EAZYPysy
Simple helper tools to interact with the output of the EAZY photometric redshift code (Brammer et al. 2008).

## Functionality
Current and planned functionality is described below.

- [x] Communicate with EAZY's binary file output
- [x] Support output with and without luminosity prior
- [x] Calculate photometric zeropoints automatically
- [x] Compute and retrieve the full PDF of an object
- [ ] Plot photometric vs. spectroscopic comparison
- [ ] Calculate useful statistics for photo-z/spec-z comparison
- [ ] Communicate with EAZY's ascii file output

## Initialising
To access basic functionality, initialise a new `eazy` object by providing the output directory and output file prefix of the files EAZY dumps after running. An example is shown below:

```python
import eazypysy as ez

outpath = "~/path/to/output/files/"
outpref = "mysurvey.output.prefix"

mysurvey = ez.eazy(outpath, outprefix)
```

## Computing the P(Z)
While the best-fit photometric redshift is an inportant quantity, the full probability density function (PDF), `P(z)`, can often be a useful tool. To retrieve the full PDF of the a galaxy with index `idx`, we can simply write the following.

```python
pdfs_one = mysurvey.getPDF(idx)
pdfs_many = mysurvey.getPDF(np.arange(2000))
```
The index `idx` can either be an integer, list or numpy array. In any case, an array will be returned with the resulting normalised PDFs inside it with shape `(len(idx),len(zgrid))`.

## Calculating Photometric Zeropoints
It is often very useful to apply small corrections to the photometry in your catalogue in order to undo systematic offsets between photometric bands, and to better match the template set you are using. These zeropoint offsets can be calculated easily using EASYPysy in the following way.

```python

inpath = "~/path/to/input/files/"
inpref = "mysurvey.input.prefix"
initial_zp_file = "initial.zeropoints"

mysurvey.calcZeropoints(inpath, inpref, tol=0.01)
```

To do this, EAZYPysy needs the input files, located in directory `inpath` and with all the same prefix `inpref`, before the .param, .translate and .zeropoint etc. This will iterate until a certain tolerance is obtained in each photometric band, such that `1 - |F_sed / F_cat| < tol`. A tolerance of `0.01` requires all photometric bands to be within 1% of the template SED fits.