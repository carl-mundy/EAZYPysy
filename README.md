# EAZYPysy
Simple helper tools to interact with the output of the EAZY photometric redshift code (Brammer et al. 2008).

## Usage
To access basic functionality, initialise a new `eazy` object by providing the output directory and output file prefix of the files EAZY dumps after running. An example is shown below:

```python
import eazypysy as ez

outpath = "~/path/to/output/files/"
outpref = "mysurvey.v2"

mysurvey = ez.eazy(outpath, outprefix)
```

