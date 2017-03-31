__author__ = "Mattias de Hollander"
__copyright__ = "Copyright 2017, Mattias de Hollander"
__email__ = "m.dehollander@nioo.knaw.nl"
__license__ = "MIT"

import textwrap
import tempfile

def report(text, prefix, **files):
    try:
        import rpy2.robjects as robjects
    except ImportError:
        raise ValueError(
            "Python 3 package rpy2 needs to be installed to use the R function.")
    rmdfile = "%s.Rmd" % prefix
    f = open(rmdfile, 'w')
    f.write(text)
    f.close()
    text = format(textwrap.dedent(text))
    #code = "rmarkdown::render('%s', c('html_document', 'pdf_document'))" % rmdfile
    code = "rmarkdown::render('%s', c('html_document'))" % rmdfile
    robjects.r(code)

