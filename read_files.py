# Reading a file list
filehandle = open('info/columns_wide_reduced.txt')
filelist = list(filehandle)

# ---------------------------------------------
# The file list looks like this:
# """
# CFHTLS_W_ugriz_141754+533431_T0007_REDUCED.FITS
# CFHTLS_W_ugriz_141155+523831_T0007_REDUCED.FITS
# CFHTLS_W_ugriz_141754+514231_T0007_REDUCED.FITS
# CFHTLS_W_ugriz_142354+523831_T0007_REDUCED.FITS
# CFHTLS_W_ugriz_141754+523831_T0007_REDUCED.FITS
# CFHTLS_W_ugriz_142401+533431_T0007_REDUCED.FITS
# CFHTLS_W_ugriz_141202+514231_T0007_REDUCED.FITS
# """
# Note: in gedit there's no empty line at the end,
# but sublime text shows and empty line for the 
# same file.
# ---------------------------------------------


