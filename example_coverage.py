import numpy
from astropy import wcs
from astropy.io import fits
from coveragemap import coordmap,coverage_multiexp

# make a test WCS
N = 4000
w = wcs.WCS(naxis=2)
w.wcs.crpix = [(N+1)/2.,(N+1)/2.]
w.wcs.cdelt = numpy.array([-1.,1.])/500.
w.wcs.crval = [40.0,-6.0]
w.wcs.ctype = ["RA---ARC", "DEC--ARC"]

ra, dec = coordmap(w,[0,N],[0,N])
deg = numpy.pi/180.
ra *= deg; dec *= deg

infile = 'tilingdata/hlwas_tiling_241117.txt.gz'

wlt = [1.1,1.5,1.9]
for i in range(3):
    map = coverage_multiexp(infile, ra, dec, 9, verbose=True, wl=wlt[i], use_t=False)
    for j in range(20):
        print('{:2d} {:7.5f}'.format(j, numpy.count_nonzero(map>=j)/numpy.size(map)))
    F = fits.PrimaryHDU(map, header=w.to_header())
    F.writeto('imageG-{:4.2f}um.fits'.format(wlt[i]), overwrite=True)
b = ['K', 'F', 'H', 'J', 'Y', 'Z', 'W']
filt = [10, 1, 2, 3, 4, 5, 0]
for i in range(6):
    map = coverage_multiexp(infile, ra, dec, filt[i], verbose=True, use_t=False)
    for j in range(20):
        print('{:2d} {:7.5f}'.format(j, numpy.count_nonzero(map>=j)/numpy.size(map)))
    F = fits.PrimaryHDU(map, header=w.to_header())
    F.writeto('image-{:s}.fits'.format(b[i]), overwrite=True)

