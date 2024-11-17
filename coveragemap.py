import numpy
from astropy import wcs
from astropy.io import fits


def coordmap(usewcs, ylims, xlims):
    """Given a WCS and the range of y and x; return the RA and Dec in degrees
    """

    x,y = numpy.meshgrid(numpy.arange(xlims[0],xlims[1]), numpy.arange(ylims[0],ylims[1]))
    sh = (ylims[1]-ylims[0], xlims[1]-xlims[0])
    pixcrd = numpy.vstack(( x.flatten(), y.flatten() )).T
    world = usewcs.wcs_pix2world(pixcrd, 0)
    return world[:,0].reshape(sh), world[:,1].reshape(sh)

def coverage_singleexp_full(exp_ra, exp_dec, exp_pa, ra, dec, wl=None):
    """
    Coverage map for a single exposure: field center at exp_ra, exp_dec, exp_pa
    The grid of points on which we compute the coverage is in ra, dec
    """

    cX = [ -22.14, -22.29, -22.44, -66.42, -66.92, -67.42,-110.70,-111.48,-112.64,
            22.14,  22.29,  22.44,  66.42,  66.92,  67.42, 110.70, 111.48, 112.64]
    cY = [  12.15, -37.03, -82.06,  20.90, -28.28, -73.06,  42.20,  -6.98, -51.06,
            12.15, -37.03, -82.06,  20.90, -28.28, -73.06,  42.20,  -6.98, -51.06]

    deg = numpy.pi/180.
    pos = numpy.stack(( numpy.cos(dec)*numpy.cos(ra), numpy.cos(dec)*numpy.sin(ra), numpy.sin(dec) ), axis=-1)

    # now rotate to the WFI coordinates
    pos = pos @ numpy.asarray([[numpy.cos(exp_ra), -numpy.sin(exp_ra), 0.], [numpy.sin(exp_ra), numpy.cos(exp_ra), 0.], [0.,0.,1.]])
    # exposure is on the 'new' Prime Meridian
    pos = pos @ numpy.asarray([[-numpy.sin(exp_dec), 0., -numpy.cos(exp_dec)], [0.,1.,0.], [numpy.cos(exp_dec), 0., -numpy.sin(exp_dec)]])
    # now Z axis is toward the observer, X is North, Y is East
    pos = pos @ numpy.asarray([[numpy.sin(exp_pa), numpy.cos(exp_pa), 0.], [-numpy.cos(exp_pa), numpy.sin(exp_pa), 0.], [0.,0.,1.]])
    # now in WFI coordinates

    # incorporate distortion
    Ropt = numpy.arctan2(numpy.sqrt(pos[:,:,0]**2+pos[:,:,1]**2), -pos[:,:,2])/deg
    phi = numpy.arctan2(pos[:,:,1],pos[:,:,0])
    XAN_opt = Ropt*numpy.cos(phi)
    YAN_opt = Ropt*numpy.sin(phi)+0.496
    Ropt = numpy.sqrt(XAN_opt**2+YAN_opt**2)
    phi = numpy.arctan2(YAN_opt,XAN_opt)
    Rpar = numpy.copy(Ropt)

    # iterate on the distortion
    for j in range(3):
        distort = -0.0050802 -0.028331*Rpar +0.088569*Rpar**2 -0.022851*Rpar**3
        Rpar = Ropt*(1+numpy.clip(distort,-.2,.2))
    EFL = 18751.5 # mm
    X = (Rpar*numpy.cos(phi))*(EFL*deg)
    Y = -(Rpar*numpy.sin(phi)-0.496)*(EFL*deg)
    # now (X,Y) are the focal plane positions in mm

    # if spectroscopic, include imaging->spectroscopic focal plane conversion
    if wl is not None:
        # build coefficient table
        xcoefs = [-5.27130647e-02,  5.91325803e-05, -5.01565305e-07, -1.72937193e-02,
                  -4.06976736e-05,  3.92229941e-07, -8.09344432e-07,  8.29889801e-11,
                   2.34328882e-03, -1.18989510e-05,  1.18622889e-08, -1.54276230e-03,
                   9.24248982e-06,  1.28885245e-08,  1.55166860e-05, -3.59608603e-07,
                   2.80664230e-05]
        ycoefs = [-5.09937114e-01, -1.60135889e-02, -3.10285074e-05,  3.26210721e-05,
                  -4.54255499e-07,  2.22885228e-09, -4.91631442e-05,  1.67082508e-07,
                   4.29940080e+00, -3.01356613e-03,  1.35029131e-05, -9.42609743e-06,
                   1.16856503e-08,  4.74517081e-06,  2.30430818e-03,  7.15718123e-05,
                  -2.80457429e-08]
        pu = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2]
        px = [0, 0, 0, 1, 1, 1, 2, 2, 0, 0, 0, 1, 1, 2, 0, 0, 1]
        py = [0, 1, 2, 0, 1, 2, 0, 1, 0, 1, 2, 0, 1, 0, 0, 1, 0]

        # position corrections
        u = (wl-1.465)/.465
        DX = numpy.zeros_like(X)
        DY = numpy.zeros_like(Y)
        for k in range(len(xcoefs)):
            DX = DX + xcoefs[k] * u**pu[k] * X**px[k] * Y**py[k]
            DY = DY + ycoefs[k] * u**pu[k] * X**px[k] * Y**py[k]
        X = X+DX
        Y = Y+DY

    map = numpy.zeros_like(ra).astype(numpy.uint8)
    for i in range(18):
        map[:,:] = map + numpy.where(numpy.logical_and.reduce((pos[:,:,2]<-.9, numpy.abs(X-cX[i])<20.44, numpy.abs(Y-cY[i])<20.44)), i+1, 0)

    return(map)

def coverage_singleexp(exp_ra, exp_dec, exp_pa, ra, dec, wl=None):
    """
    Coverage map for a single exposure: field center at exp_ra, exp_dec, exp_pa
    The grid of points on which we compute the coverage is in ra, dec
    Faster than coverage_singleexp_full
    """

    (ny,nx) = numpy.shape(ra)
    ymin = 0; ymax = ny
    xmin = 0; xmax = nx

    mu = numpy.cos(exp_dec)*numpy.cos(dec)*numpy.cos(ra-exp_ra) + numpy.sin(exp_dec)*numpy.sin(dec)
    deg=numpy.pi/180
    cut = numpy.cos(.6*deg)

    map = numpy.zeros_like(ra).astype(numpy.uint8)
    if numpy.count_nonzero(mu>cut)>0:
        mu_x = numpy.amax(mu, axis=0)
        mu_y = numpy.amax(mu, axis=1)
        ymin = numpy.amin(numpy.where(mu_y>cut))
        ymax = numpy.amax(numpy.where(mu_y>cut))+1
        xmin = numpy.amin(numpy.where(mu_x>cut))
        xmax = numpy.amax(numpy.where(mu_x>cut))+1
        map[ymin:ymax,xmin:xmax] = coverage_singleexp_full(exp_ra, exp_dec, exp_pa, ra[ymin:ymax,xmin:xmax], dec[ymin:ymax,xmin:xmax], wl)

    return(map)


def coverage_multiexp(infile, ra, dec, use_filter, wl=None, verbose=False, use_t=False):
    """Generates a coverage map from an exposure file."""

    deg = numpy.pi/180

    # first get information on the region of interest

    # the center
    x = numpy.cos(dec)*numpy.cos(ra)
    y = numpy.cos(dec)*numpy.sin(ra)
    z = numpy.sin(dec)
    xc = numpy.mean(x)
    yc = numpy.mean(y)
    zc = numpy.mean(z)
    Dec_c = numpy.arctan2(zc, numpy.sqrt(xc**2+yc**2))
    RA_c = numpy.arctan2(yc,xc)
    xc = numpy.cos(Dec_c)*numpy.cos(RA_c)
    yc = numpy.cos(Dec_c)*numpy.sin(RA_c)
    zc = numpy.sin(Dec_c)
    # this procedure always gives a unit vector
    thetamax = numpy.arccos(numpy.clip(numpy.amin(x*xc+y*yc+z*zc),-1,1))
    mu_cut = numpy.cos(numpy.clip(thetamax+0.6*deg, 0., numpy.pi))
    print(xc,yc,zc,mu_cut)
    del x,y,z # free up some memory

    # read the file, convert RA/Dec/PA from degrees to radians
    info = numpy.loadtxt(infile)
    info[:,1:4] *= deg
    ntot = numpy.shape(info)[0]
    mu = numpy.sin(info[:,1])*numpy.sin(Dec_c) + numpy.cos(info[:,1])*numpy.cos(Dec_c)*numpy.cos(info[:,2]-RA_c)

    if use_t:
        map = numpy.zeros_like(ra).astype(numpy.float32)
    else:
        map = numpy.zeros_like(ra).astype(numpy.int16)

    nuse=0
    for j in range(ntot):
        if int(info[j,4])==use_filter and mu[j]>=mu_cut:
            if use_t:
                map[:,:] = map + numpy.clip(coverage_singleexp(info[j,2],info[j,1],info[j,3],ra,dec,wl),0,1)*info[j,5]
            else:
                map[:,:] = map + numpy.clip(coverage_singleexp(info[j,2],info[j,1],info[j,3],ra,dec,wl),0,1)
            nuse+=1
            if verbose and nuse%64==0: print(nuse,j,info[j,1:4]/deg,ntot)

    return(map)
