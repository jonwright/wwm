
import numpy as np, sys
import silicon

hc = 12.3984193 # keV/Angstrom

def anglevecs2D(a, b, c, d):
    """ Angle between g=(a,b,0) and k=(c,d,0)
    g x k = |g| |k| sin(theta)    ... cross product
    g . k = |g| |k| cos(theta)    ... dot product
    tan(theta) = sin(theta) / cos(theta)
    We compute as : atan2( gxk, g.k )
    """ 
    return np.degrees(np.arctan2( a*d - b*c, a*c + b*d))


def rotmatx(r):
    """ r is in radians """
    c = np.cos(r)
    s = np.sin(r)
    return [ [ 1,0,0], [0, c, s], [0,-s,c] ]

def rotmaty(r):
    """ r is in radians """
    c = np.cos(r)
    s = np.sin(r)
    return [ [ c,0,-s], [0, 1, 0], [s,0,c] ]

def rotmatz(r):
    """ r is in radians """
    c = np.cos(r)
    s = np.sin(r)
    return [ [ c,s,0], [-s, c, 0], [0,0,1] ]

def getrxryrz(u):
    """ 
    Converts a U matrix to give rx, ry, rz angles
    Wikipedia page X1Y2Z3 with Tait-Bryan angles \
    """
    c2s1 = -u[1,2]
    c2c1 = u[2,2]
    r1   = -np.arctan2(c2s1,c2c1)
    c2c3 = u[0,0]
    c2s3 = -u[0,1]
    r3 = -np.arctan2( c2s3, c2c3 )
    s2 = u[0,2]
    if abs(np.sin(r3)) > 0.5:
        c2 = c2s3 / np.sin(r3)
    else:
        c2 = c2c3 / np.cos(r3)
    r2 = -np.arctan2( s2, c2 )
    if 1:
        utest = np.dot(np.dot(rotmatx(r1),rotmaty(r2)),rotmatz(r3))
        assert abs(utest-u).ravel().sum() < 1e-10
    return r1,r2,r3


def compute_omega( g, tilt ):
    """ Compute the diffraction angle for peaks g
    accepts g = [gx,gy,gz] with gx as scalar or vector
    tilt is approximately radians of angle between beam and axis
    """ 
    kz = g[2] # component along the rotation axis
    modg2 = (g*g).sum(axis=0)
    num = -modg2 - 2*tilt*kz
    den = 2*np.sqrt(1 - tilt*tilt)
    kx = num / den
    arg = modg2 - kx*kx - kz*kz
    mask = arg >= 0
    ky = np.sqrt( arg * mask ) # positive
    omegaplus = anglevecs2D( g[0], g[1], kx, ky )
    omegaminus= anglevecs2D( g[0], g[1], kx, -ky )
    return omegaplus, omegaminus

if __name__=="__main__":
    u = np.loadtxt( sys.argv[1] )
    energy = float( sys.argv[2] )
    rx, ry, rz = getrxryrz(u)
    print rx, ry, rz
    tilt = 0.0004
    ub = u * hc / energy / silicon.a
    x = []
    y = []
    for h in  silicon.hkls[:20]:
        print "# H", h
        for hkl in silicon.allhkls[h]: 
            g = np.dot( ub, hkl )
            o1, o2 = compute_omega( g, tilt )
            print "%4d"*3%(hkl),"%7.2f "*2%( o1, o2 )
            x.append( o1 )
            x.append( o2 )
            y.append( silicon.fhkl[h] )
            y.append( silicon.fhkl[h] )
    import pylab as pl
    pl.plot( x, y, "+")
    tth,I = np.loadtxt( "www_plot.dat").T
    pl.plot( tth, I*200/I.max(), "-")
    pl.show()


