
import numpy as np

def mat(s):
    """ Convert a string "x,y,z" to a symmop """
    return tuple([eval( "lambda x,y,z: ( %s )"%(s))(v1,v2,v3) 
        for v1,v2,v3 in [ [1,0,0], [0,1,0], [0,0,1] ] ])

def op(a,b):
    """ combines two symmops to make a third"""
    return tuple([tuple(v) for v in np.dot(a,b)])

def add( g, m ):
    """ adds a symmop to a group """
    n = [m,]
    while len(n):
        g = g+n
        n = []
        for a in g:
            for b in g:
                c = op( a, b )
                if c not in g and c not in n:
                    n.append(c)
    return g


def makegrp( symops ):
    """ Builds up group from opstring """
    g = [ ((1,0,0),(0,1,0),(0,0,1)), ]
    for s in symops.split():
        g = add( g, mat(s))
    return g

def cubic():
    return makegrp( "z,x,y -y,x,z -x,-y,-z" ) 



