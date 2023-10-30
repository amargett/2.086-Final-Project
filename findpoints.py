# -*- coding: utf-8 -*-
from  HullParameterization import Hull_Parameterization as HP

def find_ys(xs, zs, Vec): 
    Hull = HP(Vec)
    Hull.GenBowform()
    Hull.GenCrossSection()
    ys = []
    for i in range(len(xs)): 
        x = xs[i]
        z = zs[i]
        x1 = Hull.bow_profile(z)
        y = -1
        if x1< x< Hull.Lb:
            y = Hull.gen_waterline_bow(x,z)
        if Hull.Lb < x < Hull.LOA - Hull.Ls: 
            y = Hull.halfBeam_MidBody(z)
        ys.append(y)
    return ys
 

yval = find_ys(x_vals, z_vals, Vec)