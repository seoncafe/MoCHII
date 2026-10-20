"""Continuous Planck attenuation and Gaussian cell-volume integration.

Valid inside the uniform interior of the fixed-state G0 sphere. The caller
must exclude the central singularity and the voxelized outer boundary.
"""
from pathlib import Path
import sys
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import PchipInterpolator

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools/python'))
from mochii_rates import sigma_HI, sigma_HeI

PC = 3.0856776e18
EV = 1.602176634e-12
KB = 1.380649e-16


def attenuation_integral(radius_pc, order=96):
    """Gamma_HI for nH=1, xHI=0.1, yHe=0.1, xHeI=1, Tstar=1e5 K."""
    x, w = leggauss(order)
    edges = (13.598, 24.587, 54.416, 100.)
    energy = np.concatenate([(a+b)/2 + (b-a)*x/2 for a,b in zip(edges[:-1],edges[1:])])
    weight = np.concatenate([(b-a)*w/2 for a,b in zip(edges[:-1],edges[1:])])
    power = energy**3 / np.expm1(energy*EV/(KB*1.e5))
    luminosity = 1.e38 * weight*power / np.dot(weight,power)
    opacity = .1*sigma_HI(energy)+.1*sigma_HeI(energy)
    coefficient = luminosity*sigma_HI(energy)/(energy*EV)
    r = np.asarray(radius_pc)*PC
    result = np.empty(r.size)
    for start in range(0,r.size,1024):
        rr = r[start:start+1024]
        result[start:start+1024] = np.exp(-rr[:,None]*opacity) @ coefficient / (4*np.pi*rr**2)
    return result


def cell_average(xyz, size, energy_order=96, spatial_order=6, radial_nodes=8193):
    """Average Gamma_HI over each cubic cell; coordinates and sizes in pc."""
    xyz = np.asarray(xyz)
    size = np.asarray(size)
    radii = np.linalg.norm(xyz,axis=0)
    lo = np.min(radii-np.sqrt(3)*size/2)
    hi = np.max(radii+np.sqrt(3)*size/2)
    if lo <= 0:
        raise ValueError('Reference cells must exclude the source singularity')
    rg = np.linspace(lo,hi,radial_nodes)
    interp = PchipInterpolator(rg,attenuation_integral(rg,energy_order))
    nodes, weights = leggauss(spatial_order)
    value = np.zeros(size.size)
    for i,x in enumerate(nodes):
        for j,y in enumerate(nodes):
            for k,z in enumerate(nodes):
                r = np.linalg.norm(xyz+np.array([x,y,z])[:,None]*size/2,axis=0)
                value += weights[i]*weights[j]*weights[k]/8*interp(r)
    return value
