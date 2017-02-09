# -*- coding: utf-8 -*-
"""
Created on Wed May 14 15:26:28 2014
@author: pking@andersonregional.org

R. Paul King, MS MPH
Anderson Regional Medical Center
1704 23rd Avenue, 1st Floor
Meridian, MS 39301
"""
import numpy as np
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
from math import *

def leakFraction(angle, distance):
    """Returns leakage fraction
           angle; relative to isocenter->couch direction, in degrees
           distance; from isocenter, in meters""" 
    #TOMO SITE PLANNING GUIDE, T-SPG-0000 B
    #RAW VALUES -- PG 31
    dists = [1, 1.5, 2, 2.5, 3, 3.5]
    angs = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180]
    ACCURAY_LEAKAGE_RAW = np.array([
        [0,   1,       1.5,     2,       2.5,     3,       3.5    ],
        [0,   3.59E-5, 1.91E-5, 1.09E-5, 7.35E-6, 6.11E-6, 5.08E-6],
        [15,  3.35E-5, 1.99E-5, 1.36E-5, 1.04E-5, 8.23E-6, 6.62E-6],
        [30,  6.72E-5, 3.55E-5, 2.27E-5, 1.57E-5, 1.21E-5, 9.32E-6],
        [45,  5.93E-5, 3.45E-5, 2.33E-5, 1.76E-5, 1.37E-5, 1.10E-5],
        [60,  1.35E-4, 6.14E-5, 3.42E-5, 2.36E-5, 1.61E-5, 1.16E-5],
        [75,  1.94E-4, 7.96E-5, 4.25E-5, 2.59E-5, 1.73E-5, 1.24E-5],
        [90,  3.47E-4, 1.29E-4, 5.74E-5, 3.31E-5, 2.13E-5, 1.49E-5],
        [105, 3.17E-4, 1.19E-4, 6.05E-5, 3.38E-5, 2.20E-5, 1.57E-5],
        [120, 1.14E-4, 5.54E-5, 3.26E-5, 2.30E-5, 1.65E-5, 1.20E-5],
        [135, 2.24E-5, 1.49E-5, 1.11E-5, 9.05E-6, 7.61E-6, 6.28E-6],
        [150, 1.60E-5, 1.05E-5, 7.80E-6, 6.23E-6, 5.13E-6, 4.38E-6],
        [165, np.nan,  np.nan,  2.98E-6, np.nan,  np.nan,  np.nan],
        [180, np.nan,  np.nan,  2.30E-6, np.nan,  np.nan,  np.nan]])

    points = [(a,d) for a in angs for d in dists]
    values = ACCURAY_LEAKAGE_RAW[1:,1:].flatten()

    points = [points[i] for i in range(len(points)) if not np.isnan(values[i])]
    values = [values[i] for i in range(len(values)) if not np.isnan(values[i])]

    #EXTRAPOLATION FIT -- PG 32
    ACCURAY_FIT_LEAKAGE_FIT_PARAMETERS = {
        0:   (3.52E-5, -1.6058),
        15:  (3.35E-5, -1.2863),
        30:  (6.72E-5, -1.5736),
        45:  (5.93E-5, -1.3392),
        60:  (1.35E-4, -1.9428),
        75:  (1.94E-4, -2.1981),
        90:  (3.47E-4, -2.5360),
        105: (3.17E-4, -2.4154),
        120: (1.14E-4, -1.7800),
        135: (2.24E-5, -1.0010),
        150: (1.60E-5, -1.3007)}

    f = ACCURAY_FIT_LEAKAGE_FIT_PARAMETERS
    near, far, very_far = 0.5, 4, 7
    for k in f.keys():
        points.append((k, near))
        values.append(f[k][0]*near**(f[k][1]))
        points.append((k, far))
        values.append(f[k][0]*far**(f[k][1]))
        points.append((k, very_far))
        values.append(f[k][0]*very_far**(f[k][1]))
    
    #No measured falloff data at 180 deg, assume no falloff
    points.extend([(180, 0.5), (180, 1.5), (180, 2.5), (180, 3), (180, 3.5)])
    values.extend([2.30e-6, 2.30e-6, 2.30e-6, 2.30e-6, 2.30e-6])

    result = griddata(points, values, [(angle, distance)], method='cubic')
    return result[0]

def scatFraction(angle, distance):
    """Returns scatter fraction
           angle; relative to isocenter->couch direction, in degrees
           distance; from isocenter, in meters""" 
    # Ratios at r=200cm, pg 33
    ACCURAY_SCATTER_RAW = np.array([
        [0,     1.09E-5,  7.81E-5],
        [15,    1.36E-5,  8.21E-5],
        [30,    2.27E-5,  1.01E-4],
        [45,    2.33E-5,  1.14E-4],
        [60,    3.42E-5,  1.22E-4],
        [75,    4.25E-5,  1.31E-4],
        [90,    5.74E-5,  8.44E-5],
        [105,   6.05E-5,  8.43E-5],
        [120,   3.26E-5,  3.92E-5],
        [135,   1.11E-6,  2.68E-5],
        [150,   7.80E-6,  5.73E-5],
        [165,   2.98E-6,  7.59E-5],
        [180,   2.30E-5,  5.79E-5]])
        #Angle  Leakage'  Leakage&Scatter'
    s = ACCURAY_SCATTER_RAW

    scatter_per_leakage = interp1d(s[:,0], s[:,2]/s[:,1])(angle)
    leakage = leakFraction(angle, distance)
    scatFraction = scatter_per_leakage * leakage
    return scatFraction

def primFraction(angle, distance):
    """Returns primary fraction, in cGy/MU, 
           with inverse square and beamstopper attenuation.
           angle; relative to isocenter->couch direction, in degrees
           distance; from isocenter, in meters""" 
    #www.aapm.org/meetings/07ss/documents/MartinShielding_Tomotherapy.pdf
    # Max field width 5 cm at 85 cm SAD => beam angle 3.37 deg
    if angle < 90 - 1.685: 
        return 0
    elif 270 - 1.685 > angle > 90 + 1.685: 
        return 0
    elif angle > 270 + 1.685: 
        return 0
    else: 
        primary  = float(1)/16   ##for clinical use, T-SPG-0000 B, pg 33
        primary *= float(4)/100 ##beamstopper transmission, pg 33
        primary *= (float(0.85)/distance)**2
        return primary

def transmission(kind, thickness):
    """ Returns transmission through a concrete shield, thickness in cm 
        parameter kind selects the appropriate TVL at 6X """ 
    #Primary, leakage: T-SPG-0000 B, pg 35
    #Scatter, NCRP49, Fig10, at 90 deg
    TVL = {'prim':34, 'leak': 29, 'scat':17.3}
    tvl = TVL[kind]
    atten_coeff = log(10)/tvl #natural log
    return exp(-atten_coeff*thickness)

def main(): 
    assert leakFraction(45, 2) == 2.33e-05 
    assert scatFraction(45, 2) == 1.14e-04
    np.testing.assert_approx_equal(transmission('prim', 34), 0.1)
    np.testing.assert_approx_equal(transmission('leak', 29), 0.1)

    patient_dose = 200.0   #cGy/patient, average
    patients_per_day = 50  #patients per day, average
    days_per_week = 5      #treatment days per week

    cGy = patient_dose * patients_per_day * days_per_week #50000 cGy/wk
    MU  = cGy * 16   #T-SPG-0000 B, pg 32                 #800,000 MU/wk

    vault = {0: {"location": "Operator",
                 "distance": 5.94,
                 "angle": 40,
                 "thickness": 94,
                 "occupancy": 1.0,
                 "controlled": True},
            1:  {"location": "Door",
                 "distance": 6.08,
                 "angle": 10,
                 "thickness": 90,
                 "occupancy": 1.0,
                 "controlled": True},
            2:  {"location": "Hot Lab",
                 "distance": 5.14,
                 "angle": 65,
                 "thickness": 180,
                 "occupancy": 1.0,
                 "controlled": True},
            3:  {"location": "Long Room",
                 "distance": 5.18,
                 "angle": 90,
                 "thickness": 166,
                 "occupancy": 1.0,
                 "controlled": True},
            4:  {"location": "Landscape",
                 "distance": 4.14,
                 "angle": 110,
                 "thickness": 94,
                 "occupancy": 1.0/16,
                 "controlled": False},
            5: {"location": "Landscape",
                 "distance": 3.42,
                 "angle": 180,
                 "thickness": 76,
                 "occupancy": 1.0/16,
                 "controlled": False},
            6:  {"location": "Roof",
                 "distance": 2.47,
                 "angle": 90,
                 "thickness": 122,
                 "occupancy": 1.0,
                 "controlled": True},
            7:  {"location": "Roof",
                 "distance": 2.92,
                 "angle": 120,
                 "thickness": 88,
                 "occupancy": 1.0,
                 "controlled": True},
            8:  {"location": "Maze",
                 "distance": 6.67,
                 "angle": 18,
                 "thickness": 0,
                 "occupancy": 1,
                 "controlled": True}}

    print '# Location   mrem/wk   Limit    Primary  Leakage    Scatter'
    for i in range(len(vault)):
       v = vault[i]
       a = v['angle']
       d = v['distance']
       t = v['thickness']
       o = v['occupancy']
       
       primary = MU  * primFraction(a, d) * transmission('prim', t) * o * 1000
       leakage = MU  * leakFraction(a, d) * transmission('leak', t) * o * 1000
       scatter = cGy * scatFraction(a, d) * transmission('scat', t) * o * 1000

       total = (primary + leakage + scatter)

       col1 = "{0:<2}".format(i+1)
       col2 = "{0:<9}".format(v['location'])
       col3 = "{0:>9}".format('%.4f'%total)

       if v['controlled']: 
           col4 = "{0:>8}".format(100)
       else: col4 = "{0:>8}".format(2)

       col5 = "{0:>11}".format('%.4f'%primary)
       col6 = "{0:>11}".format('%.4f'%leakage)
       col7 = "{0:>11}".format('%.4f'%scatter)

       print col1 + col2 + col3 + col4 + col5 + col6 + col7

if __name__ == '__main__': main()

##  Output Results
##################
## Location   mrem/wk   Limit    Primary  Leakage    Scatter
#1 Operator    2.4485     100     0.0000     2.4438     0.0047
#2 Door        1.7190     100     0.0000     1.7136     0.0055
#3 Hot Lab     0.0026     100     0.0000     0.0026     0.0000
#4 Long Room   0.7117     100     0.7061     0.0057     0.0000
#5 Landscape   0.2868       2     0.0000     0.2867     0.0002
#6 Landscape   0.2762       2     0.0000     0.2754     0.0007
#7 Roof       62.7842     100    61.1278     1.6562     0.0002
#8 Roof       12.9396     100     0.0000    12.9309     0.0086
#9 Maze     3402.6047     100     0.0000  2506.6059   895.9988
