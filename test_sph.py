#!/usr/bin/env python2.7

import conquest_stm.plot as plot

__author__ = 'Johan G. McQuillan'
__email__ = 'johan.mcquillan.13@ucl.ac.uk'

lMax = 4

print '3-DIMENSIONAL SPHERICAL HARMONICS'
print '   Enter \'q\' to quit'

quitting = False
while not quitting:
    gotL = False
    while not gotL and not quitting:
        lString = raw_input('Enter Degree, l = ')
        if lString == 'q':
            quitting = True
        else:
            try:
                l = int(lString)
                if l > lMax or l < 0:
                    raise ValueError
                else:
                    gotL = True
            except ValueError:
                print 'l must be between an integer between 0 and {} inclusive'.format(lMax)

    if not quitting:
        if l == 0:
            m = 0
        else:
            gotM = False
            while not gotM and not quitting:
                if l == 0:
                    m = 0
                else:
                    mString = raw_input('Enter Order, m = ')
                    if mString == 'q':
                        quitting = True
                    else:
                        try:
                            m = int(mString)
                            if m > l or m < -l:
                                raise ValueError
                            else:
                                gotM = True
                        except ValueError:
                            print 'm must be and integer {} and {} inclusive'.format(-l, +l)
        if not quitting:
            print 'Showing Spherical Harmonic for l={} and m={}'.format(l, m)
            plot.plot_sph_3d(l, m)
