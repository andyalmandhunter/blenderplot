#!/bin/python


# Copy the contents of blenderplot.py here (delete its empty main()
# function), and run as described in blenderplot.py.
#
# The functions for reading data from files are provided as an example.
# Provide your own data, or generate X (1D numpy.array), Y (1D numyp.array),
# and Z (2D numpy.array) from mathematical functions.


def get_onespec(t, p):
    '''Get spectrum for temperature t (out of 13) and power p (out of 21)
    '''

    filename = '%s/EXPpow%dmask1/alpha.dat' % (get_DRPpath(t), p)
    return np.loadtxt(filename)[np.newaxis,:,1]


def get_allspec(t):
    '''Get all spectra for temperature t (out of 13)
    '''

    for p in range(1, 21):
        if p == 1:
            spec = get_onespec(t, p)
        else:
            spec = np.vstack((get_onespec(t, p), spec))
    return spec


def get_powers(t):
    '''Get array of pump powers.
    '''

    powers = np.zeros(20)
    for p in range(1, 21):
        filename = '%s/pump_powers%d.dat' % (get_DRPpath(t), p)
        f = np.loadtxt(filename)
        powers[20-p] = f[0]
    return powers


def get_Ephot(t):
    '''Get array of photon energies. Assume this is the same for all spectra
    in a given temperature, i.e., just get the first column of power 1.
    '''

    filename = '%s/EXPpow1mask1/alpha.dat' % (get_DRPpath(t))
    return np.loadtxt(filename)[:,0]


def get_DRPpath(t):
    '''Get full CET scan path.  For some reason with the dropleton studies I
    used the inconvenient convention of 'DRP%1d' for numbers 0 through
    9, and 'DR%02d' for numbers 10 and above.
    '''

    basepath = '/Users/andy/Desktop/dropleton_temp/12_06_05';
    if t < 10:
        return '%s/DRP%1d' % (basepath, t)
    else:
        return '%s/DR%02d' % (basepath, t)


def main():
    '''main
    '''
    
    # Load data
    t = 2
    X = get_Ephot(t)
    Y = get_powers(t)
    Z = get_allspec(t)

    # Get new BlenderAxis using the loaded data
    baxis = BlenderAxis(X, Y, Z)
    baxis.set_layers(2, 2, 2, 1)  # Make the contour plot visible; hide all others

    # Orthographic plot, top view
    foc = 1
    surf_alpha = 1.0
    wire_alpha = 1.0
    el = 90
    az = 0
    baxis.write_image('example00.png',
                      az, el, foc, 'O', surf_alpha, wire_alpha)

    # Perspective plot, view from bottom left
    el = 20
    az = 30
    baxis.write_image('example01.png',
                      az, el, foc, 'P', surf_alpha, wire_alpha)

    return


if __name__ == '__main__':
    main()
