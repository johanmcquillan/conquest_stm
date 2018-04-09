
import numpy as np

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mp3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import colors as clrs

from skimage import measure

from sph import sph
from vector import Vector

# Spectroscopic notation dictionary
SPECTRAL = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g', 5: 'h', 6: 'i', 7: 'j', 8: 'k'}

AXES = ('x', 'y', 'z')


def plot_2d(cell, mesh_3d, title, axis, plane_value):
    """Plots a 2D cross-section of a 3D mesh.
    Args:
        cell (Cell): Cell object
        mesh_3d (array(float)): 3D array
        title (string): Title of plot
        axis (string): Axis normal to plane to be plotted; 'x', 'y', or 'z'
        plane_value (float): Height of plane in Bohr Radii
    """

    if axis == 'x':
        index = np.where(cell.x_mesh == plane_value)[0][0]
        mesh_2d = mesh_3d[index, :, :]
        label1 = 'z'
        label2 = 'y'
    elif axis == 'y':
        index = np.where(cell.y_mesh == plane_value)[1][0]
        mesh_2d = mesh_3d[:, index, :]
        label1 = 'z'
        label2 = 'x'
    elif axis == 'z':
        index = np.where(cell.z_mesh == plane_value)[2][0]
        mesh_2d = mesh_3d[:, :, index]
        label1 = 'y'
        label2 = 'x'
    else:
        raise ValueError("Axis must be x, y, or z")

    plt.imshow(mesh_2d, interpolation='bilinear', origin='lower', cmap=cm.afmhot, extent=(0, mesh_2d.shape[0], 0, mesh_2d.shape[1]))
    plt.colorbar()
    plt.title(title)
    plt.xlabel(label1)
    plt.ylabel(label2)
    plt.show()
    plt.close()


def plot_3d_isosurface(title, mesh, fraction, alpha=1):
    """Plot isosurface of 3D mesh.

    Args:
        title (string): Title of plot
        mesh (array(float)): 3D array to be plotted
        fraction (float, opt.): Sets value of isosurface to this fraction of max charge density
        alpha (float, opt.): Transparency of plot surface
    """

    mes = measure.marching_cubes(mesh, fraction*np.max(mesh))
    verts = mes[0]
    faces = mes[1]

    # Set up plot
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    plt.title(title)

    # Set axes
    ax.view_init(elev=45, azim=-45)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2], cmap=cm.Spectral, antialiased=False, lw=0.0, alpha=alpha)

    plt.show()
    plt.close()


def plot_radials(ion, spectro=True, interpolation='cubic'):
    """Plot all radial functions of an Ion object.

    Args:
        ion (Ion): Ion object
        spectro (bool, opt.): If true, use spectroscopic notation
        interpolation (string, opt.): Method of interpolation; possible arguments are 'linear', 'quadratic', 'cubic'
    """
    fig, ax = plt.subplots(figsize=(10, 5))

    plt.title('Radial Functions for '+str(ion.ion_name)+'.ion')
    ax.set_xlabel(r'Radial Distance, $r$ / $a_0$')
    ax.set_ylabel(r'$R_{\zeta l}(r)$')
    plt.minorticks_on()
    plt.grid(b=True, which='major', alpha=0.25, linestyle='-')
    plt.grid(b=True, which='minor', alpha=0.05, linestyle='-')

    # Loop over all radial functions
    maxy = 0
    for l in ion.radials:
        for zeta in ion.radials[l]:
            # Get Radial and data from ion
            radial = ion.get_radial(l, zeta)
            my = max(radial.radial_function_values)
            if my > maxy:
                maxy = my
            n = radial.n
            r = np.linspace(0.0, radial.cutoff, 1000)
            R = np.empty_like(r)
            for i in range(0, len(r)):
                R[i] = radial.get_value(r[i], interpolation=interpolation)

            # Add radial info to legend and add to plot
            # If spectro, use spectroscopic notation for legend
            if spectro:
                label = '$\zeta ='+str(zeta)+'$; $'+str(n)+SPECTRAL[l]+'$'
            else:
                label = '$\zeta ='+str(zeta)+'$; $n='+str(n)+'$, $l='+str(l)+'$'
            ax.plot(r, R, label=label)

    ymax = 0.2 * int(maxy / 0.2 + 2)
    ax.set_ylim(ymin=0, ymax=ymax)

    plt.legend()
    plt.show()
    plt.close()


def plot_sph_2d(l, m, axis, minimum=-8.0, maximum=8.0, plane_value=0.0, step=0.1):
    """Plot cross-section of spherical harmonic.
    All lengths measured in bohr radii (a0).

    Args:
        l (int): Orbital angular momentum quantum number for spherical harmonic
        m (int): Azimuthal quantum number for spherical harmonic
        axis (string): Cartesian axis ('x', 'y', or 'z') to set to constant value given by planeValue
        minimum (int, opt.): Minimum value of coordinates measured in a0; Default is -8
        maximum (int, opt.): Maximum value of coordinates measured in a0; Default is +8
        plane_value (float, opt.): Constant value assigned to Cartesian coordinate given by axis; Default is 0.00001
        step (float, opt.): Interval between points, measured in a0; Default is 0.1
    """

    if axis not in AXES:
        raise ValueError("Axis must be x, y, or z")

    # Initialise meshes
    # 2D cartesian mesh (x, y, or z axis determined later)
    space1, space2 = np.mgrid[minimum:maximum:step, minimum:maximum:step]
    Y = np.empty_like(space1, dtype=float)  # Spherical Harmonic mesh

    maxY = 0.1  # Colour plot sets limits to -maxY and +maxY
    for i in range(0, int((maximum - minimum) / step)):
        for j in range(0, int((maximum - minimum) / step)):
            # Use axis variable to determine which axes space1 and space2 refer to
            # Evaluate spherical harmonic at mesh point
            if axis == 'z':
                r = Vector(space2[i, j], space1[i, j], plane_value)
                plt.xlabel('$x$ / $a_0$')
                plt.ylabel('$y$ / $a_0$')
            if axis == 'y':
                r = Vector(space2[i, j], plane_value, space1[i, j])
                plt.xlabel('$x$ / $a_0$')
                plt.ylabel('$z$ / $a_0$')
            if axis == 'x':
                r = Vector(plane_value, space2[i, j], space1[i, j])
                plt.xlabel('$y$ / $a_0$')
                plt.ylabel('$z$ / $a_0$')

            Y[i, j] = sph(l, m, r)
            # Update maxY
            if abs(Y[i, j]) > maxY:
                maxY = abs(Y[i, j])

    # Setup plot
    plt.imshow(Y, interpolation='bilinear', origin='center', cmap=cm.bwr, extent=(minimum, maximum, minimum, maximum), vmin=-maxY, vmax=maxY)
    plt.colorbar()
    plt.grid()
    axes = ['x', 'y', 'z']
    axes.remove(axis)
    ttl = 'Spherical Harmonic for \n \n $l='+str(l)+'$, $m_l='+str(m)+'$ in $'+axes[0]+'-'+axes[1]+'$ plane'
    plt.title(ttl)
    plt.show()
    plt.close()


def plot_sph_3d(l, m):
    """Plots 3D spherical harmonic isosurface.

    Args:
        l (int): Degree of spherical harmonic
        m (int): Order of spherical harmonic
    """

    # Get mesh of angles
    THETA, PHI = np.mgrid[0:2*np.pi:100j, 0:np.pi:100j]
    # Initialise mesh of 
    SPH = np.zeros_like(PHI, dtype=float)

    # Loop over mesh
    for i in range(0, SPH.shape[0]):
        for j in range(0, SPH.shape[1]):
            # Get cartesian point
            x = np.sin(THETA[i, j]) * np.cos(PHI[i, j])
            y = np.sin(THETA[i, j]) * np.sin(PHI[i, j])
            z = np.cos(THETA[i, j])
            r = Vector(x, y, z)
            SPH[i, j] = abs(sph(l, m, r))
    # Get cartesian mesh
    X = SPH * np.sin(THETA) * np.cos(PHI)
    Y = SPH * np.sin(THETA) * np.sin(PHI)
    Z = SPH * np.cos(THETA)

    # Plot surface
    fig, ax = plt.subplots(subplot_kw=dict(projection='3d'), figsize=(10, 10))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_xlim3d(-1.0, 1.0)
    ax.set_ylim3d(-1.0, 1.0)
    ax.set_zlim3d(-1.0, 1.0)
    title = "Real Spherical Harmonic for Degree $l="+str(l)+"$ and Order $m="+str(m)+"$"
    plt.title(title)

    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, edgecolors=None, linewidth=0)
    plt.show()
    plt.close()


def plot_ldos_2d(
        cell, min_E, max_E, T, axis, plane_value, title=True, interpolation='cubic', debug=False):
    """Plots cross-section of LDOS.

    All lengths measured in bohr radii (a0).

    Args:
        cell (Cell): Simulation cell to plot
        min_E (float): Minimum energy
        max_E (float): Maximum energy
        T (float): Absolute temperature in K
        axis (string): Cartesian axis ('x', 'y', or 'z') to set to constant value given by plane_value
        plane_value (float, opt.): Constant value assigned to Cartesian coordinate given by axis; Default is 0.0
        title (bool, opt.): If False, show no title
        interpolation (string, opt.): Method of interpolation; possible arguments are 'linear', 'quadratic', 'cubic'
        debug (bool, opt.): If true, print extra information during runtime
    """

    if axis not in AXES:
        raise ValueError("Axis must be x, y, or z")

    # Get absolute energy
    min_E_abs = min_E + cell.fermiLevel
    max_E_abs = max_E + cell.fermiLevel

    ldos_3d = cell.get_ldos_grid(min_E_abs, max_E_abs, T, interpolation=interpolation, debug=debug)

    if np.max(ldos_3d) == 0.0:
        raise ValueError("LDoS is zero at all points")

    axes = ['x', 'y', 'z']
    axes.remove(axis)
    if title:
        ttl = cell.name+' LDoS in $'+axes[0]+'-'+axes[1]+'$ plane at $'+axis+'='+str(plane_value) + 'a_0$'
    else:
        ttl = ''
    plot_2d(cell, ldos_3d, ttl, axis, plane_value)


def plot_ldos_3d(cell, min_E, max_E, T, step=0.0, fraction=0.02, title=True, recalculate=False, vectorised=True, interpolation='cubic', alpha=1, debug=False):
    """Plots charge density isosurface.

    All lengths measured in Bohr radii (a0).

    Args:
        cell (Cell): Simulation cell to plot
        min_E (float): Minimum of energy range
        max_E (float): Maximum of energy range
        T (float): Apsolute temperature in k
        step (float, opt.): Interval between Cartesian mgrid points; Default is cell.gridSpacing
        fraction (float, opt.): Sets value of isosurface to this fraction of max charge density
        title (bool, opt.): If False, show no title
        alpha (float, opt.): Transparency of plot surfaces
        recalculate (bool, opt.): Force recalculation, even if already stored
        vectorised (bool, opt.): If true, use NumPy vectorisation
        interpolation (string, opt.): Method of interpolation; possible arguments are 'linear', 'quadratic', 'cubic'
        debug (bool, opt.): If true, print extra information during runtime
    """
    x_range = (0.0, cell.vector.x)
    y_range = (0.0, cell.vector.y)
    z_range = (0.0, cell.vector.z)

    # If step not given, set to cell.grid_spacing
    if step == 0.0:
        step = cell.grid_spacing

    # Get absolute energy
    min_E_abs = min_E + cell.fermi_level
    max_E_abs = max_E + cell.fermi_level

    ldos = cell.get_ldos_grid(min_E_abs, max_E_abs, T, vectorised=vectorised, interpolation=interpolation, recalculate=recalculate, debug=debug)
    max_ldos = np.max(ldos)

    if max_ldos == 0.0:
        raise ValueError("LDoS is zero at all points")

    if title:
        ttl = (
            cell.name+' LDoS Isosurface at '+str(fraction)+' of Maximum Density for \n Energies from '
            +str(min_E)+' eV to '+str(max_E)+' eV relative to the Fermi Level')
    else:
        ttl = ''

    plot_3d_isosurface(ttl, ldos, fraction, alpha=alpha)


def plot_current_2d(cell, z, V, T, tip_work_func, tip_energy, fraction, delta_s=None, interpolation='cubic', recalculate=False, vectorised=True, title=True, debug=False):
    """Plot constant-height STM scan.

    All lengths measured in bohr radii (a0).

    Args:
        cell (Cell): Simulation cell to plot
        z (float): z-value of plane in Bohr radii; Uses nearest mesh point to given value
        V (float): Bias voltage
        T (float): Absolute temperature
        tip_work_func (float): Work function of tip
        tip_energy (float): Fermi-level of tip
        delta_s (float, opt.): Surface broadening parameter; If None, uses default value
        fraction (float, opt.): Fraction of maximum charge density to use as isovalue for isosurface
        recalculate (bool, opt.): Force recalculation, even if already stored
        vectorised (bool, opt.): If true, use NumPy vectorisation
        interpolation (string, opt.): Method of interpolation; possible arguments are 'linear', 'quadratic', 'cubic'
        title (bool, opt.): If False, show no title
        debug (bool, opt.): Print extra information during runtime
    """

    current = cell.get_current_scan(z, V, T, tip_work_func, tip_energy, delta_s, fraction=fraction, recalculate=recalculate, vectorised=vectorised, interpolation=interpolation, debug=debug)

    if title:
        ttl = cell.name + ' STM scan at $V={:.2f}V$ at $z={}a_0$ and fraction of ${}$'.format(V, z, fraction)
    else:
        ttl = ''

    plt.imshow(current, interpolation='bilinear', origin='lower', cmap=cm.afmhot)
    plt.colorbar()
    plt.title(ttl)
    plt.xlabel('y')
    plt.ylabel('x')
    plt.show()
    plt.close()


def plot_spectrum(cell, xy, min_V, max_V, sigma, T, fraction, z, delta_s=None, dE=0.005, debug=False):
    """Get spectroscopic data from a list of specific tip positions.

    Args:
        cell (Cell): Simulation cell to plot
        xy (list(list(float)): x-y points of tip in a0; Given as [[x1, y1], [x2, y2], ...]; Uses nearest mesh point
        min_V (float): Lower bound for voltage range
        max_V (float): Upper bound for voltage range
        sigma (float): State smearing parameter in eV
        T (float): Absolute temperature in K
        z (float): z-value of plane in a0; Uses nearest mesh point to given value
        fraction (float): Fraction of maximum charge density to use as isovalue for isosurface
        delta_s (float, opt.): Surface broadening parameter; If None, uses default value
        dE (float, opt.): Energy resolution of data points
        debug (bool, opt.): Print extra information during runtime
    """

    E, LDOS = cell.get_spectrum_th(xy, min_V, max_V, sigma, T, fraction, z, delta_s=delta_s, dE=dE, debug=debug)
    fig, ax = plt.subplots(1)

    for i in range(len(xy)):
        ax.plot(E, LDOS[:, i])

    ax.set_xlim(min_V, max_V)
    ax.set_xlabel(r'Sample Bias Voltage, $V / V$')
    ax.set_yticklabels(["{:.1f}".format(t) for t in ax.get_yticks()])

    ax.set_ylabel(r'Tunnelling Conductance, d$I$/d$V$ (Arbitrary Units)')

    plt.xlim(-1.5, 1.5)
    plt.show()
    plt.close()


def plot_cits(cell, V, T, fraction, sigma, title=False, delta_s=None, debug=True):
    """Plot Current Imaging Tunnelling Spectroscopy scan - UNFINISHED and UNRELIABLE"""

    scan = cell.get_cits(V, T, fraction, sigma, delta_s=delta_s, debug=debug)

    if title:
        ttl = r"{} CITS at $V = {}V$, $\sigma = {}eV$, using ".format(cell.name, V, sigma)
    else:
        ttl = ''
    plt.title(ttl)
    plt.imshow(scan, interpolation='spline16', origin='lower left', aspect='auto', cmap=cm.afmhot)
    plt.show()
    plt.close()
