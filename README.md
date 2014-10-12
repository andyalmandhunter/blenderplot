blenderplot
===============

Python class to use Blender's rendering engine to make nice-looking 3D
plots.  Drawing axes is done with LaTeX / TikZ, making it easy to
match the typographical style of the rest of your document.

I built this while writing my Ph.D. thesis, and haven't polished the
interface very much (several important setting are class variables
that have to be changed manually).  If you're interested in using
this, please send me an email!  Your enthusiasm will probably be
enough motivation for me to write a proper API / turn this into an
add-on for Blender.

A. Almand-Hunter, 2014-09-29


Help on module blenderplot:

NAME
    blenderplot

FILE
    blenderplot/blenderplot.py

DESCRIPTION
    3D plot of all spectra in a single CET scan, using Blender to render surfaces.
    Must have Blender installed.  Uses bundled version of Python.
    
Usage:

    Unless you build Blender yourself as a stand-along package, your scripts
    must be run within a Blender process:

    blender plot_scene.blend --background
            --python [your_script.py]

    Also, the Python interpreter in Blender will not be able to find
    this module if linked from another script.  The easiest thing to
    do is copy the contents of blenderplot.py to your script, and
    write your own main().

A. Almand-Hunter, 2014-09-29


CLASSES
    BlenderAxis
    
    class BlenderAxis
     |  Class for storing, manipulating, and rendering a 3D surface plot and
     |  axis in Blender.  Output uses Blender internal rendering for
     |  plotting, and LaTeX/TikZ for rendering the axes.  Output to both
     |  pdf (raster plot, vector axes and text) or pdf (all raster).  All
     |  raster output has a resolution of self.dpi
     |  
     |  Methods defined here:
     |  
     |  __init__(self, X, Y, Z)
     |      Take X (1D), Y (1D) and Z (2D) {numpy.array}s and return a new
     |      instance of BlenderAxis.  Add data, wireframe, and axes to current
     |      Blender scene.
     |  
     |  get_axis(self, X, Y, Z)
     |      Return three objects with the positions of all the X, Y, and Z tick
     |      marks (and end points, and positions of axis labels) as
     |      vertices.  This is useful so that Blender methods for object
     |      rotation and determining 2D coordinates can be used on the axis
     |      points.
     |  
     |  get_contours(self, X, Y, Z)
     |      Contour plot.  Build a solid of the plot, and use the bisect tool to
     |      make 'nslices' (instance variable) equi-Z-value filled contours
     |      of the data.  An orthographic projection of this thing viewed
     |      from directly above should look exactly like a filled contour
     |      plot.
     |  
     |  get_mesh(self, X, Y, Z)
     |      Construct filled Blender mesh from Z array, using X and Y vectors for X
     |      and Y positions of points.  Return a Blender object.
     |  
     |  get_raw(self, X, Y, Z)
     |      Construct a wireframe consisting just the rows of data, not connected to
     |      each other (useful for plotting a series of spectra, for example).
     |  
     |  get_wireframe(self, X, Y, Z)
     |      Generate wireframe as a sparse version of X, Y, Z mesh (but connected
     |      along contour of surface)
     |  
     |  regrid(self, X, Y, Z)
     |      Interpolate onto an evenly-spaced grid in x and y, spacing dx.  Return a
     |      tuple of new (X, Y, Z).  For both input and output, X and Y are
     |      1D, and Z is 2D.
     |  
     |  rotate_z(self, ob, center, angle)
     |      Rotate Blender object through angle degrees about a vector passing
     |      through center parallel to the z axis.
     |  
     |  set_layers(self, raw_layer, surface_layer, wireframe_layer, contours_layer)
     |      Set layers of objects (raw, surface, wireframe).  0 and 1 are active.
     |  
     |  update_raw(self, X, Y, Z)
     |      Replace mesh in 'CETraw' with a new one from given data.  For
     |      convenience, put 'CETraw' in an active layer before returning.
     |  
     |  write_image(self, filename, rot_angle, cam_el, foc, cam_type, surf_alpha, wire_alpha)
     |      Rotate (copy of) axis by rot_angle, elevate camera by cam_el, output
     |      .png and .tex files into temporary directory, process with pdfLaTeX,
     |      crop with pdfcrop, and convert back to .png with ghostscript.
     |      
     |      Arguments:
     |      
     |          filename -- (string) filename
     |      
     |          rot_angle -- (float) angle in degrees to rotate the plot around
     |              local z
     |      
     |          cam_el -- (float) angle to elevate the camera, still pointing at the
     |              center of the plot
     |      
     |          foc -- (float) scale the focal length of the lens, maintaining the
     |              field of view at the plot distance
     |      
     |          cam_type -- (string) 'O' for orthographic, 'P' (actually anything
     |              but 'O' for now) for perspective
     |      
     |          surf_alpha -- (float) alpha for surface-like things
     |      
     |          wire_alpha -- (float) alpha for wireframe-like things
     |              (wireframes of data, raw spectra if this is a series
     |              of spectra)
     |  
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |  
     |  DX = 1
     |  
     |  DY = 1
     |  
     |  Zscale = 20
     |  
     |  dpi = 128
     |  
     |  dx = 0.1
     |  
     |  nslices = 25.0
     |  
     |  ticklen = 0.5
     |  
     |  x_size = 2.5
     |  
     |  xlabel = r'$\hbar\omega - E_{\sf 1s}$ (meV)'
     |  
     |  xticks = [-3, 0, 3, 6]
     |  
     |  y_size = 2.5
     |  
     |  ylabel = 'Pump power (mW)'
     |  
     |  yticks = [1, 5, 9, 13]
     |  
     |  zlabel = '$1 - T$'
     |  
     |  zticks = [0.0, 0.2, 0.4]
