#!/usr/bin/python
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <http://www.gnu.org/licenses/>.


'''\
3D plot of all spectra in a single CET scan, using Blender to render surfaces.
Must have Blender installed.  Uses bundled version of Python.

Requires Blender, LaTeX, TikZ, pdfcrop, and Ghostscript

Usage:
    Unless you build Blender yourself as a stand-along package, your scripts
    must be run within a Blender process:

    blender plot_scene.blend --background
            --python [your_script.py]

    Also, the Python interpreter in Blender will not be able to find this
    module.  The easiest thing to do is copy the contents of blenderplot.py
    to your script, and write your own main().

A. Almand-Hunter, 2014-09-29
'''


import bpy
import bmesh
import bpy_extras
import mathutils
import numpy as np
import textwrap as tw
import subprocess


class BlenderAxis:
    '''Class for storing, manipulating, and rendering a 3D surface plot and
    axis in Blender.  Output uses Blender internal rendering for
    plotting, and LaTeX/TikZ for rendering the axes.  Output to both
    pdf (raster plot, vector axes and text) or png (all raster).  All
    raster output has a resolution of self.dpi
    '''

    dpi = 300  # note, 128 is macbook screen resolution
    x_size = 3  # inches
    y_size = 3

    Zscale = 20

    xticks = [-3, 0, 3, 6]  # Data units
    yticks = [1, 5, 9, 13]
    zticks = [0., 0.2, 0.4]
    ticklen = 0.5  # Blender units

    xlabel = '$\hbar\omega - E_{\sf 1s}$ (meV)'
    ylabel = 'Pump power (mW)'
    zlabel = '$1 - T$'

    dx = 0.1  # regrid spacing

    DX = 1  # wireframe spacing
    DY = 1

    nslices = 25.0  # for contour plot


    def __init__(self, X, Y, Z):
        '''Take X (1D), Y (1D) and Z (2D) {numpy.array}s and return a new
        instance of BlenderAxis.  Add data, wireframe, and axes to current
        Blender scene.
        '''

        # Get current Blender scene
        self.scn = bpy.context.scene

        # Regrid data
        self.X, self.Y, self.Z = self.regrid(X, Y, Z)

        # Save axis of rotation, and center camera vertically on plot
        self.rot_axis = mathutils.Vector((np.mean(self.X[[0,-1]]),
            np.mean(self.Y[[0,-1]]),
            np.mean(np.array([0, self.Z.max()*self.Zscale]))))
        bpy.data.objects['Camera.001'].location[2] = self.rot_axis[2]
        self.scn.cursor_location = self.rot_axis

        # Get data objects
        self.data_raw = self.get_raw(X, Y, Z)
        self.data = self.get_mesh(self.X, self.Y, self.Z)
        self.data_wireframe = self.get_wireframe(self.X, self.Y, self.Z)
        self.data_contours = self.get_contours(self.X, self.Y, self.Z)
        self.axX, self.axY, self.axZ = self.get_axis(X, Y, Z)

        # Add objects to scene
        self.scn.objects.link(self.data_raw)
        self.scn.objects.link(self.data)
        self.scn.objects.link(self.data_wireframe)
        self.scn.objects.link(self.data_contours)
        self.scn.objects.link(self.axX)
        self.scn.objects.link(self.axY)
        self.scn.objects.link(self.axZ)

        # Set smooth flag for shading
        self.data_raw.select = True
        self.data.select = True
        self.data_wireframe.select = True
        bpy.ops.object.shade_smooth()

        self.set_layers(2, 2, 2, 1)

        # Set material to Matlab-style color map (stored in .blend file)
        self.data_raw.data.materials.append(bpy.data.materials['Material.004'])
        self.data.data.materials.append(bpy.data.materials['Material.003'])
        self.data_wireframe.data.materials.append(bpy.data.materials['Material.004'])
        self.data_contours.data.materials.append(bpy.data.materials['Material.003'])

        # Set image resolution
        self.scn.render.resolution_x = self.x_size * self.dpi
        self.scn.render.resolution_y = self.y_size * self.dpi

        return


    def set_layers(self, raw_layer, surface_layer, wireframe_layer, contours_layer):
        '''Set layers of objects (raw, surface, wireframe).  0 and 1 are active,
        only 1 will cast shadows.
        '''

        self.data_raw.layers = [n == raw_layer for n in range(20)]
        self.data.layers = [n == surface_layer for n in range(20)]
        self.data_wireframe.layers = [n == wireframe_layer for n in range(20)]
        self.data_contours.layers = [n == contours_layer for n in range(20)]

        return


    def write_image(self, filename, rot_angle, cam_el, foc, cam_type,
                    surf_alpha, wire_alpha):
        '''Rotate (copy of) axis by rot_angle, elevate camera by cam_el, output
        .png and .tex files into temporary directory, process with pdfLaTeX,
        crop with pdfcrop, and convert back to .png with ghostscript.

        Arguments:

            filename -- (string) filename

            rot_angle -- (float) angle in degrees to rotate the plot around
                local z

            cam_el -- (float) angle to elevate the camera, still pointing at the
                center of the plot

            foc -- (float) scale the focal length of the lens, maintaining the
                field of view at the plot distance

            cam_type -- (string) 'O' for orthographic, 'P' (actually anything
                but 'O' for now) for perspective

            surf_alpha -- (float) alpha for surface-like things

            wire_alpha -- (float) alpha for wireframe-like things
                (wireframes of data, raw spectra if this is a series
                of spectra)
        '''

        # Save current axis
        data_raw_save = bmesh.new()
        data_raw_save.from_mesh(self.data_raw.data)
        data_save = bmesh.new()
        data_save.from_mesh(self.data.data)
        data_wireframe_save = bmesh.new()
        data_wireframe_save.from_mesh(self.data_wireframe.data)
        data_contours_save = bmesh.new()
        data_contours_save.from_mesh(self.data_contours.data)
        axX_save = bmesh.new()
        axX_save.from_mesh(self.axX.data)
        axY_save = bmesh.new()
        axY_save.from_mesh(self.axY.data)
        axZ_save = bmesh.new()
        axZ_save.from_mesh(self.axZ.data)

        # Elevate camera
        cam = bpy.data.objects['Camera.001']
        cam_loc_save = cam.location.copy()
        cam_rot_save = cam.rotation_euler.copy()
        cam_rot_matrix = mathutils.Matrix.Rotation(-1*cam_el*np.pi/180, 3, 'X')
        cam.location = cam_rot_matrix*(cam.location - self.rot_axis) + self.rot_axis
        cam.rotation_euler.rotate_axis('X', -1*cam_el*np.pi/180)

        # Deal with lens-focal-length--scale factor
        foc_rot_axis = self.rot_axis.copy()
        foc_rot_axis[0] = cam.location[0]
        cam.location = foc*(cam.location - foc_rot_axis) + foc_rot_axis
        lens_save = bpy.data.cameras['Camera.001'].lens
        bpy.data.cameras['Camera.001'].lens *= foc

        # Switch camera to orthographic projection if requested,
        # keeping the same field of view at the center of the plot
        self.scn.update()
        if cam_type == 'O':

            # Calculate ortho scale
            x0 = bpy_extras.object_utils.world_to_camera_view(self.scn,
                cam, mathutils.Vector([self.X[0],
                                       self.rot_axis[1],
                                       self.rot_axis[2]]))[0]
            x1 = bpy_extras.object_utils.world_to_camera_view(self.scn,
                cam, mathutils.Vector([self.X[-1],
                                       self.rot_axis[1],
                                       self.rot_axis[2]]))[0]
            bpy.data.cameras['Camera.001'].type = 'ORTHO'
            bpy.data.cameras['Camera.001'].ortho_scale = (self.X[-1] -
                                                          self.X[0]) / (x1-x0)

        # Rotate objects
        rot_axis = (np.mean(self.X[[0,-1]]), np.mean(self.Y[[0,-1]]), 0)
        self.rotate_z(self.data_raw, rot_axis, rot_angle)
        self.rotate_z(self.data, rot_axis, rot_angle)
        self.rotate_z(self.data_wireframe, rot_axis, rot_angle)
        self.rotate_z(self.data_contours, rot_axis, rot_angle)
        self.rotate_z(self.axX,  rot_axis, rot_angle)
        self.rotate_z(self.axY, rot_axis, rot_angle)
        self.rotate_z(self.axZ, rot_axis, rot_angle)
        
        # Set material alphas
        bpy.data.materials['Material.003'].alpha = surf_alpha
        bpy.data.materials['Material.003'].specular_alpha = surf_alpha
        bpy.data.materials['Material.004'].alpha = wire_alpha

        # Save .png
        self.scn.render.filepath = filename
        bpy.ops.render.render(write_still=True)


        # WRITE .TEX FILE
        tex_filename = 'temp.tex'
        f = open(tex_filename, 'w')

        # Beginning of .tex file
        f.write(tw.dedent('''\
        \\documentclass{beamer}
        \\setlength{\\paperwidth}{60cm}
        \\setlength{\\paperheight}{60cm}
        \\setlength{\\textwidth}{60cm}
        \\setlength{\\textheight}{60cm} 
        \\setbeamertemplate{navigation symbols}{}
        \\usepackage{tikz}

        \\begin{document}

        \\begin{frame}{}
        \\begin{tikzpicture}[line width=0.5pt]
          \\path (0,0) node{\\includegraphics[width = %.3fin]{%s}};
          \\begin{scope}[shift = {(-%.3fin,-%.3fin)}, x=1in, y=1in]
        ''') % (self.x_size, filename, self.x_size/2, self.y_size/2))

        # X axis
        axXmesh = bpy.data.meshes['CETaxisMeshX']
        axXxy = []
        for vert in axXmesh.vertices:
            axXxy.append(bpy_extras.object_utils.world_to_camera_view(self.scn,
                cam, vert.co))
        X0 = axXxy[0][0]
        X1 = axXxy[-5][0]
        Y0 = axXxy[0][1]
        Y1 = axXxy[-5][1]
        f.write('    \\draw[-latex] (%.3f,%.3f) -- (%.3f,%.3f);\n' %
                (X0*self.x_size, Y0*self.y_size,
                 X1*self.x_size+0.1, Y1*self.y_size))
        # X-axis label
        X0 = axXxy[-3][0]
        X1 = axXxy[-1][0]
        Y0 = axXxy[-3][1]
        Y1 = axXxy[-1][1]
        label_angle = 180/np.pi * np.arctan2((Y1 - Y0), (X1 - X0))
        f.write(('    \\path (%.3f,%.3f)' +
                 'node[rotate = %.2f, font = \\footnotesize]' +
                 '{%s};\n') %
                 (axXxy[-2][0]*self.x_size, axXxy[-2][1]*self.y_size,
                  label_angle, self.xlabel))
        # Tick labels
        ntick = 0
        for n in range(2, len(axXxy)-5, 2):
            tX0 = axXxy[n][0] * self.x_size
            tX1 = axXxy[n+1][0] * self.x_size
            tY0 = axXxy[n][1] * self.y_size
            tY1 = axXxy[n+1][1] * self.y_size
            anchor_angle = 180/np.pi * np.arctan2((tY0 - tY1), (tX0 - tX1))
            f.write(('    \\draw (%.3f,%.3f) -- (%.3f,%.3f) ' +
                     'node[anchor=%.2f, font=\\footnotesize] {$%.0f$};\n') %
                    (tX0, tY0, tX1, tY1,
                     anchor_angle, self.xticks[ntick]))
            ntick += 1

        # Y axis
        axYmesh = bpy.data.meshes['CETaxisMeshY']
        axYxy = []
        for vert in axYmesh.vertices:
            axYxy.append(bpy_extras.object_utils.world_to_camera_view(self.scn,
                cam, vert.co))
        X0 = axYxy[0][0]
        X1 = axYxy[-5][0]
        Y0 = axYxy[0][1]
        Y1 = axYxy[-5][1]
        f.write('    \\draw[-latex] (%.3f,%.3f) -- (%.3f,%.3f);\n' %
                (X0*self.x_size, Y0*self.y_size,
                 X1*self.x_size, Y1*self.y_size+0.1))
        # Y-axis label
        X0 = axYxy[-3][0]
        X1 = axYxy[-1][0]
        Y0 = axYxy[-3][1]
        Y1 = axYxy[-1][1]
        label_angle = 180/np.pi * np.arctan2((Y1 - Y0), (X1 - X0))
        f.write(('    \\path (%.3f,%.3f)' +
                 'node[rotate = %.2f, font = \\footnotesize]' +
                 '{%s};\n') %
                 (axYxy[-2][0]*self.x_size, axYxy[-2][1]*self.y_size,
                  label_angle, self.ylabel))
        # Tick labels
        ntick = 0
        for n in range(2, len(axYxy)-5, 2):
            tX0 = axYxy[n][0] * self.x_size
            tX1 = axYxy[n+1][0] * self.x_size
            tY0 = axYxy[n][1] * self.y_size
            tY1 = axYxy[n+1][1] * self.y_size
            anchor_angle = 180/np.pi * np.arctan2((tY0 - tY1), (tX0 - tX1))
            f.write(('    \\draw (%.3f,%.3f) -- (%.3f,%.3f) ' +
                    'node[anchor=%.2f, font=\\footnotesize] {$%.0f$};\n') %
                    (tX0, tY0, tX1, tY1,
                     anchor_angle, self.yticks[ntick]))
            ntick += 1

        # Z axis
        axZmesh = bpy.data.meshes['CETaxisMeshZ']
        axZxy = []
        for vert in axZmesh.vertices:
            axZxy.append(bpy_extras.object_utils.world_to_camera_view(self.scn,
                cam, vert.co))
        X0 = axZxy[0][0]
        X1 = axZxy[-5][0]
        Y0 = axZxy[0][1]
        Y1 = axZxy[-5][1]
        f.write('    \\draw[-latex] (%.3f,%.3f) -- (%.3f,%.3f);\n' %
                (X0*self.x_size, Y0*self.y_size,
                 X1*self.x_size, Y1*self.y_size))
        # Z-axis label
        X0 = axZxy[-3][0]
        X1 = axZxy[-1][0]
        Y0 = axZxy[-3][1]
        Y1 = axZxy[-1][1]
        label_angle = 180/np.pi * np.arctan2((Y1 - Y0), (X1 - X0))
        f.write(('    \\path (%.3f,%.3f)' +
                 'node[rotate = %.2f, font = \\footnotesize]' +
                 '{%s};\n') %
                 (axZxy[-2][0]*self.x_size, axZxy[-2][1]*self.y_size,
                  label_angle, self.zlabel))
        # Tick labels
        ntick = 0
        for n in range(2, len(axZxy)-5, 2):
            tX0 = axZxy[n][0] * self.x_size
            tX1 = axZxy[n+1][0] * self.x_size
            tY0 = axZxy[n][1] * self.y_size
            tY1 = axZxy[n+1][1] * self.y_size
            anchor_angle = 180/np.pi * np.arctan2((tY0 - tY1), (tX0 - tX1))
            f.write(('    \\draw (%.3f,%.3f) -- (%.3f,%.3f) ' +
                    'node[anchor=%.2f, font=\\footnotesize] {$%.1f$};\n') %
                    (tX0, tY0, tX1, tY1,
                     anchor_angle, self.zticks[ntick]/self.Zscale))
            ntick += 1

        # End of .tex file
        f.write(tw.dedent('''\
          \\end{scope}
        \\end{tikzpicture}
        \\end{frame}
        
        \\end{document}
        '''))
        f.close()
        print('Wrote TeX code to %s\n' % (tex_filename))
        # END WRITE .TEX FILE


        # Process .tex
        subprocess.call('pdflatex temp.tex'.split())

        # Crop .pdf, convert .pdf back to .png, and copy .tex and cropped .pdf
        # to final location
        subprocess.call('pdfcrop temp.pdf temp.pdf'.split())
        subprocess.call(['gs', '-dNOPAUSE', '-dBATCH', '-sDEVICE=pngalpha',
                         '-sOutputFile=%s' % filename, '-r%d' % self.dpi,
                         'temp.pdf'])
        subprocess.call(['mv', 'temp.tex', filename[:-3]+'tex'])
        subprocess.call(['mv', 'temp.pdf', filename[:-3]+'pdf'])

        # Restore unrotated plot and axes
        data_raw_save.to_mesh(self.data_raw.data)
        data_raw_save.free()
        data_save.to_mesh(self.data.data)
        data_save.free()
        data_wireframe_save.to_mesh(self.data_wireframe.data)
        data_wireframe_save.free()
        data_contours_save.to_mesh(self.data_contours.data)
        data_contours_save.free()
        axX_save.to_mesh(self.axX.data)
        axX_save.free()
        axY_save.to_mesh(self.axY.data)
        axY_save.free()
        axZ_save.to_mesh(self.axZ.data)
        axZ_save.free()

        # Restore camera
        cam.location = cam_loc_save
        cam.rotation_euler = cam_rot_save
        bpy.data.cameras['Camera.001'].type = 'PERSP'
        bpy.data.cameras['Camera.001'].lens = lens_save

        return


    def rotate_z(self, ob, center, angle):
        '''Rotate Blender object through 'angle' degrees about a vector passing
        through center parallel to the z axis.
        '''

        # Make a copy of the mesh from the object
        mesh = bmesh.new()
        mesh.from_mesh(ob.data)

        # Rotate the mesh
        bmesh.ops.rotate(mesh, cent=center,
                         matrix=mathutils.Matrix.Rotation(angle*np.pi/180, 3, 'Z'),
                         verts=mesh.verts)

        # Write the mesh back to the object
        mesh.to_mesh(ob.data)
        mesh.free()
        return


    def regrid(self, X, Y, Z):
        '''Interpolate onto an evenly-spaced grid in x and y, spacing dx.  Return a
        tuple of new (X, Y, Z).  For both input and output, X and Y are
        1D, and Z is 2D.
        '''

        newX = np.ogrid[X[0]:X[-1]:self.dx]
        newY = np.ogrid[Y[0]:Y[-1]:self.dx]

        # Iterate through rows of Z and interpolate using newX
        tempZ = np.zeros((len(Y), len(newX)))
        for y in range(len(Y)):
            tempZ[y,:] = np.interp(newX, X, Z[y,:])

        # Iterate through columns of Z and interpolate using newY
        newZ = np.zeros((len(newY), len(newX)))
        for x in range(len(newX)):
            newZ[:,x] = np.interp(newY, Y, tempZ[:,x])

        return newX, newY, newZ


    def get_axis(self, X, Y, Z):
        '''Return three objects with the positions of all the X, Y, and Z tick
        marks (and end points, and positions of axis labels) as
        vertices.  This is useful so that Blender methods for object
        rotation and determining 2D coordinates can be used on the axis
        points.
        '''

        self.zticks = [self.Zscale*z for z in self.zticks]
        Z = self.Zscale * Z

        axX = np.concatenate((np.array([X[0]]), self.xticks, np.array([X[-1]])))
        axY = np.concatenate((np.array([Y[0]]), self.yticks, np.array([Y[-1]])))
        axZ = np.concatenate((np.array([0]), self.zticks, np.array([Z.max()])))

        # X axis
        verts = []
        meshX = bpy.data.meshes.new('CETaxisMeshX')
        obX = bpy.data.objects.new('CETaxisX', meshX)
        obX.location = (0, 0, 0)
        obX.show_name = False
        for x in axX:
            verts.append([x, axY[0], axZ[0]])
            verts.append([x, axY[0]-self.ticklen, axZ[0]])
        verts.append([axX[0], axY[0]-8*self.ticklen, axZ[0]])
        verts.append([(axX[0]+axX[-1])/2, axY[0]-8*self.ticklen, axZ[0]])
        verts.append([axX[-1], axY[0]-8*self.ticklen, axZ[0]])
        meshX.from_pydata(verts, [], [])
        meshX.update()

        # Y axis
        verts = []
        meshY = bpy.data.meshes.new('CETaxisMeshY')
        obY = bpy.data.objects.new('CETaxisY', meshY)
        obY.location = (0, 0, 0)
        obY.show_name = False
        for y in axY:
            verts.append([axX[0], y, axZ[0]])
            verts.append([axX[0]-self.ticklen, y, axZ[0]])
        verts.append([axX[0]-8*self.ticklen, axY[0], axZ[0]])
        verts.append([axX[0]-8*self.ticklen, (axY[0]+axY[-1])/2, axZ[0]])
        verts.append([axX[0]-8*self.ticklen, axY[-1], axZ[0]])
        meshY.from_pydata(verts, [], [])
        meshY.update()

        # Z axis
        verts = []
        meshZ = bpy.data.meshes.new('CETaxisMeshZ')
        obZ = bpy.data.objects.new('CETaxisZ', meshZ)
        obZ.location = (0, 0, 0)
        obZ.show_name = False
        for z in axZ:
            verts.append([axX[-1], axY[0], z])
            verts.append([axX[-1], axY[0]-self.ticklen, z])
        verts.append([axX[-1], axY[0]-12*self.ticklen, axZ[0]])
        verts.append([axX[-1], axY[0]-12*self.ticklen, (axZ[0]+axZ[-1])/2])
        verts.append([axX[-1], axY[0]-12*self.ticklen, axZ[-1]])
        meshZ.from_pydata(verts, [], [])
        meshZ.update()

        return obX, obY, obZ


    def get_wireframe(self, X, Y, Z):
        '''Generate wireframe as a sparse version of X, Y, Z mesh (but connected
        along contour of surface)
        '''

        first = np.ceil(X[0]/self.dx)*self.dx
        last = np.floor(X[-1]/self.dx)*self.dx
        temp = np.ogrid[first:last:self.dx]
        newX = np.concatenate((np.array([X[0]]), temp, np.array([X[-1]])))

        first = np.ceil(Y[0]/self.dx)*self.dx
        last = np.floor(Y[-1]/self.dx)*self.dx
        temp = np.ogrid[first:last:self.dx]
        newY = np.concatenate((np.array([Y[0]]), temp, np.array([Y[-1]])))

        # Iterate through rows of Z and interpolate using newX
        tempZ = np.zeros((len(Y), len(newX)))
        for y in range(len(Y)):
            tempZ[y,:] = np.interp(newX, X, Z[y,:])

        # Iterate through columns of Z and interpolate using newY
        newZ = np.zeros((len(newY), len(newX)))
        for x in range(len(newX)):
            newZ[:,x] = np.interp(newY, Y, tempZ[:,x])

        # Get indices of multiples of DX in new X and Y arrays
        xI = []
        for n in range(newX.shape[0]):
            # For DX not an integer, the following is more likely to work
            # than just np.mod()
            if np.mod(np.round(np.mod(newX[n], self.DX), decimals=6), self.DX) == 0:
                xI.append(n)
        yI = []
        for n in range(newY.shape[0]):
            if np.mod(np.round(np.mod(newY[n], self.DY), decimals=6), self.DY) == 0:
                yI.append(n)

        # For new arrays, construct a list of vertices and edges
        verts = []
        nverts = -1
        edges = []

        for n in yI:
            for m in range(1, newX.shape[0]):
                if m == 1:
                    verts.append([newX[m-1], newY[n], newZ[n,m-1]])
                    nverts += 1
                verts.append([newX[m], newY[n], newZ[n,m]])
                nverts += 1
                edges.append([nverts-1, nverts])

        for m in xI:
            for n in range(1, newY.shape[0]):
                if n == 1:
                    verts.append([newX[m], newY[n-1], newZ[n-1,m]])
                    nverts += 1
                verts.append([newX[m], newY[n], newZ[n,m]])
                nverts += 1
                edges.append([nverts-1, nverts])

        # Make a new Blender mesh
        mesh = bpy.data.meshes.new('CETwireframeMesh')
        ob = bpy.data.objects.new('CETwireframe', mesh)
        ob.location = (0, 0, 0)
        ob.show_name = False
        mesh.from_pydata(verts, edges, [])
        mesh.update()

        # Scale z coordinates
        for v in mesh.vertices:
            v.co[2] *= self.Zscale

        # Add skin modifier
        mod = ob.modifiers.new('wireframeskin', 'SKIN')
        mod.use_smooth_shade = True
        for v in ob.data.skin_vertices[0].data:
            v.radius = 0.03, 0.03
            v.use_root = True

        return ob


    def get_mesh(self, X, Y, Z):
        '''Construct filled Blender mesh from Z array, using X and Y vectors for X
        and Y positions of points.  Return a Blender object.
        '''

        # Get an empty Blender mesh
        mesh = bpy.data.meshes.new('CETmesh')
        ob = bpy.data.objects.new('CET', mesh)
        ob.location = (0, 0, 0)
        ob.show_name = False

        # Iterate rows, then columns of Z, to add vertices to mesh
        verts = []
        for n in range(Z.shape[0]):
            for m in range(Z.shape[1]):
                verts.append([X[m], Y[n], Z[n,m]])

        # Iterate rows, then columns of Z, to add faces to mesh.  For
        # each point, use three nearest points (right, bottom, and
        # diagonal bottom right) to make face.  Therefore, stop one
        # from the last column and one from last row
        faces = []
        for n in range(Z.shape[0]-1):
            for m in range(Z.shape[1]-1):
                curr = m + n*Z.shape[1]
                faces.append([curr, curr+1, curr+Z.shape[1]+1, curr+Z.shape[1]])

        mesh.from_pydata(verts, [], faces)
        mesh.update()
        for v in mesh.vertices:
            v.co[2] *= self.Zscale
        return ob


    def get_raw(self, X, Y, Z):
        '''Construct a wireframe consisting just the rows of data, not connected to
        each other (useful for plotting a series of spectra, for example).
        '''

        mesh = bpy.data.meshes.new('CETrawMesh')
        ob = bpy.data.objects.new('CETraw', mesh)
        ob.location = (0, 0, 0)
        ob.show_name = False

        # Iterate rows, then columns of Z, to add vertices to mesh
        verts = []
        for n in range(Z.shape[0]):
            for m in range(Z.shape[1]):
                verts.append([X[m], Y[n], Z[n,m]])

        # Iterate rows then columns of Z, adding edges to mesh
        edges = []
        for n in range(Z.shape[0]):
            for m in range(Z.shape[1]-1):
                curr = m + n*Z.shape[1]
                edges.append([curr, curr+1])

        mesh.from_pydata(verts, edges, [])
        mesh.update()

        # Scale z coordinates
        for v in mesh.vertices:
            v.co[2] *= self.Zscale

        # Add skin modifier
        mod = ob.modifiers.new('rawdataskin', 'SKIN')
        mod.use_smooth_shade = True
        for v in ob.data.skin_vertices[0].data:
            v.radius = 0.05, 0.05
            v.use_root = True

        return ob


    def update_raw(self, X, Y, Z):
        '''Replace mesh in 'CETraw' with a new one from given data.  For
        convenience, put 'CETraw' in an active layer before returning.
        '''

        bpy.ops.object.select_all(action='DESELECT')
        self.data_raw.select = True
        bpy.ops.object.delete()
        
        self.data_raw = self.get_raw(X, Y, Z)
        self.scn.objects.link(self.data_raw)
        self.data_raw.select = True
        bpy.ops.object.shade_smooth()
        self.data_raw.data.materials.append(bpy.data.materials['Material.004'])

        self.set_layers(1, 2, 2, 2)

        return


    def get_contours(self, X, Y, Z):
        '''Contour plot.  Build a solid of the plot, and use the bisect tool to
        make 'nslices' (instance variable) equi-Z-value filled contours
        of the data.  An orthographic projection of this thing viewed
        from directly above should look exactly like a filled contour
        plot.
        '''

        print('Calculating contours ... ')

        # Make a solid mesh by connecting the data surface to a
        # projection of itself onto the z=0 plane
        solid_mesh = bpy.data.meshes.new('CETsolidMesh')
        solid_ob = bpy.data.objects.new('CETsolid', solid_mesh)
        solid_ob.location = (0, 0, 0)
        solid_ob.show_name = False
        
        # Add verts
        verts = []
        for n in range(Z.shape[0]):
            for m in range(Z.shape[1]):
                verts.append([X[m], Y[n], Z[n,m]*self.Zscale])
        # Add same verts at z=0
        for n in range(Z.shape[0]):
            for m in range(Z.shape[1]):
                verts.append([X[m], Y[n], -0.5])

        # Add data surface faces
        faces = []
        for n in range(Z.shape[0]-1):
            for m in range(Z.shape[1]-1):
                curr = m + n*Z.shape[1]
                faces.append([curr, curr+1, curr+Z.shape[1]+1, curr+Z.shape[1]])
        # Add z=0 faces
        offset = Z.size
        for n in range(Z.shape[0]-1):
            for m in range(Z.shape[1]-1):
                curr = m + n*Z.shape[1] + offset
                faces.append([curr, curr+1, curr+Z.shape[1]+1, curr+Z.shape[1]])
        # Add front and back faces
        c1 = 0
        c2 = offset-Z.shape[1]
        for n in range(Z.shape[1]-1):
            faces.append([c1, c1+1, c1+offset+1, c1+offset])
            faces.append([c2, c2+1, c2+offset+1, c2+offset])
            c1 += 1
            c2 += 1
        # Add side faces
        c1 = 0
        c2 = Z.shape[1]-1
        for n in range(Z.shape[0]-1):
            faces.append([c1, c1+Z.shape[1], c1+offset+Z.shape[1], c1+offset])
            faces.append([c2, c2+Z.shape[1], c2+offset+Z.shape[1], c2+offset])
            c1 += Z.shape[1]
            c2 += Z.shape[1]

        solid_mesh.from_pydata(verts, [], faces)
        solid_mesh.update()
        self.scn.objects.link(solid_ob)
            
        # Get slices
        slicesZ = self.Zscale * np.ogrid[0:Z.max():((Z.max()-0.0) / self.nslices)]
        verts = []
        edges = []
        faces = []
        index_offset = 0

        for n in range(1, slicesZ.size):

            solid_bmesh = bmesh.new()
            solid_bmesh.from_mesh(solid_mesh)
            bmesh.ops.bisect_plane(solid_bmesh,
                geom=solid_bmesh.verts[:]+solid_bmesh.edges[:]+solid_bmesh.faces[:],
                plane_co=mathutils.Vector([0, 0, slicesZ[n]]),
                plane_no=mathutils.Vector([0, 0, 1]),
                clear_outer=True,
                clear_inner=True)
            bmesh.ops.triangle_fill(solid_bmesh,
                                    use_beauty=True,
                                    edges=solid_bmesh.edges[:])

            # Add to verts, edges, and faces lists
            solid_bmesh.verts.index_update()
            for v in solid_bmesh.verts:
                verts.append([v.co.x, v.co.y, v.co.z])
            for e in solid_bmesh.edges:
                edges.append([e.verts[0].index+index_offset,
                              e.verts[1].index+index_offset])
            for f in solid_bmesh.faces:
                face = []
                for v in f.verts:
                    face.append(v.index+index_offset)
                faces.append(face)
            index_offset += len(solid_bmesh.verts)

            solid_bmesh.free()

        # New mesh and object for contours
        contours_mesh = bpy.data.meshes.new('CETcontourMesh')
        contours_ob = bpy.data.objects.new('CETcontour', contours_mesh)
        contours_ob.location = (0, 0, 0)
        contours_ob.show_name = False
        contours_mesh.from_pydata(verts, edges, faces)
        contours_mesh.update()

        bpy.ops.object.select_all(action='DESELECT')
        solid_ob.select = True
        bpy.ops.object.delete()
        bpy.data.meshes.remove(solid_mesh)

        return contours_ob


def main():

    return


if __name__ == '__main__':
    main()
