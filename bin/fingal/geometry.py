import numpy as np
import subprocess
from esys.finley import ReadGmsh
class MeshWithTopgraphy(object):
    def __init__(self, electrodes={}, recenter=False, basename = "tmp",
                 core_depth=40., extra_core=40., extra_padding=300.,
                 num_elements_outer_edge=10, mesh_size_core=None, rel_mesh_size_electrodes=0.1,
                 stationsFMT="s%s", GEOTOL=1e-8):
        """
        :param core_depth: depth of core relative to diameter of survey area in %
        :param extra_core: additional width of core region relative  to respective survey edge length in %
        :param extra_padding: additional padding of relative to respective core edge length in %
        :param num_elements_outer_edge: number of element on  the outer edges.
        :param
        """
        self.__electrodes=electrodes
        self.GEOTOL=GEOTOL
        self.stationsFMT = stationsFMT
        self.recenter = recenter
        self.core_depth=core_depth
        self.extra_core=extra_core
        self.extra_padding = extra_padding
        self.num_elements_outer_edge = num_elements_outer_edge
        self.mesh_size_core = mesh_size_core
        self.rel_mesh_size_electrodes = rel_mesh_size_electrodes
        
        self.positions_global = ( np.array([electrodes[s][0]  for s in electrodes]),
                                    np.array([electrodes[s][1]  for s in electrodes]) )
        if recenter:
            self.X0, self.rot_a = self.findSurveyCoordinateSystem(self.positions_global[0], self.positions_global[1], EPS=1e-10)
        else:
            self.X0, self.rot_a = np.zeros((2,)), 0.
        self.positions = self.toSurveyCoordinates(self.positions_global[0], self.positions_global[1])

        self.xmin = self.positions[0].min()
        self.xmax = self.positions[0].max()
        self.ymin = self.positions[1].min()
        self.ymax = self.positions[1].max()

        self.ebox = ( (self.xmax - self.xmin)**2 + (self.ymax - self.ymin)**2 )**0.5
        assert self.ebox > 0., "area of electrodes is zero."
        print("electrodes x range = ", self.xmin, self.xmax)
        print("           y range = ", self.ymin, self.ymax)
        print("          diameter = ", self.ebox )
        #================================
        self.x_core_extra = self.ebox  * (self.extra_core / 100.)
        self.y_core_extra = self.ebox  * (self.extra_core / 100.)
        self.zlengthCore = (self.ebox * self.core_depth) / 100.
        #if self.xmax - self.xmin < 1e-5 * self.ebox:
        #    self.x_core_extra = self.y_core_extra
        #if self.ymax - self.ymin < 1e-5 * self.ebox:
        #    self.y_core_extra = self.x_core_extra


        print("extra core space X  = ", self.x_core_extra)
        print("                 Y  = ", self.y_core_extra)
        print("                 Z  = ", self.zlengthCore)

        self.xlengthextraOuterBox = ((self.xmax - self.xmin + 2 * self.x_core_extra) * self.extra_padding) / 100.
        self.ylengthextraOuterBox = ((self.ymax - self.ymin + 2 * self.y_core_extra) * self.extra_padding) / 100.

        pz = min(self.xlengthextraOuterBox, self.ylengthextraOuterBox)
        self.xlengthextraOuterBox, self.ylengthextraOuterBox = pz, pz
        self.zlengthextraOuterBox = pz
        self.xminOuterBox = self.xmin - self.x_core_extra - self.xlengthextraOuterBox
        self.xmaxOuterBox = self.xmax + self.x_core_extra + self.xlengthextraOuterBox
        self.yminOuterBox = self.ymin - self.y_core_extra - self.ylengthextraOuterBox
        self.ymaxOuterBox = self.ymax + self.y_core_extra + self.ylengthextraOuterBox
        self.outer_mesh_size = max(self.xmaxOuterBox - self.xminOuterBox,
                                   self.ymaxOuterBox - self.yminOuterBox) / self.num_elements_outer_edge
        # ----------- make bounding box ---------------------------------------
        self.fr=1.1
        l = self.xmin - self.x_core_extra *  self.fr
        r = self.xmax + self.x_core_extra *  self.fr
        b = self.ymin - self.y_core_extra *  self.fr
        t = self.ymax + self.y_core_extra *  self.fr
        self.GlobalBoundingBox_Xmin, self.GlobalBoundingBox_Ymin = self.toGlobalCoordinates(l,b)
        self.GlobalBoundingBox_Xmax, self.GlobalBoundingBox_Ymax = self.GlobalBoundingBox_Xmin, self.GlobalBoundingBox_Ymin
        for x in [l ,r]:
            for y in [b,t]:
                XX=  self.toGlobalCoordinates(x,y)
                self.GlobalBoundingBox_Xmin = min(self.GlobalBoundingBox_Xmin, XX[0])
                self.GlobalBoundingBox_Xmax = max(self.GlobalBoundingBox_Xmax, XX[0])
                self.GlobalBoundingBox_Ymin = min(self.GlobalBoundingBox_Ymin, XX[1])
                self.GlobalBoundingBox_Ymax = max(self.GlobalBoundingBox_Ymax, XX[1])
        print(f"Global bounding box core = {self.GlobalBoundingBox_Xmin} - { self.GlobalBoundingBox_Xmax} X {self.GlobalBoundingBox_Ymin}-{self.GlobalBoundingBox_Ymax}")
        #===========================
        n = 0
        d_mean = 0
        d_min=self.ebox
        for i in range(self.positions[0].shape[0]-1):
            d = abs(self.positions[0][i]-self.positions[0][i+1:])+abs(self.positions[1][i]-self.positions[1][i+1:])
            d_mean+=1./d.min()
            d_min = min( d_min, d.min())
            n+=1
        self.d_mean = 1./( d_mean /n)
        self.d_min=d_min
        print("harmonic mean electrode distance = ", self.d_mean)
        print("minimal electrode distance = ", self.d_min)
        assert d_min >0 , "co-located electrodes found."

        if not self.mesh_size_core:
            self.mesh_size_core=self.d_mean
            self.mesh_size_electrodes = self.d_min * self.rel_mesh_size_electrodes
        else:
            self.mesh_size_electrodes = self.mesh_size_core *  self.rel_mesh_size_electrodes
        self.geofile_2D_flat = basename + '_2D_flat.geo'
        self.mshfile_2D_flat = basename + '_2D_flat.msh'
        self.mshfile_2D = basename + '_2D.msh'
        self.geofile_3D = basename + '.geo'
        self.mshfile_3D = basename + '.msh'
        # interpolation operator (x,y)->elevation
        self.interp = None

    def getElevation(self, x, y ):
        """

        """
        if self.interp:
            xg, yg = self.toGlobalCoordinates(x, y)
            return self.interp((xg, yg))
        else:
            if isinstance(x, np.nparray):
                return self.zlevel * np.ones(x.shape)
            else:
                return self.zlevel
    def setTopgraphyFromGrid(self, x=[], y=[], elevation=[], method="nearest"):
        """
        sets the interplation table for elevation over a rectangular grid with grid lines
        at x and y : elevation[i,j] at (x[i], y[j]) in global coordinates.

        locations of electrodes need to be within interval [ x[0], x[-1] ] and [y[0], y[-1]].
        locations outside this region (e.g. in the extra_padding area) are extrapolated.

        :param x: array of x coordinates of grid line in global coordinates
        :param y: array of y coordinates of grid line in global coordinates
        :param elevation: values of elevation at grid nodes.
        :param method: method of interpolation passed on to scipy.interpolate.RegularGridInterpolator
        """
        from scipy.interpolate import RegularGridInterpolator
        #if self.xmin < x[0] or x[-1] < self.xmax:
        #    raise ValueError("electrodes outside x-interpolation range")
        #if y[0] > self.ymin or y[-1] < self.ymax:
        #    raise ValueError("electrodes outside y-interpolation range")
        self.interp = RegularGridInterpolator((x, y), elevation, fill_value=None, bounds_error=False, method=method)
        self.positions_z=self.interp((self.positions_global[0], self.positions_global[1]))
        self.zlevel = np.mean(self.positions_z)
        self.zmin = self.positions_z.min()
        self.__updateGeometry()

    def setTopgraphyFromCloud(self, x=[], y=[], elevation=[], method="linear"):
        """
        sets the interplation table for elevation over point cloud: elevation[i] at (x[i], y[j]) in global coordinates.

        :param x: array of x coordinates of grid line in global coordinates
        :param y: array of y coordinates of grid line in global cordinates
        :param elevation: values of elevation at  points.
        :param method: method of interpolation
        """
        from scipy.interpolate import NearestNDInterpolator, LinearNDInterpolator
        idxX = np.logical_and(self.GlobalBoundingBox_Xmin <= x, x<=self.GlobalBoundingBox_Xmax)
        idxY = np.logical_and(self.GlobalBoundingBox_Ymin <= y, y <= self.GlobalBoundingBox_Ymax)
        idx= np.logical_and(idxX,idxY)

        zmean = np.mean(elevation)
        extraX = []
        extraZ = []
        for xx in [self.xminOuterBox, self.xmaxOuterBox]:
            for yy in [self.yminOuterBox, self.ymaxOuterBox]:
                extraX.append(self.toGlobalCoordinates(xx,yy))
                extraZ.append(zmean)
        xg = np.concatenate((x[idx], np.array([x for x,y in extraX]) ))
        yg = np.concatenate((y[idx], np.array([y for x, y in extraX]) ))
        eleg=np.concatenate((elevation[idx], np.array(extraZ)))
        if method == "linear":
            self.interp = LinearNDInterpolator((xg, yg), eleg)
        else:
            self.interp = NearestNDInterpolator((xg, yg), eleg)
        self.positions_z=self.interp((self.positions_global[0], self.positions_global[1]))
        self.zlevel = np.mean(self.positions_z)
        self.zmin = self.positions_z.min()
        self.__updateGeometry()

    def setFlatTopography(self, zlevel = 0):
        """
        set a flat geometry at level zlevel.
        overwrite any other topography seting.
        """
        self.interp = None
        self.zlevel = zlevel
        self.zmin = zlevel
        self.__updateGeometry()

    def __updateGeometry(self):

        print("elevation level at electrodes", self.zlevel)
        print("min. elevation at electrodes", self.zmin)
        # maybe it is good to use a trend plane:
        #fit, res, r, s = np.linalg.lstsq(
        #    np.stack((np.ones((len(self.positions[:, 0]),)), self.positions[:, 0], self.positions[:, 1])).T, self.positions_z,
        #    rcond=None)
        fit = np.array([self.zmin, 0,0])
        self.fit = fit

        self.zmincore = min(self.zmin,
                        self.getFlatZPosition(self.xmin - self.x_core_extra, self.ymin - self.y_core_extra),
                        self.getFlatZPosition(self.xmax + self.x_core_extra, self.ymin - self.y_core_extra),
                        self.getFlatZPosition(self.xmax + self.x_core_extra, self.ymax + self.y_core_extra),
                        self.getFlatZPosition(self.xmin - self.x_core_extra, self.ymax + self.y_core_extra))
        self.zminOuterBox = self.zmincore - self.zlengthCore - self.zlengthextraOuterBox
        self.zminCore = self.zmincore - self.zlengthCore
    def getFlatZPosition(self, x, y):
        return self.fit[0] + self.fit[1] * x + self.fit[2] * y

    def makeFlatSurfaceGeometry(self, geofile_2D_flat=None):
        """
        creates a flat 2D mesh including extra_padding region
        returns name of the 2D gmsh mesh file at self.mshfile_2D_flat
        """
        if geofile_2D_flat:
            self.geofile_2D_flat = geofile_2D_flat

        out = ""
        out += "Mesh.MshFileVersion = 2.2;\n"
        out += "// Core:\n"
        out += "xminCore = %s;\n" % (self.xmin - self.x_core_extra)
        out += "xmaxCore = %s;\n" % (self.xmax + self.x_core_extra)
        out += "yminCore = %s;\n" % (self.ymin - self.y_core_extra)
        out += "ymaxCore = %s;\n" % (self.ymax + self.y_core_extra)
        out += "zminCore = %s;\n" % self.zminCore
        out += "\n"
        out += "// big box\n"
        out += "xminOuterBox = %s;\n" % self.xminOuterBox
        out += "xmaxOuterBox = %s;\n" % self.xmaxOuterBox
        out += "yminOuterBox = %s;\n" % self.yminOuterBox
        out += "ymaxOuterBox = %s;\n" % self.ymaxOuterBox
        out += "zminOuterBox = %s;\n" % self.zminOuterBox
        out += "\n"

        out += "// element sizes\n"
        out += "meshSizeCore = %s;\n" % self.mesh_size_core
        out += "meshSizeOuterBox = %s;\n" % self.outer_mesh_size
        out += "meshSizeElectrodes = %s;\n" % self.mesh_size_electrodes
        out += "\n"
        out += "Point(1) = {xminCore, yminCore, %s, meshSizeCore};\n" % self.getFlatZPosition(self.xmin - self.x_core_extra, self.ymin - self.y_core_extra)
        out += "Point(2) = {xmaxCore, yminCore, %s, meshSizeCore};\n" % self.getFlatZPosition(self.xmax + self.x_core_extra, self.ymin - self.y_core_extra)
        out += "Point(3) = {xmaxCore, ymaxCore, %s, meshSizeCore};\n" % self.getFlatZPosition(self.xmax + self.x_core_extra, self.ymax + self.y_core_extra)
        out += "Point(4) = {xminCore, ymaxCore, %s, meshSizeCore};\n" % self.getFlatZPosition(self.xmin - self.x_core_extra, self.ymax + self.y_core_extra)
        out += "\n"
        out += "Point(5) = {xminCore, yminCore, zminCore, meshSizeCore};\n"
        out += "Point(6) = {xmaxCore, yminCore, zminCore, meshSizeCore};\n"
        out += "Point(7) = {xmaxCore, ymaxCore, zminCore, meshSizeCore};\n"
        out += "Point(8) = {xminCore, ymaxCore, zminCore, meshSizeCore};\n"
        out += "\n"
        out += "Point(9) = {xminOuterBox, ymaxOuterBox, %s, meshSizeOuterBox };\n" % self.getFlatZPosition(self.xminOuterBox, self.ymaxOuterBox)
        out += "Point(10) = {xminOuterBox, yminOuterBox, %s, meshSizeOuterBox };\n" % self.getFlatZPosition(self.xminOuterBox, self.yminOuterBox)
        out += "Point(11) = {xmaxOuterBox, yminOuterBox, %s, meshSizeOuterBox };\n" % self.getFlatZPosition(self.xmaxOuterBox, self.yminOuterBox)
        out += "Point(12) = {xmaxOuterBox, ymaxOuterBox, %s, meshSizeOuterBox };\n" % self.getFlatZPosition(self.xmaxOuterBox, self.ymaxOuterBox)
        out += "\n"
        out += "Point(17) = {xminOuterBox, ymaxOuterBox, zminOuterBox, meshSizeOuterBox };\n"
        out += "Point(18) = {xminOuterBox, yminOuterBox, zminOuterBox, meshSizeOuterBox };\n"
        out += "Point(19) = {xmaxOuterBox, yminOuterBox, zminOuterBox, meshSizeOuterBox };\n"
        out += "Point(20) = {xmaxOuterBox, ymaxOuterBox, zminOuterBox, meshSizeOuterBox };\n"
        out += "\n"
        out += "Line(1) = {10, 9};\n"
        out += "Line(3) = {10, 11};\n"
        out += "Line(4) = {11, 12};\n"
        out += "Line(5) = {12, 9};\n"
        out += "Line(6) = {4, 3};\n"
        out += "Line(7) = {3, 2};\n"
        out += "Line(8) = {2, 1};\n"
        out += "Line(9) = {1, 4};\n"
        out += "Line(10) = {4, 8};\n"
        out += "Line(11) = {9, 17};\n"
        out += "Line(12) = {12, 20};\n"
        out += "Line(13) = {3, 7};\n"
        out += "Line(14) = {2, 6};\n"
        out += "Line(15) = {6, 5};\n"
        out += "Line(16) = {1, 5};\n"
        out += "Line(17) = {5, 8};\n"
        out += "Line(18) = {17, 18};\n"
        out += "Line(19) = {18, 10};\n"
        out += "Line(20) = {18, 19};\n"
        out += "Line(21) = {19, 11};\n"
        out += "Line(22) = {20, 17};\n"
        out += "Line(23) = {20, 19};\n"
        out += "Line(24) = {8, 7};\n"
        out += "Line(25) = {7, 6};\n"
        out += "Line Loop(1) = {6, 13, -24, -10};\n"
        out += "Plane Surface(1) = {1};\n"
        out += "Line Loop(2) = {10, -17, -16, 9};\n"
        out += "Plane Surface(2) = {2};\n"
        out += "Line Loop(3) = {15, 17, 24, 25};\n"
        out += "Plane Surface(3) = {3};\n"
        out += "Line Loop(4) = {13, 25, -14, -7};\n"
        out += "Plane Surface(4) = {4};\n"
        out += "Line Loop(5) = {14, 15, -16, -8};\n"
        out += "Plane Surface(5) = {5};\n"
        out += "Line Loop(6) = {19, 1, 11, 18};\n"
        out += "Plane Surface(6) = {6};\n"
        out += "Line Loop(7) = {22, 18, 20, -23};\n"
        out += "Plane Surface(7) = {7};\n"
        out += "Line Loop(8) = {19, 3, -21, -20};\n"
        out += "Plane Surface(8) = {8};\n"
        out += "Line Loop(9) = {23, 21, 4, 12};\n"
        out += "Plane Surface(9) = {9};\n"
        out += "Line Loop(10) = {12, 22, -11, -5};\n"
        out += "Plane Surface(10) = {10};\n"
        out += "Line Loop(11) = {1, -5, -4, -3};\n"
        out += "Line Loop(12) = {9, 6, 7, 8};\n"
        out += "Plane Surface(11) = {12};\n"
        out += "Plane Surface(12) = {11, 12};\n"

        out += "// electrodes  (in the core)\n"
        out += "k=newp;\n"
        for i, s in enumerate(self.__electrodes.keys()):
            out += "Point(k+%s)={ %s, %s, %s, meshSizeElectrodes};\n" % (
                i + 1, self.positions[0][i], self.positions[1][i], self.getFlatZPosition( self.positions[0][i], self.positions[1][i],) )
            out += "Point{k+%s} In Surface{11};\n" % (i + 1)
            out += 'Physical Point("s%s")  = { k+%s } ;\n' % (s, i + 1)
        out += 'Physical Surface("innerfaces") = { 1, 2, 3, 4,5};\n'
        out += 'Physical Surface("faces") = { 6,7,8,9,10};\n'
        out += 'Physical Surface("surface") = { 11,12};\n'
        open(self.geofile_2D_flat, 'w').write(out)
        print(">> flat surface geometry was written to %s" % self.geofile_2D_flat)
        return self.geofile_2D_flat
    def makeFlatSurfaceMesh(self, geofile_2D_flat=None, mshfile_2D_flat=None):
        if geofile_2D_flat:
            self.geofile_2D_flat = geofile_2D_flat
        if mshfile_2D_flat:
            self.mshfile_2D_flat = mshfile_2D_flat

        print("generate surface mesh:")
        rp = subprocess.run(["gmsh", "-2", "-algo", "auto", "-o", self.mshfile_2D_flat, self.geofile_2D_flat ])
        rp.check_returncode()
        print(">> 2D GMSH mesh file %s generated." % self.mshfile_2D_flat)
        return self.mshfile_2D_flat
    def addTopologyTo2DMesh(self, mshfile_2D_flat=None, mshfile_2D=None ):
        """
        applies topology to flat mesh at mshfile_2D_flat and writes to mesh file mshfile_2D.
        if not given as arguments self.mshfile_2D_flat or self.mshfile_2D is used.
        the name of mshfile_2D is returned.
        """
        if mshfile_2D:
            self.mshfile_2D = mshfile_2D
        if mshfile_2D_flat:
            self.mshfile_2D_flat = mshfile_2D_flat
        fin = open(self.mshfile_2D_flat, 'r')
        fout = open(self.mshfile_2D, 'w')
        line = fin.readline()
        while line:
            if line.startswith("$NOD") or line.startswith("$Nodes"):
                process_nodes = True
                numNodes = int(fin.readline())
                if line.startswith("$NOD"):
                    fout.write("$NOD\n")
                else:
                    fout.write("$Nodes\n")
                fout.write("%d\n" % numNodes)
                print(numNodes, " nodes found.")
                cc = 0
                while cc < numNodes:
                    line = fin.readline()
                    i, x, y, z = line.split(" ")
                    level = self.fit[0]
                    zbase = self.zminOuterBox
                    #zbase = self.zminCore
                    if float(z) > zbase:
                        znew = (self.getElevation(float(x), float(y)) - zbase) * (float(z) - zbase) / (
                                    level - zbase) + zbase
                    else:
                        znew = float(z)
                    fout.write("%s %g %g %g\n" % (i, float(x), float(y), znew))
                    cc += 1
            else:
                fout.write(line)
            line = fin.readline()
    
        fin.close()
        fout.close()
        print(">> topography added to file %s and written to file %s." % (self.mshfile_2D_flat, self.mshfile_2D))
        return self.mshfile_2D

    def generate3DGeometry(self, mshfile_2D= None, geofile_3D=None):
        """

        """
        if geofile_3D:
            self.geofile_3D = geofile_3D
        if mshfile_2D:
            self.mshfile_2D = mshfile_2D
        out = ""
        out += "Mesh.MshFileVersion = 2.2;\n"
        out += "// Core:\n"
        out += "xminCore = %s;\n" % (self.xmin - self.x_core_extra)
        out += "xmaxCore = %s;\n" % (self.xmax + self.x_core_extra)
        out += "yminCore = %s;\n" % (self.ymin - self.y_core_extra)
        out += "ymaxCore = %s;\n" % (self.ymax + self.y_core_extra)
        out += "zminCore = %s;\n" % self.zminCore
        out += "// element sizes\n"
        out += "meshSizeCore = %s;\n" % self.mesh_size_core
        out += "meshSizeOuterBox = %s;\n" % self.outer_mesh_size
        out += "meshSizeElectrodes = %s;\n" % self.mesh_size_electrodes

        out += 'Merge "%s"; // merge modified msh\n' % self.mshfile_2D

        out += "\n"
        out += "Surface Loop(1) = {1,2,3,4,5,11};\n"
        out += "Volume(1) = {-1};\n"
        out += "Surface Loop(2) = {6, 8, 9, 7, 10, 2, 1, 12, 5, 4, 3};\n"
        out += "Volume(2) = {2};\n"
        out += 'Physical Volume("padding")  = { 2 } ;\n'
        out += 'Physical Volume("core")  = { 1 } ;\n'
        out += "\n"
        out += "Field[1] = Box;\n"
        out += "Field[1].VIn = meshSizeCore;\n"
        out += "Field[1].VOut = meshSizeOuterBox;\n"
        out += "Field[1].XMin = xminCore;\n"
        out += "Field[1].XMax = xmaxCore;\n"
        out += "Field[1].YMin = yminCore;\n"
        out += "Field[1].YMax = ymaxCore;\n"
        out += "Field[1].ZMin = zminCore;\n"
        out += "Field[1].ZMax = 1e55;\n"
        out += "Background Field = 1;\n"
        out += "\n"
        out += "Mesh.CharacteristicLengthExtendFromBoundary = 0;\n"
        out += "Mesh.CharacteristicLengthFromPoints = 0;\n"
        out += "Mesh.CharacteristicLengthFromCurvature = 0;\n"
        out += "Mesh.CharacteristicLengthMin=0.;\n"
        out += "Mesh.CharacteristicLengthMax=1e33;\n"

        open(self.geofile_3D, 'w').write(out)
        print(">> 3D geometry was written to %s" % self.geofile_3D)

    def generate3DMesh(self, mshfile_3D=None, geofile_3D=None):
        if mshfile_3D:
            self.mshfile_3D = mshfile_3D
        if geofile_3D:
            self.geofile_3D = geofile_3D
        rp = subprocess.run(["gmsh", "-3", "-optimize_netgen", "-algo", "auto", "-o", self.mshfile_3D, self.geofile_3D])
        rp.check_returncode()
        print(">> GMSH mesh file %s generated." % self.mshfile_3D)
        return self.mshfile_3D
    def toFlyFile(self, flyfile, mshfile_3D = None):
        if mshfile_3D:
            self.mshfile_3D = mshfile_3D
        dts = []
        dps = []
        for i, s in enumerate(self.__electrodes.keys()):
            if self.stationsFMT:
                dts.append(self.stationsFMT % s)
            else:
                dts.append(s)
            x, y = self.positions[0][i], self.positions[1][i],
            z = float(self.getElevation(x, y))
            dps.append([x, y, z])
        domain = ReadGmsh(self.mshfile_3D, 3, diracPoints=dps, diracTags=dts, optimize=True)
        domain.write(flyfile)
        print(">> Mesh written to fly file %s" % flyfile)
        return flyfile
    def makeFlyFile(self, flyfile):
        self.makeFlatSurfaceGeometry()
        self.makeFlatSurfaceMesh()
        self.addTopologyTo2DMesh()
        self.addTopologyTo2DMesh()
        self.generate3DGeometry()
        self.generate3DMesh()
        return self.toFlyFile(flyfile)
    def plotting(self, plotfile, all_domain = False, numPoints=200):
        import matplotlib.pyplot as plt
        if all_domain:
            X = np.linspace(self.xminOuterBox, self.xmaxOuterBox, num=numPoints)
            Y = np.linspace(self.yminOuterBox, self.ymaxOuterBox, num=numPoints)
        else:
            X = np.linspace(self.xmin - self.x_core_extra * self.fr, self.xmax + self.x_core_extra * self.fr, num=numPoints)
            Y = np.linspace(self.ymin - self.y_core_extra * self.fr, self.ymax + self.y_core_extra * self.fr, num=numPoints)

        X, Y = np.meshgrid(X, Y)
        Z = self.getElevation(X, Y)
        plt.figure()
        pcm = plt.pcolormesh(X, Y, Z[:, :])
        plt.scatter(self.positions[0], self.positions[1], s=4, c='r')
        plt.axis('scaled')
        plt.grid(visible=True, which='major', axis='both', color='k', linestyle='dotted', linewidth=0.5)
        plt.title("Electrodes + Elevation [m]")
        plt.xlabel("easting")
        plt.ylabel("northing")
        plt.minorticks_on()

        plt.colorbar(pcm)
        plt.savefig(plotfile)

        print("Topography pic written to %s." % plotfile)

    def findSurveyCoordinateSystem(self, xg, yg, EPS=1e-10):
        """
        find a transformation X = R * (Xg-Xg0) from the global coordinate system to an appropriate local survey coordinate system

        R = [ [cos(a), -sin(a) ], [sin(a), cos(a) ] ] where X0 is the mean of the station coordinates

        :param xg: global x coordinates of electrodes
        :type xg: `numpy.ndarray`
        :param yg: global y coordinates of electrodes
        :type yg: `numpy.ndarray`
        :return: offset Xg0, rotation angle a (in rad)
        :rtype: `tuple`, `float`
        """
        from math import atan, degrees, sin, cos
        xg0 = np.mean(xg)
        yg0 = np.mean(yg)
        print(f"survey center = {xg0}, {xg0}")

        x = xg - xg0
        y = yg - yg0

        m_xy = np.mean(x * y)
        m_x2 = np.mean(x ** 2)
        m_y2 = np.mean(y ** 2)
        print(f"m_xy = {m_xy}, m_x2={m_x2}, my_y2={m_y2}")
        if not max(m_x2, m_y2) > 0:
            a = 0.
        else:
            if abs(m_x2 - m_y2) < EPS * max(m_x2, m_y2):
                # divison by zero:
                if abs(m_xy) < EPS * max(m_x2, m_y2):
                    if m_x2 < m_y2:
                        a1 = 0.
                        a2 = np.pi / 2
                    else:
                        a1 = 0.
                        a2 = np.pi / 2
                else:
                    if (m_x2 - m_y2) * m_xy < 0:
                        a1 = - np.pi / 4
                        a2 = - np.pi / 4 + np.pi / 2
                    else:
                        a1 = np.pi / 4
                        a2 = np.pi / 4 + np.pi / 2
            else:
                a0 = atan(2 * m_xy / (m_x2 - m_y2)) / 2
                a1 = a0
                a2 = a0 + np.pi / 2
            if a1 > np.pi / 2:
                a1 = a1 - np.pi
            if a2 > np.pi / 2:
                a2 = a2 - np.pi
            print(f"rotation angle candidates are {degrees(a1)} and {degrees(a2)}")
            d1 = np.linalg.norm(x * cos(a1) + y * sin(a1))
            d2 = np.linalg.norm(x * cos(a2) + y * sin(a2))
            if d1 < d2:
                a = a2
            else:
                a = a1
        print("rotation angle = ", degrees(a))
        return np.array([xg0, yg0]), a
    def toSurveyCoordinates(self, xg, yg):
        rot_r = -self.rot_a
        x = np.cos(rot_r) * (xg - self.X0[0]) - np.sin(rot_r) * (yg - self.X0[1])
        y = np.sin(rot_r) * (xg - self.X0[0]) + np.cos(rot_r) * (yg - self.X0[1])
        return x, y
    def toGlobalCoordinates(self, x, y):
        rot_r = self.rot_a
        xg = np.cos(rot_r) * x - np.sin(rot_r) * y + self.X0[0]
        yg = np.sin(rot_r) * x + np.cos(rot_r) * y + self.X0[1]
        return xg, yg

    def writeStations(self, station_file, delimiter=','):
        """
        write station locations in survey coordinates including topograpy
        """
        f=open(station_file, 'w')
        for i, s in enumerate(self.__electrodes.keys()):
            x, y = self.positions[0][i], self.positions[1][i]
            z = float(self.getElevation(x, y))
            if abs(x) < self.ebox * self.GEOTOL:
                x=0
            if abs(y) < self.ebox * self.GEOTOL:
                y=0
            if abs(z) < self.ebox * self.GEOTOL:
                z=0

            f.write(f"{s}{delimiter} {x}{delimiter} {y}{delimiter} {z}\n")
        print("New electrodes were written to %s." % station_file)
