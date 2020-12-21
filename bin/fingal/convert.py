"""
this is class to convert VTU file into esys-escript FLY file

by Lutz Gross, Brisbane, Dec 2020.
"""
from vtk import *
import numpy as np

FINLEYELEMENTS = [ 'Point1', 'Line2', 'Line3', 'Line4', 'Tri3', 'Tri6', 'Tri9', 'Tri10', 'Rec4', 'Rec8', 'Rec9', 'Rec12', 'Rec16', 'Tet4', 'Tet10', 'Tet16', 'Hex8', 'Hex20', 'Hex27', 'Hex32', 'Line2Face', 'Line3Face', 'Line4Face', 'Tri3Face', 'Tri6Face', 'Tri9Face', 'Tri10Face', 'Rec4Face', 'Rec8Face', 'Rec9Face', 'Rec12Face', 'Rec16Face', 'Tet4Face', 'Tet10Face', 'Tet16Face', 'Hex8Face', 'Hex20Face', 'Hex27Face', 'Hex32Face', 'Point1_Contact', 'Line2_Contact', 'Line3_Contact', 'Line4_Contact', 'Tri3_Contact', 'Tri6_Contact', 'Tri9_Contact', 'Tri10_Contact', 'Rec4_Contact', 'Rec8_Contact', 'Rec9_Contact', 'Rec12_Contact', 'Rec16_Contact', 'Line2Face_Contact', 'Line3Face_Contact', 'Line4Face_Contact', 'Tri3Face_Contact', 'Tri6Face_Contact', 'Tri9Face_Contact', 'Tri10Face_Contact', 'Rec4Face_Contact', 'Rec8Face_Contact', 'Rec9Face_Contact', 'Rec12Face_Contact', 'Rec16Face_Contact', 'Tet4Face_Contact', 'Tet10Face_Contact', 'Tet16Face_Contact', 'Hex8Face_Contact', 'Hex20Face_Contact', 'Hex27Face_Contact', 'Hex32Face_Contact', 'Line3Macro', 'Tri6Macro', 'Rec9Macro', 'Tet10Macro', 'Hex27Macro' ]

    
VTK_ELEMENTS = {  VTK_VERTEX :                  {'flytype' : FINLEYELEMENTS.index('Point1'),
                                                 'dim' : 0, 
                                                 'vtkface' : None, 
                                                 'flycontact': FINLEYELEMENTS.index('Point1_Contact'), 
                                                 'tofly' : [0] }, 
                  VTK_LINE :                    {'flytype' : FINLEYELEMENTS.index('Line2'),
                                                 'dim' : 1, 
                                                 'vtkface' : VTK_VERTEX,
                                                 'flycontact': FINLEYELEMENTS.index('Point1_Contact'), 
                                                 'tofly' : [0, 1] }, 
                  VTK_TRIANGLE :                {'flytype' : FINLEYELEMENTS.index('Tri3'),
                                                 'dim' : 2, 
                                                 'vtkface' : VTK_LINE,
                                                 'flycontact': FINLEYELEMENTS.index('Line2_Contact'), 
                                                 'tofly' : [0, 1, 2] }, 
                  VTK_QUAD  :                   {'flytype' : FINLEYELEMENTS.index('Rec4'),
                                                 'dim' : 2, 
                                                 'vtkface' : VTK_LINE,
                                                 'flycontact': FINLEYELEMENTS.index('Line2_Contact'), 
                                                 'tofly' : [0, 1, 2, 3] }, 
                  VTK_TETRA  :                  {'flytype' : FINLEYELEMENTS.index('Tet4'),
                                                 'dim' : 3, 
                                                 'vtkface' : VTK_TRIANGLE,
                                                 'flycontact': FINLEYELEMENTS.index('Tri3_Contact'), 
                                                 'tofly' : [0, 1, 2, 3] }, 
                  VTK_HEXAHEDRON  :             {'flytype' : FINLEYELEMENTS.index('Hex8'),
                                                 'dim' : 4,
                                                 'vtkface' : VTK_QUAD,
                                                 'flycontact': FINLEYELEMENTS.index('Rec4_Contact'), 
                                                 'tofly' : [0, 1, 2, 3, 4, 5, 6, 7] }, 
                  VTK_QUADRATIC_EDGE  :         {'flytype' : FINLEYELEMENTS.index('Line4'),
                                                 'dim' : 1, 
                                                 'vtkface' : VTK_VERTEX,
                                                 'flycontact': FINLEYELEMENTS.index('Point1_Contact'), 
                                                 'tofly' : [0, 1, 2, 3] }, 
                  VTK_QUADRATIC_TRIANGLE  :     {'flytype' : FINLEYELEMENTS.index('Tri9'),
                                                 'dim' : 2, 
                                                 'vtkface' : VTK_QUADRATIC_EDGE,
                                                 'flycontact': FINLEYELEMENTS.index('Line3_Contact'), 
                                                 'tofly' : [0, 1, 2, 3, 4, 5, 6, 7, 8] }, 
                  VTK_QUADRATIC_QUAD  :         {'flytype' : FINLEYELEMENTS.index('Rec8'),
                                                 'dim' : 2, 
                                                 'vtkface' : VTK_QUADRATIC_EDGE,
                                                 'flycontact': FINLEYELEMENTS.index('Line3_Contact'), 
                                                 'tofly' : [0, 1, 2, 3, 4, 5, 6, 7] }, 
                  VTK_QUADRATIC_TETRA  :        {'flytype' : FINLEYELEMENTS.index('Tet10'),
                                                 'dim' : 3, 
                                                 'vtkface' : VTK_QUADRATIC_TRIANGLE, 
                                                 'flycontact': FINLEYELEMENTS.index('Tri6_Contact'), 
                                                 'tofly' : [0, 1, 2, 3, 4, 5, 6, 7, 8, 9] },
                  VTK_QUADRATIC_HEXAHEDRON  :   {'flytype' : FINLEYELEMENTS.index('Hex20') ,
                                                 'dim' : 3, 
                                                 'vtkface' : VTK_QUADRATIC_QUAD,
                                                 'flycontact': FINLEYELEMENTS.index('Rec8_Contact'), 
                                                 'tofly' : [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19] },
                  VTK_BIQUADRATIC_QUAD  :       {'flytype' : FINLEYELEMENTS.index('Rec9')  ,
                                                 'dim' : 2, 
                                                 'vtkface' : VTK_QUADRATIC_EDGE,
                                                 'flycontact': FINLEYELEMENTS.index('Line3_Contact'), 
                                                 'tofly' : [0, 1, 2, 3, 4, 5, 6, 7, 8] },
                  VTK_TRIQUADRATIC_HEXAHEDRON  :{'flytype' : FINLEYELEMENTS.index('Hex27') ,
                                                 'dim' : 3, 
                                                 'vtkface' : VTK_BIQUADRATIC_QUAD ,
                                                 'flycontact': FINLEYELEMENTS.index('Rec9_Contact'), 
                                                 'tofly' : [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25 ,26] } }

class VTK2FLY(object):
    """
    This is driver for converting vtk file to finley FLY file
    """
    def __init__(self, rocktype_tag='rocktype', node_id_tag="node number", elementid_tag="element number"):
        """
        initialize vtk->fly converter
        
        :param rocktype_tag: name of the rock type DataArray.
        :param node_id_tag: name of the node id  DataArray.
        :param elementid_tag: name of the element id DataArray.
        """
        self.__vtkfilename=None
        self.__flyfilename=None
        self.setFlipYZ()
        self.ROCKTYPE_ID_TAG_IN_VTK=rocktype_tag
        self.NODE_ID_TAG_IN_VTK=node_id_tag
        self.ELEMENT_ID_TAG_IN_VTK=elementid_tag
        
    def setVTKFileName(self, vtkfilename):
        """
        Set the file name for the vtk input file name
        """
        self.__vtkfilename=vtkfilename
        return self
        
    def getVTKFileName(self):
        """
        Get the vtk file name
        """
        return self.__vtkfilename
        
    def setFLYFileName(self, flyfilename):
        """
        Set the file name for the FLY output file name
        """
        self.__flyfilename=flyfilename
        return self
        
    def getFLYFileName(self):
        """
        Get the fly output file name
        """
        return self.__flyfilename                    

    def setFlipYZ(self, flag=False):
        """
        set flag to flip y and z direction
        """
        self.__flipyz =flag  
        return self
    
    def getFlipYZ(self, flag=False):
        """
        Get the fly output file name
        """
        return self.__flipyz 
    
    def convertVTU(self):
        """
        converts a VTU file
        """
        if not self.getFLYFileName():
            raise ValueError("No FLY output file set.")
        if not self.getVTKFileName():
            raise ValueError("No vtk input file set.")
        # Read the source file.
        reader = vtkXMLUnstructuredGridReader()
        reader.SetFileName(self.getVTKFileName())
        reader.Update() # Needed because of GetScalarRange
        inputdb = reader.GetOutput()


        numNodes=inputdb.GetNumberOfPoints()
        # make sure that there are no colocated nodes:
#        X=np.zeros((numNodes,), dtype=[('x',float), ('y',float), ('z',float) ] )
#        for i in xrange(numNodes):
#            X[i]=inputdb.GetPoint(i)
#        print X
#        ind = np.argsort(X, order=('x','y','z') )
#        ind2=[0]
#        D2=(X[ind[0]]['x']-X[ind[-1]]['x'])**2+(X[ind[0]]['y']-X[ind[-1]]['y'])**2+(X[ind[0]]['z']-X[ind[-1]]['z'])**2
#        print "Domain diameter is ",np.sqrt(D2)
#        if not D2 >0 :
#            raise ValueError("domain diameter is zero.")
#        TOL2=(1e-8)**2
#        for i in xrange(1,len(ind)):
#            r2=(X[ind2[-1]]['x']-X[ind[i]]['x'])**2+(X[ind2[-1]]['y']-X[ind[i]]['y'])**2+(X[ind2[0]]['z']-X[ind[i]]['z'])**2
#            if r2 > D2 * TOL2:
#                ind2.append(i)
#                
#        print "number of removed nodes :",numNodes-len(ind2)
        #print ind2

        
        # count elements:
        cell_counter={}
        for et in VTK_ELEMENTS:
            cell_counter[et]=0
        for i in xrange(inputdb.GetNumberOfCells()):
            c=inputdb.GetCell(i)
            if VTK_ELEMENTS.has_key(c.GetCellType()):
                cell_counter[c.GetCellType()]+=1
            else:
                raise ValueError("Unable to process cell %d of type %d."%(i,c.GetCellType()))

        # identify spatial dimension, inner and surface elements:
        DIM=-1
        VTK_INNER=None
        for et in VTK_ELEMENTS:
            if cell_counter[et]>0:
                if DIM < VTK_ELEMENTS[et]['dim']:
                    DIM = VTK_ELEMENTS[et]['dim']
                    VTK_INNER=et
        if DIM < 2:
            raise ValueError("No elements to import.")
        VTK_FACE=VTK_ELEMENTS[VTK_INNER]['vtkface']
        
        #... do we have a node id map?
        data=inputdb.GetPointData()
        datanames=[ data.GetArrayName(ia) for ia in xrange(data.GetNumberOfArrays()) ]
        if self.NODE_ID_TAG_IN_VTK in datanames:
            nodeids_dataArray=data.GetArray(datanames.index(self.NODE_ID_TAG_IN_VTK))
            if not nodeids_dataArray.GetNumberOfComponents() == 1:
                raise ValueError("Data array of '%s' is not a scalar"%NODE_ID_TAG_IN_VTK) 
        else:
            nodeids_dataArray=None
        #... do we have an element id map?
        data=inputdb.GetCellData()
        datanames=[ data.GetArrayName(ia) for ia in xrange(data.GetNumberOfArrays()) ]
        if self.ELEMENT_ID_TAG_IN_VTK in datanames:
            elementids_dataArray=data.GetArray(datanames.index(self.ELEMENT_ID_TAG_IN_VTK))
            if not elementids_dataArray.GetNumberOfComponents() == 1:
                raise ValueError("Data array of '%s' is not a scalar"%ELEMENT_ID_TAG_IN_VTK)
        else:
            elementids_dataArray=None
        #... do we have a rock type map?
        data=inputdb.GetCellData()
        datanames=[ data.GetArrayName(ia) for ia in xrange(data.GetNumberOfArrays()) ]
        if self.ROCKTYPE_ID_TAG_IN_VTK in datanames:
            rocktypes_dataArray=data.GetArray(datanames.index(self.ROCKTYPE_ID_TAG_IN_VTK))
            if not rocktypes_dataArray.GetNumberOfComponents() == 1:
                raise ValueError("Data array of '%s' is not a scalar"%ROCKTYPE_ID_TAG_IN_VTK)
        else:
            rocktypes_dataArray=None
        
        print "Spatial dimension = ",DIM
        print "Inner element type = ",VTK_INNER
        print "Inner elements found = ",cell_counter[VTK_INNER]
        print "Face element type = ",VTK_FACE
        print "Face elements found = ",cell_counter[VTK_FACE]
        print "Point sources found = ",cell_counter[VTK_VERTEX]
        if rocktypes_dataArray:
            print "rock type found!"
        if elementids_dataArray:
            print "element ids found!"
        else:
            print "WARNING: No element ids found"
        if nodeids_dataArray:
            print "node ids found!"
        else:
            print "WARNING: No node ids found"
            
        if VTK_INNER in [VTK_BIQUADRATIC_QUAD, VTK_TRIQUADRATIC_HEXAHEDRON ]:
                raise Error(">>>>>> Check node ordering for VTK element %s"%(VTK_INNER,))

        # ====== create a fly file ===============================
        fly=open(self.getFLYFileName(),'w')
        fly.write(self.getFLYFileName()+"\n")
        # ... nodes:
        fly.write("%1dD-Nodes %d\n"%(DIM,numNodes))
        nodeidmap={}
        for i in xrange(numNodes):
            if nodeids_dataArray:
                idx=int(nodeids_dataArray.GetTuple1(i)+0.5)
            else:
                idx=i
            ## @DEBUG1057642L, 1057643L
            #if i in [ 1438023, 1438022, 1057642, 1057643]:  # 1069180,  1459800
            #    print i,  idx, ":", nodeids_dataArray.GetTuple1(i-2), nodeids_dataArray.GetTuple1(i-1), nodeids_dataArray.GetTuple1(i), nodeids_dataArray.GetTuple1(i+1), nodeids_dataArray.GetTuple1(i+2)
            #if nodeidmap.values().count(idx):
            #    print("%s-th node id %s already exists (VTK value=%s)."%(i,idx, nodeids_dataArray.GetTuple1(i)))
                
            nodeidmap[i]=idx
            #if idx ==0:
            #        print "entry ",i, "ref :", idx, "@", inputdb.GetPoint(i)
            if self.getFlipYZ():
                pp=inputdb.GetPoint(i)
                fly.write(("%d %d %d "+"%13.15e "*DIM+"\n")%(idx,idx,0, pp[0], pp[2], pp[1]))
            else:
                fly.write(("%d %d %d "+"%13.15e "*DIM+"\n")%((idx,idx,0)+inputdb.GetPoint(i)[:DIM]))
        # ... inner elements:
        fly.write("%s %d\n"%(FINLEYELEMENTS[VTK_ELEMENTS[VTK_INNER]['flytype']],cell_counter[VTK_INNER]))
        nodeidx=VTK_ELEMENTS[VTK_INNER]['tofly']
        numNodes=len(nodeidx)
        domId=0 
        tags_found=set()
        for i in xrange(inputdb.GetNumberOfCells()):
            c=inputdb.GetCell(i)
            if c.GetCellType() == VTK_INNER:
                if rocktypes_dataArray:
                    domId=int(rocktypes_dataArray.GetTuple1(i)+0.5)
                    tags_found.add(domId)
                if elementids_dataArray:
                    idx=int(elementids_dataArray.GetTuple1(i)+0.5)
                else:
                    idx=i
                ### @DEBUG
                if nodeidx[3] == nodeidx[2]:
                    print("nodeidx failed in ", idx,domId, nodeidx)
                if c.GetPointId(nodeidx[3]) == c.GetPointId(nodeidx[2]):
                    print("GetPointId failed in ", idx,domId, str([ c.GetPointId(n) for n in nodeidx ]))
                if nodeidmap[c.GetPointId(nodeidx[3])] == nodeidmap[c.GetPointId(nodeidx[2])]:
                    print("nodeidmap failed in ", idx,domId, ":", str([ c.GetPointId(n) for n in nodeidx ]), "->", str([ nodeidmap[c.GetPointId(n)] for n in nodeidx ]))
                ### END DEBUG
                fly.write(("%d %d "+"%d "*numNodes+"\n")%((idx,domId)+tuple([ nodeidmap[c.GetPointId(n)] for n in nodeidx ]) ))
        print "Inner element tags =",  tags_found
        # face elemets:
        fly.write("%s %d\n"%(FINLEYELEMENTS[VTK_ELEMENTS[VTK_FACE]['flytype']],cell_counter[VTK_FACE]))
        nodeidx=VTK_ELEMENTS[VTK_FACE]['tofly']
        numNodes=len(nodeidx)
        domId+=1
        tags_found=set()
        for i in xrange(inputdb.GetNumberOfCells()):
            c=inputdb.GetCell(i)
            if c.GetCellType() == VTK_FACE:
                if rocktypes_dataArray:
                    domId=int(rocktypes_dataArray.GetTuple1(i)+0.5)
                    tags_found.add(domId)
                if elementids_dataArray:
                    idx=int(elementids_dataArray.GetTuple1(i)+0.5)
                else:
                    idx=i
                fly.write(("%d %d "+"%d "*numNodes+"\n")%((idx,domId)+tuple([ nodeidmap[c.GetPointId(n)] for n in nodeidx ]) ))
        print "Face element tags = ",tags_found
        # contact elements:
        fly.write("%s %d\n"%(FINLEYELEMENTS[VTK_ELEMENTS[VTK_INNER]['flycontact']],0))

        # nodal elements:
        fly.write("%s %d\n"%(FINLEYELEMENTS[VTK_ELEMENTS[VTK_VERTEX]['flytype']],cell_counter[VTK_VERTEX]))
        nodeidx=VTK_ELEMENTS[VTK_VERTEX]['tofly']
        numNodes=len(nodeidx)
        domId+=1
        tags_found=set()
        for i in xrange(inputdb.GetNumberOfCells()):
            c=inputdb.GetCell(i)
            if c.GetCellType() == VTK_VERTEX:
                if rocktypes_dataArray:
                    domId=int(rocktypes_dataArray.GetTuple1(i)+0.5)
                    tags_found.add(domId)
                if elementids_dataArray:
                    idx=int(elementids_dataArray.GetTuple1(i)+0.5)
                else:
                    idx=i
                fly.write(("%d %d "+"%d "*numNodes+"\n")%((idx,domId)+tuple([ nodeidmap[c.GetPointId(n)] for n in nodeidx ]) ))
        print "Nodal element tags = ",tags_found
        # write tags:
        fly.write("Tags\n")
        fly.close()
        return self
        
class ExtractFromVTK(object):
    """
    a reader to extract a data set (Point or Cell) from a VTK (VTU) file_name
    to generate a corresponding escript data object
    """
    NODE_ID_TAG_IN_VTK="node number"
    elementid_TAG_IN_VTK="element number"
    def __init__(self, domain, node_id_tag="node number", elementid_tag="element number"):
        """
        opens the reader for a domain. Typically the domain is generated via ``VTK2FLY``
        """
        self.NODE_ID_TAG_IN_VTK=node_id_tag
        self.ELEMENT_ID_TAG_IN_VTK=elementid_tag
        self.__domain=domain
        self.__vtkfn=None
        self.__vtkdb=None
    def getDomain(self):
        """
        return target domain
        """
        return self.__domain
    def setVTKFileName(self, filename):
        """
        set the vtk file name and opens the vtk data file for reading
        """
        self.__vtkfn=filename
        self.__nodeidmap=None
        self.__elementidmap=None
        reader = vtkXMLUnstructuredGridReader()
        reader.SetFileName(filename)
        reader.Update() # Needed because of GetScalarRange
        self.__vtkdb = reader.GetOutput()
        print "VTK file %s ready for reading."%(filename,)
        # do we have node ids:
        data=self.__vtkdb.GetPointData()
        datanames=[ data.GetArrayName(ia) for ia in xrange(data.GetNumberOfArrays()) ]
        if self.NODE_ID_TAG_IN_VTK in datanames:
            print "point ids found"
            nodeids_dataArray=data.GetArray(datanames.index(self.NODE_ID_TAG_IN_VTK))
            if not nodeids_dataArray.GetNumberOfComponents() == 1:
                raise ValueError("Data array of node ids is not a scalar") 
            nodeidmap={}
            for i in xrange(self.__vtkdb.GetNumberOfPoints()):
                idx=int(nodeids_dataArray.GetTuple1(i)+0.5)
                #if idx ==0:
                #    print "mapping entry ",i, "ref :", idx, "@", self.__vtkdb.GetPoint(i)
                nodeidmap[idx] = i
            self.__nodeidmap=nodeidmap
        else:
            print "WARNING: no node ids found"
        # do we have element ids:
        data=self.__vtkdb.GetCellData()
        datanames=[ data.GetArrayName(ia) for ia in xrange(data.GetNumberOfArrays()) ]
        if self.elementid_TAG_IN_VTK in datanames:
            print "element ids found"
            elementid_dataArray=data.GetArray(datanames.index(self.elementid_TAG_IN_VTK))
            if not elementid_dataArray.GetNumberOfComponents() == 1:
                raise ValueError("Data array of element ids is not a scalar") 
            elementidmap={}
            for i in xrange(self.__vtkdb.GetNumberOfCells()):
                idx=int(elementid_dataArray.GetTuple1(i)+0.5)
                elementidmap[idx] = i
            self.__elementidmap=elementidmap
        else:
             print "WARNING: no element ids found"
        return self
    def getVTKFileName(self):
        """
        returns the VTK file name
        """
        return self.__vtkfn
    def getData(self, dname):
        """
        return the data with given name dname from the VTK file.
        The function returns a ``Scalar``, ``Vector`` or ``Tensor`` escript Data object with
        the appropriate ``FunctionSpace``.
        """
        from esys.escript import Scalar, Vector, ReducedFunction, ContinuousFunction
        print "searching for data set `%s`"%(dname,)
        edata=None
        # check cell data:
        if edata is None:
            data=self.__vtkdb.GetCellData()
            for ia in xrange(data.GetNumberOfArrays()):
                if data.GetArrayName(ia) == dname:
                    a=data.GetArray(ia)
                    if a.GetNumberOfComponents() == 1:

                        edata=Scalar(0.,ReducedFunction(self.getDomain()))
                    elif a.GetNumberOfComponents() == 3:
                        edata=Vector(0., ReducedFunction(self.getDomain()))
                    else:
                        raise ValueError("unable to load %s components."%(a.GetNumberOfComponents(),))
                    print  "'%s' imported as cell data with shape %s."%(dname,str(edata.getShape()))
                    indxmap=self.__elementidmap
                    ndata=self.__vtkdb.GetNumberOfPoints()
        # check point data:
        if edata is None:
            data=self.__vtkdb.GetPointData()
            for ia in xrange(data.GetNumberOfArrays()):
                if data.GetArrayName(ia) == dname:
                    a=data.GetArray(ia)
                    if a.GetNumberOfComponents() == 1:
                        edata=Scalar(0.,ContinuousFunction(self.getDomain()))
                    elif a.GetNumberOfComponents() == 3:
                        edata=Vector(0., ContinuousFunction(self.getDomain()))
                    else:
                        raise ValueError("unable to load %s components."%(a.GetNumberOfComponents(),))
                    print  "'%s' imported as point data with shape %s."%(dname,str(edata.getShape()))
                    indxmap=self.__nodeidmap
                    ndata=self.__vtkdb.GetNumberOfCells()
        # all went wrong
        if edata is None:
            raise ValueError("unable to find data for '%s' in file %s."%(dname, self.getVTKFileName()))
    
        #edata.expand()
        if indxmap:
            if a.GetNumberOfComponents() == 1:
                for j in xrange(edata.getNumberOfDataPoints()):
                    idx=edata.getFunctionSpace().getReferenceIDFromDataPointNo(j)
                    if idx in indxmap:
                        edata.setValueOfDataPoint(j, a.GetTuple1(indxmap[idx]))
            else:
                for j in xrange(edata.getNumberOfDataPoints()):
                    idx=edata.getFunctionSpace().getReferenceIDFromDataPointNo(i)
                    if idx in indxmap:
                        edata.setValueOfDataPoint(j, a.GetTuple(indxmap[idx]))
        else:
            if a.GetNumberOfComponents() == 1:
                for j in xrange(min(edata.getNumberOfDataPoints(), ndata)):
                   edata.setValueOfDataPoint(j, a.GetTuple1(j))
            else:
                for j in xrange(min(edata.getNumberOfDataPoints(), ndata)):
                   edata.setValueOfDataPoint(j, a.GetTuple(j))
            
        return edata
