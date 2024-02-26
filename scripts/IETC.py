# state file generated using paraview version 5.12.0-RC1
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 12

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [3058, 1188]
renderView1.AxesGrid = 'Grid Axes 3D Actor'
renderView1.CenterOfRotation = [0.0, 0.0, 13.0]
renderView1.UseAmbientOcclusion = 1
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [6.979554540305161, -62.6840105805603, -54.04532959065491]
renderView1.CameraFocalPoint = [6.979548109885778, -62.68393323785537, -54.045254376087904]
renderView1.CameraViewUp = [-0.00808094212739432, -0.6974974366701516, 0.7165417114257212]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 30.907118921051183
renderView1.LegendGrid = 'Legend Grid Actor'
renderView1.UseColorPaletteForBackground = 0
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.BackEnd = 'OSPRay raycaster'

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(3058, 1188)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

hex_layout_skeletons = []

for i in range(1,41):
    hex_layout_skeletons.append( XMLUnstructuredGridReader(registrationName=f"is{i}", FileName=[f"/Users/caleb/sweeps/build/hex_layout_skeleton{i}.vtu"]) )
    # Show( hex_layout_skeletons[-1], renderView1, 'UnstructuredGridRepresentation')

for i in range(0, 16):
    hex_layout_skeletons.append( XMLUnstructuredGridReader(registrationName=f"os{i}", FileName=[f"/Users/caleb/sweeps/build/outer_surface{i}.vtu"]) )
    # Show( hex_layout_skeletons[-1], renderView1, 'UnstructuredGridRepresentation')

vol_ii = 0
def build_solid( surf_id_1 : int, surf_id_2 : int, surf_id_3 : int, surf_id_4 : int, func : str, vol_ii : int ):
    ad = AppendDatasets(registrationName=f"AD{vol_ii}", Input=[hex_layout_skeletons[surf_id_1-1], hex_layout_skeletons[surf_id_2-1], hex_layout_skeletons[surf_id_3-1], hex_layout_skeletons[surf_id_4-1]])
    delaunay3D = Delaunay3D(registrationName=f"D3D{vol_ii}", Input=ad)
    delaunay3D.Alpha = 2.0
    delaunay3D.AlphaTris = 0
    groupDatasets = GroupDatasets(registrationName=f"GD{vol_ii}", Input=[ad, delaunay3D])
    calculator = Calculator(registrationName=f"C{vol_ii}", Input=groupDatasets)
    calculator.Function = func
    warpByVector = WarpByVector(registrationName=f"WBV{vol_ii}", Input=calculator)
    warpByVector.Vectors = ['POINTS', 'Result']
    Show( warpByVector, renderView1, 'UnstructuredGridRepresentation' )

build_solid( 13, 14, 41, 22, '1.5*jHat', vol_ii )

build_solid( 21, 22, 53, 40, '1.5*(0.5*(coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))+sqrt(3)/2*jHat)', vol_ii )
build_solid( 32, 40, 44, 31, '1.5*(sqrt(3)/2*(coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))+0.5*jHat)', vol_ii )

build_solid( 30, 31, 46, 29, '1.5*(coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))', vol_ii )

build_solid( 28, 29, 48, 33, '1.5*(sqrt(3)/2*(coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))-0.5*jHat)', vol_ii )
build_solid( 35, 34, 33, 49, '1.5*(0.5*(coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))-sqrt(3)/2*jHat)', vol_ii )

build_solid( 35, 36, 37, 43, '-1.5*jHat', vol_ii )

build_solid( 55, 39, 38, 37, '1.5*(-0.5*(coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))-sqrt(3)/2*jHat)', vol_ii )
build_solid( 8, 9, 39, 56, '1.5*(-sqrt(3)/2*(coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))-0.5*jHat)', vol_ii )

build_solid( 9, 17, 18, 50, '-1.5*(coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))', vol_ii )

build_solid( 18, 54, 2, 1, '1.5*(-sqrt(3)/2*(coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))+0.5*jHat)', vol_ii )
build_solid( 2, 51, 14, 3, '1.5*(-0.5*(coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))+sqrt(3)/2*jHat)', vol_ii )

build_solid( 12, 13, 15, 20, 'jHat', vol_ii )
build_solid( 32, 21, 20, 24, '1/sqrt(2)*((coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))+jHat)', vol_ii )
build_solid( 24, 30, 26, 19, '(coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))', vol_ii )
build_solid( 34, 25, 26, 28, '1/sqrt(2)*((coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))-jHat)', vol_ii )
build_solid( 10, 25, 36, 5, '-jHat', vol_ii )
build_solid( 5, 38, 8, 7, '1/sqrt(2)*(-(coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))-jHat)', vol_ii )
build_solid( 7, 4, 16, 17, '-(coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))', vol_ii )
build_solid( 1, 3, 12, 16, '1/sqrt(2)*(-(coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))+jHat)', vol_ii )

build_solid( 11, 15, 23, 47, '0.5*jHat', vol_ii )
build_solid( 23, 19, 27, 42, '0.5*(coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))', vol_ii )
build_solid( 6, 10, 27, 45, '-0.5*jHat', vol_ii )
build_solid( 6, 4, 11, 52, '-0.5*(coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))', vol_ii )




# # create a new 'XML Unstructured Grid Reader'
# outer_surface1vtu = XMLUnstructuredGridReader(registrationName='outer_surface1.vtu', FileName=['/Users/caleb/sweeps/build/outer_surface1.vtu'])
# outer_surface1vtu.TimeArray = 'None'

# # create a new 'Append Datasets'
# appendDatasets10 = AppendDatasets(registrationName='AppendDatasets10', Input=[hex_layout_skeleton19vtu_1, hex_layout_skeleton23vtu, hex_layout_skeleton27vtu, outer_surface1vtu])

# # create a new 'Delaunay 3D'
# delaunay3D10 = Delaunay3D(registrationName='Delaunay3D10', Input=appendDatasets10)
# delaunay3D10.Alpha = 2.0
# delaunay3D10.AlphaTris = 0

# # create a new 'XML Unstructured Grid Reader'
# hex_layout_skeleton8vtu = XMLUnstructuredGridReader(registrationName='hex_layout_skeleton8.vtu', FileName=['/Users/caleb/sweeps/build/hex_layout_skeleton8.vtu'])
# hex_layout_skeleton8vtu.TimeArray = 'None'

# # create a new 'Append Datasets'
# appendDatasets6 = AppendDatasets(registrationName='AppendDatasets6', Input=[hex_layout_skeleton38vtu, hex_layout_skeleton5vtu, hex_layout_skeleton7vtu, hex_layout_skeleton8vtu])

# # create a new 'Clip'
# clip7 = Clip(registrationName='Clip7', Input=appendDatasets6)
# clip7.ClipType = 'Plane'
# clip7.HyperTreeGridClipper = 'Plane'
# clip7.Scalars = ['POINTS', '']

# # init the 'Plane' selected for 'ClipType'
# clip7.ClipType.Origin = [0.014100074768066406, -0.36231493949890137, 0.01]
# clip7.ClipType.Normal = [0.0, 0.0, 1.0]

# # init the 'Plane' selected for 'HyperTreeGridClipper'
# clip7.HyperTreeGridClipper.Origin = [0.014100074768066406, -0.36231493949890137, 11.956950187683105]

# # create a new 'Delaunay 3D'
# delaunay3D6 = Delaunay3D(registrationName='Delaunay3D6', Input=clip7)
# delaunay3D6.Alpha = 2.0
# delaunay3D6.Tolerance = 0.01

# # create a new 'Clip'
# clip8 = Clip(registrationName='Clip8', Input=appendDatasets7)
# clip8.ClipType = 'Plane'
# clip8.HyperTreeGridClipper = 'Plane'
# clip8.Scalars = ['POINTS', '']

# # init the 'Plane' selected for 'ClipType'
# clip8.ClipType.Origin = [0.014100074768066406, -0.36231493949890137, 0.01]
# clip8.ClipType.Normal = [0.0, 0.0, 1.0]

# # init the 'Plane' selected for 'HyperTreeGridClipper'
# clip8.HyperTreeGridClipper.Origin = [0.014100074768066406, -0.36231493949890137, 11.956950187683105]

# # create a new 'Delaunay 3D'
# delaunay3D7 = Delaunay3D(registrationName='Delaunay3D7', Input=clip8)
# delaunay3D7.Alpha = 2.0
# delaunay3D7.Tolerance = 0.01

# # create a new 'Group Datasets'
# groupDatasets7 = GroupDatasets(registrationName='GroupDatasets7', Input=[delaunay3D7, appendDatasets7])
# groupDatasets7.BlockNames = ['Delaunay3D7', 'AppendDatasets7']

# # create a new 'Calculator'
# calculator2 = Calculator(registrationName='Calculator2', Input=groupDatasets7)
# calculator2.Function = '-jHat'

# # create a new 'Warp By Vector'
# warpByVector2 = WarpByVector(registrationName='WarpByVector2', Input=calculator2)
# warpByVector2.Vectors = ['POINTS', 'Result']

# # create a new 'XML Unstructured Grid Reader'
# hex_layout_skeleton28vtu = XMLUnstructuredGridReader(registrationName='hex_layout_skeleton28.vtu', FileName=['/Users/caleb/sweeps/build/hex_layout_skeleton28.vtu'])
# hex_layout_skeleton28vtu.TimeArray = 'None'

# # create a new 'Append Datasets'
# appendDatasets8 = AppendDatasets(registrationName='AppendDatasets8', Input=[hex_layout_skeleton34vtu, hex_layout_skeleton28vtu, hex_layout_skeleton25vtu, hex_layout_skeleton26vtu])

# # create a new 'Clip'
# clip9 = Clip(registrationName='Clip9', Input=appendDatasets8)
# clip9.ClipType = 'Plane'
# clip9.HyperTreeGridClipper = 'Plane'
# clip9.Scalars = ['POINTS', '']

# # init the 'Plane' selected for 'ClipType'
# clip9.ClipType.Origin = [0.014100074768066406, -0.36231493949890137, 0.01]
# clip9.ClipType.Normal = [0.0, 0.0, 1.0]

# # init the 'Plane' selected for 'HyperTreeGridClipper'
# clip9.HyperTreeGridClipper.Origin = [0.014100074768066406, -0.36231493949890137, 11.956950187683105]

# # create a new 'Delaunay 3D'
# delaunay3D8 = Delaunay3D(registrationName='Delaunay3D8', Input=clip9)
# delaunay3D8.Alpha = 2.0
# delaunay3D8.Tolerance = 0.01

# # create a new 'Group Datasets'
# groupDatasets8 = GroupDatasets(registrationName='GroupDatasets8', Input=[delaunay3D8, appendDatasets8])
# groupDatasets8.BlockNames = ['Delaunay3D8', 'AppendDatasets8']

# # create a new 'Calculator'
# calculator1 = Calculator(registrationName='Calculator1', Input=groupDatasets8)
# calculator1.Function = '(coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))-jHat'

# # create a new 'Warp By Vector'
# warpByVector1 = WarpByVector(registrationName='WarpByVector1', Input=calculator1)
# warpByVector1.Vectors = ['POINTS', 'Result']
# warpByVector1.ScaleFactor = 0.7

# # create a new 'Delaunay 3D'
# delaunay3D9 = Delaunay3D(registrationName='Delaunay3D9', Input=appendDatasets9)
# delaunay3D9.Alpha = 2.0
# delaunay3D9.AlphaTris = 0

# # create a new 'Calculator'
# calculator9 = Calculator(registrationName='Calculator9', Input=delaunay3D9)
# calculator9.Function = 'jHat'

# # create a new 'Warp By Vector'
# warpByVector9 = WarpByVector(registrationName='WarpByVector9', Input=calculator9)
# warpByVector9.Vectors = ['POINTS', 'Result']
# warpByVector9.ScaleFactor = 1.5

# # create a new 'XML Unstructured Grid Reader'
# outer_surface6vtu = XMLUnstructuredGridReader(registrationName='outer_surface6.vtu', FileName=['/Users/caleb/sweeps/build/outer_surface6.vtu'])
# outer_surface6vtu.TimeArray = 'None'

# # create a new 'Append Datasets'
# appendDatasets11 = AppendDatasets(registrationName='AppendDatasets11', Input=[outer_surface6vtu, hex_layout_skeleton11vtu, hex_layout_skeleton15vtu_1, hex_layout_skeleton23vtu])

# # create a new 'Delaunay 3D'
# delaunay3D11 = Delaunay3D(registrationName='Delaunay3D11', Input=appendDatasets11)
# delaunay3D11.Alpha = 2.0
# delaunay3D11.AlphaTris = 0

# # create a new 'Calculator'
# calculator11 = Calculator(registrationName='Calculator11', Input=delaunay3D11)
# calculator11.Function = 'jHat'

# # create a new 'Warp By Vector'
# warpByVector11 = WarpByVector(registrationName='WarpByVector11', Input=calculator11)
# warpByVector11.Vectors = ['POINTS', 'Result']
# warpByVector11.ScaleFactor = 0.5

# # create a new 'Calculator'
# calculator10 = Calculator(registrationName='Calculator10', Input=delaunay3D10)
# calculator10.Function = '(coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))'

# # create a new 'Warp By Vector'
# warpByVector10 = WarpByVector(registrationName='WarpByVector10', Input=calculator10)
# warpByVector10.Vectors = ['POINTS', 'Result']
# warpByVector10.ScaleFactor = 0.5

# # create a new 'Clip'
# clip4 = Clip(registrationName='Clip4', Input=appendDatasets3)
# clip4.ClipType = 'Plane'
# clip4.HyperTreeGridClipper = 'Plane'
# clip4.Scalars = ['POINTS', '']

# # init the 'Plane' selected for 'ClipType'
# clip4.ClipType.Origin = [0.014100074768066406, -0.36231493949890137, 0.01]
# clip4.ClipType.Normal = [0.0, 0.0, 1.0]

# # init the 'Plane' selected for 'HyperTreeGridClipper'
# clip4.HyperTreeGridClipper.Origin = [0.014100074768066406, -0.36231493949890137, 11.956950187683105]

# # create a new 'Extract Edges'
# extractEdges2 = ExtractEdges(registrationName='ExtractEdges2', Input=clip4)

# # create a new 'Delaunay 3D'
# delaunay3D3 = Delaunay3D(registrationName='Delaunay3D3', Input=extractEdges2)
# delaunay3D3.Alpha = 2.0

# # create a new 'Group Datasets'
# groupDatasets3 = GroupDatasets(registrationName='GroupDatasets3', Input=[appendDatasets3, delaunay3D3])
# groupDatasets3.BlockNames = ['AppendDatasets3', 'Delaunay3D3']

# # create a new 'Calculator'
# calculator5 = Calculator(registrationName='Calculator5', Input=groupDatasets3)
# calculator5.Function = 'jHat'

# # create a new 'Warp By Vector'
# warpByVector5 = WarpByVector(registrationName='WarpByVector5', Input=calculator5)
# warpByVector5.Vectors = ['POINTS', 'Result']

# # create a new 'Clip'
# clip3 = Clip(registrationName='Clip3', Input=appendDatasets2)
# clip3.ClipType = 'Plane'
# clip3.HyperTreeGridClipper = 'Plane'
# clip3.Scalars = ['POINTS', '']

# # init the 'Plane' selected for 'ClipType'
# clip3.ClipType.Origin = [0.014100074768066406, -0.36231493949890137, 0.01]
# clip3.ClipType.Normal = [0.0, 0.0, 1.0]

# # init the 'Plane' selected for 'HyperTreeGridClipper'
# clip3.HyperTreeGridClipper.Origin = [0.014100074768066406, -0.36231493949890137, 11.956950187683105]

# # create a new 'Delaunay 3D'
# delaunay3D2 = Delaunay3D(registrationName='Delaunay3D2', Input=clip3)
# delaunay3D2.Alpha = 2.0

# # create a new 'Group Datasets'
# groupDatasets1 = GroupDatasets(registrationName='GroupDatasets1', Input=[appendDatasets2, delaunay3D2])
# groupDatasets1.BlockNames = ['AppendDatasets2', 'Delaunay3D2']

# # create a new 'Calculator'
# calculator7 = Calculator(registrationName='Calculator7', Input=groupDatasets1)
# calculator7.Function = '(coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))+jHat'

# # create a new 'Warp By Vector'
# warpByVector7 = WarpByVector(registrationName='WarpByVector7', Input=calculator7)
# warpByVector7.Vectors = ['POINTS', 'Result']
# warpByVector7.ScaleFactor = 0.7

# # create a new 'Clip'
# clip5 = Clip(registrationName='Clip5', Input=appendDatasets4)
# clip5.ClipType = 'Plane'
# clip5.HyperTreeGridClipper = 'Plane'
# clip5.Scalars = ['POINTS', '']

# # init the 'Plane' selected for 'ClipType'
# clip5.ClipType.Origin = [0.014100074768066406, -0.36231493949890137, 0.01]
# clip5.ClipType.Normal = [0.0, 0.0, 1.0]

# # init the 'Plane' selected for 'HyperTreeGridClipper'
# clip5.HyperTreeGridClipper.Origin = [0.014100074768066406, -0.36231493949890137, 11.956950187683105]

# # create a new 'Delaunay 3D'
# delaunay3D4 = Delaunay3D(registrationName='Delaunay3D4', Input=clip5)
# delaunay3D4.Alpha = 2.0
# delaunay3D4.Tolerance = 0.01

# # create a new 'Group Datasets'
# groupDatasets4 = GroupDatasets(registrationName='GroupDatasets4', Input=[appendDatasets4, delaunay3D4])
# groupDatasets4.BlockNames = ['AppendDatasets4', 'Delaunay3D4']

# # create a new 'Calculator'
# calculator8 = Calculator(registrationName='Calculator8', Input=groupDatasets4)
# calculator8.Function = '-(coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))+jHat'

# # create a new 'Warp By Vector'
# warpByVector8 = WarpByVector(registrationName='WarpByVector8', Input=calculator8)
# warpByVector8.Vectors = ['POINTS', 'Result']
# warpByVector8.ScaleFactor = 0.7

# # create a new 'Group Datasets'
# groupDatasets6 = GroupDatasets(registrationName='GroupDatasets6', Input=[delaunay3D6, appendDatasets6])
# groupDatasets6.BlockNames = ['Delaunay3D6', 'AppendDatasets6']

# # create a new 'Calculator'
# calculator3 = Calculator(registrationName='Calculator3', Input=groupDatasets6)
# calculator3.Function = '-(coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))-jHat'

# # create a new 'Warp By Vector'
# warpByVector3 = WarpByVector(registrationName='WarpByVector3', Input=calculator3)
# warpByVector3.Vectors = ['POINTS', 'Result']
# warpByVector3.ScaleFactor = 0.7

# # create a new 'Clip'
# clip6 = Clip(registrationName='Clip6', Input=appendDatasets5)
# clip6.ClipType = 'Plane'
# clip6.HyperTreeGridClipper = 'Plane'
# clip6.Scalars = ['POINTS', '']

# # init the 'Plane' selected for 'ClipType'
# clip6.ClipType.Origin = [0.014100074768066406, -0.36231493949890137, 0.01]
# clip6.ClipType.Normal = [0.0, 0.0, 1.0]

# # init the 'Plane' selected for 'HyperTreeGridClipper'
# clip6.HyperTreeGridClipper.Origin = [0.014100074768066406, -0.36231493949890137, 11.956950187683105]

# # create a new 'Delaunay 3D'
# delaunay3D5 = Delaunay3D(registrationName='Delaunay3D5', Input=clip6)
# delaunay3D5.Alpha = 2.0
# delaunay3D5.Tolerance = 0.01

# # create a new 'Group Datasets'
# groupDatasets5 = GroupDatasets(registrationName='GroupDatasets5', Input=[delaunay3D5, appendDatasets5])
# groupDatasets5.BlockNames = ['Delaunay3D5', 'AppendDatasets5']

# # create a new 'Calculator'
# calculator4 = Calculator(registrationName='Calculator4', Input=groupDatasets5)
# calculator4.Function = '-(coordsX*iHat+coordsZ*kHat)/(sqrt(coordsX^2+coordsZ^2))'

# # create a new 'Warp By Vector'
# warpByVector4 = WarpByVector(registrationName='WarpByVector4', Input=calculator4)
# warpByVector4.Vectors = ['POINTS', 'Result']

# # ----------------------------------------------------------------
# # setup the visualization in view 'renderView1'
# # ----------------------------------------------------------------

# # show data from warpByVector1
# warpByVector1Display = Show(warpByVector1, renderView1, 'UnstructuredGridRepresentation')

# # trace defaults for the display properties.
# warpByVector1Display.Representation = 'Surface'
# warpByVector1Display.ColorArrayName = ['POINTS', '']
# warpByVector1Display.SelectTCoordArray = 'None'
# warpByVector1Display.SelectNormalArray = 'None'
# warpByVector1Display.SelectTangentArray = 'None'
# warpByVector1Display.OSPRayScaleArray = 'Result'
# warpByVector1Display.OSPRayScaleFunction = 'Piecewise Function'
# warpByVector1Display.Assembly = 'Hierarchy'
# warpByVector1Display.SelectOrientationVectors = 'Result'
# warpByVector1Display.ScaleFactor = 4.597971343994141
# warpByVector1Display.SelectScaleArray = 'None'
# warpByVector1Display.GlyphType = 'Arrow'
# warpByVector1Display.GlyphTableIndexArray = 'None'
# warpByVector1Display.GaussianRadius = 0.22989856719970703
# warpByVector1Display.SetScaleArray = ['POINTS', 'Result']
# warpByVector1Display.ScaleTransferFunction = 'Piecewise Function'
# warpByVector1Display.OpacityArray = ['POINTS', 'Result']
# warpByVector1Display.OpacityTransferFunction = 'Piecewise Function'
# warpByVector1Display.DataAxesGrid = 'Grid Axes Representation'
# warpByVector1Display.PolarAxes = 'Polar Axes Representation'
# warpByVector1Display.ScalarOpacityUnitDistance = 1.7838278401781058
# warpByVector1Display.OpacityArrayName = ['POINTS', 'Result']
# warpByVector1Display.SelectInputVectors = ['POINTS', 'Result']
# warpByVector1Display.WriteLog = ''

# # init the 'Piecewise Function' selected for 'ScaleTransferFunction'
# warpByVector1Display.ScaleTransferFunction.Points = [-0.7071067811865475, 0.0, 0.5, 0.0, 0.7071067811865475, 1.0, 0.5, 0.0]

# # init the 'Piecewise Function' selected for 'OpacityTransferFunction'
# warpByVector1Display.OpacityTransferFunction.Points = [-0.7071067811865475, 0.0, 0.5, 0.0, 0.7071067811865475, 1.0, 0.5, 0.0]

# # show data from warpByVector2
# warpByVector2Display = Show(warpByVector2, renderView1, 'UnstructuredGridRepresentation')

# # trace defaults for the display properties.
# warpByVector2Display.Representation = 'Surface'
# warpByVector2Display.ColorArrayName = ['POINTS', '']
# warpByVector2Display.SelectTCoordArray = 'None'
# warpByVector2Display.SelectNormalArray = 'None'
# warpByVector2Display.SelectTangentArray = 'None'
# warpByVector2Display.OSPRayScaleArray = 'Result'
# warpByVector2Display.OSPRayScaleFunction = 'Piecewise Function'
# warpByVector2Display.Assembly = 'Hierarchy'
# warpByVector2Display.SelectOrientationVectors = 'Result'
# warpByVector2Display.ScaleFactor = 4.046080017089844
# warpByVector2Display.SelectScaleArray = 'None'
# warpByVector2Display.GlyphType = 'Arrow'
# warpByVector2Display.GlyphTableIndexArray = 'None'
# warpByVector2Display.GaussianRadius = 0.20230400085449218
# warpByVector2Display.SetScaleArray = ['POINTS', 'Result']
# warpByVector2Display.ScaleTransferFunction = 'Piecewise Function'
# warpByVector2Display.OpacityArray = ['POINTS', 'Result']
# warpByVector2Display.OpacityTransferFunction = 'Piecewise Function'
# warpByVector2Display.DataAxesGrid = 'Grid Axes Representation'
# warpByVector2Display.PolarAxes = 'Polar Axes Representation'
# warpByVector2Display.ScalarOpacityUnitDistance = 1.1559300961987278
# warpByVector2Display.OpacityArrayName = ['POINTS', 'Result']
# warpByVector2Display.SelectInputVectors = ['POINTS', 'Result']
# warpByVector2Display.WriteLog = ''

# # init the 'Piecewise Function' selected for 'ScaleTransferFunction'
# warpByVector2Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.402408285460702e+38, 1.0, 0.5, 0.0]

# # init the 'Piecewise Function' selected for 'OpacityTransferFunction'
# warpByVector2Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.402408285460702e+38, 1.0, 0.5, 0.0]

# # show data from warpByVector3
# warpByVector3Display = Show(warpByVector3, renderView1, 'UnstructuredGridRepresentation')

# # trace defaults for the display properties.
# warpByVector3Display.Representation = 'Surface'
# warpByVector3Display.ColorArrayName = ['POINTS', '']
# warpByVector3Display.SelectTCoordArray = 'None'
# warpByVector3Display.SelectNormalArray = 'None'
# warpByVector3Display.SelectTangentArray = 'None'
# warpByVector3Display.OSPRayScaleArray = 'Result'
# warpByVector3Display.OSPRayScaleFunction = 'Piecewise Function'
# warpByVector3Display.Assembly = 'Hierarchy'
# warpByVector3Display.SelectOrientationVectors = 'Result'
# warpByVector3Display.ScaleFactor = 1.9226402282714845
# warpByVector3Display.SelectScaleArray = 'None'
# warpByVector3Display.GlyphType = 'Arrow'
# warpByVector3Display.GlyphTableIndexArray = 'None'
# warpByVector3Display.GaussianRadius = 0.09613201141357422
# warpByVector3Display.SetScaleArray = ['POINTS', 'Result']
# warpByVector3Display.ScaleTransferFunction = 'Piecewise Function'
# warpByVector3Display.OpacityArray = ['POINTS', 'Result']
# warpByVector3Display.OpacityTransferFunction = 'Piecewise Function'
# warpByVector3Display.DataAxesGrid = 'Grid Axes Representation'
# warpByVector3Display.PolarAxes = 'Polar Axes Representation'
# warpByVector3Display.ScalarOpacityUnitDistance = 0.9917697446492034
# warpByVector3Display.OpacityArrayName = ['POINTS', 'Result']
# warpByVector3Display.SelectInputVectors = ['POINTS', 'Result']
# warpByVector3Display.WriteLog = ''

# # init the 'Piecewise Function' selected for 'ScaleTransferFunction'
# warpByVector3Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# # init the 'Piecewise Function' selected for 'OpacityTransferFunction'
# warpByVector3Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# # show data from warpByVector4
# warpByVector4Display = Show(warpByVector4, renderView1, 'UnstructuredGridRepresentation')

# # trace defaults for the display properties.
# warpByVector4Display.Representation = 'Surface'
# warpByVector4Display.ColorArrayName = ['POINTS', '']
# warpByVector4Display.SelectTCoordArray = 'None'
# warpByVector4Display.SelectNormalArray = 'None'
# warpByVector4Display.SelectTangentArray = 'None'
# warpByVector4Display.OSPRayScaleArray = 'Result'
# warpByVector4Display.OSPRayScaleFunction = 'Piecewise Function'
# warpByVector4Display.Assembly = 'Hierarchy'
# warpByVector4Display.SelectOrientationVectors = 'Result'
# warpByVector4Display.ScaleFactor = 2.2215499877929688
# warpByVector4Display.SelectScaleArray = 'None'
# warpByVector4Display.GlyphType = 'Arrow'
# warpByVector4Display.GlyphTableIndexArray = 'None'
# warpByVector4Display.GaussianRadius = 0.11107749938964843
# warpByVector4Display.SetScaleArray = ['POINTS', 'Result']
# warpByVector4Display.ScaleTransferFunction = 'Piecewise Function'
# warpByVector4Display.OpacityArray = ['POINTS', 'Result']
# warpByVector4Display.OpacityTransferFunction = 'Piecewise Function'
# warpByVector4Display.DataAxesGrid = 'Grid Axes Representation'
# warpByVector4Display.PolarAxes = 'Polar Axes Representation'
# warpByVector4Display.ScalarOpacityUnitDistance = 0.829808008030397
# warpByVector4Display.OpacityArrayName = ['POINTS', 'Result']
# warpByVector4Display.SelectInputVectors = ['POINTS', 'Result']
# warpByVector4Display.WriteLog = ''

# # init the 'Piecewise Function' selected for 'ScaleTransferFunction'
# warpByVector4Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]

# # init the 'Piecewise Function' selected for 'OpacityTransferFunction'
# warpByVector4Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]

# # show data from warpByVector5
# warpByVector5Display = Show(warpByVector5, renderView1, 'UnstructuredGridRepresentation')

# # trace defaults for the display properties.
# warpByVector5Display.Representation = 'Surface'
# warpByVector5Display.ColorArrayName = ['POINTS', '']
# warpByVector5Display.SelectTCoordArray = 'None'
# warpByVector5Display.SelectNormalArray = 'None'
# warpByVector5Display.SelectTangentArray = 'None'
# warpByVector5Display.OSPRayScaleArray = 'Result'
# warpByVector5Display.OSPRayScaleFunction = 'Piecewise Function'
# warpByVector5Display.Assembly = 'Hierarchy'
# warpByVector5Display.SelectOrientationVectors = 'Result'
# warpByVector5Display.ScaleFactor = 4.2211198806762695
# warpByVector5Display.SelectScaleArray = 'None'
# warpByVector5Display.GlyphType = 'Arrow'
# warpByVector5Display.GlyphTableIndexArray = 'None'
# warpByVector5Display.GaussianRadius = 0.21105599403381348
# warpByVector5Display.SetScaleArray = ['POINTS', 'Result']
# warpByVector5Display.ScaleTransferFunction = 'Piecewise Function'
# warpByVector5Display.OpacityArray = ['POINTS', 'Result']
# warpByVector5Display.OpacityTransferFunction = 'Piecewise Function'
# warpByVector5Display.DataAxesGrid = 'Grid Axes Representation'
# warpByVector5Display.PolarAxes = 'Polar Axes Representation'
# warpByVector5Display.ScalarOpacityUnitDistance = 1.2056904582287655
# warpByVector5Display.OpacityArrayName = ['POINTS', 'Result']
# warpByVector5Display.SelectInputVectors = ['POINTS', 'Result']
# warpByVector5Display.WriteLog = ''

# # init the 'Piecewise Function' selected for 'ScaleTransferFunction'
# warpByVector5Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# # init the 'Piecewise Function' selected for 'OpacityTransferFunction'
# warpByVector5Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# # show data from warpByVector6
# warpByVector6Display = Show(warpByVector6, renderView1, 'UnstructuredGridRepresentation')

# # trace defaults for the display properties.
# warpByVector6Display.Representation = 'Surface'
# warpByVector6Display.ColorArrayName = ['POINTS', '']
# warpByVector6Display.SelectTCoordArray = 'None'
# warpByVector6Display.SelectNormalArray = 'None'
# warpByVector6Display.SelectTangentArray = 'None'
# warpByVector6Display.OSPRayScaleArray = 'Result'
# warpByVector6Display.OSPRayScaleFunction = 'Piecewise Function'
# warpByVector6Display.Assembly = 'Hierarchy'
# warpByVector6Display.SelectOrientationVectors = 'Result'
# warpByVector6Display.ScaleFactor = 4.788080024719238
# warpByVector6Display.SelectScaleArray = 'None'
# warpByVector6Display.GlyphType = 'Arrow'
# warpByVector6Display.GlyphTableIndexArray = 'None'
# warpByVector6Display.GaussianRadius = 0.23940400123596192
# warpByVector6Display.SetScaleArray = ['POINTS', 'Result']
# warpByVector6Display.ScaleTransferFunction = 'Piecewise Function'
# warpByVector6Display.OpacityArray = ['POINTS', 'Result']
# warpByVector6Display.OpacityTransferFunction = 'Piecewise Function'
# warpByVector6Display.DataAxesGrid = 'Grid Axes Representation'
# warpByVector6Display.PolarAxes = 'Polar Axes Representation'
# warpByVector6Display.ScalarOpacityUnitDistance = 1.2505994185950335
# warpByVector6Display.OpacityArrayName = ['POINTS', 'Result']
# warpByVector6Display.SelectInputVectors = ['POINTS', 'Result']
# warpByVector6Display.WriteLog = ''

# # init the 'Piecewise Function' selected for 'ScaleTransferFunction'
# warpByVector6Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]

# # init the 'Piecewise Function' selected for 'OpacityTransferFunction'
# warpByVector6Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]

# # show data from warpByVector7
# warpByVector7Display = Show(warpByVector7, renderView1, 'UnstructuredGridRepresentation')

# # trace defaults for the display properties.
# warpByVector7Display.Representation = 'Surface'
# warpByVector7Display.ColorArrayName = ['POINTS', '']
# warpByVector7Display.SelectTCoordArray = 'None'
# warpByVector7Display.SelectNormalArray = 'None'
# warpByVector7Display.SelectTangentArray = 'None'
# warpByVector7Display.OSPRayScaleArray = 'Result'
# warpByVector7Display.OSPRayScaleFunction = 'Piecewise Function'
# warpByVector7Display.Assembly = 'Hierarchy'
# warpByVector7Display.SelectOrientationVectors = 'Result'
# warpByVector7Display.ScaleFactor = 4.776149940490723
# warpByVector7Display.SelectScaleArray = 'None'
# warpByVector7Display.GlyphType = 'Arrow'
# warpByVector7Display.GlyphTableIndexArray = 'None'
# warpByVector7Display.GaussianRadius = 0.23880749702453613
# warpByVector7Display.SetScaleArray = ['POINTS', 'Result']
# warpByVector7Display.ScaleTransferFunction = 'Piecewise Function'
# warpByVector7Display.OpacityArray = ['POINTS', 'Result']
# warpByVector7Display.OpacityTransferFunction = 'Piecewise Function'
# warpByVector7Display.DataAxesGrid = 'Grid Axes Representation'
# warpByVector7Display.PolarAxes = 'Polar Axes Representation'
# warpByVector7Display.ScalarOpacityUnitDistance = 1.8637723983424388
# warpByVector7Display.OpacityArrayName = ['POINTS', 'Result']
# warpByVector7Display.SelectInputVectors = ['POINTS', 'Result']
# warpByVector7Display.WriteLog = ''

# # init the 'Piecewise Function' selected for 'ScaleTransferFunction'
# warpByVector7Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# # init the 'Piecewise Function' selected for 'OpacityTransferFunction'
# warpByVector7Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# # show data from warpByVector8
# warpByVector8Display = Show(warpByVector8, renderView1, 'UnstructuredGridRepresentation')

# # trace defaults for the display properties.
# warpByVector8Display.Representation = 'Surface'
# warpByVector8Display.ColorArrayName = ['POINTS', '']
# warpByVector8Display.SelectTCoordArray = 'None'
# warpByVector8Display.SelectNormalArray = 'None'
# warpByVector8Display.SelectTangentArray = 'None'
# warpByVector8Display.OSPRayScaleArray = 'Result'
# warpByVector8Display.OSPRayScaleFunction = 'Piecewise Function'
# warpByVector8Display.Assembly = 'Hierarchy'
# warpByVector8Display.SelectOrientationVectors = 'Result'
# warpByVector8Display.ScaleFactor = 2.4427024841308596
# warpByVector8Display.SelectScaleArray = 'None'
# warpByVector8Display.GlyphType = 'Arrow'
# warpByVector8Display.GlyphTableIndexArray = 'None'
# warpByVector8Display.GaussianRadius = 0.12213512420654297
# warpByVector8Display.SetScaleArray = ['POINTS', 'Result']
# warpByVector8Display.ScaleTransferFunction = 'Piecewise Function'
# warpByVector8Display.OpacityArray = ['POINTS', 'Result']
# warpByVector8Display.OpacityTransferFunction = 'Piecewise Function'
# warpByVector8Display.DataAxesGrid = 'Grid Axes Representation'
# warpByVector8Display.PolarAxes = 'Polar Axes Representation'
# warpByVector8Display.ScalarOpacityUnitDistance = 1.2092921500663771
# warpByVector8Display.OpacityArrayName = ['POINTS', 'Result']
# warpByVector8Display.SelectInputVectors = ['POINTS', 'Result']
# warpByVector8Display.WriteLog = ''

# # init the 'Piecewise Function' selected for 'ScaleTransferFunction'
# warpByVector8Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# # init the 'Piecewise Function' selected for 'OpacityTransferFunction'
# warpByVector8Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# # show data from warpByVector9
# warpByVector9Display = Show(warpByVector9, renderView1, 'UnstructuredGridRepresentation')

# # trace defaults for the display properties.
# warpByVector9Display.Representation = 'Surface'
# warpByVector9Display.ColorArrayName = [None, '']
# warpByVector9Display.SelectTCoordArray = 'None'
# warpByVector9Display.SelectNormalArray = 'None'
# warpByVector9Display.SelectTangentArray = 'None'
# warpByVector9Display.OSPRayScaleArray = 'Result'
# warpByVector9Display.OSPRayScaleFunction = 'Piecewise Function'
# warpByVector9Display.Assembly = ''
# warpByVector9Display.SelectOrientationVectors = 'Result'
# warpByVector9Display.ScaleFactor = 4.394599914550781
# warpByVector9Display.SelectScaleArray = 'None'
# warpByVector9Display.GlyphType = 'Arrow'
# warpByVector9Display.GlyphTableIndexArray = 'None'
# warpByVector9Display.GaussianRadius = 0.21972999572753907
# warpByVector9Display.SetScaleArray = ['POINTS', 'Result']
# warpByVector9Display.ScaleTransferFunction = 'Piecewise Function'
# warpByVector9Display.OpacityArray = ['POINTS', 'Result']
# warpByVector9Display.OpacityTransferFunction = 'Piecewise Function'
# warpByVector9Display.DataAxesGrid = 'Grid Axes Representation'
# warpByVector9Display.PolarAxes = 'Polar Axes Representation'
# warpByVector9Display.ScalarOpacityUnitDistance = 1.1499400392976236
# warpByVector9Display.OpacityArrayName = ['POINTS', 'Result']
# warpByVector9Display.SelectInputVectors = ['POINTS', 'Result']
# warpByVector9Display.WriteLog = ''

# # init the 'Piecewise Function' selected for 'ScaleTransferFunction'
# warpByVector9Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# # init the 'Piecewise Function' selected for 'OpacityTransferFunction'
# warpByVector9Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# # show data from warpByVector10
# warpByVector10Display = Show(warpByVector10, renderView1, 'UnstructuredGridRepresentation')

# # trace defaults for the display properties.
# warpByVector10Display.Representation = 'Surface'
# warpByVector10Display.ColorArrayName = [None, '']
# warpByVector10Display.SelectTCoordArray = 'None'
# warpByVector10Display.SelectNormalArray = 'None'
# warpByVector10Display.SelectTangentArray = 'None'
# warpByVector10Display.OSPRayScaleArray = 'Result'
# warpByVector10Display.OSPRayScaleFunction = 'Piecewise Function'
# warpByVector10Display.Assembly = ''
# warpByVector10Display.SelectOrientationVectors = 'Result'
# warpByVector10Display.ScaleFactor = 4.388019943237305
# warpByVector10Display.SelectScaleArray = 'None'
# warpByVector10Display.GlyphType = 'Arrow'
# warpByVector10Display.GlyphTableIndexArray = 'None'
# warpByVector10Display.GaussianRadius = 0.21940099716186523
# warpByVector10Display.SetScaleArray = ['POINTS', 'Result']
# warpByVector10Display.ScaleTransferFunction = 'Piecewise Function'
# warpByVector10Display.OpacityArray = ['POINTS', 'Result']
# warpByVector10Display.OpacityTransferFunction = 'Piecewise Function'
# warpByVector10Display.DataAxesGrid = 'Grid Axes Representation'
# warpByVector10Display.PolarAxes = 'Polar Axes Representation'
# warpByVector10Display.ScalarOpacityUnitDistance = 1.1726719022924894
# warpByVector10Display.OpacityArrayName = ['POINTS', 'Result']
# warpByVector10Display.SelectInputVectors = ['POINTS', 'Result']
# warpByVector10Display.WriteLog = ''

# # init the 'Piecewise Function' selected for 'ScaleTransferFunction'
# warpByVector10Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]

# # init the 'Piecewise Function' selected for 'OpacityTransferFunction'
# warpByVector10Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]

# # show data from warpByVector11
# warpByVector11Display = Show(warpByVector11, renderView1, 'UnstructuredGridRepresentation')

# # trace defaults for the display properties.
# warpByVector11Display.Representation = 'Surface'
# warpByVector11Display.ColorArrayName = [None, '']
# warpByVector11Display.SelectTCoordArray = 'None'
# warpByVector11Display.SelectNormalArray = 'None'
# warpByVector11Display.SelectTangentArray = 'None'
# warpByVector11Display.OSPRayScaleArray = 'Result'
# warpByVector11Display.OSPRayScaleFunction = 'Piecewise Function'
# warpByVector11Display.Assembly = ''
# warpByVector11Display.SelectOrientationVectors = 'Result'
# warpByVector11Display.ScaleFactor = 4.134900093078613
# warpByVector11Display.SelectScaleArray = 'None'
# warpByVector11Display.GlyphType = 'Arrow'
# warpByVector11Display.GlyphTableIndexArray = 'None'
# warpByVector11Display.GaussianRadius = 0.20674500465393067
# warpByVector11Display.SetScaleArray = ['POINTS', 'Result']
# warpByVector11Display.ScaleTransferFunction = 'Piecewise Function'
# warpByVector11Display.OpacityArray = ['POINTS', 'Result']
# warpByVector11Display.OpacityTransferFunction = 'Piecewise Function'
# warpByVector11Display.DataAxesGrid = 'Grid Axes Representation'
# warpByVector11Display.PolarAxes = 'Polar Axes Representation'
# warpByVector11Display.ScalarOpacityUnitDistance = 1.2065586382914517
# warpByVector11Display.OpacityArrayName = ['POINTS', 'Result']
# warpByVector11Display.SelectInputVectors = ['POINTS', 'Result']
# warpByVector11Display.WriteLog = ''

# # init the 'Piecewise Function' selected for 'ScaleTransferFunction'
# warpByVector11Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.402408285460702e+38, 1.0, 0.5, 0.0]

# # init the 'Piecewise Function' selected for 'OpacityTransferFunction'
# warpByVector11Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.402408285460702e+38, 1.0, 0.5, 0.0]

# # ----------------------------------------------------------------
# # setup animation scene, tracks and keyframes
# # note: the Get..() functions create a new object, if needed
# # ----------------------------------------------------------------

# # get time animation track
# timeAnimationCue1 = GetTimeTrack()

# # initialize the animation scene

# # get the time-keeper
# timeKeeper1 = GetTimeKeeper()

# # initialize the timekeeper

# # initialize the animation track

# # get animation scene
# animationScene1 = GetAnimationScene()

# # initialize the animation scene
# animationScene1.ViewModules = renderView1
# animationScene1.Cues = timeAnimationCue1
# animationScene1.AnimationTime = 6.0
# animationScene1.EndTime = 38.0
# animationScene1.PlayMode = 'Snap To TimeSteps'

# # ----------------------------------------------------------------
# # restore active source
# SetActiveSource(warpByVector9)
# # ----------------------------------------------------------------


# ##--------------------------------------------
# ## You may need to add some code at the end of this python script depending on your usage, eg:
# #
# ## Render all views to see them appears
# # RenderAllViews()
# #
# ## Interact with the view, usefull when running from pvpython
# # Interact()
# #
# ## Save a screenshot of the active view
# # SaveScreenshot("path/to/screenshot.png")
# #
# ## Save a screenshot of a layout (multiple splitted view)
# # SaveScreenshot("path/to/screenshot.png", GetLayout())
# #
# ## Save all "Extractors" from the pipeline browser
# # SaveExtracts()
# #
# ## Save a animation of the current active view
# # SaveAnimation()
# #
# ## Please refer to the documentation of paraview.simple
# ## https://kitware.github.io/paraview-docs/latest/python/paraview.simple.html
# ##--------------------------------------------