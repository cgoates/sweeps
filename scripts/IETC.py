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

