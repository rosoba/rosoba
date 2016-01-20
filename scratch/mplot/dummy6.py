from mayavi import mlab as p3d
from mayavi.sources.vtk_data_source import VTKDataSource
from tvtk.api import tvtk
from vtk.util import colors
from vtk.util.colors import peacock, tomato


e1 = p3d.get_engine()
# First start by reading a cow model. We also generate surface normals for
# prettier rendering.

cow = tvtk.BYUReader()
cow.file_name = r"/usr/share/VTKData/Data/Viewpoint/cow.g"

#cow = tvtk.STLReader()
# cow.file_name=r"C:\Downloads\junk_adjusted11.stl"
# cow=tvtk.OBJReader()
# cow.file_name=r"D:\SciPylot2\diamlrge.obj"

cowNormals = tvtk.PolyDataNormals()
cowNormals.set_input(cow.get_output())

# We clip with an implicit function. Here we use a plane positioned near
# the center of the cow model and oriented at an arbitrary angle.
plane = tvtk.Plane()

# for diamlarge.obj
##plane.origin=[0.0, 0.0, 0.0]
##plane.normal=[0.0, 0.0, 1.0]
##
# clip_value=20.0

# For cow.g

plane.origin = [0.25, 0.0, 0.0]
plane.normal = [1.0, 1.0, 0.0]

clip_value = .50

# vtkClipPolyData requires an implicit function to define what it is to
# clip with. Any implicit function, including complex boolean combinations
# can be used. Notice that we can specify the value of the implicit function
# with the SetValue method.
clipper = tvtk.ClipPolyData()
clipper.set_input(cowNormals.get_output())
clipper.clip_function = plane
clipper.generate_clip_scalars = 1
clipper.generate_clipped_output = 1
clipper.value = clip_value
#clipM.apper = tvtk.PolyDataMapper()
# clipMapper.set_input_connection(clipper.output_port)
# clipMapper.scalar_visibility=0
backProp = tvtk.Property()
backProp.diffuse_color = tomato
#clipActor = tvtk.Actor()
# clipActor.mapper=clipMapper
# clipActor.property.color=peacock
# clipActor.backface_property=backProp
##bmp1 = tvtk.JPEGReader()
##bmp1.file_name=r"C:\Program Files\Cut3D Trial\Textures\1_Materials\SeamlessPine_rotated.jpg"
# my_texture=tvtk.Texture()
# my_texture.interpolate=1
# my_texture.set_input(0,bmp1.get_output())
my_scene = e1.new_scene()
# my_scene.scene.add_actor(clipActor)

src1 = VTKDataSource(data=clipper.get_output())
e1.add_source(src1)
asurf = p3d.pipeline.surface(
    src1, opacity=1.0, name="main_cow", color=colors.white)
asurf.actor.actor.backface_property = backProp

p3d.show_pipeline()
# asurf.actor.actor.texture=my_texture

# Here we are cutting the cow. Cutting creates lines where the cut
# function intersects the model. (Clipping removes a portion of the
# model but the dimension of the data does not change.)
#
# The reason we are cutting is to generate a closed polygon at the
# boundary of the clipping process. The cutter generates line
# segments, the stripper then puts them together into polylines. We
# then pull a trick and define polygons using the closed line
# segements that the stripper created.
cutEdges = tvtk.Cutter()
cutEdges.set_input(cowNormals.get_output())
cutEdges.cut_function = plane
cutEdges.generate_cut_scalars = 1
cutEdges.set_value(0, clip_value)
cutStrips = tvtk.Stripper()
cutStrips.set_input_connection(cutEdges.output_port)
cutStrips.update()
cutPoly = tvtk.PolyData()
cutPoly.points = cutStrips.get_output().points
cutPoly.polys = cutStrips.get_output().lines

# Triangle filter is robust enough to ignore the duplicate point at
# the beginning and end of the polygons and triangulate them.
cutTriangles = tvtk.TriangleFilter()
cutTriangles.input = cutPoly
cutMapper = tvtk.PolyDataMapper()
cutMapper.input = cutPoly
cutMapper.set_input_connection(cutTriangles.output_port)
cutActor = tvtk.Actor()
cutActor.mapper = cutMapper
cutActor.property.color = peacock

src2 = VTKDataSource(data=cutTriangles.get_output())
e1.add_source(src2)
asurf = p3d.pipeline.surface(
    src2, opacity=.50, name="cutting_plane", color=colors.peacock)


# my_scene.scene.add_actor(cutActor)
##
# The clipped part of the cow is rendered wireframe.
##restMapper = tvtk.PolyDataMapper()
# restMapper.set_input_connection(clipper.clipped_output_port)
# restMapper.scalar_visibility=0
##restActor = tvtk.Actor()
# restActor.mapper=restMapper
# restActor.property.representation='wireframe'
# my_scene.scene.add_actor(restActor)

src3 = VTKDataSource(data=clipper._get_clipped_output())
e1.add_source(src3)
asurf = p3d.pipeline.surface(
    src3, opacity=.50, name="clipped_output", color=colors.peacock)

##
# Create graphics stuff
##ren = vtk.vtkRenderer()
##renWin = vtk.vtkRenderWindow()
# renWin.AddRenderer(ren)
##iren = vtk.vtkRenderWindowInteractor()
# iren.SetRenderWindow(renWin)
##
# Add the actors to the renderer, set the background and size
# ren.AddActor(clipActor)
# ren.AddActor(cutActor)
# ren.AddActor(restActor)
##ren.SetBackground(1, 1, 1)
# ren.ResetCamera()
# ren.GetActiveCamera().Azimuth(30)
# ren.GetActiveCamera().Elevation(30)
# ren.GetActiveCamera().Dolly(1.5)
# ren.ResetCameraClippingRange()
##
##renWin.SetSize(300, 300)
# iren.Initialize()
##
# Lets you move the cut plane back and forth by invoking the function
# Cut with the appropriate plane value (essentially a distance from
# the original plane).  This is not used in this code but should give
# you an idea of how to define a function to do this.


def Cut(v):
    clipper.value = v
    cutEdges.set_value(0, v)
    cutStrips.update()
    cutPoly.points = cutStrips.get_output().points
    cutPoly.polys = cutStrips.get_output().lines
    cutMapper.update()
    my_scene.scene.render()
##
# renWin.Render()
# iren.Start()
