# Code to Generate a Local Coordinate System on an STL file from a Vertebrae using the VTK Library
# Andres Gonzalez 2024

import vtk
import numpy as np
from sklearn.decomposition import PCA
import os


def read_stl(file_path):
    """ Function to read the STL file """
    reader = vtk.vtkSTLReader()
    reader.SetFileName(file_path)
    reader.Update()
    return reader.GetOutput()

def create_bounding_box(polydata):

    """ Function to create a bounding box around the polydata """

    bounds = polydata.GetBounds()
    cube = vtk.vtkCubeSource()
    cube.SetBounds(bounds)
    cube_mapper = vtk.vtkPolyDataMapper()
    cube_mapper.SetInputConnection(cube.GetOutputPort())
    
    cube_actor = vtk.vtkActor()
    cube_actor.SetMapper(cube_mapper)
    cube_actor.GetProperty().SetRepresentationToWireframe()
    cube_actor.GetProperty().SetColor(1.0, 0, 0)  # Red color

    centerOfMass = vtk.vtkCenterOfMass()
    centerOfMass.SetInputConnection(cube.GetOutputPort())
    centerOfMass.Update()
    localCoordinateCenter = centerOfMass.GetCenter()

    sphereSource = vtk.vtkSphereSource()
    sphereSource.SetCenter(localCoordinateCenter)
    sphereSource.SetRadius(2)  # Adjust the radius as needed
        
    localCoordinateCenter_mapper = vtk.vtkPolyDataMapper()
    localCoordinateCenter_mapper.SetInputConnection(sphereSource.GetOutputPort())
        
    localCoordinateCenter_actor = vtk.vtkActor()
    localCoordinateCenter_actor.SetMapper(localCoordinateCenter_mapper)
    localCoordinateCenter_actor.GetProperty().SetColor(1,1,1)

    return cube_actor, localCoordinateCenter_actor, localCoordinateCenter


def calculate_centroid(bounds):
    """  Function to calculate the center of the bounding box (centroid) using the bounds  """
    x_center = (bounds[0] + bounds[1]) / 2
    y_center = (bounds[2] + bounds[3]) / 2
    z_center = (bounds[4] + bounds[5]) / 2
    return x_center, y_center, z_center

def find_extreme_points(polydata):
    """  Function to calculate the extreme point located on the surface of the STL (Polydata)  """
    bounds = polydata.GetBounds()
    most_posterior = bounds[0]  # Minimum X
    most_anterior = bounds[1]  # Maximum X
    most_inferior = bounds[4]  # Minimum Z
    most_superior = bounds[5]  # Maximum Z
    most_lateral_left = bounds[2]  # Minimum Y
    most_lateral_right = bounds[3]  # Maximum Y

    # Initialize points to the bounds as a starting reference
    posterior_point = [most_posterior, 0, 0]
    anterior_point = [most_anterior, 0, 0]
    inferior_point = [0, 0, most_inferior]
    superior_point = [0, 0, most_superior]
    lateral_left_point = [0, most_lateral_left, 0]
    lateral_right_point = [0, most_lateral_right, 0]

    points = polydata.GetPoints()
    
    for i in range(points.GetNumberOfPoints()):
        point = points.GetPoint(i)
        if point[0] == most_posterior:
            posterior_point = point
        if point[0] == most_anterior:
            anterior_point = point
        if point[2] == most_inferior:
            inferior_point = point
        if point[2] == most_superior:
            superior_point = point
        if point[1] == most_lateral_left:
            lateral_left_point = point
        if point[1] == most_lateral_right:
            lateral_right_point = point

    return posterior_point, anterior_point, inferior_point, superior_point, lateral_left_point, lateral_right_point

def visualize_extreme_points(renderer, points):
    """ Function to create an Actor and add it to Renderer for each calculated extreme point """

    colors = [
        (1, 0, 0),  # Red for the most posterior
        (1, 1, 0),  # Yellow for the most anterior
        (0, 1, 0),  # Green for the most inferior
        (0, 0, 1),  # Blue for the most superior
        (1, 0, 1),  # Magenta for the most lateral left
        (0, 1, 1)   # Cyan for the most lateral right
    ]
    for i, point in enumerate(points):
        sphereSource = vtk.vtkSphereSource()
        sphereSource.SetCenter(point)
        sphereSource.SetRadius(1)  # Adjust the radius as needed
        
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphereSource.GetOutputPort())
        
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(colors[i])
        
        renderer.AddActor(actor)

# 
def segment_anterior(polydata, max_point, min_point):
    """  Function to segment the anterior part of the STL based on two opposite laying extreme points (MIN/MAX points)  
        IMPORTANT: Adjust according to your coordinate system.
    """
    # Calculate a threshold for segmentation based on two opposite laying extreme points
    segmentation_threshold = max_point[0] - (max_point[0] - min_point[0]) * 0.25

    # Create a mask for points that are in the anterior segment (Quarter)
    ids = vtk.vtkIdTypeArray()
    ids.SetNumberOfComponents(1)
    points = polydata.GetPoints()
    for i in range(polydata.GetNumberOfPoints()):
        # IMPORTANT: Adapt the for x,y or z in order to properly segment the BODY of the vertebrae STL 
        _, y, _ = points.GetPoint(i) 
        if y > segmentation_threshold:
            ids.InsertNextValue(i)

    selection_node = vtk.vtkSelectionNode()
    selection_node.SetFieldType(vtk.vtkSelectionNode.POINT)
    selection_node.SetContentType(vtk.vtkSelectionNode.INDICES)
    selection_node.SetSelectionList(ids)

    selection = vtk.vtkSelection()
    selection.AddNode(selection_node)

    extract_selection = vtk.vtkExtractSelection()
    extract_selection.SetInputData(0, polydata)
    extract_selection.SetInputData(1, selection)
    extract_selection.Update()

    # Convert the extracted selection to polydata
    geometry_filter = vtk.vtkGeometryFilter()
    geometry_filter.SetInputConnection(extract_selection.GetOutputPort())
    geometry_filter.Update()

    return geometry_filter.GetOutput()

def extract_surface_points(polydata):
    """Extract points from the surface of the polydata.
    IMPORTANT: Assuming Z is up-direction. Adjust according to your coordinate system.
    """
    
    bounds = polydata.GetBounds()
    top_z = bounds[5]  # Upper bound of Z-coordinate
    points = []
    for i in range(polydata.GetNumberOfPoints()):
        x, y, z = polydata.GetPoint(i)
        # If the full Polydata is to be considered remove the if below to add all points
        if z > top_z - 10:  # Adjust the threshold according to your needs
            points.append([x, y, z])
    return np.array(points)

def PCA_normal(points):
    """Calculate normal using PCA of the points."""
    pca = PCA(n_components=3)
    pca.fit(points)
    normal = pca.components_[-1]  # The component with the smallest variance
    centroid = np.mean(points, axis=0)
    return normal, centroid

def create_line_actor(start_point, end_point):
    """Create a VTK actor for a line given start and end points."""
    # Create line
    line_source = vtk.vtkLineSource()
    line_source.SetPoint1(start_point)
    line_source.SetPoint2(end_point)
    
    # Mapper
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(line_source.GetOutputPort())
    
    # Actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(1, 1, 0)  # Yellow color
    actor.GetProperty().SetLineWidth(2)  # Adjust line width as needed
    
    return actor

def calculate_perpendicular_axes(normal):
    """Calculate two orthogonal axes perpendicular to the given normal vector.
    this is used to generate the Local Coordinate System
    """
    # Generate an arbitrary vector. Change this if it might be parallel to your normal.
    arbitrary_vector = np.array([1, 0, 0]) if not np.allclose(normal, [1, 0, 0]) else np.array([0, 1, 0])
    
    # X-axis: Cross product of normal and arbitrary vector
    x_axis = np.cross(normal, arbitrary_vector)
    x_axis /= np.linalg.norm(x_axis)  # Normalize
    
    # Y-axis: Cross product of normal and x_axis
    y_axis = np.cross(normal, x_axis)
    y_axis /= np.linalg.norm(y_axis)  # Normalize
    
    return x_axis, y_axis

def set_camera_for_view(renderer, view_name, zoom_factor=2):
    """Adjust the camera settings based on the desired view and zoom out."""
    camera = renderer.GetActiveCamera()
    position = None
    focal_point = [0, 0, 0]
    view_up = [0, 0, 0]
    if view_name == "top":
        position = [0, 0, 100 * zoom_factor]
        view_up = [0, 1, 0]
    elif view_name == "side":
        position = [100 * zoom_factor, 0, 0]
        view_up = [0, 0, 1]
    elif view_name == "back":
        position = [0, -100 * zoom_factor, 0]
        view_up = [0, 0, 1]

    camera.SetPosition(position)
    camera.SetFocalPoint(focal_point)
    camera.SetViewUp(view_up)
    renderer.ResetCameraClippingRange()


def capture_screenshot(render_window, filename):
    """Capture a screenshot of the current VTK render window."""
    window_to_image_filter = vtk.vtkWindowToImageFilter()
    window_to_image_filter.SetInput(render_window)
    window_to_image_filter.SetScale(1)  # Adjust scale to increase resolution
    window_to_image_filter.ReadFrontBufferOff()  # Use the back buffer
    window_to_image_filter.Update()

    writer = vtk.vtkPNGWriter()
    writer.SetFileName(filename)
    writer.SetInputConnection(window_to_image_filter.GetOutputPort())
    writer.Write()

# Main function to execute the visualization
if __name__ == "__main__":
    file_path = 'Bvretebc2.stl'  # Update this with your actual file path
    polydata = read_stl(file_path)
    bounding_box_actor = create_bounding_box(polydata)[0]
    extreme_points = find_extreme_points(polydata)
    anterior_point = extreme_points[3]  # Assuming this is the most anterior point
    posterior_point = extreme_points[2]  # Assuming this is the most posterior point

    # Now segment the anterior part using these points
    anterior_segment = segment_anterior(polydata, anterior_point, posterior_point)
    anterior_bounding_box_actor, localCoordinatePoint_actor, localCoordinatePoint  = create_bounding_box(anterior_segment)
    anterior_bounds = anterior_segment.GetBounds()
    anterior_extreme_points = find_extreme_points(anterior_segment)
    
    # Extract the surface points  
    surfpoints = extract_surface_points(anterior_segment)
    normal = PCA_normal(surfpoints)[0]

    # Calculate end point for the normal vector line (adjust length as needed)
    line_length = 20  # Adjust length as needed

    normal_actor = create_line_actor(localCoordinatePoint, localCoordinatePoint + normal * line_length)
    x_axis, y_axis = calculate_perpendicular_axes(normal)
    x_actor = create_line_actor(localCoordinatePoint, localCoordinatePoint + x_axis  * line_length)
    y_actor = create_line_actor(localCoordinatePoint, localCoordinatePoint + y_axis  * line_length)

    # Setup visualization
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polydata)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetOpacity(0.3)
    
    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    
    renderer.AddActor(actor)
    renderer.AddActor(bounding_box_actor)
    renderer.AddActor(anterior_bounding_box_actor)
    renderer.AddActor(localCoordinatePoint_actor)
    renderer.AddActor(normal_actor)
    renderer.AddActor(x_actor)
    renderer.AddActor(y_actor)
    visualize_extreme_points(renderer, extreme_points)  # Add extreme points to the renderer
    visualize_extreme_points(renderer, anterior_extreme_points)  # Add extreme points to the renderer
    

    renderer.SetBackground(.3, .2, .1)
    
    # Generate Pictures of the TOP, SIDE and BACK Views of the STL with everything generated
    views = ["top", "side", "back"]
    for view in views:
        set_camera_for_view(renderer, view)
        renderWindow.Render()  # Ensure the scene is rendered before capturing
        screenshot_filename = f"{os.path.splitext(file_path)[0]}_{view}.png"
        capture_screenshot(renderWindow, screenshot_filename)
        print(f"Screenshot saved: {screenshot_filename}")

    renderWindow.Render()
    renderWindowInteractor.Start()

