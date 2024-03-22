# Vertebrae local coordinate system generator
Python Code to Generate a Local Coordinate System on an STL file from a Vertebrae using the VTK Library

Here's a step-by-step guide for your Python code, formatted in Markdown. This guide explains how the code works for generating a local coordinate system on an STL file from a vertebrae using the VTK library.

## Main Execution Flow

1. **Read STL File**: The STL file is loaded into polydata.
2. **Bounding Box Creation**: A bounding box is created around the entire STL polydata.
3. **Extreme Points Identification**: Identifies the most anterior, posterior, superior, inferior, lateral left, and right points. These points are then color coded and visualized using the respective function. Here it is important to note that the points represent the closest points from the STL to the Bounds from the bounding box. 
4. **Anterior Segmentation**: Segments the anterior part of the STL for further analysis. This is done using opposing extreme points calculated in the aforementioned step. The distance between these points is calculated and then a threshold (in this case 25%) is taken to essentially cut and segment the STL to the desired zone. 
5. **Surface Point Extraction**: Extracts points from the surface of the anterior segment.
6. **Normal Calculation Using PCA**: Determines the normal to the surface best fitting the extracted points.
7. **Line Actors Creation**: Creates line actors to represent the local coordinate system's axes. It first uses the Normal calculated in step 6 and then the other two axis are generated using the calculate_perpendicular_axis function to be visualized in the same way. 
8. **Visualization**: Visualizes the original STL, the bounding boxes, extreme points, and the local coordinate system in a VTK render window. 
9. Pictures are also taken to see the Top, Side and Back view of the generated render window.

---

# Guide to Generating a Local Coordinate System on a Vertebrae STL File using VTK

This Python script utilizes the Visualization Toolkit (VTK) and other libraries to generate a local coordinate system for STL files, particularly those representing vertebrae. The process involves reading the STL file, creating a bounding box, finding extreme points, segmenting the anterior part of the STL, extracting surface points, and ultimately determining the local coordinate system.

## Dependencies
- `vtk`: For 3D computer graphics, image processing, and visualization.
- `numpy`: For numerical operations.
- `sklearn.decomposition`: Specifically, PCA (Principal Component Analysis) for dimensionality reduction tasks.

## Functions Overview

### `read_stl(file_path)`
- **Purpose**: Reads an STL file and returns its content as polydata.
- **Parameters**: `file_path` - Path to the STL file.
- **Returns**: Polydata of the STL file.

### `create_bounding_box(polydata)`
- **Purpose**: Creates a wireframe bounding box around the given polydata.
- **Parameters**: `polydata` - The polydata for which to create a bounding box.
- **Returns**: Actors for both the bounding box and its center point, along with the local coordinate center.

### `calculate_centroid(bounds)`
- **Purpose**: Calculates the centroid of the bounding box.
- **Parameters**: `bounds` - Bounds of the polydata.
- **Returns**: Coordinates of the centroid.

### `find_extreme_points(polydata)`
- **Purpose**: Identifies extreme points (most anterior, posterior, superior, inferior, lateral left, and lateral right) on the polydata surface.
- **Parameters**: `polydata` - The polydata to analyze.
- **Returns**: Coordinates of each extreme point.

### `visualize_extreme_points(renderer, points)`
- **Purpose**: Visualizes extreme points on the renderer.
- **Parameters**: `renderer` - The VTK renderer, `points` - Extreme points to visualize.
- **Effects**: Adds visual representation of extreme points to the renderer.

### `segment_anterior(polydata, max_point, min_point)`
- **Purpose**: Segments the anterior part of the polydata based on extreme points.
- **Parameters**: `polydata`, `max_point`, `min_point` - Points to define the segmentation.
- **Returns**: Polydata of the segmented anterior part.

### `extract_surface_points(polydata)`
- **Purpose**: Extracts points from the surface of the polydata.
- **Parameters**: `polydata` - Polydata from which to extract points.
- **Returns**: Array of surface points.

### `PCA_normal(points)`
- **Purpose**: Calculates the normal to the plane best fitting the given points using PCA.
- **Parameters**: `points` - Points to analyze.
- **Returns**: Normal vector and centroid of the points.

### `create_line_actor(start_point, end_point)`
- **Purpose**: Creates a line actor between two points.
- **Parameters**: `start_point`, `end_point` - Coordinates of the line's endpoints.
- **Returns**: A line actor for visualization.

### `calculate_perpendicular_axes(normal)`
- **Purpose**: Calculates two orthogonal axes perpendicular to a given normal vector.
- **Parameters**: `normal` - The normal vector.
- **Returns**: Two orthogonal axes as numpy arrays.

### `set_camera_for_view(renderer, view_name, zoom_factor)`
- **Purpose**: Adjusts camera settings for a specific view.
- **Parameters**: `renderer`, `view_name`, `zoom_factor` - Determines the view and zoom level.
- **Effects**: Modifies camera position and orientation.

### `capture_screenshot(render_window, filename)`
- **Purpose**: Captures a screenshot of the current render window.
- **Parameters**: `render_window`, `filename` - The render window and filename for the screenshot.
- **Effects**: Saves a screenshot to the specified file.

