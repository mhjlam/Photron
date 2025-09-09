# Architecture Refactoring

There are two main parts of the application: the Simulator and the Renderer. There is some overlap between these.

## Renderer

Renderer visualizes a scene, handles user input, the graphical user interface.

- Scene contains a list of Models, a Camera, and Light sources.
- Model represents a 3D object in the scene. Each Model has at least one a Mesh and a Shader to render it.
- Camera represents the viewpoint from which the scene is rendered. It handles user input to move around the scene and change the view.
- Light represents a light source in the scene. It can be a point light, directional light, or spotlight.
- Shader represents an OpenGL shader program that is used to render objects in the scene. It handles loading, compiling, and linking vertex and fragment shaders.
- Mesh represents a geometric shape that can be rendered as polygonal 3D object. It contains vertex data (positions, normals, texture coordinates) and index data for drawing the model using OpenGL.
- Volume represents a voxelized geometric shape. It contains a list of Voxels that make up the volume. It might be useful to have a Mesh associated with a Volume for rendering purposes, or have these classes inherit from a common base class.
- Voxel represents a single voxel in a Volume. It is an axis-aligned bounding box (AABB) that contains a position, size, and color.

## Simulator

Simulator controls and maintains the Monte Carlo simulation, which propogates photons through a (multi-layered) medium.

- Medium represents the (multi-layered) medium that photons are traced through. It contains one or a number of Layers and keeps track of overall radiance metrics, like total absorption, total reflection, total emittance, etc. Medium is associated with a Model for rendering purposes.
- Layer represents a single layer in the Medium and has certain optical properties (Material). Layer is associated with a Mesh for both simulation and rendering purposes. The simulator needs to be able to determine where a photon is inside the layer geometry (for surface interactions). The layer geometry is voxelized into a Volume for simulation purposes (interaction inside the layer and recording of radiance).
- Material represents the optical properties of a Layer. It contains properties like refractive index, absorption coefficient, scattering coefficient, anisotropy factor, etc.
- Cell represents a single cell in the Volume. It is an axis-aligned bounding box (AABB) that contains a position, size, and optical properties like absorption and emittance. A Cell is associated with a Voxel for rendering purposes.
- Photon represents a single photon in the simulation. It contains properties like position, direction, weight, and state (alive or dead), methods for propagating the photon through the medium, handling interactions with layers, and recording radiance in the medium. Each Photon consists of a Path with Nodes that represent a history of interactions like scattering events, reflections, absorptions, etc.
