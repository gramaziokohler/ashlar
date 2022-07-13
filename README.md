# Ashlar:  Mesh Properties for Stone Datasets

This repository is provided as a supplement to the paper _Multimodal Robotic Construction with On-Site Materials_.  
For closed triangular meshes (i.e. representing stones), we provide methods for computing meta-faces using iterative variational shape approximation, form and mass properties, and the ashlarness metric.

## Dependencies
This project uses the geometry processing library [libigl](http://libigl.github.io/libigl/), 
and the only dependencies are STL, Eigen, libigl and the dependencies
of the `igl::opengl::glfw::Viewer` (OpenGL, glad and GLFW).
The CMake build system will automatically download libigl and its dependencies using
[CMake FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html),
thus requiring no setup on your part.

## Compile

Compile this project using the standard cmake routine:

    mkdir build && cd build
    cmake ..
    make

This should find and build the dependencies and create the ``ashlar`` binary.

## Run

For computing the properties of a single stone, execute with the additional argument of the 3D triangle mesh file location.
For example, to use the provided sample mesh `923.obj`, execute:

```
./ashlar  ../mesh/923.obj
```

For processing many mesh files at once (i.e. to assess the volume, dimensions, and form properties of a dataset of stone models),
simply launch with a wildcard variable for all files in a specified directory:

```
./ashlar  ../mesh/*.obj
```

In both cases, the computed properties should print to the console, and a viewer will launch displaying the latest stone, with
colored meta-faces, edge midpoints, and bounding box (the mesh is reoriented using PCA).

![ashlar](images/preview.png)

## Output
At each run, the software should output a log file `logger.csv` in the main project directory which includes salient properties of each stone.
This can be useful, for example, for plotting or assessing the distribution of properties in a dataset of stones using external tools.
Note that this log file is overwritten at each run as implemented.

## File formats
Input files should be in one of the following formats:  obj, off, stl, wrl, ply, mesh
Models should be clean, closed triangle meshes of stone-like objects:  disjoint pieces, many duplicate vertices, etc. will likely cause issues.

## Local libigl
To use a local copy of libigl rather than downloading the repository via FetchContent, you can use
the CMake cache variable `FETCHCONTENT_SOURCE_DIR_LIBIGL` when configuring your CMake project for
the first time:
```
cmake -DFETCHCONTENT_SOURCE_DIR_LIBIGL=<path-to-libigl> ..
```
When changing this value, do not forget to clear your `CMakeCache.txt`, or to update the cache variable
via `cmake-gui` or `ccmake`.

## References
The Variational Shape Approximation approach for finding meta-faces is documented in

```
@inproceedings{
}
```

Mass properties are computed using the method described in

```
@inproceedings{
}
```

For the visualization of meta-faces, with non-overlappling colors, we use the four color implementation of [okaydemir](https://github.com/okaydemir/4-color-theorem), modified for a simplified eigen interface.

