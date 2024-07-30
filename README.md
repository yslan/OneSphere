# OneSphere

All-hex Nek mesh generator for a single pebble inside a box (or cylinder) container.

### Features
- sphere in box
- sphere in tube
- boundary layers, wall refinement
- dump vtk files for checking
- (WIP) high order curves in Nek's sphere and cylinder format.

### Usage:
See drivers 
- `driver_box.m`: box container
- `driver_cyl.m`: cylinder container


### (WIP) Code Structure

- Basic Workflow
  1. Make the surface mesh for sphere following cubed sphere.
  2. track the surface, extrude a fer sph boundary layers
  3. extrude the surface to a box to get "sph in box
  4. (cyl) extrude to make a fab of cylinder, extrude again for inlet/outlet.
  4. (box) append box meshes.

- Structure: fbox 
  6 faces of the cubed sphere, nbx * nbx Quad elements each.
  We track the element id for this "face front" for extrusions

  ```
  Hfront: nbx x nbx x 6, stores the elements id that is on the frontier of the fbox.
  ```

- Structure: Hex
  ```
  X:     npts x 3, coordinates of points
  Hexes: nel x 8, storing the vertices id of Hex elemenets, not watertight
  Hbc:   A template for Nek's CBC, but it's not watertight

  Hcurves:  6 (type + 5 curve) x 12 (nface or nedge) x nel
    type 1 = 'm', (edge) midpt for H20 mesh, 2-4: coordinates
    type 2 = 's', (face) sphere where 2-4: center, 5: radius
    type 3 = 'c', (not used, use midpt instead) (edge) cylinder, centered at origin
                  2: radius (positive=inside, negative=outside)
                  Nek only reads edge 1-8, won't work for face 5-6
    type 4 = 'e', (not supprted) (face) some parametric surface
  ```

- Mesh Actions
  ```
  fbox_extrude_general:  extrude fbot from geom1 to geom2
      geom type 1 = box
      geom type 2 = sphere
      geom type 3 = cylinder (z-axis)
  ```

- Vtk files

- Curves

- Mesh Checkes

- Boundary Conditiones

- IO for Nek
  ```
  X, Hexes, CBC
  ```

### TODOs
- clean up
- logfile
- some documentationn, and data structure
- some figures
- Nek5000 / NekRS case files and quick verification for fluid results
- verbose, write control
- cht for solid of sphere
- co2 (this will take tons of conversion)
- octave
- vectorize the computation for E>>100k
- Extend Nek5000's "cyl" curve to face 5-6


