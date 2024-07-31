# OneSphere

All-hex Nek mesh generator for a single pebble inside a box (or cylinder) container.

### Features
- sphere in box
- sphere in tube
- boundary layers, wall refinement
- dump vtk files for checking
- High order curves in Nek's sphere (ok) and cylinder (2nd order via midpt) format.

### Usage:
See drivers 
- `driver_box.m`: box container
- `driver_cyl.m`: cylinder container

| box | cyl |
|:---:|:---:|
| [](demo_figs/mesh_box.png) | [](demo_figs/mesh_cyl.png) |


### (WIP) Code Structure

- Basic Workflow
  1. Make the surface mesh for sphere following cubed sphere.
  2. track the surface, extrude a few sph boundary layers
  3. extrude the surface to a box to get "sph in box"
  4. Fill the remaining space
     - (cyl) extrude to make a fab of cylinder, extrude again for inlet/outlet.
     - (box) append box meshes.

- Structure: fbox 
  6 faces of the cubed sphere, `nbx * nbx` Quad elements each.
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
  ```
  Hexes_<name>       dump hex elements
  Hexes_<name>_f     dump 6 faces of hex, so we can check the sideset at each stage
  Hexes_<name>_cbc   same as the _f one, but the values are from (final) CBC
  ```

- Curves    
  Sphere curves is easy but Cylinder is tricky.
  Nek's cyl curve format requires face 5-6 parallel to bottom and top of the cyl. 
  However, after extrusion, can be on face 6. 
  Until we extend the Nek's cyl format, the workaround is to use midpt (2nd order, hex20). 
  Here, we regroup those sidesets and manally reproject the midpt onto the cyl.

- Mesh Checkes

- Boundary Conditiones     
  At every mesh stage, we also print the setdiff of the existing sideset id.     
  Then, we manually group those bcid in Hbc and assign the boudaryID via the function `set_cbc`    

- IO for Nek     
  Basic mesh should contain `(X,Hexes,CBC)`     
  We map the `boundaryID` to a dump BC set `W  ,W01,W02, ...` that can be recovered in `usdat2`.
  Hcurve is optional with `ifcurve=1` for curved sides.     

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


