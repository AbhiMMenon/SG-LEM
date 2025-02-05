# 2D MRB test case

- Mesh created by `blockMesh`
- To change cell size, modify line 24 of `system/blockMeshDict.m4` and then
    
    m4 -P system/blockMeshDict.m4 > system/blockMeshDict
    blockMesh
    setFields
    decomposePar -force
    ./allRun

- Files `scalarZC.dat` and `scalarZC_flag.dat` are the conditioned mass
  fractions and "holes" generated using the 1D LEM solver
  `../spliceFlameInitZ/` for the mechanism specified in
  `constant/chemistryProperties` using the `zBinSize` and `cBinSize` in
  `constant/RILEMProperties`; to use a different mechanism or bin sizes, these
  two files will have to be regenerated using the 1D solver.

- solver settings in `constant/combustionProperties` and `constant/chemistryProperties` are needed for OpenFOAM but ignored in the solver.
 
- Processor decomposition using the `structured` option leads to flame stability as described in the MRB paper.
- `0.org/` contains initial values, i.e.,  before `setFields` is used.
