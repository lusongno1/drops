set(HOME surfactant)

libs(ifacetransp)
libs(sfpde)
libs(femP3)

target_link_libraries(surfactant-ifacetransp levelset-levelsetmapper)

exec_ser(surfactant geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-unknowns geom-topo num-fe misc-problem levelset-levelset levelset-marking_strategy levelset-adaptriang misc-scopetimer misc-progressaccu misc-utils out-output num-discretize misc-params num-interfacePatch levelset-fastmarch surfactant-ifacetransp num-fe out-ensightOut levelset-surfacetension out-vtkOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-dynamicload misc-funcmap misc-scopetimer num-gradient_recovery levelset-levelsetmapper -fopenmp -lphg -lmpi surfactant-sfpde surfactant-femP3)

exec_ser(surfacestokes geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-unknowns geom-topo num-fe misc-problem levelset-levelset levelset-marking_strategy levelset-adaptriang misc-scopetimer misc-progressaccu misc-utils out-output num-discretize misc-params num-interfacePatch levelset-fastmarch surfactant-ifacetransp num-fe out-ensightOut levelset-surfacetension out-vtkOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-dynamicload misc-funcmap misc-scopetimer -fopenmp -lphg -lmpi surfactant-sfpde femP3)

exec_ser(test_tmp geom-boundary geom-builder geom-simplex geom-multigrid geom-deformation num-unknowns geom-topo num-fe misc-problem levelset-levelset levelset-marking_strategy levelset-adaptriang misc-scopetimer misc-progressaccu misc-utils out-output num-discretize misc-params num-interfacePatch levelset-fastmarch surfactant-ifacetransp num-fe out-ensightOut levelset-surfacetension out-vtkOut geom-principallattice geom-reftetracut geom-subtriangulation num-quadrature misc-dynamicload misc-funcmap misc-scopetimer -fopenmp -lphg -lmpi surfactant-sfpde femP3)

if(NOT MPI)
add_dependencies(surfactant misc-scalarFunctions misc-vectorFunctions)
add_dependencies(surfacestokes misc-scalarFunctions misc-vectorFunctions levelset-twophaseCoeff)
endif(NOT MPI)

add_my_custom_targets(surfactant surfacestokes)
