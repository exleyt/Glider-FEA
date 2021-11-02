model = createpde();
importGeometry(model, 'template-model\files\CORNER.STL');
generateMesh(model);
pdeplot3D(model)


