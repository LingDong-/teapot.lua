require "raytrace"

MODE = 2

if MODE == 1 then
  local bvh = newBVH(OBJToTriangles(readOBJFile("models/box.obj"),{generateNormals=false}),4)

  local cam = newCamera({1.5,2,2}, {0,0,0}, 60, {640,480})
  setCameraLookAt(cam,{0,0,0})

  local frame = renderNDotL(cam, {1,0.8,0.5}, bvh)
  writePPM('output/render.ppm',frame)

elseif MODE == 2 then

  local triangles = OBJToTriangles(readOBJFile("models/box.obj"),{generateNormals=false})

  local cam = newCamera({2,2,2}, {0,0,0}, 60, {320,240})
  setCameraLookAt(cam,{0,0,0})

  local light = newDirectionalLight({1,1,1},{-1,-1,-1})

  local material = newMaterial(newDiffuseBSDF({1,0,1}),nil)

  local scene = newScene()

  sceneAddMaterial(scene,material)
  sceneAddPrimitives(scene,triangles,{materialIndex=1})

  sceneAddLight(scene,light)

  local frame = render(cam,scene)

  writePPM('output/render.ppm',frame)
end