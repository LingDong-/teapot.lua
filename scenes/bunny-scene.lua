require "raytrace"

NUM_SAMPLES_AREA_LIGHT = 4
NUM_SAMPLES_ANTIALIAS = 64
MAX_RAY_DEPTH = 4

MODE = 2

if MODE == 1 then

  local bvh = newBVH(OBJToTriangles(readOBJFile("models/bunny.obj")),4)

  local cam = newCamera({0.17,0.2,0.17}, {0,0,0}, 60, {640,480})
  setCameraLookAt(cam,boundingBoxCentroid(bvh.root.boundingBox))

  local frame = renderNDotL(cam, {1,0.8,0.5}, bvh)
  writePPM('output/render.ppm',frame)

elseif MODE == 2 then


  local triangles = OBJToTriangles(readOBJFile("models/bunny.obj"),{generateNormals=true})

  local material = newMaterial(newGlossyBSDF({0.8,0.5,1}),nil)

  local scene = newScene()

  sceneAddMaterial(scene,material)
  sceneAddPrimitives(scene,triangles,{materialIndex=1})

  sceneWrapCornellBox(scene,{scale=0.6})

  sceneCompileBVH(scene)

  local cam = newCamera({0,0.5,1.45}, {0,0,0}, 60, {720,720})

  local frame = render(cam,scene)

  writePPM('output/render.ppm',frame)


end