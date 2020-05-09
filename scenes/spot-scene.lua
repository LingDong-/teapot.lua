require "raytrace"

NUM_SAMPLES_AREA_LIGHT = 8
NUM_SAMPLES_ANTIALIAS = 32
MAX_RAY_DEPTH = 4

MODE = 2

if MODE == 1 then

  local texture = readPPM("models/spot_texture.ppm")

  local bvh = newBVH(OBJToTriangles(readOBJFile("models/spot.obj")),4)

  local cam = newCamera({1.4,1.4,-1.4}, {-math.pi/6,3*math.pi/4,0}, 60, {640,480})

  local frame = renderNDotL(cam, {1,1,-0.5}, bvh, texture)
  writePPM('output/render.ppm',frame)

elseif MODE == 2 then

  local texture = readPPM("models/spot_texture.ppm")

  local triangles = OBJToTriangles(readOBJFile("models/spot.obj"),{generateNormals=true})

  local material = newMaterial(newDiffuseBSDF({1,1,1}),texture)

  local scene = newScene()

  sceneAddMaterial(scene,material)
  sceneAddPrimitives(scene,triangles,{materialIndex=1,transform=newTransformationMatrix("rotateY",math.pi*3/4)})

  sceneWrapCornellBox(scene,{scale=0.7})

  sceneCompileBVH(scene)

  local cam = newCamera({0,0.5,1.45}, {0,0,0}, 60, {200,200})

  local frame = render(cam,scene)

  writePPM('output/render.ppm',frame)

end