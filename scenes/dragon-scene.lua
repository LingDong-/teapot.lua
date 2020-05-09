require "raytrace"

NUM_SAMPLES_AREA_LIGHT = 8
NUM_SAMPLES_ANTIALIAS = 256
MAX_RAY_DEPTH = 6

MODE = 4

sceneWrapInvisibleFloor = function(scene, options)
  if options == nil then options = {} end
  options.scale = options.scale or 0.9
  options.spectrum = options.spectrum or {1,1,1}
  if options.normalize == nil then options.normalize = true end

  if options.normalize then
    sceneNormalizePrimitives(scene, options)
  end

  local v1={-100, 0,  100}
  local v2={-100, 0, -100}
  local v5={ 100, 0,  100}
  local v6={ 100, 0, -100}  

  local material = newMaterial(newDiffuseBSDF(options.spectrum))
  material.invisible = true

  sceneAddMaterial(scene, material)
  sceneAddPrimitives(scene,{newTriangle(v2,v6,v5),newTriangle(v2,v5,v1)},{materialIndex=#scene.materials})

end

if MODE == 1 then

  local bvh = newBVH(OBJToTriangles(readOBJFile("models/dragon-s1.obj"),{generateNormals=false}),4)

  dumpToFile("models/dragon-simp-bvh.lua",bvh)

  -- local bvh = loadDumpFile("models/dragon-s1-bvh.lua")

  local cam = newCamera({-13,9,13}, {0,0,0}, 60, {640,480})
  setCameraLookAt(cam,boundingBoxCentroid(bvh.root.boundingBox))

  local frame = renderNDotL(cam, {-0.8,1,0.5}, bvh)
  writePPM('output/render.ppm',frame)




elseif MODE == 2 then

  local triangles = OBJToTriangles(readOBJFile("models/dragon-s2nb.obj"),{generateNormals=false})

  local cam = newCamera({-10,9,10}, {0,0,0}, 60, {120,120})


  local material = newMaterial(newSSSBSDF({0.96,1,0.85}),nil)
  material.bsdf.radius = 1

  local scene = newScene()

  sceneAddMaterial(scene,material)
  sceneAddPrimitives(scene,triangles,{materialIndex=1})

  -- sceneAddLight(scene,newPointLight({20,20,17},{0,1,-4}))
  -- sceneAddLight(scene,newDirectionalLight({5,5,4},{-0.5,0,-0.3}))


  -- sceneAddLight(scene,newPointLight({1,1,1},{0,5,5}))

  -- sceneAddLight(scene,newAmbientLight({0.9,0.9,0.9}))

  sceneAddLight(scene,newEnvironmentLight(tintTexture(readPPM("models/pisa.ppm"),{0.5,0.5,0.5})))

  sceneCompileBVH(scene)

  --setCameraLookAt(cam,boundingBoxCentroid(scene.bvh.root.boundingBox))
  setCameraLookAt(cam,{-1.5,4,0})

  local frame = render(cam,scene)

  writePPM('output/render.ppm',frame)


elseif MODE == 3 then

  local triangles = OBJToTriangles(readOBJFile("models/dragon-s2nb.obj"),{generateNormals=false})

  local material = newMaterial(newMirrorBSDF({1,0.75,0.2}),nil)

  local scene = newScene()

  sceneAddMaterial(scene,material)
  sceneAddPrimitives(scene,triangles,{materialIndex=1})


  local sphere = newSphere({0,3,4},3)

  sceneAddMaterial(scene,newMaterial(newGlassBSDF({1,1,1}),nil))
  sceneAddPrimitives(scene,{sphere},{materialIndex=2})

  sceneWrapCornellBox(scene,{scale=0.9,leftSpectrum={0.2,0.6,0.9},rightSpectrum={1,0,0}})

  sceneCompileBVH(scene)


  local cam = newCamera({0,0.5,1.45}, {0,0,0}, 60, {100,100})

  local frame = render(cam,scene)

  writePPM('output/render.ppm',frame)


elseif MODE == 4 then
  local triangles = OBJToTriangles(readOBJFile("models/dragon-s2nb.obj"),{generateNormals=false})

  local cam = newCamera({-10,9,10}, {0,0,0}, 60, {120,120})

  local material = newMaterial(newSSSBSDF({0.96,1,0.85}),nil)
  material.bsdf.radius = 1

  local scene = newScene()

  sceneAddMaterial(scene,material)
  sceneAddPrimitives(scene,triangles,{materialIndex=1})

  sceneWrapInvisibleFloor(scene,{normalize=false})

  sceneAddLight(scene,newDirectionalLight({20,19.9,19.5},{-0.8,-0.1,0.5}))
  sceneAddLight(scene,newAmbientLight({0.1,0.1,0.1}))

  sceneCompileBVH(scene)

  setCameraLookAt(cam,{-1.5,4,0})

  local frame = render(cam,scene)

  writePPM('output/render.ppm',frame)

end