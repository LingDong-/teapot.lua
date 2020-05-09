require "raytrace"

NUM_SAMPLES_AREA_LIGHT = 8
NUM_SAMPLES_ANTIALIAS = 256
MAX_RAY_DEPTH = 4


local scene = newScene()

local obj = readOBJFile("models/silverpot.obj")

local triangles = OBJToTriangles(obj,{generateNormals=true})

-- local material = newMaterial(newDiffuseBSDF({1,1,1}))
local material = newMaterial(newGlossyBSDF({0.5,0.5,0.5}))
material.bsdf.mirrorEffect=0.75
material.bsdf.specular=0.6

sceneAddMaterial(scene, material)

sceneAddPrimitives(scene, triangles, {materialIndex=1,transform=newTransformationMatrix("rotateY",math.pi*0.2)})---math.pi*1.2

sceneNormalizePrimitives(scene)

local env = readPPM("models/ennis.ppm")

-- sceneAddLight(scene,newAmbientLight({1,1,1}))

sceneAddLight(scene,newEnvironmentLight(tintTexture(env,{3.6,3.6,3.6})))

sceneAddEnvironmentMap(scene, env, {spectrum={4,4,4}, invisible=true})

local cam = newCamera({0,1.1,1.2}, {-math.pi*0.15,0,0}, 60, {720,720})

local frame = render(cam,scene)

writePPM('output/render.ppm',frame)