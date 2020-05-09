require "raytrace"

NUM_SAMPLES_AREA_LIGHT = 8
NUM_SAMPLES_ANTIALIAS = 256
MAX_RAY_DEPTH = 3


local scene = newScene()

local triangles = OBJToTriangles(readOBJFile("models/happy.obj"),{generateNormals=false})

local material = newMaterial(newGlossyBSDF({0.5,0.18,0.1}),nil)
-- material.bsdf.specular=0.9
-- material.bsdf.diffuse=0.0001
-- material.bsdf.shininess=1000

sceneAddMaterial(scene,material)
sceneAddPrimitives(scene,triangles,{materialIndex=1})

local env = readPPM("models/pisa.ppm")

sceneAddLight(scene,newEnvironmentLight(tintTexture(env,{0.9,0.9,0.9})))

-- sceneAddEnvironmentMap(scene, env, {blur=10,spectrum={0.1,0.1,0.1}})

sceneNormalizePrimitives(scene)

local cam = newCamera({0,0.5,1.05}, {0,0,0}, 60, {640,640})

local frame = render(cam,scene)

writePPM('output/render.ppm',frame)