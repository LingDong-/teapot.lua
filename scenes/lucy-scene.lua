require "raytrace"

NUM_SAMPLES_AREA_LIGHT = 8
NUM_SAMPLES_ANTIALIAS = 256
MAX_RAY_DEPTH = 3


local scene = newScene()

local triangles = OBJToTriangles(readOBJFile("models/lucy.obj"),{generateNormals=false})

sceneAddPrimitives(scene,triangles,{transform=newTransformationMatrix("rotateY",math.pi*2.85)})

local env = readPPM("models/doge.ppm")

sceneAddLight(scene,newEnvironmentLight(tintTexture(env,{1,1,1})))

sceneAddEnvironmentMap(scene, env, {blur=8,spectrum={1,1,1}})

sceneNormalizePrimitives(scene)

local cam = newCamera({-0.8,0.3,0.6}, {math.pi*0.05,math.pi*1.7,0}, 60, {720,720})

local frame = render(cam,scene)

writePPM('output/render.ppm',frame)