require "raytrace"

NUM_SAMPLES_AREA_LIGHT = 8
NUM_SAMPLES_ANTIALIAS = 64
MAX_RAY_DEPTH = 4


local scene = newScene()

local sphere1 = newSphere({0,3,0},2)
local sphere2 = newSphere({5,3,0},2)

sceneAddMaterial(scene,newMaterial(newGlassBSDF({1,1,1}),nil))
sceneAddMaterial(scene,newMaterial(newGlossyBSDF({0.8,0.5,1}),nil))

sceneAddPrimitives(scene,{sphere1},{materialIndex=1})
sceneAddPrimitives(scene,{sphere2},{materialIndex=2})

sceneAddLight(scene,newEnvironmentLight(tintTexture(readPPM("models/pisa.ppm"),{0.5,0.5,0.5})))


-- sceneWrapCornellBox(scene,{scale=0.85})

-- local cam = newCamera({0,0.5,1.45}, {0,0,0}, 60, {100,100})

local cam = newCamera({0,3,10}, {0,0,0}, 60, {100,100})

local frame = render(cam,scene)

writePPM('output/render.ppm',frame)