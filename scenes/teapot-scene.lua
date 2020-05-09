require "teapot"

NUM_SAMPLES_AREA_LIGHT = 4
NUM_SAMPLES_ANTIALIAS = 32
MAX_RAY_DEPTH = 4


local triangles = OBJToTriangles(readOBJFile("models/teapot.obj"),{generateNormals=true})

local material = newMaterial(newDiffuseBSDF({1,1,1}),nil)

local scene = newScene()

sceneAddMaterial(scene,material)
sceneAddPrimitives(scene,triangles,{materialIndex=1})

sceneAddLight(scene,newEnvironmentLight(readPPM("models/pisa.ppm")))

sceneWrapCornellBox(scene)

-- pprint(scene)

sceneCompileBVH(scene)

-- dumpSceneToFile("tmp/teapot-compiled2.lua",scene)



-- scene = loadDumpFile("tmp/teapot-compiled2.lua")

local cam = newCamera({0,0.5,1.45}, {0,0,0}, 60, {200,200})

local frame = render(cam,scene)

writePPM('output/render.ppm',frame)