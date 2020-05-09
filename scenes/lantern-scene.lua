require "raytrace"

NUM_SAMPLES_AREA_LIGHT = 8
NUM_SAMPLES_ANTIALIAS = 256
MAX_RAY_DEPTH = 3


local scene = newScene()

local obj = readOBJFile("models/lantern.obj")

local triangles = OBJToTriangles(obj,{generateNormals=true})

sceneAddMaterial(scene, newMaterial(newEmissionBSDF({0.5,0.45,0.4})))  -- bottom paper
sceneAddMaterial(scene, newMaterial(newGlossyBSDF({0.3,0.17,0.1})))  -- wood
sceneAddMaterial(scene, newMaterial(newEmissionBSDF({0.4,0.35,0.3})))  -- top paper
sceneAddMaterial(scene, newMaterial(newMirrorBSDF({1,0.6,0.1})))  -- gold foil
sceneAddMaterial(scene, newMaterial(newDiffuseBSDF({0.4,0.08,0.08})))  -- tassel

sceneAddPrimitives(scene, slice(triangles, obj.groups[1].iBegin, obj.groups[1].iEnd), {materialIndex=1})
sceneAddPrimitives(scene, slice(triangles, obj.groups[2].iBegin, obj.groups[2].iEnd), {materialIndex=2})
sceneAddPrimitives(scene, slice(triangles, obj.groups[3].iBegin, obj.groups[3].iEnd), {materialIndex=3})
sceneAddPrimitives(scene, slice(triangles, obj.groups[4].iBegin, obj.groups[4].iEnd), {materialIndex=4})
sceneAddPrimitives(scene, slice(triangles, obj.groups[5].iBegin, obj.groups[5].iEnd), {materialIndex=5})

sceneNormalizePrimitives(scene)

sceneAddPrimitives(scene,scene.primitives,{transform=newTransformationMatrix("translate",{0,0,-1})})
sceneAddPrimitives(scene,scene.primitives,{transform=newTransformationMatrix("translate",{0,0,-2})})


local env = readPPM("models/pisa.ppm")

sceneAddLight(scene,newEnvironmentLight(tintTexture(env,{0.5,0.45,0.35})))

sceneAddLight(scene,newDirectionalLight({5,5,5},{0.5,0,0.1}))

for i=0,3 do
  sceneAddLight(scene,newPointLight({18,17,10},{0,0.2,0.05-i}))
  sceneAddLight(scene,newPointLight({4,3,0},{0,0.05,0.05-i}))
end

sceneAddEnvironmentMap(scene, env, {spectrum={0.8,0.8,0.8}, invisible=true})


local cam = newCamera({-0.15,0,0.4}, {math.pi*0.16,-math.pi*0.08,0}, 60, {640,640})

local frame = render(cam,scene)

writePPM('output/render.ppm',frame)