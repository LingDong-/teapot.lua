require "raytrace"

NUM_SAMPLES_AREA_LIGHT = 8
NUM_SAMPLES_ANTIALIAS = 256
MAX_RAY_DEPTH = 4

local newJadeBSDF = function(spectrum, refractiveIndex)
  return {type="jadeBSDF",
          emission={0,0,0},
          spectrum=spectrum,
          isDelta=false,
          transmitRatio=0.5,
          translucency=0.1,
          clarity=0.8,
          radius=0.1}
end
evaluateBSDF.jadeBSDF = function(bsdf,wo,wi)
  return scale(bsdf.spectrum , (1 / math.pi));
end
computeBSDFIncident.jadeBSDF = function(bsdf,wo,side)
  local wi = uniformHemisphereSampler()
  local dir = 1

  if (math.random() < bsdf.transmitRatio and side > 0) 
  or (side < 0) then
    dir = -1
    wi = scale(lerp(wi,wo,bsdf.clarity),dir)
  else
    wi = lerp(wi,reflect(wo),bsdf.clarity)
  end

  local tMax = INFINITY
  if side > 0 and dir == -1 then
    tMax = bsdf.radius
  end
  return {
    wi = wi,
    pdf = math.abs(wi[3]),
    tMax = tMax,
    spectrum = bsdf.spectrum,
  }
end

local floorSize = 100
local floor = {newTriangle({-floorSize,0,-floorSize},{floorSize,0,-floorSize},{floorSize,0,floorSize}),
               newTriangle({-floorSize,0,-floorSize},{floorSize,0,floorSize},{-floorSize,0,floorSize})}

local scene = newScene()

local obj = readOBJFile("models/bowl.obj")

local triangles = OBJToTriangles(obj,{generateNormals=true})

-- local jadeBSDF = newSSSBSDF({0.96,0.96,0.64})
-- jadeBSDF.radius = 0.08

local jadeBSDF = newJadeBSDF({0.96,0.92,0.64})
local floorMaterial = newMaterial(newDiffuseBSDF({1,1,1}))
floorMaterial.invisible = true

sceneAddMaterial(scene, newMaterial(newGlossyBSDF({0.3,0.28,0.22})))  -- base
sceneAddMaterial(scene, newMaterial(jadeBSDF))  -- bowl
sceneAddMaterial(scene, floorMaterial)

sceneAddPrimitives(scene, slice(triangles, obj.groups[1].iBegin, obj.groups[1].iEnd), {materialIndex=1,transform=newTransformationMatrix("rotateY",-math.pi*0.1)})
sceneAddPrimitives(scene, slice(triangles, obj.groups[2].iBegin, obj.groups[2].iEnd), {materialIndex=2,transform=newTransformationMatrix("rotateY",-math.pi*0.1)})

-- sceneWrapCornellBox(scene,{scale=0.7,leftSpectrum={1,1,1},rightSpectrum={1,1,1}})

sceneNormalizePrimitives(scene)

sceneAddPrimitives(scene, floor, {materialIndex=3})

local env = readPPM("models/pisa.ppm")

sceneAddLight(scene,newEnvironmentLight(tintTexture(env,{0.45,0.45,0.45})))

sceneAddLight(scene,newDirectionalLight({3,3,3},{0.2,-0.2,0.5}))


-- sceneAddEnvironmentMap(scene, env, {spectrum={0.8,0.8,0.8}, invisible=true})

local cam = newCamera({0,1.1,1.2}, {-math.pi*0.15,0,0}, 60, {640,640})

local frame = render(cam,scene)

writePPM('output/render.ppm',frame)