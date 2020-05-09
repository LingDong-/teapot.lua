INFINITY = 1/0
EPSILON = 0.0001
NUM_SAMPLES_AREA_LIGHT = 8
NUM_SAMPLES_ANTIALIAS = 256
MAX_RAY_DEPTH = 4

------- FUNCTIONAL -------

newPolymorphicFunction = function(index)
  index = index or 1
  local f = {}
  setmetatable(f,{
    __index = function(_) return f.default end,
    __call = function(...) local a = {...} return f[a[1+index].type](select(2,...)) end
  })
  return f
end

------- VECTORS -------

dot = function(u, v)
  local n = 0
  for i = 1,#v do
    n = n + u[i] * v[i]
  end
  return n
end

norm2 = function(v)
  return dot(v, v)
end

norm = function(v)
  return math.sqrt(norm2(v))
end

normalize = function(v)
  return scale(v, 1 / (norm(v)))
end

scale = function(v, x)
  local u = {}
  for i = 1,#v do
    u[i] = v[i]*x
  end
  return u
end

sub = function(u, v)
  return add(u, scale(v, -1))
end

add = function(u, v)
  local w = {}
  for i = 1,#v do
    w[i] = u[i]+v[i]
  end
  return w
end

mul = function(u,v)
  local w = {}
  for i = 1,#v do
    w[i] = u[i]*v[i]
  end
  return w
end

cross = function(u, v)
  return {u[2] * v[3] - u[3] * v[2], u[3] * v[1] - u[1] * v[3], u[1] * v[2] - u[2] * v[1]}
end

lerp = function(u,v,t)
  return add(scale(u,1-t),scale(v,t))
end

tabulate = function(n, f)
  local v = {}
  local g
  if (type(f)~='function') then
    g = function() return f end
  else
    g = f
  end
  for i = 1,n do
    v[i]=g(i)
  end
  return v
end

map = function(v, f)
  local w = {}
  for i = 1,#v do
    w[i] = f(v[i])
  end
  return w
end

mergeSort = function(arr, i0, i1, f)
  local merge = function(arr, i0, m, i1, f)
    
    local n1 = m - i0 + 1
    local n2 = i1 - m
    
    local lhs = {}
    local rhs = {}
    for i = 1 , n1 do
      lhs[i] = arr[i0+i-1]
    end
    for i = 1 , n2 do
      rhs[i] = arr[m+i]
    end
    
    local i = 1
    local j = 1
    local k = i0
    while i <= n1 and j <= n2 do
      local ai = lhs[i]
      local aj = rhs[j]
      if f(ai) <= f(aj) then
        arr[k] = ai
        i = i + 1
      else
        arr[k] = aj
        j = j + 1
      end
      k = k + 1
    end
    
    while i <= n1 do
      arr[k] = lhs[i]
      i = i + 1
      k = k + 1
    end
    while j <= n2 do
      arr[k] = rhs[j]
      j = j + 1
      k = k + 1
    end
  end

  if i0 < i1 then
    local m = math.floor((i1 + i0) / 2)
    mergeSort(arr, i0, m, f)
    mergeSort(arr, m + 1, i1, f)
    merge(arr, i0, m, i1, f)
  end
end

slice = function(arr, i0, i1)
  if i0 == nil then
    i0 = 1
  end
  if i1 == nil then
    i1 = #arr
  end
  local brr = {}
  for i = i0, math.min(#arr, i1) do
    brr[i-i0+1] = arr[i]
  end
  return brr
end

clone = function(a)
  if type(a) ~= "table" then
    return a
  end
  local b = {}
  for k,v in pairs(a) do
    b[k] = clone(v)
  end
  return b
end

concat = function(arr,brr)
  local crr = {}
  for i = 1, #arr do
    crr[i] = arr[i]
  end
  for i = 1, #brr do
    crr[#arr+i] = brr[i]
  end
  return crr
end

contains = function(arr,a)
  for i = 1, #arr do
    if arr[i] == a then
      return true
    end
  end
  return false
end

------- MATRICES -------

newTransformationMatrix = function(typ, a)
  if typ == "rotateX" then
    return {{1, 0, 0, 0}, {0, math.cos(a), -(math.sin(a)), 0}, {0, math.sin(a), math.cos(a), 0}, {0, 0, 0, 1}}
  elseif typ == "rotateY" then
    return {{math.cos(a), 0, math.sin(a), 0}, {0, 1, 0, 0}, {-(math.sin(a)), 0, math.cos(a), 0}, {0, 0, 0, 1}}
  elseif typ == "rotateZ" then
    return {{math.cos(a), -(math.sin(a)), 0, 0}, {math.sin(a), math.cos(a), 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}}
  elseif typ == "translate" then
    return {{1, 0, 0, a[1]}, {0, 1, 0, a[2]}, {0, 0, 1, a[3]}, {0, 0, 0, 1}}
  elseif typ == "scale" then
    if type(a) ~= "table" then
      a = {a, a, a}
    end
    return {{a[1], 0, 0, 0}, {0, a[2], 0, 0}, {0, 0, a[3], 0}, {0, 0, 0, 1}}
  elseif typ == "id" then
    return {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}}
  else
    error "Unsupported matrix type"
  end
  
end

multiplyMatrices = function(m1, m2)
  local mat = {}
  for i = 1 , #m1 do
    mat[i] = {}
    for j = 1, #(m2[1]) do
      local sum = 0
      for k = 1, #(m1[1]) do
        sum = sum + m1[i][k] * m2[k][j]
      end
      mat[i][j] = sum
    end
  end
  return mat
end

transpose = function(m)
  local mat = {}
  for i = 1 , #(m[1]) do
    mat[i] = {}
    for j = 1 , #m do
      mat[i][j] = m[j][i]
    end
  end
  return mat
end

transform = newPolymorphicFunction(2)

transform.default = function(m, v)
  local vh = transpose({ concat(v,{1})})
  local wh = multiplyMatrices(m, vh)
  local w = {}
  for i = 1 , (#wh - 1) do
    w[i] = wh[i][1] / wh[#wh][1]
  end
  return w
end

transform.triangle = function(m,p)
  local transformVertexInPlace = function(m,v)
    local u = transform(m,v)
    if v.normal then
      local nx = transform(m, add(v,v.normal))
      local n = normalize(sub(nx,u))
      v.normal = n
    end
    v[1] = u[1]
    v[2] = u[2]
    v[3] = u[3]
  end
  local q = clone(p)
  transformVertexInPlace(m,q[1])
  transformVertexInPlace(m,q[2])
  transformVertexInPlace(m,q[3])
  return q
end

transform.sphere = function(m,p)
  local q = clone(p)
  local o = transform(m,q.center)
  local r = transform(m, add(q.center,{0,0,q.radius}))
  q.radius = norm(sub(r,o))
  q.center = o
  return q
end

hasNaN = function(m)
  if type(m) == 'number' then
    if (m==m) then
      return false
    else
      return true
    end
  end
  for i = 1, #m do
    if hasNaN(m[i]) then
      return true
    end
  end
  return false
end

------- RAYS -------

newRay = function(o, d)
  return {type="ray", origin=o, direction=normalize(d), 
          tMin=0, tMax=INFINITY, depth=0}
end

rayIntersect = {}

rayIntersect.triangle = function(r, tri)
  local p1 = tri[1]
  local p2 = tri[2]
  local p3 = tri[3]
  local d = r.direction
  local o = r.origin
  local e1 = sub(p2, p1)
  local e2 = sub(p3, p1)
  local s = sub(o, p1)
  local _d = scale(d, -1)
  local normal = cross(e1, e2)
  local denom = dot(normal, _d)
  if denom == 0 then
    return nil
  end
  local u = dot(cross(s, e2), _d) / denom
  local v = dot(cross(e1, s), _d) / denom
  local t = dot(normal, s) / denom
  local side = 1

  if t > (r.tMax) or t < (r.tMin) then
    return nil
  end

  if u < 0 or v < 0 or (1 - u - v) < 0 then
    return nil
  end
  normal = normalize(normal)
  if p1.normal and p2.normal and p3.normal then
    local n1 = scale(p1.normal, 1-u-v)
    local n2 = scale(p2.normal, u)
    local n3 = scale(p3.normal, v)
    local n = add(add(n1,n2),n3)
    if not (dot(d,n) > 0 and dot(d,normal) < 0) then
      normal = n
    end
  end
  if (dot(d,normal) > 0) then
    normal =scale(normal,-1)
    side = -1
  end
  
  local uv
  if p1.uv and p2.uv and p3.uv then
    local uv1 = scale(p1.uv, 1-u-v)
    local uv2 = scale(p2.uv, u)
    local uv3 = scale(p3.uv, v)
    uv = add(add(uv1,uv2),uv3)
  end

  r.tMax = t

  return {uv=uv, t=t, normal=normal, primitive=tri, side=side}
end

rayIntersect.sphere = function(r,sph)
  local d = sub(r.origin,sph.center)

  local a = dot(r.direction,r.direction)
  local b = 2*dot(r.direction, d)
  local c = dot(d,d) - sph.radius*sph.radius
  local discr = b*b-4*a*c
  if discr < 0 then
    return nil
  end
  local t1 = (-b-math.sqrt(discr))/(2*a)
  local t2 = (-b+math.sqrt(discr))/(2*a)
  if t1 > r.tMax or t2 < r.tMin then
    return nil
  end
  local t
  local side
  if t1 < r.tMin then
    t = t2
    side = -1
  else
    t = t1
    side = 1
  end
  local x = add(r.origin, scale(r.direction, t1))
  local n = normalize(sub(x, sph.center))

  r.tMax = t2

  return {t=t, normal=n, primitive=sph, side=side}
end

rayIntersect.boundingBox = function(r, bbox)
  local tx1 = ((bbox.min[1]) - (r.origin[1])) / (r.direction[1])
  local tx2 = ((bbox.max[1]) - (r.origin[1])) / (r.direction[1])
  local ty1 = ((bbox.min[2]) - (r.origin[2])) / (r.direction[2])
  local ty2 = ((bbox.max[2]) - (r.origin[2])) / (r.direction[2])
  local tz1 = ((bbox.min[3]) - (r.origin[3])) / (r.direction[3])
  local tz2 = ((bbox.max[3]) - (r.origin[3])) / (r.direction[3])
  
  local t1 = math.max(math.min(tx1, tx2), math.min(ty1, ty2), math.min(tz1, tz2))
  local t2 = math.min(math.max(tx1, tx2), math.max(ty1, ty2), math.max(tz1, tz2))
  
  if t2 - t1 < 0 then
    return nil
  end
  if t1 > (r.tMax) or t2 < (r.tMin) then
    return nil
  end
  return {t1=t1,t2=t2}
  
end


------- TRIANGLES -------

newTriangle = function(v1, v2, v3)
  return {type="triangle", [1]=v1, [2]=v2, [3]=v3}
end

triangleNormal = function(tri)
  local p1 = tri[1]
  local p2 = tri[2]
  local p3 = tri[3]
  local e1 = sub(p2, p1)
  local e2 = sub(p3, p1)
  return normalize(cross(e1, e2))
end

triangleArea = function(tri)
  local a = norm(sub(tri[1],tri[2]))
  local b = norm(sub(tri[2],tri[3]))
  local c = norm(sub(tri[3],tri[1]))
  local s = (a+b+c)/2
  return math.sqrt(s*(s-a)*(s-b)*(s-c))
end

------- SPHERES -------

newSphere = function(o,r)
  return {type="sphere", center=o, radius=math.abs(r)}
end

------- CAMERAS -------

setCameraTranslation = function(camera, translation)
  camera.translationMatrix = newTransformationMatrix("translate", translation)
end

setCameraRotation = function(camera, euler)
  camera.rotationMatrix = multiplyMatrices(multiplyMatrices(
    newTransformationMatrix("rotateZ", euler[3]), 
    newTransformationMatrix("rotateY", euler[2])), 
    newTransformationMatrix("rotateX", euler[1]))
end

setCameraLookAt = function(camera,target)
  -- https://github.com/mrdoob/three.js/blob/dev/src/math/Matrix4.js
  local eye = transform(camera.translationMatrix, {0, 0, 0})
  local up = {0,1,0}
  local z = sub(eye,target)
  if norm2(z) == 0 then
    z[3] = 1
  end
  z = normalize(z)
  local x = cross(up,z)
  if norm2(x) == 0 then
    if math.abs(up[3]) == 1 then
      z[1] = z[1] + 0.0001
    else
      z[3] = z[3] + 0.0001
    end
    z = normalize(z)
    x = cross(up,z)
  end
  normalize(x)
  local y = cross(z,x)
  camera.rotationMatrix = {{x[1],y[1],z[1],0},{x[2],y[2],z[2],0},{x[3],y[3],z[3],0},{0,0,0,1}}
end

newCamera = function(position, euler, fov, screenSize)
  local camera = {translationMatrix=nil, rotationMatrix=nil, screenSize=screenSize, fov=fov}
  setCameraTranslation(camera, position)
  setCameraRotation(camera, euler)
  return camera
end

generateCameraRay = function(camera, x, y)
  x = x / camera.screenSize[1]
  y = 1 - y / camera.screenSize[2]
  local hFov = camera.fov * math.pi / 180
  local vFov = hFov * (camera.screenSize[2]) / (camera.screenSize[1])
  local sx = math.tan(hFov / 2) * 2
  local sy = math.tan(vFov / 2) * 2
  local o = transform(camera.translationMatrix, {0, 0, 0})
  local d = transform(camera.rotationMatrix, {(x - 0.5) * sx, (y - 0.5) * sy, -1})
  return newRay(o, d)
end

------- BOUNDING BOXES -------

newEmptyBoundingBox = function(dimension)
  local min = tabulate(dimension, INFINITY)
  local max = tabulate(dimension, -INFINITY)
  return {type="boundingBox", min=min, max=max}
end

newBoundingBox = function(v1, v2)
  return getBoundingBox({v1, v2})
end

getBoundingBox = newPolymorphicFunction()

getBoundingBox.default = function(vs)
  local n = #(vs[1])
  local bb = newEmptyBoundingBox(n)
  for i = 1 , #vs do
    local v = vs[i]
    for j = 1 , #v do
      if v[j] < bb.min[j] then
        bb.min[j] = v[j]
      end
      if v[j] > bb.max[j] then
        bb.max[j] = v[j]
      end
    end
  end
  return bb
end

getBoundingBox.triangle = function(tri)
  return getBoundingBox({tri[1], tri[2], tri[3]})
end

getBoundingBox.sphere = function(sph)
  local o = sph.center
  local r = sph.radius
  return newBoundingBox(sub(o,{r,r,r}), add(o,{r,r,r}))
end

expandBoundingBox = function(b, v)
  local bb = newEmptyBoundingBox(#v)
  bb.min = slice(b.min, 1, INFINITY)
  bb.max = slice(b.max, 1, INFINITY)
  for j = 1 , #v do
    if v[j] < bb.min[j] then
      bb.min[j] = v[j]
    end
    if v[j] > bb.max[j] then
      bb.max[j] = v[j]
    end
  end
  return bb
end

mergeBoundingBoxes = function(b1, b2)

  local bb = newEmptyBoundingBox(#(b1.min))
  bb.min = slice(b1.min)
  bb.max = slice(b1.max)
  for j = 1 , #(b2.min) do
    if b2.min[j] < bb.min[j] then
      bb.min[j] = b2.min[j]
    end
  end
  for j = 1 , #(b2.max) do
    if b2.max[j] > bb.max[j] then
      bb.max[j] = b2.max[j]
    end
  end
  return bb
end

inBoundingBox = function(b, v)
  for i = 1 , #(b.min) do
    if b.min[i] > v[i] or b.max[i] < v[i] then
      return false
    end
  end
  return true
end

boundingBoxCentroid = function(b)
  return scale(add(b.min, b.max), 0.5)
end

boundingBoxSurfaceArea = function(b)
  local extent = sub(b.max, b.min)
  local x = math.max(0,extent[1])
  local y = math.max(0,extent[2])
  local z = math.max(0,extent[3])
  return 2 * (x * z + x * y + y * z)
end

------- BVH -------

newBVHNode = function(boundingBox, i0, i1)
  return {type="BVHNode", boundingBox=boundingBox, iBegin=i0, iEnd=i1, leftChild=null, rightChild=null}
end


newBVH = function(primitives, maxLeafSize)
  local buildBVH
  buildBVH = function(primitives, i0, i1, maxLeafSize, bucketCount)
    if i1 - i0 > 500 then
      print("Building BVH from " .. tostring(i0) .. " to " .. tostring(i1))
    end
    local bb = newEmptyBoundingBox(3)
    for i = i0 , i1 do
      bb = mergeBoundingBoxes(bb, getBoundingBox(primitives[i]))
    end
    if ((i1 - i0 + 1) <= maxLeafSize) then
      return newBVHNode(bb, i0, i1)
    end

    local parts = {}
    for ax = 1 , 3 do
      local buckets = {}
      local lo = bb.min[ax]
      local hi = bb.max[ax]
      for i = 1 , bucketCount do
        local b = {type="BVHBucket", min=nil, max=nil, primitiveCount=0, area=0, boundingBox=newEmptyBoundingBox(3)}
        b.min = lo + (i - 1) / bucketCount * (hi - lo)
        b.max = (b.min) + (hi - lo) / bucketCount
        buckets[i] = b
      end
      
      for i = i0, i1 do
        local pbb = getBoundingBox(primitives[i])
        local c = boundingBoxCentroid(pbb)
        for j = 1 , bucketCount do
          if buckets[j].min <= c[ax] and c[ax] <= buckets[j].max then
            buckets[j].primitiveCount = buckets[j].primitiveCount + 1
            buckets[j].boundingBox = mergeBoundingBoxes(buckets[j].boundingBox, pbb)
            buckets[j].area = buckets[j].area + boundingBoxSurfaceArea(pbb)
            break
          end
        end
      end
    
      for i = 1 , bucketCount - 1 do
        local part = {type="BVHPartition", planeIndex=i, axis=ax, 
          leftPrimitiveCount=0, rightPrimitiveCount=0, 
          leftArea=0, rightArea=0, 
          leftBoundingBox=newEmptyBoundingBox(3), 
          rightBoundingBox=newEmptyBoundingBox(3), SAH=0}
        for j = 1 , i do
          part.leftPrimitiveCount = (part.leftPrimitiveCount) + (buckets[j].primitiveCount)
          part.leftArea = (part.leftArea) + buckets[j].area
          part.leftBoundingBox = mergeBoundingBoxes(part.leftBoundingBox, buckets[j].boundingBox)
        end
        for j = i + 1 , bucketCount do
          part.rightPrimitiveCount = (part.rightPrimitiveCount) + buckets[j].primitiveCount
          part.rightArea = (part.rightArea) + buckets[j].area
          part.rightBoundingBox = mergeBoundingBoxes(part.rightBoundingBox, buckets[j].boundingBox)
        end
        if (part.leftPrimitiveCount) > 0 and (part.rightPrimitiveCount) > 0 then
         
          part.SAH = boundingBoxSurfaceArea(part.leftBoundingBox) / (part.leftPrimitiveCount) + 
                     boundingBoxSurfaceArea(part.rightBoundingBox) / (part.rightPrimitiveCount)
          -- part.SAH = math.abs(part.leftPrimitiveCount-part.rightPrimitiveCount)
          parts[#parts+1]=part
        end
      end
    
    end
    
    if #parts == 0 then
      return newBVHNode(bb, i0, i1)
    end
    
    local minSAH = INFINITY
    local minPart = nil
    
    for i = 1 , #parts do
      if parts[i].SAH < minSAH then
        minSAH = parts[i].SAH
        minPart = parts[i]
      end
    end
    -- pprint(minPart)
    

    local comp = function(x)
      return boundingBoxCentroid(getBoundingBox(x))[minPart.axis]
    end

    mergeSort(primitives, i0, i1, comp)
    local node = newBVHNode(bb, i0, i1)
    
    local m = i0 + minPart.leftPrimitiveCount
    node.leftChild = buildBVH(primitives, i0, m - 1, maxLeafSize, bucketCount)
    node.rightChild = buildBVH(primitives, m, i1, maxLeafSize, bucketCount)
    return node
  end


  print "Initializing BVH..."
  local BVH = {type="BVH", primitives=slice(primitives), root=nil}
  BVH.root = buildBVH(BVH.primitives, 1, #primitives, maxLeafSize, 4)
  return BVH
end

isBVHLeaf = function(node)
  if node.type ~= "BVHNode" then
    return false
  end
  if node.leftChild == nil and node.rightChild == nil then
    return true
   end
  return false
end

rayIntersect.BVH = function(ray, bvh)
  local findClosestHit
  findClosestHit = function(ray, node)
    local closest = nil
    if isBVHLeaf(node) then
      local closestT = INFINITY
      for i = node.iBegin, node.iEnd do
        local p = bvh.primitives[i]
        local isect = rayIntersect[p.type](ray,p);
        if isect then
          if closestT > isect.t then
            closest = isect
          end
        end
      end
      return closest
    end
    
    local hitL = rayIntersect.boundingBox(ray, node.leftChild.boundingBox);
    local hitR = rayIntersect.boundingBox(ray, node.rightChild.boundingBox);
    
    if hitL and (not hitR) then
      return findClosestHit(ray, node.leftChild)
    elseif (not hitL) and hitR then
      return findClosestHit(ray, node.rightChild)
    elseif (not hitL) and (not hitR) then
      return nil
    end
    
    local first, second
    if hitL.t1 < hitR.t1 then
      first = node.leftChild
      second = node.rightChild
    else
      first = node.rightChild
      second = node.leftChild
    end

    closest = findClosestHit(ray, first)
    if (not closest) or (closest.t >= math.max(hitL.t1,hitR.t1)) then
      local closest2 = findClosestHit(ray, second)
      if closest2 then
        if (not closest) or closest2.t < closest.t then
          return closest2
        end
      end
    end
    return closest
  end

  return findClosestHit(ray, bvh.root)

end



------- IO -------

split = function(str, delim)
  local t={}
  for s in string.gmatch(str, "([^"..delim.."]+)") do
    table.insert(t, s)
  end
  return t
end

isdir = function(path)
  --https://stackoverflow.com/questions/1340230/check-if-directory-exists-in-lua
  local exists = function(filename)
    local ok, err, code = os.rename(filename, filename)
    if not ok then
      if code == 13 then
        -- Permission denied, but it exists
        return true
      end
    end
    return ok
  end
  return exists(path.."/")
end

dump = function(o,indent,printer)
  if indent == nil then indent = INFINITY end
  local _dump
  _dump = function(o,lvl)
    if type(o) == 'table' then
      printer('{')
      if indent and lvl <= indent then
        printer('\n')
      end
      for k,v in pairs(o) do
        if type(k) == 'number' then
           k = '['..k..']'
        end
        if indent and lvl <= indent then
          for i = 1, lvl do
            printer('  ')
          end
        end
        printer(k..'=')
        _dump(v,lvl+1)
        printer(', ')
        if indent and lvl <= indent then
          printer('\n')
        end
      end
      if indent and lvl <= indent then
        for i = 1, lvl-1 do
          printer('  ')
        end
      end
      printer('}')
      return s
    elseif type(o) == 'string' then
      printer('"'..o..'"')
    else
      printer(tostring(o))
    end
  end
  _dump(o,1)
end

pprint = function(o)
  dump(o,3,io.write)
  io.write("\n")
end

dumpToFile = function(filename,o)
  local file = io.open(filename,'w')
  io.output(file)
  io.write("return ")
  dump(o,3,io.write)
  io.close(file)
end

dumpToString = function(o,indent)
  local s = ""
  dump(o,indent,function(x) s = s..x end)
  return s
end

loadDumpFile = function(filename)
  return dofile(filename)
end



readOBJFile = function(filename)
  local vertices = {}
  local faces = {}
  local normals = nil
  local uvs = nil
  local groups = nil

  local vts = {}
  local vns = {}
  local objLines = {}
  print("Loading OBJ file ...")
  for line in io.lines(filename) do
    objLines[#objLines+1] = line;
  end
  for i = 1, #objLines do 
    if i % 2000 == 0 then
      print("Reading OBJ file " .. tostring(i) .. "/" .. tostring(#objLines))
    end
    repeat
      if objLines[i][1] == "#" then
        break -- continue
      end
      local lineStr = split(objLines[i], "#")[1]

      if lineStr == nil or #lineStr == 0 then
        break -- continue
      end

      local cmd = split(lineStr, " ")
      if cmd[1] == "v" then
        local v = {}
        for j = 2 , #cmd do
          v[#v+1] = tonumber(cmd[j])
        end
        vertices[#vertices+1] = v
      elseif cmd[1] == "vt" then
        local v = {}
        for j = 2, #cmd do
          v[#v+1] = tonumber(cmd[j])
        end
        vts[#vts+1]=v
      elseif cmd[1] == "vn" then
        local v = {}
        for j = 2, #cmd do
          v[#v+1] = tonumber(cmd[j])
        end
        vns[#vns+1]=v
      elseif cmd[1] == "f" then
        local f = {}
        local fuv = {}
        for j = 2 , #cmd do
          local iii = split(cmd[j],"/")
          local vi = tonumber(iii[1])
          f[#f+1] = vi
          if iii[2] and #(iii[2]) then
            fuv[#fuv+1] = tonumber(iii[2])
          end
          if iii[3] and #(iii[3]) then
            if normals == nil then
              normals = {}
            end
            normals[vi] = tonumber(iii[3])
          end
        end
        faces[#faces+1] = f
        if #fuv then
          if uvs == nil then
            uvs = {}
          end
          uvs[#uvs+1]=fuv
        end

      elseif cmd[1] == "g" then
        if groups then
          groups[#groups].iEnd = #faces
        else
          groups = {}
        end
        groups[#groups+1] = {name=cmd[2],iBegin=#faces+1,iEnd=-1}

      else
        print("Unsupported command "..cmd[1])
      end
    until true
  end
  if groups then
    groups[#groups].iEnd = #faces
  end
  if normals then
    for i=1,#normals do
      normals[i] = vns[normals[i]]
    end
  end
  if uvs then
    for i=1,#uvs do
      for j=1,#(uvs[i]) do
        uvs[i][j] = vts[uvs[i][j]]
      end
    end
  end
  print("Done loading OBJ file.")

  return {type="OBJ", vertices=vertices, faces=faces, normals=normals, uvs=uvs, groups=groups}
end

generateOBJNormals = function(obj)
  print("Generating OBJ normals...")
  local vertices = obj.vertices
  local faces = obj.faces
  local normals = {}
  for i = 1, #vertices do
    if i % 100 == 0 then
      print("Generating normal "..tostring(i).."/"..tostring(#vertices))
    end
    local ns = {}
    for j = 1, #faces do
      if contains(faces[j], i) then
         ns[#ns+1] = triangleNormal(newTriangle(
           vertices[faces[j][1]],
           vertices[faces[j][2]],
           vertices[faces[j][3]]
         ))
      end
    end
    local n = {0,0,0}
    for j = 1, #ns do
      n = add(n, ns[j])
    end
    n = scale(n,1/#ns)
    normals[i] = n
  end
  print("Done generating normals.");
  return normals
end

OBJToTriangles = function(obj, options)
  if options == nil then options = {} end
  if options.generateNormals == nil then 
    options.generateNormals = (obj.normals == nil)
  end

  local vertices = obj.vertices
  local faces = obj.faces
  local normals = nil
  if options.generateNormals then
    normals = generateOBJNormals(obj)
  end

  local triangles = {}
  for i = 1 , #faces do
    if i % 2000 == 0 then
      print ("Compiling face " .. tostring(i) .. "/" .. tostring(#faces))
    end
    local n = #faces[i]
    for j = 2 , n - 1 do
      local i1 = faces[i][1]
      local i2 = faces[i][j]
      local i3 = faces[i][j+1]
      local v1 = slice(vertices[i1])
      local v2 = slice(vertices[i2])
      local v3 = slice(vertices[i3])

      if normals then
        v1.normal = normals[i1]
        v2.normal = normals[i2]
        v3.normal = normals[i3]
      end
      
      if obj.uvs then
        v1.uv = obj.uvs[i][1]
        v2.uv = obj.uvs[i][j]
        v3.uv = obj.uvs[i][j+1]
      end
      triangles[#triangles+1] = newTriangle(v1,v2,v3);
      
    end
  end
  print "Done converting to triangles."
  return triangles
end

readPPM = function(filename)
  local idx = 1
  local buffer = {}
  local w
  local h
  local cmax
  for line in io.lines(filename) do
    if line[1] ~= "#" then
      local tokens = split(line," ")
      for i = 1, #tokens do
        if idx == 1 then
          if tokens[i] ~= "P3" then
            error "Only P3 PPM are supported"
          end
        elseif idx == 2 then
          w = tonumber(tokens[i])
        elseif idx == 3 then
          h = tonumber(tokens[i])
          for j=1,h do
            buffer[j] = {}
            for k=1,w do
              buffer[j][k] = {}
            end
          end
        elseif idx == 4 then
          cmax = tonumber(tokens[i])
          
        else
          local y = math.floor((idx-5) / (w*3))+1
          local x = (math.floor((idx-5)/3) % w) + 1
          local c = (idx-5)%3+1
          --pprint({y,x,c})
          buffer[y][x][c] = tonumber(tokens[i])/cmax
        end
        idx = idx + 1
      end
    end
  end
  return buffer
end


writePPM = function(filename, buffer)
  local overwrite = io.open(filename,'w')
  io.output(overwrite)
  io.write('')
  io.close(overwrite)


  local file = io.open(filename,'a')
  io.output(file)

  io.write("P3\n")
  io.write((#(buffer[1])) .. " " .. #(buffer) .. "\n")
  io.write("255\n")
  for i=1,#buffer do
    for j=1,#(buffer[i]) do
      local c = clampSpectrum(buffer[i][j])
      io.write( math.floor(c[1]*255) .. " "
             .. math.floor(c[2]*255) .. " "
             .. math.floor(c[3]*255) .. " ")
    end
    io.write("\n")
  end
  
  io.close(file)
end

------- SAMPLER ------- 

uniformHemisphereSampler = function()
  local Xi1 = math.random()
  local Xi2 = math.random()
  local theta = math.acos(Xi1)
  local phi = 2*math.pi*Xi2
  local xs = math.sin(theta)*math.cos(phi)
  local ys = math.sin(theta)*math.sin(phi)
  local zs = math.cos(theta)
  return normalize({xs,ys,zs})
end

textureSampler = function(texture,u,v)
  u = math.min(math.max(u,0),1)
  v = math.min(math.max(1-v,0),1)
  local h = #texture
  local w = #(texture[1])
  local i = v*(h-1)+1
  local j = u*(w-1)+1
  local i1 = math.floor(i)
  local i2 = math.ceil(i)
  local j1 = math.floor(j)
  local j2 = math.ceil(j)
  local ti = i-i1
  local tj = j-j1
  local s11 = texture[i1][j1]
  local s12 = texture[i1][j2]
  local s21 = texture[i2][j1]
  local s22 = texture[i2][j2]
  return lerp(lerp(s11,s12,tj),lerp(s21,s22,tj),ti)

  --return texture[math.floor(v*(h-1))+1][math.floor(u*(w-1))+1]

end

environmentMapSampler = function(texture, direction)
  local theta = math.acos(  dot(direction,{0,1,0})/norm(direction) );
  local phi = math.atan2(direction[3], direction[1]) + math.pi;
  
  theta = math.min(math.max(theta,0),math.pi)
  phi = math.min(math.max(phi,0),math.pi*2)
  local x = phi / (math.pi*2)
  local y = 1-theta / math.pi
  return textureSampler(texture, x, y)
end

------- BSDF -------

newBSDFSpaceMatrix = function(n)
  local z = {n[1],n[2],n[3]}
  local h = {n[1],n[2],n[3]}
  if (math.abs(h[1]) < math.abs(h[2])) and (math.abs(h[1]) < math.abs(h[3])) then
    h[1] = 1
  elseif (math.abs(h[2]) < math.abs(h[1])) and (math.abs(h[2]) < math.abs(h[3])) then
    h[2] = 1
  else
    h[3] = 1
  end
  z = normalize(z)
  local y = cross(h,z)
  if norm2(y) == 0 then
    if math.abs(h[3]) == 1 then
      z[1] = z[1] + 0.0001
    else
      z[3] = z[3] + 0.0001
    end
    z = normalize(z)
    y = cross(h,z)
  end
  y = normalize(y)
  local x = cross(z,y)
  x = normalize(x)
  return {{x[1],y[1],z[1],0},{x[2],y[2],z[2],0},{x[3],y[3],z[3],0},{0,0,0,1}}
end

newDiffuseBSDF = function(spectrum)
  return {type="diffuseBSDF",
          emission={0,0,0},
          spectrum=spectrum,
          isDelta=false}
end

newEmissionBSDF = function(spectrum)
  return {type="emissionBSDF",
          emission=slice(spectrum),
          spectrum=spectrum,
          isDelta=false}
end

newMirrorBSDF = function(spectrum)
  return {type="mirrorBSDF",
          emission={0,0,0},
          spectrum=spectrum,
          isDelta=true}
end

newGlassBSDF = function(spectrum, refractiveIndex)
  return {type="glassBSDF",
          emission={0,0,0},
          spectrum=spectrum,
          refractiveIndex = refractiveIndex or 1.45,
          isDelta=true,
          translucency=0.2}
end

newGlossyBSDF = function(spectrum)
  return {type="glossyBSDF",
          emission={0,0,0},
          spectrum=spectrum,
          specular=0.5,
          diffuse=0.05,
          shininess=20,
          mirrorEffect=0.3,
          isDelta=false}
end


newSSSBSDF = function(spectrum)
  return {type="SSSBSDF",
          emission={0,0,0},
          spectrum=spectrum,
          isDelta=false,
          transmitRatio=0.5,
          translucency=0.1,
          radius=0.1}
end

evaluateBSDF = newPolymorphicFunction()

evaluateBSDF.default = function(bsdf,wo,wi)
  return bsdf.spectrum
end

evaluateBSDF.diffuseBSDF = function(bsdf,wo,wi)
  return scale(bsdf.spectrum , (1 / math.pi));
end

evaluateBSDF.glossyBSDF = function(bsdf,wo,wi)
  local L = wi
  local V = wo
  local LN = wi[3]
  local kdLN = LN*bsdf.diffuse
  local R = reflect(wi)
  local RV =dot(normalize(R), normalize(V))
  local RVa = math.max(0,math.pow(RV, bsdf.shininess))

  local ksRVa = RVa * bsdf.specular
  local I = scale(bsdf.spectrum, kdLN+ksRVa)
  return scale(I , (0.5+bsdf.shininess*0.5) / (wi[3]*math.pi));

end

evaluateBSDF.SSSBSDF = function(bsdf,wo,wi)
  return scale(bsdf.spectrum , (1 / math.pi));
end


computeBSDFIncident = newPolymorphicFunction()

computeBSDFIncident.diffuseBSDF = function(bsdf,wo)
  return {
    wi = uniformHemisphereSampler(),
    pdf = 1/(2*math.pi)
  }
end

computeBSDFIncident.emissionBSDF = function(bsdf,wo)
  return {
    wi = uniformHemisphereSampler(),
    pdf = 1/(2*math.pi)
  }
end

computeBSDFIncident.mirrorBSDF = function(bsdf,wo)
  return {
    wi = reflect(wo),
    pdf = dot(wo,{0,0,1})
  }
end

computeBSDFIncident.glassBSDF = function(bsdf,wo)
  local ior = 1/bsdf.refractiveIndex

  local r0 = math.pow((ior-1)/(ior+1),2)
  local fr = r0 + (1-r0) * math.pow(1-dot(wo,{0,0,1}),5)
  local wi, pdf

  if math.random() < fr then
    wi = reflect(wo)
    pdf = dot(wo,{0,0,1})
  else
    wi = refract(wo, ior)
    pdf = dot(wi,{0,0,1})
    if not wi then
      wi = reflect(wo)
    end
    
  end
  return {
    wi = wi,
    pdf = pdf
  }
end

computeBSDFIncident.glossyBSDF = function(bsdf,wo)
  local wi = lerp(uniformHemisphereSampler(), reflect(wo), bsdf.mirrorEffect)
  return {
    wi = wi,
    pdf = wi[3],
    spectrum = bsdf.spectrum,
  }
end

computeBSDFIncident.SSSBSDF = function(bsdf,wo,side)
  local wi = uniformHemisphereSampler()
  local dir = 1

  if (math.random() < bsdf.transmitRatio and side > 0) 
  or (side < 0) then
    dir = -1
  end
  wi = scale(wi,dir)

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


reflect = function(wo)
  return normalize(sub(scale({0,0,1}, dot(wo,{0,0,1})*2) , wo))

end

refract = function(wo, ior)
  -- do return scale(wo,-1) end
  local n = {0,0,1}
  local test = dot(wo,n)
  local r = ior
  local cosTheta = wo[3]
  local ct = 1 - r*r*(1-cosTheta*cosTheta)
  if ct < 0 then
    return nil
  end
  return add(scale(wo,-r), scale(n, r*cosTheta-math.sqrt(ct)))
end


------- LIGHTS -------

newDirectionalLight = function(radiance,direction)
  return {type="directionalLight",
          isDelta=true,
          direction=normalize(direction),
          radiance=radiance}
end

newAreaLight = function(radiance,v1,v2,v3,v4)
  local tri1 = newTriangle(v1,v2,v3)
  local tri2 = newTriangle(v1,v3,v4)
  local direction = triangleNormal(tri1)
  local area = triangleArea(tri1)+triangleArea(tri2)
  return {type="areaLight",
          isDelta=false,
          direction=direction,
          vertices = {v1,v2,v3,v4},
          area=area,
          radiance=radiance}

end

newHemisphereLight = function(radiance)
  return {type="hemisphereLight",
          isDelta=false,
          radiance=radiance}
end

newAmbientLight = function(radiance)
  return {type="ambientLight",
          isDelta=false,
          radiance=radiance}
end


newPointLight = function(radiance,position)
  return {type="pointLight",
          isDelta=true,
          position=position,
          radiance=radiance}
end


newEnvironmentLight = function(texture)
  return {type="environmentLight",
          isDelta=false,
          texture=texture}
end

sampleLight = newPolymorphicFunction()

sampleLight.directionalLight = function(light,p)
  return {
    wi = scale(light.direction,-1),
    distance = INFINITY,
    pdf = 1.0,
    spectrum = light.radiance,
  }
end

sampleLight.areaLight = function(light,p)
  local sample = {math.random(),math.random()}
  local q = lerp(
    lerp(light.vertices[1],light.vertices[2],sample[1]),
    lerp(light.vertices[4],light.vertices[3],sample[1]),
    sample[2]
  )
  local d = sub(q,p)
  local cosTheta = dot(d, light.direction)
  local sqDist = norm2(d)
  local distance = math.sqrt(sqDist)
  local wi = normalize(d)
  local pdf = sqDist / (light.area * math.abs(cosTheta))
  local spectrum = {0,0,0}
  if cosTheta < 0 then
    spectrum = light.radiance
  end
  return {wi=wi, distance=distance, pdf=pdf, spectrum=spectrum}
end

sampleLight.hemisphereLight = function(light,p)
  local dir = uniformHemisphereSampler()
  local wi = transform({{1,0,0,0},{0,0,1,0},{0,-1,0,0},{0,0,0,1}}, dir)
  return {
    wi = wi,
    distance = INFINITY,
    pdf = 1.0/(2*math.pi),
    spectrum = light.radiance,
  }
end

sampleLight.ambientLight = function(light,p)
  local dir = uniformHemisphereSampler()
  if math.random() < 0.5 then
    dir = scale(dir,-1)
  end
  return {
    wi = normalize(dir),
    distance = INFINITY,
    pdf = 1.0/(4*math.pi),
    spectrum = light.radiance,
  }
end

sampleLight.pointLight = function(light,p)
  local d = sub(light.position,p)
  return {
    wi = normalize(d),
    distance = norm(d),
    pdf = 1.0,
    spectrum = light.radiance,
  }
end

sampleLight.environmentLight = function(light,p)
  local dir = uniformHemisphereSampler();
  if math.random() < 0.5 then
    dir = scale(dir,-1)
  end
  return {
    wi = dir,
    distance = INFINITY,
    pdf = 1/(math.pi*4),
    spectrum = environmentMapSampler(light.texture,dir),
  }
end

luminance = function(spectrum)
  return 0.2126 * spectrum[1] + 0.7152 * spectrum[2] + 0.0722 * spectrum[3];
end


tintTexture = function(texture, spectrum)
  return map(texture, function(x) return map(x,function(y) return mul(y,spectrum) end) end)
end

blurTexture = function(texture, halfKernelSize)
  local hks = halfKernelSize
  local ks = hks*2+1

  local gaussian2d = function(x,y)
    local sig = 1/5
    return math.exp(-math.pow(x-0.5,2)/(2*sig*sig)-math.pow(y-0.5,2)/(2*sig*sig))
  end
  local kernel = {}
  local s = 0
  for i=1,ks do
    kernel[i] = {}
    for j=1,ks do
      z = gaussian2d((j-0.5)/ks,(i-0.5)/ks)
      kernel[i][j] = z
      s = s+z
    end
  end
  kernel = map(kernel, function(x) return scale(x,1/s) end)

  local m = #texture
  local n = #(texture[1])
  local blurred = {}
  for i=1,m do
    if i % 20 == 0 then
      print("Blurring "..i.."/"..m)
    end
    blurred[i] = {}
    for j=1,n do
      local s = {0,0,0}
      for ki=-hks,hks do
        for kj=-hks,hks do
          local ii = (i+ki+m-1)%m+1
          local jj = (j+kj+n-1)%n+1
          s = add(s, scale(texture[ii][jj], kernel[hks+ki+1][hks+kj+1]))
        end
      end
      blurred[i][j] = s
    end
  end
  return blurred
end

------- SCENE -------

newScene = function()
  return {lights={}, materials={}, primitives={}, bvh=nil}
end

sceneAddPrimitives = function(scene,primitives,options)
  if options == nil then options = {} end
  for i=1,#primitives do
    if options.materialIndex then
      primitives[i].materialIndex = options.materialIndex
    end
    if options.transform then
      scene.primitives[#scene.primitives+1] = transform(options.transform,primitives[i])
    else
      scene.primitives[#scene.primitives+1] = primitives[i]
    end
  end
end

newMaterial = function(bsdf,texture)
  return {type="Material",bsdf=bsdf,texture=texture}
end

sceneAddMaterial = function(scene,material)
  scene.materials[#scene.materials+1] = material
end

sceneAddLight = function(scene,light)
  scene.lights[#scene.lights+1] = light
end

sceneAddEnvironmentMap = function(scene, texture, options)
  if options == nil then options = {} end
  if options.spectrum then
    texture = tintTexture(texture, options.spectrum)
  end
  if options.blur then
    texture = blurTexture(texture, options.blur)
  end
  scene.environmentMap = texture
  scene.environmentMap.invisible = options.invisible
end

sceneCompileBVH = function(scene)
  scene.bvh = newBVH(scene.primitives,4);
  scene.primitives = nil
end

sceneNormalizePrimitives = function(scene, options)
  if options == nil then options = {} end
  options.scale = options.scale or 1.0

  local bb = newEmptyBoundingBox(3)
  for i=1,#scene.primitives do
    bb = mergeBoundingBoxes(bb, getBoundingBox(scene.primitives[i]))
  end
  local centroid = boundingBoxCentroid(bb)
  local extent = sub(bb.max, bb.min)
  -- pprint(bb)
  local scl = options.scale/math.max(extent[1],extent[2],extent[3])
  local offset = sub({0,0,0},centroid)
  local y = extent[2]*scl/2

  local trans = multiplyMatrices(multiplyMatrices(
    newTransformationMatrix("translate",{0,y,0}),
    newTransformationMatrix("scale",scl)),
    newTransformationMatrix("translate",offset))
  for i=1,#scene.primitives do
    scene.primitives[i] = transform(trans,scene.primitives[i])
  end

end

sceneWrapCornellBox = function(scene, options)
  if options == nil then options = {} end
  options.scale = options.scale or 0.9
  options.leftSpectrum = options.leftSpectrum or {1,0,0}
  options.rightSpectrum = options.rightSpectrum or {0,1,1}

  sceneNormalizePrimitives(scene, options)

  local v1={-0.5, 0,  0.5}
  local v2={-0.5, 0, -0.5}
  local v3={-0.5, 1, -0.5}
  local v4={-0.5, 1,  0.5}
  local v5={ 0.5, 0,  0.5}
  local v6={ 0.5, 0, -0.5}
  local v7={ 0.5, 1, -0.5}
  local v8={ 0.5, 1,  0.5}

  sceneAddMaterial(scene, newMaterial(newDiffuseBSDF(options.leftSpectrum),nil))
  sceneAddPrimitives(scene,{newTriangle(v4,v3,v2),newTriangle(v4,v2,v1)},{materialIndex=#scene.materials})

  sceneAddMaterial(scene, newMaterial(newDiffuseBSDF(options.rightSpectrum),nil))
  sceneAddPrimitives(scene,{newTriangle(v6,v7,v8),newTriangle(v6,v8,v5)},{materialIndex=#scene.materials})

  sceneAddMaterial(scene, newMaterial(newDiffuseBSDF({1,1,1}),nil))
  sceneAddPrimitives(scene,{newTriangle(v2,v6,v5),newTriangle(v2,v5,v1)},{materialIndex=#scene.materials})
  sceneAddPrimitives(scene,{newTriangle(v3,v7,v6),newTriangle(v3,v6,v2)},{materialIndex=#scene.materials})
  sceneAddPrimitives(scene,{newTriangle(v8,v7,v3),newTriangle(v8,v3,v4)},{materialIndex=#scene.materials})
  
  local lw = 0.25
  local l1={-lw,0.99,-lw}
  local l2={lw,0.99,-lw}
  local l3={lw,0.99,lw}
  local l4={-lw,0.99,lw}

  sceneAddMaterial(scene, newMaterial(newEmissionBSDF({100,100,100}),nil))
  sceneAddPrimitives(scene,{newTriangle(l1,l2,l3),newTriangle(l1,l3,l4)},
    {materialIndex=#scene.materials})

  local ep = {0,-0.01,0}
  local light = newAreaLight({10,10,10},add(l1,ep),add(l2,ep),add(l3,ep),add(l4,ep))
  -- local light = newDirectionalLight({1,1,1},{0,0,-1})
  sceneAddLight(scene,light)

end



------- RENDERING -------

renderNDotL = function(cam, lightDirection, bvh, texture)
  local l = normalize(lightDirection)
  local buffer = {}
  for y=1,cam.screenSize[2] do
    if y % 20 == 0 then
      print("Rendering "..y.."/"..cam.screenSize[2])
    end
    buffer[y] = {}
    for x=1,cam.screenSize[1] do
      local ray = generateCameraRay(cam,x,y)
      local isect = rayIntersect.BVH(ray,bvh)
      if isect then
        local n = isect.normal
        local ndl = dot(n,l)

        --local p = add(ray.origin,scale(ray.direction, isect.t))
        local c = math.max(0,math.max(0,ndl)*0.9)+0.1
        
        local spec = {1,1,1}
        if isect.uv and texture then
          spec = textureSampler(texture,isect.uv[1],isect.uv[2])
          --pprint(spec)
        end
        buffer[y][x] = {c*spec[1],c*spec[2],c*spec[3]}
      else
        buffer[y][x] = {0,0,0}
      end
    end
  end
  return buffer
end

clampSpectrum = function(s)
  return {
    math.min(math.max(math.abs(s[1]),0),1),
    math.min(math.max(math.abs(s[2]),0),1),
    math.min(math.max(math.abs(s[3]),0),1),
  }
end

traceRay = function(ray,scene)

  local getMaterial = function(isect)
    if isect == nil then
      return nil
    end
    if isect.primitive.materialIndex then
      return scene.materials[isect.primitive.materialIndex]
    else
      return newMaterial(newDiffuseBSDF({1,1,1}),nil)
    end
  end

  if not scene.bvh then
    sceneCompileBVH(scene)
  end

  local isect = rayIntersect.BVH(ray,scene.bvh)
  local material = getMaterial(isect)
  local void = {0,0,0}

  if isect == nil or (material.invisible and ray.depth == 0) then
    if scene.environmentMap then
      if (scene.environmentMap.invisible and ray.depth == 0) then
        return void
      else
        return environmentMapSampler(scene.environmentMap, ray.direction)
      end
    end
    return void
  end
  

  local fragColor = {1,1,1}
  if material.texture and isect.uv then
    fragColor = textureSampler(material.texture,isect.uv[1],isect.uv[2])
  end

  local Lout = material.bsdf.emission;

  local hitP = add(ray.origin, scale(ray.direction, isect.t))
  local hitN = isect.normal;

  local o2w = newBSDFSpaceMatrix(slice(hitN));
  local w2o = transpose(o2w)

  local wOut = transform(w2o, sub(ray.origin, hitP))
  wOut = normalize(wOut)

  if not material.bsdf.isDelta then
    for _,light in pairs(scene.lights) do
      local numLightSamples = 1
      if not light.isDelta then
        numLightSamples = NUM_SAMPLES_AREA_LIGHT
      end
      for i = 1, numLightSamples do
        local s = sampleLight(light, hitP)
        local dirToLight = s.wi
        local distToLight = s.distance
        local pr = s.pdf
        local lightL = s.spectrum

        local wIn = transform(w2o, dirToLight)

        if wIn[3] >= 0 then
          local shadowRay = newRay(add(hitP,scale(dirToLight,EPSILON)),dirToLight)
          shadowRay.tMax = distToLight
          local inShadow = rayIntersect.BVH(shadowRay,scene.bvh)
          local m = 1
          if inShadow then
            m = getMaterial(inShadow).bsdf.translucency or 0
          end
          if m > 0 then
            local f = mul(evaluateBSDF(material.bsdf,wOut,wIn), fragColor)
            Lout = add( Lout, scale( mul(f,lightL), m*wIn[3]/(numLightSamples*pr) ))
          end
        end
      end
    end

  end
  local c = computeBSDFIncident(material.bsdf,wOut,isect.side)
  local f = mul(c.spectrum or evaluateBSDF(material.bsdf, wOut, c.wi), fragColor)
  local terminateProbability = 1 - luminance(f)

  if math.random() < terminateProbability then
    return Lout--clampSpectrum(Lout)
  end

  local rd = transform(o2w, c.wi)
  local r = newRay(add(hitP,scale(rd,EPSILON)),rd)
  r.depth = ray.depth + 1
  r.tMax = c.tMax or INFINITY
  if r.depth > MAX_RAY_DEPTH then
    return Lout--clampSpectrum(Lout)
  end

  local s = mul(f,traceRay(r, scene))
  Lout = add(Lout, scale(s, math.abs(c.wi[3])/(c.pdf*(1-terminateProbability))))
  -- Lout = scale(s, c.wi[3]/(c.pdf*(1-terminateProbability)))
  return Lout--clampSpectrum(Lout)

end

render = function(cam,scene)
  local numSamples = NUM_SAMPLES_ANTIALIAS
  

  local tracePixel = function(x,y)
    local d
    if numSamples == 1 then
      d = {0.5,0.5}
    else
      d = {math.random(),math.random()}
    end
    local ray = generateCameraRay(cam,x+d[1],y+d[2])
    local s = traceRay(ray,scene)
    return s
  end

  local buffer = {}
  for y=1,cam.screenSize[2] do
    buffer[y] = {}
    for x=1,cam.screenSize[1] do
      buffer[y][x] = {0,0,0}
    end
  end
  for i=1,numSamples do
    for y=1,cam.screenSize[2] do
      if y % 10 == 0 then
        print("Rendering: pass "..i.."/"..numSamples..", row "..y.."/"..cam.screenSize[2])
      end
      for x=1,cam.screenSize[1] do
        local s = clampSpectrum(tracePixel(x,y))
        local t = 1/i
        buffer[y][x] = lerp(buffer[y][x],s,t)
      end
    end
    if isdir("tmp") and i%4==1 then
      writePPM("tmp/render-tmp-"..i..".ppm", buffer)
    end
  end

  return buffer

end
