import cv2
import numpy as np
import sys

def fromPNG(path,w=0,h=0):
  im = cv2.imread(path)
  if w == 0:
    h,w,c = im.shape
  else:
    im = cv2.resize(im,(w,h))
  result = "P3\n"
  result += str(w)+" "+str(h)+"\n255\n"
  for i in range(0,h):
    for j in range(0,w):
      result += str(im[i][j][2])+" "+str(im[i][j][1])+" "+str(im[i][j][0])+" "
    result += "\n"
  return result

def fromEXR(path,w=0,h=0,precision=3):
  p = 0.1**precision
  im = cv2.imread(path,cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)
  if w == 0:
    h,w,c = im.shape
  else:
    im = cv2.resize(im,(w,h))
  m = np.max(im)
  result = "P3\n"
  result += str(w)+" "+str(h)+"\n"+str(int(round(1/p)))+"\n"
  for i in range(0,h):
    for j in range(0,w):
      result += str(int(im[i][j][2]/p))+" "+str(int(im[i][j][1]/p))+" "+str(int(im[i][j][0]/p))+" "
    result += "\n"
  return result

if __name__ == "__main__":
  if len(sys.argv) == 2:
    pth = sys.argv[1]
    w = 0
    h = 0
  elif len(sys.argv) == 3:
    pth = sys.argv[1]
    w = int(sys.argv[2].split("x")[0])
    h = int(sys.argv[2].split("x")[1])
  else:
    print("Usage: python convert_ppm.py path/to/image.png 512x512 > path/to/output.ppm")
    exit()

  if pth.split(".")[-1] in ["jpg","png"]:
    print(fromPNG(pth,w=w,h=h))
  elif pth.split(".")[-1] in ["exr"]:
    print(fromEXR(pth,w=w,h=h))
  else:
    print("Image format not supported")



