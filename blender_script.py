#!/usr/bin/python2.5

import bpy #Blender Python
import bmesh #Blender mesh creation pkg
import glob
import os

#change working dir.
os.chdir("/home/alan/Documents/Weber n-body sim/")

#function def:
def fileNameToData(filepath):
     #import first data file:
     file = open(filepath,'r')
     lines = file.readlines()[17:] #read 18th lines onward (used Python slice)
                                                #(ignore data files header info)
     data=[]
     for i in lines:
         data.append(list(map(float,i.split())))
     return data

#number of objects simulated
n = len(fileNameToData("000000.txt")) 

#Read in all data files.
#(ibid. pp. 8-6 – 8-7)
files = glob.glob('*.txt')
files.sort()
for f in files:
    csvin = fileNameToData(f)
    data = [row for row in csvin]

#add circle:
bpy.ops.object.mode_set(mode='OBJECT')
bpy.ops.mesh.primitive_circle_add()

obj = bpy.data.objects['Circle']
mesh = obj.data
bpy.ops.object.mode_set(mode='EDIT')
bm = bmesh.from_edit_mesh(obj.data)
rowcount = 0

#delete all vertices
for v in bm.verts:
    bm.verts.remove(v)

#add back correct number vertices
for i in range(n):
    bm.verts.new((0,0,0))

#Create the mesh object for animating changes
# cf. «Scientific Vis. w/ Blender» p. 8-6 (PDF p. 65)
#Create initial Basis Shape Key
bpy.ops.object.mode_set(mode='OBJECT')
bpy.ops.object.shape_key_add()

#add in keyframes
#Key a single frame per data file
frames = len(files)
framecount = 0
for j in range(1, frames+1):
    for keyname in bpy.data.shape_keys[0].key_blocks.keys():
        if keyname != os.path.basename(files[framecount-1]):
            bpy.data.shape_keys[0].key_blocks[keyname].value = 0
            bpy.data.shape_keys[0].key_blocks[keyname].keyframe_insert("value", frame = framecount)
            for i in range(0,len(bm.verts)-1):
                try:
                     bm.verts[i].co.x = data[rowcount][1]
                     bm.verts[i].co.y = data[rowcount][2]
                     bm.verts[i].co.z = data[rowcount][3]
                     rowcount += 1
                     bmesh.update_edit_mesh(obj.data)
                except:
                    pass
    framecount += 1
