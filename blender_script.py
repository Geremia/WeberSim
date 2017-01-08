#!/usr/bin/python3.5

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
data=[]
for f in files:
    for i in fileNameToData(f):
        data.append(i)

#add circle:
bpy.ops.object.mode_set(mode='OBJECT')
bpy.ops.mesh.primitive_circle_add()

bpy.ops.object.mode_set(mode='EDIT')
obj = bpy.data.objects['Circle']
bm = bmesh.from_edit_mesh(obj.data)

#delete all vertices
for v in bm.verts:
    bm.verts.remove(v)

#add back correct number vertices
for i in range(n):
    bm.verts.new((0,0,0))

#Switch to object mode
bpy.ops.object.mode_set(mode = 'OBJECT')
#Create initial Basis Shape Key
bpy.ops.object.shape_key_add()
#Create initial Shape Key
bpy.ops.object.mode_set(mode = 'OBJECT')
bpy.ops.object.shape_key_add()
bpy.data.shape_keys[0].key_blocks["Key 1"].name = 'Key 1'
bpy.data.shape_keys[0].key_blocks["Key 1"].value = 1
bpy.data.shape_keys[0].key_blocks["Key 1"].keyframe_insert("value", frame=0)

#add in keyframes
#Key a single frame per data file
frames = len(files)
framecount = 0
rowcount = 0
obj = bpy.data.objects['Circle']
bpy.ops.object.mode_set(mode='EDIT')
mesh = obj.data
bm = bmesh.from_edit_mesh(obj.data)
for j in range(1, frames):
    for keyname in bpy.data.shape_keys[0].key_blocks.keys():
        bpy.data.shape_keys[0].key_blocks[keyname].value = 0
        bpy.data.shape_keys[0].key_blocks[keyname].keyframe_insert("value", frame = framecount)
        for i in range(0,len(bm.verts)-1):
            try:
                bpy.ops.object.mode_set(mode = 'OBJECT')
                bm.verts[i].co.x = data[rowcount][1]
                bm.verts[i].co.y = data[rowcount][2]
                bm.verts[i].co.z = data[rowcount][3]
                rowcount += 1
                bmesh.update_edit_mesh(obj.data)
            except:
                pass
    framecount += 1
