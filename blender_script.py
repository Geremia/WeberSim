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
         data.append(list(map(float,i.split()[1:4])))
     return data

#number of objects / vertices
n = len(fileNameToData("000000.txt")) 

#Read in all data files.
#(ibid. pp. 8-6 – 8-7)
files = glob.glob('*.txt')
files.sort()
data = [] #data is of the format ((frame 1 verts: (2,3,3), (3,6,2), …), (frame 2 verts: (0,5,2), (3,3,5), …), …)
for f in files:
    data.append(fileNameToData(f))

n_frames = len(files)
bpy.context.scene.frame_end = n_frames

#add mesh object:
bpy.ops.mesh.primitive_circle_add(vertices=n)
obj = bpy.context.active_object

for i_frame in range(n_frames):
    block = obj.shape_key_add(name=str(i_frame), from_mix=False)  # returns a key_blocks member
    block.value = 1.0
    block.mute = True
    vertex_idx = 0
    for (vert, co) in zip(block.data, data[i_frame]):
        vert.co = co
    block.mute = True
    block.keyframe_insert(data_path='mute', frame=0, index=-1)
    block.mute = False
    block.keyframe_insert(data_path='mute', frame=i_frame + 1, index=-1)
    block.mute = True
    block.keyframe_insert(data_path='mute', frame=i_frame + 2, index=-1)
