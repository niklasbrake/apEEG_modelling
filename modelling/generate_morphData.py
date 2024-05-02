"""Converts the morphology information of a neuron into a matlab file.

To run this script, navigate to the folder of the neuron model, e.g.
    cd $baseFolder/data_files/neuron_models/L6_BPC_cADpyr231_3

Ensure that the mechansim files are compiled using the neuron funciton
    nrnivmodl ./mechanisms

Then run this function as a script
    python $baseFolder/modelling/unitary_AP_response/analysis/generate_morphData.py

Output: out.swc and morphData.mat
"""

import sys
import os
import neuron
import LFPy
import numpy as np
from scipy.io import savemat
from hoc2swc import neuron2swc

def morph2Segs(cellID):
    data = np.loadtxt(cellID)
    arbours = list()
    somaPoints = np.argwhere(data[:,1]==1)
    for i in somaPoints:
        temp = get_child(int(data[i,0]),data).tolist()
        if(len(temp)>0):
            arbours.append((i,temp))
    segs = list()
    cons = list()
    for i in range(len(arbours)):
        for j in arbours[i][1]:
            segs,cons = get_segments(np.array([j]),arbours[i][0],-1,segs,cons,data)

    endPoints = [[x[i] for i in (0,-1)] for x in segs]
    c0 = np.array([data[i[0]-1,2:6] for i in endPoints])
    c1 = np.array([data[i[1]-1,2:6] for i in endPoints])
    compType0 = np.array([data[i[0]-1,1,None] for i in endPoints])
    compType1 = np.array([data[i[1]-1,1,None] for i in endPoints])
    pts3d = np.concatenate((c0,c1,compType0,compType1),1)
    return pts3d,cons,segs,data

def get_segments(node,parentNode,parentSegment,segs,cons,data):
    i = len(segs)
    segs.append(np.append(parentNode,advance_segment(node,data)))
    cons.append([parentSegment,i])
    child = get_child(segs[-1][-1],data).tolist()
    if(len(child)>0):
        segs,cons = get_segments([child[0]],segs[i][-1],i,segs,cons,data)
        for j in range(1,len(child)):
            segs,cons = get_segments([child[j]],segs[i][-1],i,segs,cons,data)
    return segs,cons

def advance_segment(branch,data):
    node = branch[-1]
    childIdcs = get_child(node,data)
    if(len(childIdcs)==1):
        if(len(branch)==1):
            branch = np.append(branch,childIdcs)
            branch = advance_segment(branch,data)
        else:
            L = np.sum(np.linalg.norm(np.diff(data[branch-1,2:5],1,0),2,1)) # DEFINE MAX SEGMENT LENGTH TO BE 20 UM
            if(L<=20):
                branch = np.append(branch,childIdcs)
                branch = advance_segment(branch,data)
    elif(len(branch)>1):
        branch = branch[:-1]
    return branch

def get_child(node,data):
    isChild = data[node-1:,-1]==node
    isDendrite = (data[node-1:,1]==3)+(data[node-1:,1]==4)
    return node+np.argwhere(isDendrite*isChild).flatten()

if __name__ == '__main__':
    path0 = os.getcwd();
    path1 = os.path.join(path0,'morphology.hoc')

    # Find template name within the file template.hoc
    templatefile = os.path.join(path0,'template.hoc')
    with open(templatefile, newline='') as file:
         lines = file.readlines()
         for line in lines:
            if(line.find('begintemplate')>-1):
                temp = line.split(' ')
                templatename = temp[1][:-1]
                print(templatename)

    # LFPy stuff to generate python object for neuron model
    cellParameters = {
        'morphology' : path1,
        'templatefile' :  templatefile,
        'templatename' :  templatename,
        'templateargs' :  0,
        'Ra': 100,
        'passive' : False,
        'celsius': 34,
        'v_init': -65
    }
    neuron.h.load_file("biophysics.hoc")
    cell = LFPy.TemplateCell(**cellParameters)

    # Use neuron2sec to convert morphology of current neuron model to swc
    neuron2swc("out.swc")

    # Use function morph2Segs to read swc into python and then save as mat file
    pts3d,connections,segs,morphData = morph2Segs("out.swc")
    S = np.array(segs,dtype=object)
    savemat('morphData.mat', {'connections':connections,'data':morphData,'segs':S})
