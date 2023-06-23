 # -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 10:11:53 2021

@author: admin
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage.morphology as ni
from skimage.filters import gaussian
from skimage import io
import os
import nrrd
from scipy import ndimage

from allensdk.api.queries.ontologies_api import OntologiesApi
from allensdk.core.structure_tree import StructureTree
from allensdk.core.reference_space import ReferenceSpace
from allensdk.api.queries.mouse_connectivity_api import MouseConnectivityApi
from allensdk.config.manifest import Manifest

def rotate_3d_matrix(array_3d,angle = 0):   
    this_slice = array_3d[0,:,:]
    this_rotated = ndimage.rotate(this_slice, angle, reshape=True)
    rotated_matrix = np.ndarray([array_3d.shape[0],this_rotated.shape[0],this_rotated.shape[1]],dtype=np.int8)
    for idx in range(array_3d.shape[0]):
        this_slice = array_3d[idx,:,:]
        if len(np.where(this_slice==1)[0])>0:
            this_rotated = ndimage.rotate(this_slice, angle, reshape=True)
            this_rotated[this_rotated>=.2]=1
            this_rotated[this_rotated<.2]=0
            rotated_matrix[idx,:,:] = this_rotated
    return rotated_matrix

def smooth_masks(masks):
    for i in range(masks.shape[0]):
        img = masks[i,:,:]
        dst = gaussian(1-img,sigma=(1))
        dst[dst<=.01]=1
        dst[dst<.1]=0
        masks[i,:,:]=dst
    return masks

def reduce_mask_overlap(masks):
    for i in range(masks.shape[0]):
        this_mask = masks[i,:,:]
        other_masks = masks[np.where([np.arange(masks.shape[0])!=i])[1],:,:]
        other_masks_merge = np.max(other_masks,axis=0)
        overlap = np.where((this_mask+other_masks_merge)>1)
        this_mask[overlap[0],overlap[1]]=0
        masks[i,:,:]=this_mask
    return masks

def masks_to_outlines(masks):
    outlines = np.zeros(masks.shape,dtype=np.int8)
    for i in range(outlines.shape[0]):
        mask = masks[i,:,:]
        outlines[i,:,:] = mask - ni.binary_erosion(mask)
    return outlines

def surface_projection(volume):
    shell_volume = np.zeros(volume.shape,np.int8)
    for ii in range(volume.shape[0]):
        for jj in range(volume.shape[2]):
            zplane = np.where(volume[ii,:,jj]>0)[0]
            if len(zplane)>0:
                shell_volume[ii,zplane[0],jj] = 1
    return shell_volume

### access allen data 
# the annotation download writes a file, so we will need somwhere to put it
annotation_dir = 'annotation'
Manifest.safe_mkdir(annotation_dir)

annotation_path = os.path.join(annotation_dir, 'annotation.nrrd')

# this is a string which contains the name of the latest ccf version
annotation_version = MouseConnectivityApi.CCF_VERSION_DEFAULT

mcapi = MouseConnectivityApi()
mcapi.download_annotation_volume(annotation_version, 10, annotation_path)

annotation, meta = nrrd.read(annotation_path)

oapi = OntologiesApi()
structure_graph = oapi.get_structures_with_sets([1])  # 1 is the id of the adult mouse structure graph

# This removes some unused fields returned by the query
structure_graph = StructureTree.clean_structures(structure_graph)  

tree = StructureTree(structure_graph)

# build a reference space from a StructureTree and annotation volume, the third argument is 
# the resolution of the space in microns
rsp = ReferenceSpace(tree, annotation, [10, 10, 10])

# make list of cortical areas
area_list = [
#    tree.get_structures_by_acronym(['Isocortex'])[0]['id'],
    tree.get_structures_by_acronym(['VISp'])[0]['id'],
    tree.get_structures_by_acronym(['VISpl'])[0]['id'],
    tree.get_structures_by_acronym(['VISpor'])[0]['id'],
    tree.get_structures_by_acronym(['VISl'])[0]['id'],
    tree.get_structures_by_acronym(['VISli'])[0]['id'],
    tree.get_structures_by_acronym(['VISal'])[0]['id'],
    tree.get_structures_by_acronym(['VISrl'])[0]['id'],
    tree.get_structures_by_acronym(['VISa'])[0]['id'],
    tree.get_structures_by_acronym(['VISam'])[0]['id'],
    tree.get_structures_by_acronym(['VISpm'])[0]['id'],
    tree.get_structures_by_acronym(['FRP'])[0]['id'],
    tree.get_structures_by_acronym(['PL'])[0]['id'],
    tree.get_structures_by_acronym(['ACAd'])[0]['id'],
    tree.get_structures_by_acronym(['MOs'])[0]['id'],
    tree.get_structures_by_acronym(['MOp'])[0]['id'],
    tree.get_structures_by_acronym(['RSPv'])[0]['id'],
    tree.get_structures_by_acronym(['RSPd'])[0]['id'],
    tree.get_structures_by_acronym(['RSPagl'])[0]['id'],
    tree.get_structures_by_acronym(['SSp-tr'])[0]['id'],
    tree.get_structures_by_acronym(['SSp-ll'])[0]['id'],
    tree.get_structures_by_acronym(['SSp-ul'])[0]['id'],
    tree.get_structures_by_acronym(['SSp-m'])[0]['id'],
    tree.get_structures_by_acronym(['SSp-n'])[0]['id'],
    tree.get_structures_by_acronym(['SSp-un'])[0]['id'],
    tree.get_structures_by_acronym(['SSp-bfd'])[0]['id'],
    tree.get_structures_by_acronym(['SSs'])[0]['id'],
    tree.get_structures_by_acronym(['AUDp'])[0]['id'],
    tree.get_structures_by_acronym(['AUDv'])[0]['id'],
    tree.get_structures_by_acronym(['AUDd'])[0]['id'],
    tree.get_structures_by_acronym(['AUDpo'])[0]['id'],
    tree.get_structures_by_acronym(['TEa'])[0]['id'],
]

areas = ['VISp','VISpl','VISpor','VISl','VISli','VISal','VISrl','VISa','VISam','VISpm','FRP',
         'PL','ACAd','MOs','MOp','RSPv','RSPd','RSPagl','SSp-tr','SSp-ll','SSp-ul','SSp-m',
         'SSp-n','SSp-un','SSp-bfd','SSs','AUDp','AUDv','AUDd','AUDpo','TEa']

# make surface mask for isocortex. This is to avoid sub-surface borders appearing on surface border map.
isoctx = rsp.make_structure_mask([315])

# copy and split whole isoctx to allow R and L identities
right_isoctx = np.array(isoctx,dtype=np.int8)
left_isoctx = np.array(isoctx,dtype=np.int8)
right_isoctx[:,:,:int(isoctx.shape[2]/2)] = 0
left_isoctx[:,:,int(isoctx.shape[2]/2):] = 0

# rotate L and R hemispheres
rotated_r_ctx = rotate_3d_matrix(right_isoctx,angle = 0)
rotated_l_ctx = rotate_3d_matrix(left_isoctx,angle = 0)

## get a surface erosion from rotated isoctx
rotated_rctx_surface = surface_projection(rotated_r_ctx)
rotated_lctx_surface = surface_projection(rotated_l_ctx)
    
# plt.figure()
# plt.imshow(rotated_rctx_surface.max(axis=1))
# plt.figure()
# plt.imshow(rotated_lctx_surface.max(axis=1))


#### rotate each area and multiply by rotated surface shell
r_masks = np.zeros([len(area_list),rotated_rctx_surface.shape[0],rotated_rctx_surface.shape[2]],dtype=np.int8)
l_masks = np.zeros([len(area_list),rotated_lctx_surface.shape[0],rotated_lctx_surface.shape[2]],dtype=np.int8)
for ii in range(len(area_list)):
    print(str(ii) + ' of ' + str(len(area_list)))
    # get target area and rotat
    mask = rsp.make_structure_mask([area_list[ii]])
    rotate_mask = rotate_3d_matrix(mask,angle = 0)
    #multiply by l and r shell masks to get surface of each. Flatten along top axis
    r_masks[ii,:,:] = np.max(rotate_mask * rotated_rctx_surface,axis=1)
    l_masks[ii,:,:] = np.max(rotate_mask * rotated_lctx_surface,axis=1)
    
    
#combine L and R masks. Since R will occlude parts of L, place R after for reducing overlap code to work properly
all_masks = np.concatenate((l_masks,r_masks))
all_masks = smooth_masks(all_masks)
all_masks = reduce_mask_overlap(all_masks)
outlines = masks_to_outlines(all_masks)

# plt.imshow(sum(all_masks),cmap='Greys') # show all masks
# plt.imshow(sum(all_masks[:31,:,:]),cmap='Greys') # show left hemispherel masks
# plt.imshow(sum(all_masks[31:,:,:]),cmap='Greys') # show right hemispherel masks

all_areas=[]
for i in range(0,len(areas*2)):
    if i >= len(areas):
        direction = 'r_'
        i=i-len(areas)
    else:
        direction = 'l_'
    area = areas[i]
    all_areas.append(direction+area)

#flatten all masks
flat_ccf = outlines.max(axis=0)
plt.imshow(flat_ccf, cmap='Greys')
plt.imshow(all_masks[31], cmap='Greys')

all_masks = np.array(all_masks,dtype=np.int8)
outlines = np.array(outlines,dtype=np.int8)


if 0:
    svfolder = "C:\\Users\\McCormick Lab\\Documents\\Python\\Hulsey_A1V1M2_CCF_affine_Jan1122\\0deg\\"
    os.mkdir(svfolder)
    np.save(svfolder+'outlines',outlines)
    np.save(svfolder+'masks',all_masks)
    np.save(svfolder+'areas',all_areas)
    
    img = np.zeros([flat_ccf.shape[0],flat_ccf.shape[1],3])
    img[:,:,0]=(flat_ccf-1)*-255
    img[:,:,1]=(flat_ccf-1)*-255
    img[:,:,2]=(flat_ccf-1)*-255
    #img=img.astype(uint8)
    io.imsave(svfolder+'CCF.png', img)
        
    # for idx in range(outlines.shape[0]):
    #     outline = outlines[idx,:,:]
    #     img = np.zeros([outline.shape[0],outline.shape[1],3])
    #     img[:,:,0]=(outline-1)*-255
    #     img[:,:,1]=(outline-1)*-255
    #     img[:,:,2]=(outline-1)*-255
    #     cv2.imwrite(svfolder+'outlines//'+areas[idx]+'.png', img)
    
    # for idx in range(masks.shape[0]):
    #     mask = masks[idx,:,:]
    #     img = np.zeros([mask.shape[0],mask.shape[1],3])
    #     img[:,:,0]=(mask-1)*-255
    #     img[:,:,1]=(mask-1)*-255
    #     img[:,:,2]=(mask-1)*-255
    #     cv2.imwrite(svfolder+'masks//'+areas[idx]+'.png', img)