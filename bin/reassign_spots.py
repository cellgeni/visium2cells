import pandas as pd
import anndata
import numpy as np
import cv2 as cv
import fire
import yaml
import json



def ReadConfFile(FilePath):
    with open(FilePath, 'r') as file:
        data = yaml.safe_load(file)
    path_adata = data['path_adata']
    path_tmat = data['path_tmat_json']
    downscale = data['downscale']
    out_file_path = data['out_file_path']
    cx = data['center_image_x']
    cy = data['center_image_y']
    return path_adata, path_tmat,downscale, out_file_path, cx, cy

def adjust_transformation_for_new_center(T_center, cx, cy):
    # Translation matrix to move points to the origin (0, 0)
    T_translate = np.array([
        [1, 0, -cx],
        [0, 1, -cy],
        [0, 0, 1]
    ])
    # Translation matrix to move points back after transformation
    T_translate_back = np.array([
        [1, 0, cx],
        [0, 1, cy],
        [0, 0, 1]
    ])
    # New transformation matrix with rotation around (0, 0)
    T_new = T_translate @ T_center @ T_translate_back

    return T_new

def update_spot_pos(adata, tr_mat_dict, downscale=2, Cx = 15000, Cy = 15000):
    list_of_section_names = list(adata.uns['spatial_affine']['alignment_metadata'].keys())
    idx0 = adata.obs.index.str.contains(list_of_section_names[0])
    list_of_section_names = list_of_section_names[1:]
    i=0; spot_pos_upd = np.zeros_like(adata.obsm['spatial_affine'])
    spot_pos_upd[idx0] = adata.obsm['spatial_affine'][idx0]
    TrMat = np.array([[1,0,0],[0,1,0], [0,0,1]])
    for section_name in list_of_section_names:
        print(section_name)
        idx = adata.obs.index.str.contains(section_name)
        spot_pos_section = adata.obsm['spatial_affine'][idx]
        points = np.expand_dims(spot_pos_section, axis=1) # Convert points to shape (N, 1, 2) as expected by OpenCV
        tr_mat = np.array(tr_mat_dict[str(i)]);
        tr_mat = adjust_transformation_for_new_center(tr_mat, Cx, Cy)
        
        tr_mat[0,2]*=downscale; tr_mat[1,2]*=downscale
        
        tr_mat = np.linalg.inv(tr_mat)

        spot_pos_section_upd = cv.transform(points, tr_mat[:2])
        spot_pos_upd[idx] = spot_pos_section_upd[:,0,:]
        i+=1
    return spot_pos_upd


def main(ConfFilePath):
    path_adata, path_tmat, downscale, out_file_path, cx, cy = ReadConfFile(ConfFilePath)
    with open(path_tmat) as json_data:
        tmat_all = json.load(json_data)
    adata = anndata.read_h5ad(path_adata)
    spot_pos_upd = update_spot_pos(adata, tmat_all, downscale, cx, cy)
    adata.obsm['spatial_affine_postreg'] = spot_pos_upd
    adata.write_h5ad(out_file_path)
    
if __name__ == "__main__":
    fire.Fire(main) 
