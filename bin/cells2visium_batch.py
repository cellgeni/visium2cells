import yaml
import pandas as pd
from cells2visium import main as C2V
import fire

def ReadConfFile(FilePath):
    with open(FilePath, 'r') as file:
        data = yaml.safe_load(file)
    csv_file_path = data['table_with_paths']
    out_folder = data['output_folder']
    background_thresh = data['background_threshold_intensity']
    skip_failed = data['skip_failed_samples']#
    prob_thresh = data['segmentation_prob_thresh']
    nms_thresh = data['segmentation_nms_thresh']
    pmin = data['norm_pmin']
    pmax = data['norm_pmax']
    scale_factor = data['scale_factor']
    save_h5ad = data['save_h5ad']
    save_csv = data['save_csv']
    save_polygons = data['save_segmentation_polygons']
    save_norm_img = data['save_normalised_image']
    return csv_file_path, out_folder, background_thresh, save_csv, save_h5ad, skip_failed, prob_thresh, nms_thresh, pmin, pmax, scale_factor, save_polygons, save_norm_img


def main(ConfFilePath):
    csv_file_path, out_folder, background_thresh, save_csv, save_h5ad, skip_failed, prob_thresh, nms_thresh, pmin, pmax, scale_factor, save_polygons, save_norm_img = ReadConfFile(ConfFilePath)
    csv_table = pd.read_csv(csv_file_path)
    list_of_failed_sections = []
    for i in range(csv_table.shape[0]):
        img_path = csv_table['image_path'][i]
        spaceranger_path = csv_table['spaceranger_path'][i]
        sample_name = csv_table['sample_name'][i]
        
        print(sample_name)
        
        if skip_failed:
            try:
                C2V(img_path, spaceranger_path, sample_name, out_folder, background_thresh = background_thresh, save_csv = save_csv, save_h5ad = save_h5ad, prob_thresh = prob_thresh, nms_thresh = nms_thresh, pmin = pmin, pmax = pmax, scale_factor = scale_factor, save_segm_polygons = save_polygons, save_normalised_img = save_norm_img)
            except:
                list_of_failed_sections.append(sample_name)
                
        else:
            C2V(img_path, spaceranger_path, sample_name, out_folder, background_thresh = background_thresh, save_csv = save_csv, save_h5ad = save_h5ad, prob_thresh = prob_thresh, nms_thresh = nms_thresh, pmin = pmin, pmax = pmax, scale_factor = scale_factor, save_segm_polygons = save_polygons, save_normalised_img = save_norm_img)
    
    if len(list_of_failed_sections)>0:
        txt_path = out_folder + '/failed_samples.txt'
        with open(txt_path, 'w') as f:
            for line in list_of_failed_sections:
                f.write(f"{line}\n")
        
if __name__ == "__main__":
    fire.Fire(main)
