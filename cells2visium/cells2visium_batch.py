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
    skip_failed = data['skip_failed_samples']
    save_h5ad = data['save_h5ad']
    save_csv = data['save_csv']
    return csv_file_path, out_folder, background_thresh, save_csv, save_h5ad, skip_failed


def main(ConfFilePath):
    csv_file_path, out_folder, background_thresh, save_csv, save_h5ad, skip_failed = ReadConfFile(ConfFilePath)
    csv_table = pd.read_csv(csv_file_path)
    list_of_failed_sections = []
    for i in range(csv_table.shape[0]):
        img_path = csv_table['image_path'][i]
        spaceranger_path = csv_table['spaceranger_path'][i]
        sample_name = csv_table['sample_name'][i]
        
        print(sample_name)
        
        if skip_failed:
            try:
                C2V(img_path, spaceranger_path, sample_name, out_folder, background_thresh = background_thresh, save_csv = save_csv, save_h5ad = save_h5ad)
            except:
                list_of_failed_sections.append(sample_name)
                
        else:
            C2V(img_path, spaceranger_path, sample_name, out_folder, background_thresh = background_thresh, save_csv = save_csv, save_h5ad = save_h5ad)
    
    if len(list_of_failed_sections)>0:
        txt_path = out_folder + '/failed_samples.txt'
        with open(txt_path, 'w') as f:
            for line in list_of_failed_sections:
                f.write(f"{line}\n")
        
if __name__ == "__main__":
    fire.Fire(main)
