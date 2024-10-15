## Preparation
Upload whole folder. Create conda environment using **environment.yaml** file:

`conda env create -f environment.yml`

Alternatively you an use **requirements.txt** to create python environment.
For each dataset you will need to prepare: (i) it's name, (ii) path to the full resolution H&E image and (iii) path to the spaceranger output folder. It is important that image was the same as used in spacerranger analysis.
In case you want to run a batch of the samples you will need to save this information for all samples in csv file with the same column names as in **table_with_paths.csv**
Also, you will need to prepare configuration file (see **conf.yaml**)

## Running
Activate conda environment

`conda activate visium_segm`

To run single sample:

`python cells2visium.py /path/to/img path/to/spaceranger/folder path/to/output/folder`

you can also add flags:
 - *background_thresh* [0-255] - intensity in R and B channels, below which pixels will be assigned to tissue (only for "occupancy_tissue" parameter)
 - *save_csv* [boolean] - whether to save result as csv file or not
 - *save_h5ad* [boolean] - whether to save anndata from spaceranger output with additional folder (both flags can be *True*)
 - *prob_thresh* - probability threshold defined for segmentation (segmented cells with detection probability less than this value will be filtered out)
 - *nms_thresh* - overlapping threshold. Higher value will lead for larger area fraction of neighbouring segmented cells being overlapped
 - *pmin* - min boundary of percentile-based image normalisation *pmin*=[0,1] (everything below this value will be 0)
 - *pmax* - max boundary of percentile-based image normalisation *pmax*=[0,1]>*pmin* (everything above this value will be 1)
 - *scale_factor* - a multiplier to be used for spot position and sizes. SHould be used if spot positions were defined in image with different resolution to the image you use.
 - *save_segm_polygons* - save segmented polygons in json format. the json file will be saved as dictionary with fields `['coord']` (actual corrdinates of polygons), `['points']` (center positions of each cell/polygon), and `['prob']` (detection probability for each cell)
To run batch of samples (sequentially, no parallelisation is used):

`python cells2visium_batch.py conf.yaml`

Make sure to fill configuration file as well as csv table which is points to

## Output parameters
- **y**, **x** - positions of visium spot inside the full-resolution image
- **n_cell** - number of cell at the visium spot (even if cell only partially occupies visium spot it will be counted)
- **occupancy** - percentage of visium spot occupied by segmented cells (their nuclei)
- **R**, **G**, **B** - sum of all pixels inside cellular nuclei in each of the channels respectively (red, green, blue)
- **area_cell_px** - average area of one cell in pixels. It is set to 0, if no cells were detected for the visium spot
- **seg_prob** - average probability of segmentation (kind of QC). It is set to 0, if no cells were detected
- **elongation** - average nuclei elongation (each cell nuclei is fitted to ellipse, and aspect ratio is computed). elongation >= 1 (elongation = 1 for perfect circle). It is set to 0, if no cells were detected
- **occupancy_tissue** - percentage of visium spot occupied by the tissue. Tissue is defined as an area where red and blue intensity are less than background_thresh

## Notes
- processing of one sample requires quite a lot of memory - I usaully ask for 300-400 Gb
- in case of working with batch of files, I suggest firstly to use flag *skip_failed_samples: True*, then all failed samples will not stop the whole process, and you will have an utput txt with the list of failed samples
  

