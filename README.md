## Preparation
Upload whole folder. Create conda environment using environment yaml file:
`conda env create -f environment.yml`
For each dataset you will need to prepare: (i) it's name, (ii) path to the full resolution H&E image and (iii) path to the spaceranger output folder. It is important that image was the same as used in spacerranger analysis.
In case you want to run a batch of the samples you will need to save this information for all samples in csv file with the same column names as in **table_with_paths.csv**
Also, you will need to prepare configuration file (see **conf.yaml**)

## Running
To run single sample:

`python cells2visium.py /path/to/img path/to/spaceranger/folder path/to/output/folder`

you can also add flags:
 - background_thresh [0-255] - intensity in R and B channels, below which pixels will be assigned to tissue (only for "occupancy_tissue" parameter)
 - save_csv [boolean] - whether to save result as csv file or not
 - save_h5ad [boolean] - whether to save anndata from spaceranger output with additional folder (both flags can be *True*)

To run batch of samples (sequentially, no parallelisation is used):

`python cells2visium_batch.py conf.yaml`

Make sure to fill configuration file as well as csv table which is points to

## Notes
- processing of one sample requires quite a lot of memory - I usaully ask for 200-300 Gb
- in case of working with batch of files, I suggest firstly to use flag *skip_failed_samples: True*, then all failed samples will not stop the whole process, and you will have an utput txt with the list of failed samples

