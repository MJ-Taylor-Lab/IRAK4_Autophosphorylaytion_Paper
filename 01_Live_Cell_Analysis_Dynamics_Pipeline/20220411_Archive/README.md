# Image Analysis Pipeline
## Input
The input data goes into ~/Input/parameter_tables. There are seven files, including:
* constants.csv
* dark_frames.csv
* directories.csv
* exclusion_channels.csv
* flat_fields.csv
* images.csv
* ligand.csv

### constants.csv
Numbers which will be constant throughout the analysis
| parameter              | value | comments       |
|------------------------|-------|----------------|
| tiff_compression_level | 5     | out of 10      |
| cell_diameter          | 25    | px, odd number |
| puncta_diameter        | 5     | px, odd number |

### dark_frames.csv
 The dark frame is the camera noise (https://en.wikipedia.org/wiki/Dark-frame_subtraction). This typically is 1000 frames averaged, though 50 frames could do, so long as the standard deviation does not change with more images added. It should be at the same exposure as the images using the same camera as the microscopy images. Thus, one image could be used for multiple channels.
 
 The table contains the image names of the dark frame average.

| image                               | exposure |
|-------------------------------------|----------|
| 20201026 Darkfield 200ms binned.tif | 200 ms   |
| 20201026 Darkfield 50ms binned.tif  | 50 ms    |
| 20201026 Darkfield 100ms binned.tif | 100 ms   |

### directories.csv
| contains    | path                      |
|-------------|---------------------------|
| input       | ~/Input                   |
| processing  | ~/Processing              |
| output      | ~/Output                  |
| dark_frames | ~/dark_frames             |
| flat_fields | ~/flat_fields             |
| ImageJ      | ~/Fiji.app/ImageJ-linux64 |

### exclusion_channels.csv
Channels to exclude from the pipeline analysis.

| value       |
|-------------|
| IL-1        |
| Brightfield |
| WideField   |

### flat_fields.csv
Image of the illumination function
| image                       | channel   | power | exposure |
|-----------------------------|-----------|-------|----------|
| 20200326 IlluminationFx.tif | T Cy5     | 30    | 200      |
| 20200326 IlluminationFx.tif | T GFP     | 30    | 200      |
| 20200326 IlluminationFx.tif | T RFP Cy3 | 30    | 200      |

### images.csv
### ligand.csv
