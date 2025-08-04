# Author -  Rafael Deliz-Aguirre


# Installation
## Requirements
* Linux cluster
* R 4.0.2

## Process
Connect to the cluster computer:
    
        ssh username@raven.mpcdf.mpg.de
     
Create Python Packages list. Open the terminal text editor, then type `nano` and paste:

        aiohttp
        aiohttp-cors
        aioredis
        appdirs
        async-timeout
        attrs
        blessings
        boto3
        botocore
        cachetools
        certifi
        chardet
        click
        colorama
        colorful
        cycler
        decorator
        distlib
        et-xmlfile
        filelock
        future
        google-api-core
        google-auth
        googleapis-common-protos
        gpustat
        grpcio
        hiredis
        idna
        imageio
        imglyb
        jgo
        jmespath
        JPype1
        jsonschema
        kiwisolver
        matplotlib
        msgpack
        multidict
        nd2reader
        networkx
        numpy
        nvidia-ml-py3
        opencensus
        opencensus-context
        opencv-python
        openpyxl
        packaging
        pandas
        Pillow
        PIMS
        pims-nd2
        pipenv
        prometheus-client
        protobuf
        psutil
        py-spy
        pyasn1
        pyasn1-modules
        pyparsing
        pyrsistent
        python-dateutil
        python-javabridge
        pytz
        PyWavelets
        PyYAML
        ray
        requests
        rsa
        s3transfer
        scikit-image
        scipy
        scyjava
        six
        slicerator
        tifffile
        typing-extensions
        urllib3
        virtualenv
        virtualenv-clone
        xarray
        xlrd
        xmltodict
        yarl
        
Close nano:
* Press [CTRL] + [X] to close
* Press [Y] to save
* Save as `'python_requirements.txt'`

> **Optional**: You may use screen to let process run on the background
>    * Type `screen`
>    * Wait 5s to load
>    * Press [CTRL] + [A] and let go
>    * Press [D]
>    * To resume, type `screen -r`
>> _If there's more than one screen, type `screen -r` . to get the index of screens and then replace 00000 with the index (`screen -r 00000`)_

### Python packages installation

Paste the following in the terminal (**change `username` with your user name**):

        # Install ImageJ
        wget https://downloads.imagej.net/fiji/latest/fiji-linux64.zip
        unzip fiji-linux64.zip
        
        # Install conda
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
        chmod +x Miniconda3-latest-Linux-x86_64.sh
        ./Miniconda3-latest-Linux-x86_64.sh
        export PATH=~/miniconda/bin:$PATH
        source ~/miniconda3/bin/activate
        export PATH="/miniconda3/bin":$PATH

        # Follow Terminal instructions to install conda
        
        # Create conda environment
        conda create -n dynamics_pipeline python=3.8 anaconda
        conda activate dynamics_pipeline
        
        # Install conda packages
        conda install -c conda-forge opencv
        conda install -c conda-forge libglib
        conda install -c conda-forge r-xml
        
        # Install libtiff
        mkdir libtiff
        cd libtiff
        wget https://download.osgeo.org/libtiff/tiff-4.3.0.tar.gz
        tar -xvf tiff-4.3.0.tar.gz
        mkdir install
        cd tiff-4.3.0
        mkdir compile
        cd compile
        # Change username here
        ../configure --prefix=/u/username/libtiff/install
        make
        make install
        # Change username here
        export PKG_CONFIG_PATH=/u/username/libtiff/install/lib/pkgconfig/
        
        # Install Python packages
        python -m pip install --user -r python_requirements.txt

### R packages installation

Load R

        cd ~
        module purge
        module load jdk/8.265 gcc/10 impi/2021.2 fftw-mpi R/4.0.2
        echo 'modules loaded'
        conda activate dynamics_pipeline
        echo 'conda activated'
        R

Paste the following in the R console inside the terminal:

        if("pacman" %in% rownames(installed.packages()) == FALSE)
        {install.packages("pacman")}

When prompted, select which repository you would like to use. Enter the number and then press [Enter]

Then, paste the following:

        pacman::p_load(dplyr, stringr, parallel, tidyr, data.table, ff, dtplyr, compiler, changepoint, R.utils, ijtiff, XML, lemon, ggquiver, ggplot2, ggdark, scales, ggforce, viridis, RcppRoll, metR)

Exit R by typing `q()` and then `N` to not save

### Git installation

Paste the following and rename the last element to yours:

        git config --global user.name "username"
        git config --global user.email email@mpcdf.mpg.de
        ssh-keygen -t rsa -b 4096 -C "name@raven.mpcdf.mpg.de"
        
Press [Enter] to skip some steps and get the ssh key. Then, copy the entire block of string from start to finish, including the _ssh-rsa_ and the ending name

Sign-in to GitHub, https://github.com/login

Paste the key into GitHub, https://github.com/settings/keys

* New SSH key
* Type in the cluster name as **Title** and paste the key under **Key**

Create a directory to save the pipeline scripts

        mkdir dynamics_pipeline
        cd dynamics_pipeline

Paste the following in the terminal to clone the pipeline:

        git init
        git remote add origin git@github.com:MJ-Taylor-Lab/DynamicsPipeline.git
        git remote set-url origin git@github.com:MJ-Taylor-Lab/DynamicsPipeline.git
        git fetch --all
        git reset --hard origin/master
        git pull origin master

---

# Image Analysis Pipeline
## Input
The input data goes into ~/new_pipeline/pending_processing/batch_date/Input/parameter_tables. There are five files, including:
* constants.csv
* dark_frames.csv
* directories.csv
* exclusion_channels.csv
* images.csv

### constants.csv
Numbers which will be constant throughout the analysis
| parameter              | value | comments       |
|------------------------|-------|----------------|
| tiff_compression_level | 5     | out of 10      |
| cell_diameter          | 25    | px, odd number |
| puncta_diameter        | 5     | px, odd number |

### dark_frames.csv
 The dark frame is the camera noise (https://en.wikipedia.org/wiki/Dark-frame_subtraction). This typically is 1000 frames averaged, though 50 frames could do, so long as the standard deviation does not change with more images added. It should be at the same exposure as the images using the same camera as the microscopy images. Thus, one image could be used for multiple channels.
 
 The table contains the image names of the dark frame average and their exposures **with units**.

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

### images.csv
| image | cohort | segment_with | ligand | ligand_density | trackmate_max_link_distance | trackmate_threshold | trackmate_frame_gap | T Cy5 protein_name | T GFP protein_name | T RFP protein_name | WideField protein_name |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 20211218 0p8nM 069-1R_TRAF6_MyD88   Grid_1um_11mol 001.nd2 | MyD88 TRAF6 1um_grid | MyD88 | 0.8 nM IL-1 | 11 | 5 | 1.5 | 5 | IL-1 | MyD88 | TRAF6 | Brightfield |
| 20211218 GFP calibration_10pct_60ms   005.nd2 | Calibrations | GFP |  |  | 2.5 | 1.5 | 5 | IL-1 | GFP | mScarlet | Brightfield |
| 20211218 mScarlet calibration_10pct_60ms   001.nd2 | Calibrations | mScarlet |  |  | 2.5 | 1.5 | 5 | IL-1 | GFP | mScarlet | Brightfield |

## Run

Connect to the cluster computer:
    
        ssh username@raven.mpcdf.mpg.de
> If you need the latest scripts, paste in the Terminal:
>
>       cd dynamics_pipeline
>       git pull origin master
>

### SLURM Instructions

Create SLURM instructions file. Open the terminal text editor, then type `nano` and paste:

        #!/bin/bash -l
        
        #SBATCH -o ./job.out.%j
        #SBATCH -e ./job.err.%j
        #SBATCH -D ./
        #SBATCH -J 20211218
        #SBATCH --mail-type=ALL
        #SBATCH --mail-user=email@mpcdf.mpg.de
        #SBATCH --nodes=1
        #SBATCH --ntasks-per-node=1
        #SBATCH --cpus-per-task=72
        #SBATCH --time=24:00:00
        
        # Load all needed packages
        module purge
        module load jdk/8.265 gcc/10 impi/2021.2 fftw-mpi R/4.0.2
        echo 'modules loaded'
        conda activate dynamics_pipeline
        echo 'conda activated'
        
        # Specify parameters
        ## Path of parameters table
        ## Change username to your cluster user name
        path=$'/raven/u/username/new_pipeline/pending_processing/batch_date/Input/parameter_tables'
        
        ## Scripts folder
        ## Change username to your cluster user name
        cd /raven/u/username/dynamics_pipeline
        
        ## Cores for parallel processing in R
        export OMP_NUM_THREDS=144
        
        # Run scripts
        ## Python scripts
        python mission_control.py $path 12
        
        ## Run R Scripts
        Rscript --vanilla --verbose r_scripts/extract_intensity.R $path
        Rscript --vanilla --verbose r_scripts/colocalization.R $path
        Rscript --vanilla --verbose r_scripts/compile_tables.R $path
        #Rscript --vanilla --verbose r_scripts/compress_everything.R $path
        
        sleep 10

Press [CTRL] + [X] to close, then [Y] to save and type `submit_node.sh` to save it under that name

Paste `sbatch submit_node` to submit to SLURM



---


# Pipeline Details
## 1) nd2 to TIFF

Nikon microscopy images of different channels and meta data are captured and stored in nd2 format. The first step converts each channel of image captured to a tiff file. Simultaneously meta data is extracted (Eg: laser wavelength used, exposure time, laser power, real time of each frame etc.) 

## 2) Dark-current noise subtraction

Camera noise was estimated by averaging several pictures with the
camera shutter closed. The script identifies the most appropriate dark frame based on the exposure time and subtracts this dark frame image from the one to be analysed.


## 3) Images for tracking protein paths

A edge preserving technique called median filter subtraction is used to remove blur and create sharp puncta. A user supplied median value (25 pixels) is used by each pixel to create a box the specified size centered around itself. The median value in this box is calculated and assigned as a new value to the pixel creating a blurred imaged. The median blurred image is then subtracted from the dark frame subtracted image. This reduces the cytoplasmic noise accquired during imaging. This step is done for images of each channel. In order to counteract random pixel fluctations two filters were used 1) Then a time-independent median filter (3 pixels) was used to sharpen the puncta similar to the median filter above 2) A time-depndent process using moving averages was used.  The moving average filter is achived by grouping a set of frames (3 frames) together and finding the average value for each pixel. For multiple channels of images the images processed up to this point are added up. The image is now suitable for tracking puncta using trackmate in the later steps and is called tracking reference image. 


## 4) Images for extracting protein intensities

Similar to the previous step, the original image is dark-frame subtracted followed by subtraction of the median blurred image (25 pixels). This image is used for intenisty extraction and referred to as intensities reference image. This step is done for images of each channel.


## 5) Cell Segemntation

Marker calculation - A maximum instensity projection image is calculalted by selecting maximum brightness value at each pixel position using the intensities reference image. Contrast is adjusted by scaling intensities. Dimmest pixel value being set to 0 and highest value set to max value based on the bit rate of the image. Using this 2 images were created. One with a gaussian blur (15 pixels) and other one with a median blur (5 pixels). The median blurred image is subtracted from the gaussian blur to generate an image with the brightest spots. This image was then log transformed to generate a more uniform intensity. Puncta are made clearer by applying a median blur (3 pixels). Lastly, the image was thresholded using connected components labeling - a method for identifying blobs by looking at neighbouring pixels. Lastly image is made binary.
Mask Caluculation - A mean instensity projection image is calculalted by selecting average brightness value at each pixel position using the intensities reference image. Contrast is adjusted by scaling intensities as above. Similarly gaussian blur (15 pixels) and median blur (5 pixels) pixels were created. However, the gaussian blur image dialtes when filter applied on a mean instensity projection image makiing dim regions brighter and vice versa. The median blurred image is subtracted from the gaussian blur image. A median blur (5 pixels) is applied. The image is log transformed and made binary.
Next an algorithm called marker controlled watershed segmentation was used. The mask and marker images are used in this algorithm to determine the location of each cell. Each cell is then cropped and stored as a separate tiff file. Each cell folder contains an image of the cropped tracking reference image and intensity reference image.


## 6) Tracking Puncta

TrackMate is used to identify spots and track their position over time. Trackmate uses a Laplacian of Gausian (LoG) filter to identify spots in the tracking reference image of a single cell.  An alogrithm called Linear Assignment Problem (LAP) is used to calculated the trajectories of each puncta by doing a cost analysis. The coordinates are then stored as an xml file for each cell. Note: currently tracking reference is a combination of images from multiple channels because expect spots to be colocalised.  


## 7) Intenisty extraction

The coordinates containing xml file for each cell determined from the tracking reference image are extracted. Missing frame information for each track was approximated using the mean distance between the points prceeding and occuring after the missing event. With all the cooridnates available, around each spot a circular mask coordinates are calculated. This mask is then applied to the intensity reference image for the cell and intensity value is extracted for each spot at time point. The value extracted for each puncta for each time point is the sum of all pixels in the mask. This is done for all channels. Simultaneously meta data extracted in step 1 is combined with the intensity data and coordinate data to generate a analysis table. 

## 8) Intenisty Normalisation

Indiviual flurophores which are imaged as controls to estimate intensity of single flurophores undergo all steps of analysis above except step 5. The median intenisty is then extracted and for each image and set as the value of a single flurophore. Using parameters such as Power, Exposure, Angle (in radians), Direction, and Date the flurophore images are matched to the cell images that fit the these parameters the best. The intensity values of each spot of each channel at each time point is divded by the appropriate median flurophore intensity value to estimate the number of molecules within a puncta. 


## 9) Parameter extraction

Using the coordinates, time and intensity value of each puncta in each channel at each time for every track various parameters are caluculated like Max Normalised Intenisty, Lifetime of track, Start to Max Intensity etc. All of this information is saved as a csv.gz file for further analysis. 
