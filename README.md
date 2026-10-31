# Scalable Spatial Stream Modelling (S3N)

This project implements a scalable spatial stream model (S3N) using Environment Agency (EA) water quality data (NCEA) on the UK Open Rivers Network (ORN) riverline dataset.

This project was built using Google Antigravity and the Gemini Argon 4 model, and reviewed by a human :)

It is based on this paper: Kunke, J.P., Olden, J.D. and McCormick, T.H., 2025. Scalable Spatial Stream Network (S3N) Models. arXiv preprint arXiv:2512.12398. [arXiv preprint arXiv:2512.12398.](https://arxiv.org/abs/2512.12398)

It uses the R `S3N` package, available from: [Github](https://github.com/jpierkunke/S3N) (tested with version 1.0.5), with a custom `C++` script to compute stream distances, weights and nearest neighbours.

River lines/reaches are preprocessed using a custom Python script.

River monitoring data comes from the EA NCEA project, available from: [Data.gov.uk](https://environment.data.gov.uk/natural-capital-ecosystem-assessment/search/21763290-f395-4ca1-bcf6-884c2e7548d4?q=analysis+ready+data#general)

The riverline data comes from the Open Rivers Network, available from: [Open Rivers](https://openrivers.net/)

<img width="2600" height="1450" alt="fig1_england_stream_0117_map" src="https://github.com/user-attachments/assets/12cdd4ab-f316-463c-bfac-f6855b946db7" />

## Usage

Clone this repository to get the scripts. 

Install the S3N package from Github:

```
# Install remotes (if needed)
install.packages("remotes")
# Install S3N
remotes::install_github("jpierkunke/S3N")
```

You'll also need to install the S3N dependencies (this should happen automatically).

Within the cloned project folder, create two additional folders; `data` and `gis`. 

In the `data` folder, download and unzip the `NCEA_River_Surveillance_Network_ARD_Open.zip` into it which you can get from: [Data.gov.uk](https://environment.data.gov.uk/natural-capital-ecosystem-assessment/search/21763290-f395-4ca1-bcf6-884c2e7548d4?q=analysis+ready+data#general), and in the `gis` folder, download the ORN GeoPackage file and unzip. [Open Rivers](https://openrivers.net/download/).

If using the terminal... in the root project folder run:

```
mkdir data gis
```

...in the data folder run:

```
wget -c -O NCEA_River_Surveillance_Network_ARD_Open.zip \ "https://environment.data.gov.uk/file-management-open/data-sets/68444e44-461e-45f8-a7f8-d5296eaf2da9/files/NCEA_River_Surveillance_Network_ARD_Open.zip"
unzip NCEA_River_Surveillance_Network_ARD_Open.zip
```

...in the gis folder run:

```
wget -c -O ORN_v2_GeoPackage.zip \ "https://openrivers.net/download/ORN_v2_GeoPackage.zip"
unzip ORN_v2_GeoPackage.zip
```

The final folder structure should look like this:

```
spatial_stream_modelling
├── data
│   └── NCEA_River_Surveillance_Network_ARD_Open        <- Download from Data.gov.uk
├── docs
├── gis         <- Download the ORN GeoPackage and place in here        
├── readings 
├── scripts
└── src
```

The scripts are run from the command line (rather than from within R/Python).

Refer to the `usage_guide.md` file, found in the `docs` folder for full details on how to use these scripts.

## Requirements

Refer to the `preprocessing_pipeline.md` file for a full list of requirements.

This has been tested on Fedora v44 (Linux), R v4.5+, and Python 3.14.7+.