# Crossmatching and Spectral Analysis

This repository contains a Python script for crossmatching radio sources and performing spectral analysis using data from the TGSS (TIFR GMRT Sky Survey) and NVSS (NRAO VLA Sky Survey) catalogs. The script calculates spectral indices for matched sources and identifies potential pulsar candidates. The code is designed to run in a specific directory or region file created for bookkeeping purposes, which may be removed in the future.

## Requirements

To run the script, make sure you have the following prerequisites installed:

- Python 
- Required Python packages: pandas, astropy, numpy, matplotlib, scipy, psrqpy, astroquery, bdsf
  
## Download Data

Before running the script, you need to download the following data files:

1. Download TGSS Catalog (TGSSADR1_7sigma_catalog.tsv):
   [Download TGSS Catalog](http://tgssadr.strw.leidenuniv.nl/dokuwiki/doku.php?id=tgssadr)

2. Download Spectral Index Catalog (spidxcat_v1.1b.fits):
   [Download Spectral Index Catalog](http://svo2.cab.inta-csic.es/vocats/tools/spidxcat/)

## Usage

1. Clone this repository to your local machine:

   ```
   git clone https://github.com/viren1408/SIBPF_new.git
   cd SIBPF_new
   ```

2. Prepare a configuration JSON file (e.g., `config.json`) with the following structure:

   ```json
   {
     "dir": "/path/to/data",
     "region": "your_region",
     "image": "/path/to/image.fits",
     "spidx": "/path/to/spidx.fits",
     "tgss": "/path/to/tgss.tsv",
     "show_matches": true,
     "get_spectral_index": true,
     "get_candidates": true,
     "get_pulsars": true
   }
   ```
    please note that in the current version save the image file as "region.fits" eg . "G033.0-5.0.FITS" and save it in "/path/to/data/region" or "dir/region" this has been done for my personal         
    bookeeping    and will be edited out soon , The output files will be saved in the same "dir/region" with proper names as given below. 
  
3. Execute the script using the command:

   ```
   python parser.py --config config.json
   ```

## Functionality

The script provides the following functionalities:

### Crossmatching

- Matches sources within a specified circular region in the provided FITS image (`image.fits`) using the TGSS and NVSS catalogs.
- Generates plots showing matched sources on an RA-Dec plot, with marker sizes proportional to resolution.

### Spectral Analysis

- Calculates spectral indices for the matched sources using TGSS and NVSS flux data.
- Produces individual plots of spectral index analysis for each matched source, along with a summary histogram of spectral indices.
- Identifies potential pulsar candidates based on a spectral index threshold (default: -0.9).

### Output

The script generates the following outputs:

- Matched sources information in CSV format (`Matched_sources_your_region.csv`).
- Individual spectral index analysis plots for each source (`spectral_index_plot_i.png`).
- Summary plots of spectral index distribution and SPIDX in the image field (`spectral_index_hist.png` and `spectral_index_hist_spidx.png`).
- Pulsar candidates list in CSV format (`Pulsar_candidates_your_region.csv`).



---

# Note 
Make sure you have properly referenced the TGSS  and the Spidx.fits catalogs and provided accurate paths to the required files in your configuration.The files are present for download
