# Spectral-Index Based Pulsar Finder 

This repository contains a Python script for crossmatching radio sources and performing spectral analysis using data from the TGSS (TIFR GMRT Sky Survey) and NVSS (NRAO VLA Sky Survey) catalogs. The script calculates spectral indices for matched sources and identifies potential pulsar candidates. 

## Requirements

To run the script, make sure you have the following prerequisites installed:

- Python 
- Required Python packages: pandas, astropy, numpy, matplotlib, scipy, psrqpy, astroquery, bdsf
  
## Download Data

Before running the script, you need to download the following data files:

1. Download TGSS Catalog (TGSSADR1_7sigma_catalog.tsv):
   [Download TGSS Catalog](http://tgssadr.strw.leidenuniv.nl/catalogs/TGSSADR1_7sigma_catalog.tsv)

2. Download Spectral Index Catalog (spidxcat_v1.1b.fits):
   [Details about the catalog](https://tgssadr.strw.leidenuniv.nl/doku.php?id=spidx)
   [Download Spectral Index Catalog](http://tgssadr.strw.leidenuniv.nl/spidx/spidxcat_v1.1b.fits)

## Usage

1. Clone this repository to your local machine:

   ```
   git clone https://github.com/viren1408/SIBPF_new.git
   
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
    please note that in the current version save the image file as "region.fits" eg. "G033.0-5.0.FITS" and save it in "/path/to/data/region" or "dir/region" this has been done for my personal         
    bookeeping  and will be edited out soon , The output files will be saved in the same "dir/region" with proper names as given below. 
  
3. Execute the script using the command:

   ```
   python parser.py --config config.json
   ```

## Functionality

The script provides the following functionalities:

### Crossmatching

- Matches the sources within a specified circular region in the provided FITS image (`image.fits`) obtained using AIPS or CASA with the TGSS and NVSS catalogs.
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


# Acknowledgments

## TGSS Catalog

I acknowledge the use of data from the TGSS (TIFR GMRT Sky Survey) catalog in this project. The TGSS survey has provided invaluable radio astronomy data that has contributed to the analysis. For more information about the TGSS catalog and data access.

## PyBDSF

I acknowledge the developers of PyBDSF (Python Blob Detection and Source Finder), a powerful tool that has been instrumental in our source extraction and analysis. PyBDSF's functionality and capabilities have greatly facilitated our radio astronomy research. For more information about PyBDSF and its features, please refer to their GitHub repository: [PyBDSF GitHub Repository](https://github.com/lofar-astron/PyBDSF).

I acknowledge the use of various other libraries and tools that have contributed to the success of this project.

## Citation
If you use this code or the results obtained from it in your research or work , kindly consider citing relevent publications or resources associated with the TGSS catalog , PyBDSF and NVSS catalog including any other relevent tool that have contributed to your analysis

---

# Note 
Make sure you have properly referenced the TGSS  and the Spidx.fits catalogs and provided accurate paths to the required files in your configuration.The files are present for download
