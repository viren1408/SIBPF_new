import numpy as np
# Set up matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
from astropy.wcs import WCS
from astropy.io import fits
import bdsf
import pandas as pd
# Replace 'your_file.fits' with the actual file path to your FITS file
region = 'G033.0-5.0'
image_file = '/home/viren/work/DATA/'+str(region)+'.FITS'


# Open the original FITS file
hdulist = fits.open(image_file)
hdulist.info()
# Access the image data
image_data = hdulist[0].data[0,0]

print(image_data)


# Invert the pixel values
inverted_image_data = np.max(image_data) - image_data

# Create a new FITS header by copying the original header
new_header = hdulist[0].header.copy()

# Create a new FITS HDU with the inverted image and the copied header
new_hdulist = fits.HDUList([fits.PrimaryHDU(inverted_image_data, header=new_header)])

# Save the new FITS file
new_hdulist.writeto('/home/viren/work/DATA/'+str(region)+'_inverted_image.FITS', overwrite=True)

#Add Pybdsf block to get the artifact locations 
# Specify the file path and name
#file_path = str(dir)+'/'+str(region)+'/'+str(region)+'.sav'
file_path = '/home/viren/work/DATA/'+str(region)+'_inverted_image.sav' #create sav file
# Open the file in write mode
open(file_path, mode='w', newline='')
# Load the input FITS file
input_image = ('/home/viren/work/DATA/'+str(region)+'_inverted_image.FITS')
save_file = '/home/viren/work/DATA/'+str(region)+'_inverted_image.sav' #sav file path 
img = bdsf.process_image(save_file, filename=input_image, quiet=True)
img.write_catalog(format='csv', catalog_type='srl',clobber = True)
img.export_image(img_format = 'fits',img_type = 'gaus_model',clobber = True)
data = pd.read_csv('/home/viren/work/DATA/'+str(region)+'_inverted_image.pybdsm.srl',index_col = None,skiprows=5) #change 1

# Close both FITS files
hdulist.close()
new_hdulist.close()


