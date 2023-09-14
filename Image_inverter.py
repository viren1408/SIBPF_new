import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
from astropy.wcs import WCS
from astropy.io import fits
import bdsf
import pandas as pd
# Replace 'your_file.fits' with the actual file path to your FITS file
region = 'G024.0-5.6'
dir = '/home/viren/work/DATA/inverter'
image_file = str(dir)+'/'+str(region)+'.FITS'


# Open the original FITS file
hdulist = fits.open(image_file)
hdulist.info()
# Access the image data
image_data = hdulist[0].data[0,0]

#print(image_data)


# Invert the pixel values
inverted_image_data =  -1.0*(image_data)

# Create a new FITS header by copying the original header
new_header = hdulist[0].header.copy()

# Create a new FITS HDU with the inverted image and the copied header
new_hdulist = fits.HDUList([fits.PrimaryHDU(inverted_image_data, header=new_header)])

# Save the new FITS file
new_hdulist.writeto(str(dir)+'/'+str(region)+'_inverted_image.FITS', overwrite=True)

#Pybdsf block to get the artifact locations 
# Specify the file path and name
file_path = str(dir)+'/'+str(region)+'_inverted_image.sav' #create sav file inverted
file_path_1 = str(dir)+'/'+str(region)+'.sav' #create sav file
# Open the file in write mode
open(file_path, mode='w', newline='') #inverted
open(file_path_1, mode='w', newline='')
# Load the input FITS file
input_image = (str(dir)+'/'+str(region)+'_inverted_image.FITS') #inverted
input_image_1 = (str(dir)+'/'+str(region)+'.FITS') #original

save_file = str(dir)+'/'+str(region)+'_inverted_image.sav' #sav file path 
save_file_1 = str(dir)+'/'+str(region)+'.sav' #sav file path original

img = bdsf.process_image(save_file, filename=input_image, quiet=True) #inverted
img_1 = bdsf.process_image(save_file_1, filename=input_image_1, quiet=True) #original

img.write_catalog(format='csv', catalog_type='gaul',clobber = True) #inverted
img_1.write_catalog(format='csv', catalog_type='gaul',clobber = True) #original

img.export_image(img_format = 'fits',img_type = 'gaus_model',clobber = True) #inverted
img_1.export_image(img_format = 'fits',img_type = 'gaus_model',clobber = True) #original

data_inverted= pd.read_csv(str(dir)+'/' +str(region)+'_inverted_image.pybdsm.gaul',index_col = None,skiprows=5) # inverted
data = pd.read_csv(str(dir)+'/' +str(region)+'.pybdsm.gaul',index_col = None,skiprows=5) 



print(data_inverted)

#Fraction of Sources 
Inverted_fraction = len(data_inverted)/len(data)

print(f"False detection fraction: {Inverted_fraction}")


# Close both FITS files
hdulist.close()
new_hdulist.close()


