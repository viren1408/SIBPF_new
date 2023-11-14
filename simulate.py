from astropy.io import fits
import numpy as np

# Load the original FITS file
original_fits_file = 'G033.0-5.0.FITS'
original_hdul = fits.open(original_fits_file)
original_data = original_hdul[0].data
original_header = original_hdul[0].header
original_hdul.close()

# Create and save 5 new FITS files with varied pixel values
for i in range(1, 6):
    # Generate simulated data with variations
    simulated_data = i*(original_data)
    
    # Create a new HDU (Header Data Unit)
    hdu = fits.PrimaryHDU(simulated_data, header=original_header)
    # Save the new FITS file with a different name for each iteration
    new_fits_file = f'G033.0-5.0_mod{i}.FITS'
    hdu.writeto(new_fits_file, overwrite=True)
    
    print(f'Saved simulated FITS file: {new_fits_file}')
