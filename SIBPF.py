import pandas as pd
import warnings
warnings.simplefilter(action='ignore')
from astroquery.vizier import Vizier
import astropy.units as u
import astropy.coordinates as coord
import numpy as np
import os
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
from astropy.coordinates import SkyCoord
from astropy import units as u
from psrqpy import QueryATNF
from astropy.table import Table

def Crossmatching(dir,region,file_path_image,file_path_spidx,file_path_TGSS,show_matches =True,get_spectral_index =True,get_candidates=True,get_pulsars =True ):


  # Open the FITS file and get the header

  Image_file = fits.open(file_path_image)
  header = Image_file[0].header
  file_name, file_ext = os.path.splitext(file_path_image)

  # Open the FITS image file and get the header
  with fits.open(file_path_image) as hdul:
    header = hdul[0].header

  # pixel scale from the header
  cdelt1 = header['CDELT1']
  cdelt2 = header['CDELT2']
  pixel_scale = np.sqrt(cdelt1**2 + cdelt2**2)

  # size of the image in pixels
  nx = header['NAXIS1']
  ny = header['NAXIS2']

  # Calculate the radius of the image in degrees
  radius_deg = (np.sqrt((nx/2)**2 + (ny/2)**2) * pixel_scale)/2 # to be understood why we are getting double value

  # Get the reference RA and Dec from the header
  ref_ra = header['CRVAL1']
  ref_dec = header['CRVAL2']

  #Creating Source Files

  #GMRT sources in the processed image (sources extracted using BDSF)

  import bdsf

  # Specify the file path and name
  file_path = str(dir)+'/'+str(region)+'/'+str(region)+'.sav'

  # Open the file in write mode
  open(file_path, mode='w', newline='')
  # Load the input FITS file
  input_image = (file_path_image)
  save_file = str(dir)+'/'+str(region)+'/'+str(region)+'.sav' #create .sav file
  img = bdsf.process_image(save_file, filename=input_image, quiet=True)
  img.write_catalog(format='csv', catalog_type='srl',clobber = True)
  img.export_image(img_format = 'fits',img_type = 'gaus_model',clobber = True)
  data = pd.read_csv(str(dir)+'/'+str(region)+'/'+str(region)+'.pybdsm.srl',index_col = None,skiprows=5) #change 1
  observed_sources_main = pd.DataFrame(data)
  observed_sources = pd.DataFrame(data)
  observed_sources.set_index("# Source_id", inplace = True)


  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
  #extract the header information from the Pybdsf fits file.
  Image_file = fits.open(str(dir)+'/'+str(region)+'/'+str(region)+'.pybdsm_gaus_model.fits') #change 2
  header = Image_file[0].header
  bmaj = header['BMAJ']
  bmin = header['BMIN']
  theta_obs = 0.5*(header['BMAJ']+header['BMIN'])




  # Initialize lists for matched Spidx and E_Spidx values

  dat = Table.read(file_path_spidx, format='fits')
  spix_main = dat.to_pandas()
  spix = spix_main[(spix_main['S_code'] == b'S') | (spix_main['S_code'] == b'M')] #Get the matched and double matched sources from SPIDX as given in the text


  # Create a boolean mask for sources within the circular region
  distances = np.sqrt((spix['RA'] - ref_ra)**2 + (spix['DEC'] - (ref_dec))**2)
  mask = distances <= radius_deg

  # Extract the sources within the circular region
  spix  = spix[mask]

  # Reset the index of the sources DataFrame
  spix  = spix.reset_index(drop=True)



  # Loop through each source in spix
  for i in range(len(spix)):
      ra_1 = spix['RA'][i]
      dec_1 = spix['DEC'][i]

      # Initialize variables to store the closest matched source information
      min_theta = np.inf
      matched_source_index = None

      # Loop through each matched source in Matched_sources
      for j in range(len(observed_sources)):
          ra_2 = observed_sources[' RA'][j]
          dec_2 = observed_sources[' DEC'][j]

          # Calculate the angular separation between the two sources
          theta = np.degrees(np.arccos(np.sin(np.radians(dec_1)) * np.sin(np.radians(dec_2)) + np.cos(np.radians(dec_1)) * np.cos(np.radians(dec_2)) * np.cos(np.radians(ra_1 - ra_2))))

          # Check if the angular separation is within the constraint
          if abs(theta) <= theta_obs:
              # Update the closest matched source information if this is a better match
              if theta < min_theta:
                  matched_source_index = j
                  min_theta = theta

      # If a match was found, update the 'Spidx' and 'E_Spidx' columns in matched_sources
      if matched_source_index is not None:
          matched_spidx = spix['Spidx'][i]
          matched_e_spidx = spix['E_Spidx'][i]
          matched_sep = min_theta
          ra_spidx = spix['RA'][i]
          dec_spidx = spix['DEC'][i]
          flux_NVSS = spix['Total_flux_NVSS'][i]
          e_flux_NVSS = spix['E_Total_flux_NVSS'][i]
          flux_TGSS = spix['Total_flux_TGSS'][i]
          e_flux_TGSS = spix['E_Total_flux_TGSS'][i]
          observed_sources.at[matched_source_index, 'RA_SPIDX'] = ra_spidx
          observed_sources.at[matched_source_index, 'DEC_SPIDX'] = dec_spidx
          observed_sources.at[matched_source_index, 'sep_arcsecods'] = matched_sep*3600
          observed_sources.at[matched_source_index, 'Total_flux_NVSS'] = flux_NVSS
          observed_sources.at[matched_source_index, 'E_Total_flux_NVSS'] = e_flux_NVSS
          observed_sources.at[matched_source_index, 'Total_flux_TGSS'] = flux_TGSS
          observed_sources.at[matched_source_index, 'E_Total_flux_TGSS'] = e_flux_TGSS
          observed_sources.at[matched_source_index, 'Spidx'] = matched_spidx
          observed_sources.at[matched_source_index, 'E_Spidx'] = matched_e_spidx



  # If you want to remove rows with no matches, you can do the following:
  observed_sources.dropna(subset=['Spidx', 'E_Spidx'], inplace=True)

  # Reset the index of the matched_sources DataFrame
  observed_sources.reset_index(drop=True, inplace=True)
  filename = str(dir)+'/'+str(region)+'/Matched_sources_'+str(region)+'.csv'
  matched_sources = observed_sources.copy()
  matched_sources.to_csv(filename, index=False)

  if show_matches == True :
    Vizier.ROW_LIMIT = -1
    result = Vizier.query_region(coord.SkyCoord(ra=ref_ra, dec=ref_dec,unit=(u.deg, u.deg),frame='icrs'),radius= radius_deg*u.degree,catalog=["NVSS"])
    table = result[0]
    NVSS = table.to_pandas()
    # Convert RA and Dec from sexagesimal units to decimal degrees

    coords = SkyCoord(NVSS['RAJ2000'], NVSS['DEJ2000'], unit=(u.hourangle, u.deg))
    NVSS['RAJ2000'] = coords.ra.degree
    NVSS['DEJ2000'] = coords.dec.degree

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

    #TGSS sources in the region of the target image:
    main_TGSS = pd.read_csv(file_path_TGSS,sep='\t')


    # Create a boolean mask for sources within the circular region
    distances = np.sqrt((main_TGSS['RA'] - ref_ra)**2 + (main_TGSS['DEC'] - (ref_dec))**2)
    mask = distances <= radius_deg

    # Extract the sources within the circular region
    TGSS = main_TGSS[mask]

    # Reset the index of the sources DataFrame
    TGSS = TGSS.reset_index(drop=True)

   # Define the resolutions for each telescope
    resolution_nvss = 45
    resolution_tgss = 25
    resolution_observed = theta_obs*3600

# Calculate the marker sizes based on the resolution (you can adjust the scale factor as desired) in arcseconds
    marker_size_nvss = 5*resolution_nvss
    marker_size_tgss = 5*resolution_tgss
    marker_size_gmrt = 5*resolution_observed

# Plot the matched sources on a RA-Dec plot, with marker sizes proportional to resolution
    fig, ax = plt.subplots(figsize=(15,8))
    ax.scatter(NVSS['RAJ2000'], NVSS['DEJ2000'], color="black", label='NVSS', s=marker_size_nvss,alpha=0.2)
    ax.scatter(TGSS['RA'], TGSS['DEC'], color="black", label='TGSS', s=marker_size_tgss,alpha=0.3)
    ax.scatter(observed_sources_main[' RA'], observed_sources_main[' DEC'], color="black", label='Observed_sources', s=marker_size_gmrt,alpha=0.4)

    #ax.scatter(spix['RA'], spix['DEC'], color="blue",label='SPIX',s=marker_size_gmrt,alpha=1)
    ax.scatter(matched_sources[' RA'], matched_sources[' DEC'], color="none", edgecolor="red", label='Matched_sources',s=marker_size_gmrt)
    plt.gca().invert_xaxis() #orientation according to ds9 and other sofwares
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('Dec (deg)')
    ax.set_title('Matched sources_' + str(region))
    ax.legend()
    ax.grid(True)

    file_loc = str(dir) + '/' + str(region) + '/matched_sources.png'
    plt.savefig(file_loc)
    plt.show()

  if get_spectral_index == True:
    frequencies = np.array([1.50*1e08,3.2270*1e08,1.4*1e09]) #in Hz
    spectral_index = matched_sources.copy()
    def f(x):
      return slope*x + intercept

    x = np.linspace(np.log10(1.0*1e08),np.log10(frequencies[-1]*2), 100)

    #####Initialization of Loop#####

    for i in range(len(spectral_index)):
      #title = matched_sources['TGSS_id'][i]
      Ra = matched_sources[' RA'][i]
      Dec = matched_sources[' DEC'][i]
      log_of_flux = np.array([np.log10(matched_sources['Total_flux_TGSS'][i]),np.log10(matched_sources[' Total_flux'][i]),np.log10(matched_sources['Total_flux_NVSS'][i])])
      flux_1 = np.array([matched_sources['Total_flux_TGSS'][i],matched_sources[' Total_flux'][i],matched_sources['Total_flux_NVSS'][i]])
      e_flux = np.array([matched_sources['E_Total_flux_TGSS'][i],matched_sources[' E_Total_flux'][i],matched_sources['E_Total_flux_NVSS'][i]])
      spidx = matched_sources['Spidx'][i]
      e_spidx = matched_sources['E_Spidx'][i]
      slope, intercept, r_value, p_value, std_err = linregress(np.log10(frequencies),log_of_flux)
      error = std_err / abs(slope)
      # Update values in the DataFrame using .at[]
      spectral_index.at[i, 'Sr.no'] = i
      spectral_index.at[i, 'Spectral_index'] = slope
      spectral_index.at[i, 'E_Spectral_index'] = error
     


      #print("The spectral index",slope)

      fig, ax = plt.subplots()
      plt.grid()
      plt.scatter(np.log10(frequencies),log_of_flux,marker='.')
      plt.errorbar(np.log10(frequencies),log_of_flux,yerr = 0.434*(e_flux/flux_1),alpha=0.7,ecolor='black',capsize=2,ls = 'none')
      plt.title('RA,DEC :'+ str(Ra.round(4))+','+ str(Dec.round(4)))
      plt.xlabel('Log of Frequency')
      plt.ylabel('Log of Flux')
      plt.tight_layout


      ax.plot(x,f(x),'--g',label = "α = " + str(slope.round(3)) + "±" + str(error.round(3)) + ' ' + "α_spidx" + str(spidx.round(3)) + "±" + str(e_spidx.round(3)))
      plt.legend(loc='best')

      filename = str(dir)+'/'+str(region)+'/spectral_index_plot_{}.png'.format(i)

      plt.savefig(filename)
      #plt.show()

      #Close plot
      plt.close()

 
    filename =  str(dir)+'/'+str(region)+'/spectral_index_'+str(region)+'.csv'
    spectral_index.to_csv(filename, index=False)

  if get_candidates == True:
    pulsar_candidates = spectral_index[spectral_index['Spectral_index'] < -0.9]
    filename = str(dir)+'/'+str(region)+'/Pulsar_candidates_'+str(region)+'.csv'
    pulsar_candidates.to_csv(filename, index=False)
    fig, ax = plt.subplots()
    #plotting the candidates
    plt.scatter(spectral_index['Sr.no'],spectral_index['Spectral_index'],marker = '.')
    plt.scatter(pulsar_candidates['Sr.no'],pulsar_candidates['Spectral_index'],s=50, linewidth=0.5, alpha=0.7,color='r',label = 'Selected_candidates')
    plt.errorbar(spectral_index['Sr.no'],spectral_index['Spectral_index'],yerr =spectral_index['E_Spectral_index'],alpha=0.7,ecolor='black',capsize=2,ls = 'none')
    plt.axhline(y=-0.9, color='b', linestyle='--', label='Spectral index = -0.9')
    plt.axhspan(ymin=spectral_index['Spectral_index'].min(), ymax=-0.9, alpha=0.3, color='g', label='Spectral index < -0.9')
    plt.legend()
    plt.grid()

    plt.title('Pulsar_candidates')
    plt.ylabel('Spectral Index')
    plt.xlabel('Sr.no')
    loc = str(dir)+'/'+str(region)+'/Pulsar_candidates.png'
    plt.savefig(loc)
    plt.show()

    #Get spectral_index histogram
    ax_1 = spectral_index.hist(column='Spectral_index', bins=25, grid=True, figsize=(12,8), color='#86bf91')
    # add a vertical line at spectral index = -0.9
    plt.axvline(x=-0.9, color='b', linestyle='--', label='Spectral index = -0.9')

    # add a translucent shade to the region with spectral index steeper than -0.9
    plt.axvspan(xmin=spectral_index['Spectral_index'].min(), xmax=-0.9, alpha=0.3, color='r', label='Spectral index < -0.9')

    # set the plot title and labels
    plt.title('Spectral Index Distribution')
    plt.xlabel('Spectral Index')
    plt.ylabel('Counts')

    # add a legend
    plt.legend()

    loc = str(dir)+'/'+str(region)+'/spectral_index_hist.png'
    plt.savefig(loc)
    plt.show


    ax_2 = spectral_index.hist(column='Spidx', bins=25, grid=True, figsize=(12,8), color='#86bf91')
    # add a vertical line at spectral index = -0.9
    plt.axvline(x=-0.9, color='b', linestyle='--', label='Spectral index = -0.9')

    # add a translucent shade to the region with spectral index steeper than -0.9
    plt.axvspan(xmin=spectral_index['Spidx'].min(), xmax=-0.9, alpha=0.3, color='r', label='Spectral index < -0.9')

    # set the plot title and labels
    plt.title('Spectral Index Distribution of SPIDX')
    plt.xlabel('Spectral Index')
    plt.ylabel('Counts')

    # add a legend
    plt.legend()

    loc = str(dir)+'/'+str(region)+'/spectral_index_hist_spidx.png'
    plt.savefig(loc)
    plt.show

  if get_pulsars == True :
    #Pulsars = pd.read_csv(str(dir)+'/'+str(region)+'/Pulsar_candidates_'+str(region)+'.csv')
    Image_file = fits.open(file_path_image)
    header = Image_file[0].header
    file_name, file_ext = os.path.splitext(file_path_image)

    # Open the FITS image file and get the header
    with fits.open(file_path_image) as hdul:
      header = hdul[0].header

    # pixel scale from the header
    cdelt1 = header['CDELT1']
    cdelt2 = header['CDELT2']
    pixel_scale = np.sqrt(cdelt1**2 + cdelt2**2)

    # size of the image in pixels
    nx = header['NAXIS1']
    ny = header['NAXIS2']

    # Calculate the radius of the image in degrees
    radius_deg = (np.sqrt((nx/2)**2 + (ny/2)**2) * pixel_scale)/2 # to be understood why we are getting double value

    # Get the reference RA and Dec from the header
    ref_ra = header['CRVAL1']
    ref_dec = header['CRVAL2']



    # convert the input values to sexagesimal format
    c = SkyCoord(ra=ref_ra*u.deg, dec=ref_dec*u.deg, frame='icrs')
    ra_sex = c.ra.to_string(unit=u.hour, sep=':')
    dec_sex = c.dec.to_string(unit=u.degree, sep=':')

    # set the circular boundary centre (RAJ then DECJ) and radius in sexagesimal format
    circular_boundary = [ra_sex, dec_sex, radius_deg]

    # query the ATNF database
    query = QueryATNF(params=['JNAME', 'RAJ', 'DECJ', 'S400', 'S1400'], circular_boundary=circular_boundary)

    # create NumPy arrays for each column in the query
    jname = np.array(query['JNAME'])
    raj = np.array(query['RAJ'])
    decj = np.array(query['DECJ'])
    s400 = np.array(query['S400'])
    s1400 = np.array(query['S1400'])

    # calculate spectral index using NumPy operations
    spec_idx = np.where((np.isnan(s400)) | (np.isnan(s1400)), np.nan, (np.log10(s400) - np.log10(s1400)) / (np.log10(400) - np.log10(1400)))

    # combine the arrays into a Pandas DataFrame
    pulsar_spectral_index = pd.DataFrame({
      'JNAME': jname,
      'RAJ': raj,
      'DECJ': decj,
      'S400' : s400,
      's1400' : s1400,
      'spectral_index': spec_idx
    })
    Pulsars = pd.DataFrame(pulsar_spectral_index)
    filename =  str(dir)+'/'+str(region)+'/Pulsars_'+str(region)+'.csv'
    Pulsars.to_csv(filename, index=False)
    print(Pulsars)

  return spectral_index, pulsar_candidates, matched_sources
