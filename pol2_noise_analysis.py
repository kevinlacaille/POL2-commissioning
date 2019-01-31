import star
import os
import logging
from datetime import datetime
import numpy as np
import astropy
from astropy.table import Table, Column
from  scipy.signal import convolve2d
from astropy.io import fits

"""Script to run noise analysis on various POL-2 observations.

Requires a list of files representing maps (full path).

All maps must have pixel position 0,0 at the centre of the map, this
will not work otherwise.

All maps must have sensible FITS header values for:
 - TELAPSE (or TOTEXP)
 - EFFBOLO
 - TAUST and TAUEN

(Note that if you use wcsmosaic, makemos or pol2stack to combine
images *ALL THESE HEADERS WILL BE WRONG*. You will have to manually
fix this up with 'magic' P2 prefixed fits headers -- see code to see
what is required.)


"""

### values to set (move to command line arguments?)
# list of files
# filelist = '/export/data/sgraves/POL2/analysis_scratch/all_3c_and_uranus_maps.lis'
filelist = '/export/data/sgraves/POL2/analysis_scratch/subscans_20160112_00056.lis'



# Central radius to block out due to having a source.
source_radius=[0,15,30] # arcseconds

# Central circular region to analyse noise within
noise_radius = [2*60.0,3*60.0, 4*60.0 ] # arcseconds

# name of output catalog (without extension)
# outputcat = 'single_observations_multiple_radii'
outputcat = 'subscans_20160112_0056.lis'

# Logging.
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# Create ard mask for source and analysis area.
def make_ard_mask(inputfile, noise_radius, source_radius,
                  ardname = 'ardmask.ARD'):
    """Create an annuluar ardmask blocking out a source.

    noise_radius and source_radius in arcseconds.

    assumes FPIXSCALE first and second values from ndftrace give
    pixelscale in arcseconds.

    returns name of output ardmask file.

    """

    # First get pixel size. from data file
    star.kappa('ndftrace', inputfile)
    pixelscales = [float(i) for i in star.read_starval('ndftrace', 'FPIXSCALE')]
    pixelsize = np.sqrt(pixelscales[0] *  pixelscales[1])
    pixel_noise_radius = int(noise_radius / pixelsize)
    pixel_source_radius = int(source_radius / pixelsize)

    # Create ARD MASK
    mask ="\n".join(
        ["COFRAME(PIXEL)",
         "CIRCLE( 0, 0, " + str(pixel_noise_radius) + " )",
         ".AND. .NOT. ",
         "CIRCLE( 0, 0, " + str(pixel_source_radius) + " )",
         ])
    f = open(ardname, 'w')
    f.writelines(mask)
    f.close()

    return ardname, pixelsize


# Get the fitsheaders we care about for the noise.
def get_noise_fitsvalues(inputfile):
    """Get noise related values.

    Need values for nboleff, tau225, wvmtau, elevation, elapsed
    time, airmass. (averaged over observation).

    First of all check for 'magic' P2 values:
    P2TAU225
    P2WVMTAU
    P2EL
    P2ELAPS

    Return dictionary of values.

    """


    star.convert('ndf2fits', inputfile, '!tempfile.fits')
    hdr = fits.getheader('tempfile.fits')
    results = {}

    results['effective_bolometers'] = hdr['NBOLOEFF']
    results['object'] = hdr['OBJECT']
    results['utdate'] = hdr.get('UTDATE', '')
    results['obsnum'] = hdr.get('OBSNUM', np.nan)
    results['filter'] = hdr['FILTER']
    results['mapheight'] = hdr.get('MAP_HGHT', np.nan)
    results['mapwidth'] = hdr.get('MAP_WDTH', np.nan)

    for key, value in results.items():
        if isinstance(value, astropy.io.fits.card.Undefined):
            results[key] = np.nan

    # Try to check the special P2 values first; if they are missing
    # assume its a single observation and calculate from the ST/EN
    # headers.
    results['average_225GHz_tau'] = hdr.get('P2TAU225',
                                            default=(hdr['TAU225ST'] + hdr['TAU225EN'])/2.0)
    if results['average_225GHz_tau'] == 0:
        results['average_225GHz_tau'] = np.nan

    results['average_wvm_tau'] = hdr.get('P2WVMTAU',
                                         default=(hdr['WVMTAUST'] + hdr['WVMTAUEN'])/2.0)
    if results['average_wvm_tau'] == 0:
        results['average_wvm_tau'] = np.nan

    results['average_elevation'] = hdr.get('P2EL',
                                     default=(hdr['ELSTART'] + hdr['ELEND'])/2.0)
    results['average_airmass'] = hdr.get('P2AM',
                                     default=(hdr['AMSTART'] + hdr['AMEND'])/2.0)

    # Date obs elapsedtime
    dateobs = datetime.strptime(hdr['DATE-OBS'].split('.')[0], '%Y-%m-%dT%H:%M:%S')

    dateend = datetime.strptime(hdr['DATE-END'].split('.')[0], '%Y-%m-%dT%H:%M:%S')
    dateobs_duration = (dateend - dateobs).total_seconds()
    results['date-obs-duration'] = dateobs_duration

    # Check for elapsed time -- either P2ELAP, ELAPTIME or TOTEXPT in that order
    if 'P2ELAP' in hdr.keys():
        results['elapsed_time'] = hdr['P2ELAP']
    elif 'ELAPTIME' in hdr.keys():
        results['elapsed_time'] = hdr['ELAPTIME']
    elif 'TOTEXP' in hdr.keys():
        results['elapsed_time'] = hdr['TOTEXP']
    else:

        logger.warning('%s has no elapsed time header: DATE-END - DATE-OBS gives %f ' % (inputfile,
                                                                             dateobs_duration))
        results['elapsed_time'] = dateobs_duration

    # Calculate the transmisison
    results['850transmission_from_225'] = np.exp( - results['average_airmass'] * 4.6 *
                                                      (results['average_225GHz_tau'] - 0.00435))

    results['850transmission_from_wvm'] = np.exp( - results['average_airmass'] * 4.6 *
                                                      (results['average_wvm_tau'] - 0.00435))



    return results



def create_pers_noise_maps(inputfile, output_noise_ext = '_per_noise',
                           output_expt_ext = '_per_expt'
                       ):
    """
    Create pers noise and exposure time maps.

    Return tuple of output (sdf) base file names.
    (noise, exp_time)
    """
    #file names
    basename = os.path.splitext(os.path.split(inputfile)[1])[0]
    noiseout = basename + output_noise_ext
    exptout = basename + output_expt_ext


    # Turn input into fits files
    exp_time_fits = basename + 'exp_time.fits'
    star.convert('ndf2fits', inputfile,  '!'+basename + '.fits')
    star.convert('ndf2fits', inputfile + '.more.SMURF.EXP_TIME', '!' +exp_time_fits)

    # data and exp_time
    datahdulist = fits.open(basename + '.fits')
    logger.debug('Opening ' + basename +'.fits for creation of pers noise map')
    data = datahdulist[0].data
    if len(data.shape) == 3:
        data = data[0,:,:]

    expthdulist = fits.open(exp_time_fits)
    expt_data = expthdulist[0].data
    if len(expt_data.shape) == 3:
        expt_data = expt_data[0,:,:]

    # Carry out the conversion
    box = np.ones([9,9])
    N = np.sum(box)
    mean = convolve2d(data, box, mode='same', boundary='symm') / N
    x2 = convolve2d(data * data, box, mode='same', boundary='symm')
    var = (x2 - N * mean * mean) / (N - 1)
    sigma = np.sqrt(var)
    exp_time = convolve2d(expt_data, box, mode='same', boundary='symm') / N

    # Set the data arrays and write out the files
    datahdulist[0].data = sigma
    datahdulist.writeto(noiseout + '.fits', clobber=True)
    expthdulist[0].data = exp_time
    expthdulist.writeto(exptout + '.fits', clobber=True)

    # convert back to sdf
    star.convert('fits2ndf', noiseout+'.fits', noiseout)
    star.convert('fits2ndf', exptout + '.fits', exptout)

    return noiseout, exptout






def get_map_information(filepath, noise_radius, source_radius):

    """
    Get various bits of map information from a file.

    filepath: string, full path to .sdf file on disk

    returns a dictionary of values.
    """

    basename = os.path.splitext(os.path.split(filepath)[1])[0]
    masked = basename + '_masked'
    masked_exptime = basename + '_exptime_masked'
    persnoise = basename + '_per_noise'
    persexpt = basename + '_per_exptime'


    # Get noise related values
    noisevalues = get_noise_fitsvalues(filepath)

    # Add in filepath and name
    noisevalues['filepath'] = filepath

    # Create Pers noise maps
    persnoise, persexpt = create_pers_noise_maps(filepath)
    logger.debug('persnoise is ' + persnoise + ' and persexpt is ' + persexpt)

    # Mask out central and other areas, both data and exp_time
    ardmask, pixelsize = make_ard_mask(filepath, noise_radius, source_radius)
    star.kappa('ardmask', in_=filepath, out=masked, ardfile=ardmask, inside='NO', comp='ALL' )
    star.kappa('ardmask', in_=filepath + '.more.SMURF.EXP_TIME', out= masked_exptime, ardfile=ardmask, inside='NO')
    star.kappa('ardmask', in_=persnoise, out=persnoise + '_masked', ardfile=ardmask, inside='NO', comp='DATA' )
    star.kappa('ardmask', in_=persexpt, out=persexpt + '_masked', ardfile=ardmask, inside='NO', comp='DATA' )

    noisevalues['pixelsize'] = pixelsize

    # Now do stats on masked map.
    star.kappa('stats', ndf=masked, comp='DATA')
    rms = star.read_starval('stats', 'sigma')[0]
    numpix = star.read_starval('stats', 'numpix')[0]
    mean = star.read_starval('stats', 'mean')[0]
    noisevalues['rms_from_map'] = rms
    noisevalues['rms_from_map_numpix'] = numpix
    noisevalues['mean_from_map'] = mean

    # Stats on variance array
    star.kappa('stats', ndf=masked, comp='ERR')

    # Note that I read the ERROR component, so do not need to take sqrt of the mean.
    rms_from_variance = star.read_starval('stats', 'mean')[0]
    numpix_from_variance = star.read_starval('stats', 'numpix')[0]
    rms_sd_from_variance = star.read_starval('stats', 'sigma')[0]
    noisevalues['rms_from_variance'] = rms_from_variance
    noisevalues['rms_from_variance_numpix'] = numpix_from_variance
    noisevalues['rms_sd_from_variance'] = rms_sd_from_variance

    # Stats on Pers noise map
    star.kappa('stats', ndf=persnoise+'_masked', comp='DATA')
    rms_from_per = star.read_starval('stats', 'mean')[0]
    numpix_from_per = star.read_starval('stats', 'numpix')[0]
    rms_sd_from_per = star.read_starval('stats', 'sigma')[0]
    noisevalues['rms_from_persmap'] = rms_from_per
    noisevalues['rms_sd_from_pers'] = rms_sd_from_per


    # Now do stats on exp_time array.
    star.kappa('stats', ndf=masked_exptime)
    expt_sigma = star.read_starval('stats', 'sigma')[0]
    expt_numpix = star.read_starval('stats', 'numpix')[0]
    expt_mean = star.read_starval('stats', 'mean')[0]
    noisevalues['expt_mean' ] = expt_mean
    noisevalues['expt_sigma'] = expt_sigma
    noisevalues['expt_numpix'] = expt_numpix

    # Stats on Pers exposure time map
    star.kappa('stats', ndf=persexpt+'_masked', comp='DATA')
    persexpt_sigma = star.read_starval('stats', 'sigma')[0]
    persexpt_numpix = star.read_starval('stats', 'numpix')[0]
    persexpt_mean = star.read_starval('stats', 'mean')[0]

    noisevalues['expt_from_pers_mean' ] = persexpt_mean
    noisevalues['expt_from_pers_sigma' ] = persexpt_sigma
    noisevalues['expt_from_pers_numpix' ] = persexpt_numpix

    #Calculate various nefd
    noisevalues['nefd_map'] = 725 * noisevalues['rms_from_map'] * np.sqrt(noisevalues['expt_mean'])
    noisevalues['nefd_pers'] = 725 * noisevalues['rms_from_persmap'] * np.sqrt(noisevalues['expt_from_pers_mean'])
    noisevalues['nefd_variance'] = 725 * noisevalues['rms_from_variance']* np.sqrt(noisevalues['expt_mean'])
    noisevalues['nefd_map_effective'] = noisevalues['nefd_map'] / noisevalues['effective_bolometers']
    noisevalues['nefd_pers_effective'] = noisevalues['nefd_pers'] / noisevalues['effective_bolometers']
    noisevalues['nefd_variance_effective'] = noisevalues['nefd_variance'] / noisevalues['effective_bolometers']

    # Add on the source and noise masking radii for completeenss
    noisevalues['source_radius'] = source_radius
    noisevalues['noise_radius'] = noise_radius

    return noisevalues

# Create output table
from astropy.units import arcsec, Jy, mJy, second, pW, micron, beam

columns = [Column(name='filepath',dtype='|S100'),
               Column(name='name',dtype='|S50'),
               Column(name='object',dtype='|S20'),
               Column(name='utdate',dtype='int64'),
               Column(name='obsnum',dtype='int64'),
               Column(name='pixelsize', unit=arcsec),
               Column(name='mapheight', unit=arcsec),
               Column(name='mapwidth', unit=arcsec),

               Column(name='effective_bolometers'),
               Column(name='elapsed_time', unit=second, description='from ELAPTIME, or DATE-OBS hdr if not set'),

               Column(name='average_elevation'),
               Column(name='average_airmass'),
               Column(name='average_225GHz_tau'),
               Column(name='average_wvm_tau'),
               Column(name='850transmission_from_225'),
               Column(name='850transmission_from_wvm'),

               Column(name='nefd_map',unit=(Jy/beam)*(second**0.5) ),
               Column(name='nefd_pers', unit=(Jy/beam)*(second**0.5) ),
               Column(name='nefd_variance', unit=(Jy/beam)*(second**0.5)  ),

               Column(name='nefd_map_effective', unit=(Jy/beam)*(second**0.5)),
               Column(name='nefd_pers_effective', unit=(Jy/beam)*(second**0.5)),
               Column(name='nefd_variance_effective', unit=(Jy/beam)*(second**0.5)),

               Column(name='rms_from_map', unit=pW),
               Column(name='rms_from_variance', unit=pW),
               Column(name='rms_from_persmap', unit=pW),

               Column(name='rms_sd_from_pers', unit=pW),
               Column(name= 'rms_sd_from_variance', unit=pW),

               Column(name='expt_mean', unit=second),
               Column(name='expt_from_pers_mean', unit=second),

               Column(name='mean_from_map', unit=pW),
               Column(name='rms_from_map_numpix'),
               Column(name='rms_from_variance_numpix'),

               Column(name='expt_numpix', unit=second),
               Column(name='expt_sigma', unit=second),

               Column(name='expt_from_pers_numpix', unit=second),
               Column(name='expt_from_pers_sigma',unit=second),

               Column(name='filter', unit=micron),
               Column(name='noise_radius', unit=arcsec),
               Column(name='source_radius', unit=arcsec),
               Column(name='date-obs-duration', unit=second, description='DATE-END - DATE-OBS')
           ]




outtab = Table(columns)

# Go through maps, set at various radiss
fh = open(filelist, 'r')
maps = fh.readlines()
fh.close()


for m  in maps:
    fullpath = m.strip()
    logger.info('Working on %s' %fullpath)

    for sr in source_radius:
        for nr in noise_radius:
            # Get values
            noisevalues = get_map_information(fullpath, nr, sr)

            # Update catalogue
            outtab.add_row(noisevalues)


outtab.write(outputcat + '.fits', format='fits')
outtab.write(outputcat + '.csv', format='ascii.csv')


