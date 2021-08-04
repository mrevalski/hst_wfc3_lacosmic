# Define a function that calls lacosmic on each file.
# https://docs.astropy.org/en/stable/stats/index.html
# https://docs.scipy.org/doc/numpy-1.14.0/reference/routines.random.html
# https://learning.oreilly.com/library/view/python-for-data/9781449323592/ch04.html

def run_lacosmic_parallel(file, modify_flux, my_sigclip, fix_neg_pix, neg_pix_sig, random_seed, save_cr_msk, show_result):

    '''
    Parameters
    ----------

    files : list
            The list of *flc.fits filenames that will be run through the cleaning procedure.

    modify_flux : bool
            The user may choose to clean the cosmic rays from the data and allow lacosmic to
            replace the flux values with interpolated values, or to only use the results to
            flag cosmic rays as bad pixels in the DQ arrays of the files. Setting this to
            False is recommended for > 4-6 exposures, where there is sufficient clean data
            at each pixel location in other exposures for astrodrizzle to determine a proper
            flux value. In cases with < 2-3 exposures, setting this to True may be desired
            to produce significantly cleaner images, possibly at the cost of affecting the
            photometry for a small number of sources as the interpolation is approximate.

    my_sigclip : float
            The float value that will be passed to lacosmic's 'sigclip' input parameter. The
            recommended value is between 5 and 7, with 6 the most ideal. Setting lower may flag
            the centers of bright stars as cosmic rays, while higher values will not clean the
            images sufficiently.

    fix_neg_pix : bool
            Setting this option to True will replace negative outlier pixels so they are not
            accidently flagged as cosmic rays by lacosmic. The recommended value it True, and
            the threshold for pixel replacement is set by the 'neg_pix_sig' parameter.

    neg_pix_sig : float
            The float value that will be multiplied by the negative, sigma-clipped median flux
            (the background level) to set the negative pixel replacement threshold. Specifically,
            if 'fix_neg_pix' is True, then pixels below (-1 * neg_pix_sig * sigma-clipped median)
            will be replaced. The recommended value is between 3 and 8, with 5 the most ideal.
            Setting this lower will replace real noise in the data that would not be flagged by
            lacosmic, while higher values will not replace the negative pixels that lacomic flags
            as cosmic rays and then mistakenly replaces with a larger group of negative pixels.

    random_seed : bool
            The user may choose to "seed" the random number generator so the negative pixel
            replacement values drawn by numpy.random.random_sample() are exactly repeatable.
            This is required to generate cleaned files that are absolutely identical between
            different runs of the lacosmic cleaning function, so the recommended value is True.

    save_cr_msk : bool
            The user may choose to save the cosmic ray masks that are determined by lacosmic.
            This is useful for determining if bright sources are being accidently flagged, and
            for additional processing with other codes. Note that these files tend to have larger
            filesizes than the input images, so is not recommended when cleaning many files.

    show_result : bool
            The user may choose whether to display the original and cleaned images in the notebook.
            This is recommended to see the overall affect of the cleaning. Note that when 'modify_flux'
            is False, the input and output images should look identical as only the DQ arrays are
            modified.

    Notes:
            The cosmic ray detection thresholds, other than sigclip set by the user through 'my_sigclip',
            are hard-coded into the function below. In general, these settings should be robust and
            conservative for WFC3 UVIS images. However, if the centers of bright stars are still being
            flagged as cosmic rays, then increase the 'objlim' value above 18. See the end of this
            notebook for the original lacosmic documentation and explanation of each parameter. The
            defaults used in this function are copied below for reference:

    ccdproc.cosmicray_lacosmic(chip1, sigclip=my_sigclip, sigfrac=0.3, objlim=18.0, gain=1.0,
                                      readnoise=3.2, satlevel=65000.0, pssl=0.0, niter=4, sepmed=True,
                                      cleantype='meanmask', fsmode='median', psfmodel='gauss', psffwhm=1.82,
                                      psfsize=7, psfk=None, psfbeta=4.765, verbose=False, gain_apply=True)
    '''

    # Import packages.
    import ccdproc
    from astropy.io import fits
    from astropy.stats import sigma_clipped_stats
    import numpy
    import matplotlib.pyplot

    if (modify_flux == False):

        print('\033[1m\nThe modify_flux flag is set to FALSE. Only the DQ arrays will be modified...\033[0m')

    # Loop over each file, clean both extensions, and export the results.
    #for file in files:

    # Import the data file.
    image = fits.open(file)
    print('\033[1m\nWorking on file '+str(file) + '\033[0m:\n')

    chip1 = image[1].data # SCI
    chip2 = image[4].data # SCI



    # Replace very negative pixels with sigma-clipped median.
    if (fix_neg_pix == True):

        print('The negative pixel outliers flag is enabled.')

        # Determine the standard deviation of negative pixels.
        # sigma_clipped_stats() returns the mean, median, and standard deviation.
        neg_pix_bkg_chip1 = sigma_clipped_stats(chip1, sigma=1)
        neg_pix_bkg_chip2 = sigma_clipped_stats(chip2, sigma=1)

        print('The sigma-clipped median = '+str('{:.1f}'.format(neg_pix_bkg_chip1[1])+str('+/-')+\
                                    str('{:.1f}'.format(neg_pix_bkg_chip1[2])))+' counts.')
        print('The sigma-clipped median = '+str('{:.1f}'.format(neg_pix_bkg_chip2[1])+str('+/-')+\
                                    str('{:.1f}'.format(neg_pix_bkg_chip2[2])))+' counts.')

        neg_pix_bkg = -1.0*((abs(neg_pix_bkg_chip1[1]) + abs(neg_pix_bkg_chip2[1]))/2.0)

        if (neg_pix_bkg >= 0.0):

            print('\n\033[1mWARNING: The negative background value is positive!\033[0m')
            print('\033[1mThis will have undesired clipping effects on the data!\033[0m\n')

        print('Using a sigma-clipped thresh of '+str('{:.1f}'.format(neg_pix_bkg)+' counts.'))
        print('Replacing values < '+str(neg_pix_sig)+'*'+str('{:.1f}'.format(neg_pix_bkg))+' = '+\
              str('{:.1f}'.format(neg_pix_bkg*neg_pix_sig))+' in data.')

        if (random_seed == True):

            numpy.random.seed(16246)



        # Loop over each pixel in each chip so the random replacement values are unique to each pixel.
        num_orig = num_fix = a = b = 0

        for row in range(chip1.shape[0]):
            b = 0
            for col in range(chip1.shape[1]):
                if (chip1[a][b] < neg_pix_bkg*neg_pix_sig):
                    # Replace with random negative 1-sigma value.
                    chip1[a][b] = (numpy.random.random_sample()*1.0*neg_pix_bkg)
                    num_fix +=1
                else:
                    num_orig +=1
                b += 1
            a += 1

        print('Replaced '+str('{:.2f}'.format(round(((num_fix/(num_fix+num_orig))*100.0),2)))+\
              '% of the pixels in extension 1.')

        num_orig = num_fix = a = b = 0

        for row in range(chip2.shape[0]):
            b = 0
            for col in range(chip2.shape[1]):
                if (chip2[a][b] < neg_pix_bkg*neg_pix_sig):
                    # Replace with random negative 1-sigma value.
                    chip2[a][b] = (numpy.random.random_sample()*1.0*neg_pix_bkg)
                    num_fix +=1
                else:
                    num_orig +=1
                b += 1
            a += 1

        print('Replaced '+str('{:.2f}'.format(round(((num_fix/(num_fix+num_orig))*100.0),2)))+\
              '% of the pixels in extension 2.\n')

        # Replace data with the cleaned arrays.
        if (modify_flux == True):

            image[1].data = chip1 # SCI
            image[4].data = chip2 # SCI



    # Call LACOSMIC on the first channel.
    print('Cleaning the cosmic rays in extension 1...')
    chip1_clean = ccdproc.cosmicray_lacosmic(chip1, sigclip=my_sigclip, sigfrac=0.3, objlim=18.0, gain=1.0,
                                          readnoise=3.2, satlevel=65000.0, pssl=0.0, niter=4, sepmed=True,
                                          cleantype='meanmask', fsmode='median', psfmodel='gauss', psffwhm=1.82,
                                          psfsize=7, psfk=None, psfbeta=4.765, verbose=False, gain_apply=True)

    # Print the percentage of pixels flagged as cosmic rays.
    num_flag = numpy.sum(chip1_clean[1].astype(int))
    num_pix = (chip1_clean[1].size)
    print("Flagged "+str('{:.1f}'.format(round(((num_flag/num_pix)*100.0),2)))+'% of the pixels as cosmic rays.')

    # Call LACOSMIC on the second channel.
    print('Cleaning the cosmic rays in extension 2...')
    chip2_clean = ccdproc.cosmicray_lacosmic(chip2, sigclip=my_sigclip, sigfrac=0.3, objlim=18.0, gain=1.0,
                                          readnoise=3.2, satlevel=65000.0, pssl=0.0, niter=4, sepmed=True,
                                          cleantype='meanmask', fsmode='median', psfmodel='gauss', psffwhm=1.82,
                                          psfsize=7, psfk=None, psfbeta=4.765, verbose=False, gain_apply=True)

    # Print the percentage of pixels flagged as cosmic rays.
    num_flag = numpy.sum(chip2_clean[1].astype(int))
    num_pix = (chip2_clean[1].size)
    print("Flagged "+str('{:.1f}'.format(round(((num_flag/num_pix)*100.0),2)))+'% of the pixels as cosmic rays.')



    if (modify_flux == True):

        # Replace data with the cleaned arrays.
        image[1].data = chip1_clean[0] # SCI
        image[4].data = chip2_clean[0] # SCI

    if (modify_flux == False):

        # Replace only the DQ arrays with a bad value of "4" where cosmic rays were flagged.
        dq1 = image[3].data # DQ
        dq2 = image[6].data # DQ
        print('\nUpdating the DQ arrays...')
        dq1[numpy.where(chip1_clean[1].astype(int) == 1)] = 4
        dq2[numpy.where(chip2_clean[1].astype(int) == 1)] = 4

    # Write out the results with identical headers.
    print('\nExporting cleaned file...')
    image.writeto(file.split('_flc.fits', 1)[0]+'_lacos_flc.fits', overwrite=True)

    if (save_cr_msk == True):

        # Replace data with the mask arrays. Convert T/F -> 1/0.
        image[1].data = chip1_clean[1].astype(int)
        image[4].data = chip2_clean[1].astype(int)

        # Write out the results with identical headers.
        print('Exporting CR mask file...')
        image.writeto(file.split('_flc.fits', 1)[0]+'_lacos_flc_crmask.fits', overwrite=True)

    if (show_result == True):

        image_input = fits.open(file)[1].data
        image_clean = fits.open(file.split('_flc.fits', 1)[0]+'_lacos_flc.fits')[1].data
        minmax = (numpy.amax(image_input) / 1000.) # Change the denominator for scaling.
        print('\n                    Original Image                                         Cleaned Image')
        print('Chip 1:')
        figure, mysubplot = matplotlib.pyplot.subplots(1, 2, figsize=(15, 15))
        mysubplot[0].imshow(image_input, cmap='gray', vmin=0.0, vmax=minmax, origin='lower')
        mysubplot[1].imshow(image_clean, cmap='gray', vmin=0.0, vmax=minmax, origin='lower')
        matplotlib.pyplot.show()

        image_input = fits.open(file)[4].data
        image_clean = fits.open(file.split('_flc.fits', 1)[0]+'_lacos_flc.fits')[4].data
        minmax = (numpy.amax(image_input) / 1000.) # Change the denominator for scaling.
        print('Chip 2:')
        figure, mysubplot = matplotlib.pyplot.subplots(1, 2, figsize=(15, 15))
        mysubplot[0].imshow(image_input, cmap='gray', vmin=0.0, vmax=minmax, origin='lower')
        mysubplot[1].imshow(image_clean, cmap='gray', vmin=0.0, vmax=minmax, origin='lower')
        matplotlib.pyplot.show()
