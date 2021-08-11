# hst_wfc3_lacosmic

Jupyter Notebook that cleans cosmic rays from HST WFC3/UVIS exposures

This notebook was created by Dr. Mitchell Revalski to clean cosmic rays from Hubble Space Telescope (HST) Wide Field Camera 3 (WFC3) UVIS single-visit images (FLCs). The core algorith uses Astropy's implementation of the [ccdproc.cosmicray_lacosmic](https://ccdproc.readthedocs.io/en/latest/api/ccdproc.cosmicray_lacosmic.html) routine, which is based on the original Laplacian Cosmic Ray Identification ([L.A.Cosmic](http://www.astro.yale.edu/dokkum/lacosmic/)) algorithm developed by  Pieter G. van Dokkum. Please see the software notes below for the appropriate references.

This Jupyter Notebook includes two refinements for users to fine-tune the L.A.Cosmic settings for their data. First, the user can choose to replace negative outlier pixels with local random noise so they are not incorrectly flagged by lacosmic. Second, the user can choose to have the flux values interpolated and replaced (the default behavior of L.A.Cosmic), or to only flag these pixels in the WFC3/UVIS Data Quality (DQ) arrays. This can result in cleaner drizzles with more accurate photometry when multiple exposures are available, as these pixels are ignored by [DrizzlePac](https://www.stsci.edu/scientific-community/software/drizzlepac.html) when combining exposures with [AstroDrizzle](https://drizzlepac.readthedocs.io/en/latest/astrodrizzle.html).

Please send questions, comments, and suggestions to [Mitchell Revalski](https://www.mitchellrevalski.com). Thank you, and have a nice day!

Mitchell Revalski gratefully acknowledges Marc Rafelski and Ben Sunnquist for helping to improve this code.
