# Extracting cutouts

[...]


# Building segmentation map

You can build a segmentation map with the ```make_detect``` tool. It is used as follows:
```bash
./make_detect files=[wfc3-f105w,wfc3-f125w.fits,wfc3-f160w.fits] flux_image=wfc3-f160w.fits
```

This will try to detect galaxies in each image listed in ```files=[...]```, as well as on the stacked image (to increase S/N and possibly find fainter sources).

* Before looking for detections, it will convolve the images by a Gaussian profile of width 0.2 arcseconds (can be changed with ```conv_fwhm=...```). Increasing the size of the convolution radius will increase sensitivity to extended sources, but it will reduce sensitivity to point sources and it can also make it harder to separate neighboring galaxies. The default value of 0.2 arcseconds is optimal for point sources, because this is the width of the Hubble PSF in the H band.

* Once galaxies are identified in the images, the code will discard any galaxy with a segmentation area less then 10 pixels (can be changed with ```min_area=...```). This threshold avoids including image artifacts as "sources" (i.e., single pixel noise spikes). Increasing this value will improve the purity of the detection catalog (make sure that all detections are real), at the expense of completeness (will miss smaller and fainter galaxies).

* Then, it computes fluxes on the image specified in ```flux_image=...```, and discard any source with a S/N lower than 10 (can be changed with ```min_snr=...```). This is a double check that the detections are real.

* This is the end of the automatic detection algorithm. Then, the code does the manual detection pass. There are two steps in this. First, it will open the DS9 region file ```merged_segments.reg``` (if it exists). For each circular region listed in this file, the code will find every segment (or source) that is included, even partially, within the circle, and it will merge all these segments into a single segment. This is useful to deal with cases where the code has over-deblended the images, for example by splitting a large clumpy galaxy into multiple sources.

* In the second step of the manual detection pass, the code will open the DS9 region file ```det_manual.reg``` (if it exists). For each circular region listed in this file, it will create a new segment (source) at the position of the circle, with the size specified by the circle radius. This is useful if you have a source that exists in other images (e.g., Spitzer IRAC) but is undetected in the HST images.

* Once all of this is done, it will create a file called ```det_seg.fits``` which contains the segmentation map. This is an image where each pixel is given the value of the index of the source it belongs to (and zero if there is no source there). It will also save a FITS catalog called ```det_cat.fits``` which contains the IDs, positions, areas, and detection flux for each source in the segmentation map. By construction, the sources in the catalog will be roughly sorted by decreasing flux, so the first sources (with the first IDs) will be the brightest.

* The recommended way to run this is to do a first pass without any of the manual detection steps (i.e., do not immediately create the ```merged_segments.reg``` and ```det_manual.reg``` files), inspect the resulting segmentation map, and only create the manual detection region files if needed.



# Creating models

Once the segmentation map is built, you need to create a model for each source that was detected. For each source with ID ```{ID}```, there must exist a model called ```{ID}.fits``` in a sub-directory called ```imfit/models```. These models must be FITS images with valid WCS astrometry, they must be un-convolved (not convolved by any PSF), and they can be stored with any pixel size as you wish (although smaller pixels will always give better results). The position of the model galaxy in the FITS image is not important, as long as the WCS information in the FITS header is specified such that the center of the model corresponds to the real physical coordinates of the corresponding galaxy as found in the detection image.

If you want to, you can build these models manually with any method you prefer. However there is an automatic method provided here that can produce these models for you. The method fits Sersic profiles to the galaxies in a reference image (the stacked detection image). This is done with the program ```go_imfit```. It has no command line parameters or configurable options, so all you have to do is to run it:
```bash
./go_imfit
```

It will automatically read the reference image as ```det_stack.fits```, its PSF as ```det_stack_psf.fits``` (which you can create, for example, as a copy of the HST F160W PSF), the segmentation map as ```det_seg.fits```, the detection catalog as ```det_cat.fits```. Any image artifact in ```det_stack.fits``` can (optionally) be masked with DS9 circular regions listed in the file ```det_stack-bad.reg```.

When you run it, the program does the following.

* It starts with the first source in the detection catalog. Based on the area found in the detection image, it extracts a cutout of size about twice larger (to include the outskirts of the galaxy) on both the reference image and on the segmentation map.

* If the cutout is large, it may contain other galaxies too. If there are other galaxies that are only partly included in the cutout (i.e., they extend beyond the edges of the cutout), they will be automatically masked and ignored. If there are other galaxies that are fully contained in the cutout, they will be kept and modeled together with the galaxy.

* For each galaxy kept in the cutout, the code will mask all the other galaxies, and use IMFIT to find a first estimate of the best-fit Sersic profile of the galaxy.

* Once it has a first estimate for each galaxy, it will run a last pass of IMFIT by fitting all the galaxies of the cutout together, using the first estimate of their profiles as starting guess. This step is skipped if there is only one galaxy in the cutout.

* Then for each galaxy in the cutout, it will use IMFIT to create the best-fit un-convolved model, and save it in the ```imfit/models``` folder.

* It then subtracts all the galaxies in the cutout from the reference image, marks them as "extracted", and goes on to the next source in the catalog that is not yet extracted.

* This is repeated until there is no source left in the catalog.

The whole process can take some time, and you can monitor the progress of the code by opening the ```imfit/det_stack-sub.fits``` file. This will show you the reference image, where all the galaxies that have been extracted have been subtracted from the image, and it will be updated each time the program has finished extracting a galaxy.


# Adjusting the IMFIT models

In some cases IMFIT can crash, or it can produce a bad fit and you wanted to stop it. To avoid loosing time, the ```go_imfit``` program will automatically re-read the models found in ```imfit/models``` and it will re-use them without re-computing the Sersic profiles. This means that, if you think a model is not correct, you can stop the program, delete the corresponding model file from the ```imfit/models``` directory, make some adjustments, and run ```go_imfit``` to resume fitting.

If a model is incorrect, you should therefore stop ```go_imfit``` and adjust the IMFIT configuration. This can be a bit tedious, so only do this if you think the model is very wrong, or if the galaxy with a suspicious model is very close to the galaxy you're interested in. The way to do this depends on whether the galaxy you want to adjust was alone it its cutout, or if it was fitted simultaneously with other galaxies. In all cases, when you need to run one of the ```*.sh``` scripts, you *have* to run the script from the directory with the images (where ```det_stack.fits``` is located). For example, if you need to run the script ```imfit/7/fit.sh```, do not go to the directory ```imfit/7``` and run ```fit.sh```! Stay in the base directory, and run ```imfit/7/fit.sh``` directly from there.

**For a galaxy alone it its cutout**. Find the directory ```imfit/{ID}```, where ```{ID}``` is the ID of the source you want to fix. If there is no directory with this name, this means that the galaxy was included in another galaxy's cutout. Likewise, if this directory contains a script file called ```fit.sh```, this means there are other galaxies in the cutout. See further below for instructions for such cases.

Open the file ```imfit/{ID}/indiv/imfit_config_{ID}.dat```. Adjust the starting parameters to your liking (see below for advice), save the file, then run ```imfit/{ID}/indiv/fit_{ID}.sh``` to let IMFIT do the fit again. This will update the best-fit model in ```imfit/{ID}/indiv/imfit_bestfit_{ID}.dat```. When IMFIT is finished, look at the residual image ```imfit/{ID}/indiv/imfit_bestfit_res_{ID}.fits``` to see if you have improved the model.

Once you are happy with the model, copy ```imfit/{ID}/indiv/imfit_bestfit_{ID}.dat``` into ```imfit/{ID}/indiv/imfit_bestfit_full_{ID}.dat```, then execute the scripts ```imfit/{ID}/indiv/make_{ID}.sh``` and ```imfit/{ID}/indiv/make_conv_{ID}.sh```. These scripts will produce the updated model.

**For a galaxy with other galaxies in the cutout**. Find the directory ```imfit/{ID}```, where ```{ID}``` is the ID of the source you want to fix. If there is no directory with this name, this means that the galaxy was included in another galaxy's cutout. If this is the case, look inside each directory ```imfit/{OTHER_ID}/indiv``` for a file called ```imfit_config_{ID}.dat```. This will tell you that your source was included in the cutout of ```{OTHER_ID}```.

In the directory ```imfit/{MAIN_ID}``` (where ```{MAIN_ID}``` is either ```{ID}``` or ```{OTHER_ID}```, see above), open the file ```imfit/{MAIN_ID}/indiv/imfit_config_{ID}.dat```. Adjust the starting parameters to your liking (see below for advice), save the file, then run ```imfit/{MAIN_ID}/indiv/fit_{ID}.sh``` to let IMFIT do the fit again. This will update the best-fit model in ```imfit/{MAIN_ID}/indiv/imfit_bestfit_{ID}.dat```. When IMFIT is finished, look at the residual image ```imfit/{MAIN_ID}/indiv/imfit_bestfit_res_{ID}.fits``` to see if you have improved the model.

Once you are happy with the model, you have to update the full fit. To do so, open the file ```imfit/{MAIN_ID}/imfit_config.dat``` and look for the line that starts with ```X0 ...``` and ends with ```# ID={ID}``` (again, where ```{ID}``` is your galaxy ID). Then, modify the value of each parameter to match the best-fit values from ```imfit/{MAIN_ID}/indiv/imfit_bestfit_{ID}.dat```. When doing so, make sure that the best-fit value is still within the parameter bounds specified in ```imfit/{MAIN_ID}/imfit_config.dat```, and if not, adjust the bounds accordingly. When you are done, run the script ```imfit/{MAIN_ID}/fit.sh``` to let IMFIT do the full fit. It will save the best-fit values in ```imfit/{MAIN_ID}/imfit_bestfit.dat```, and produce a residual image in ```imfit/{MAIN_ID}/imfit_bestfit_res.fits```. You can iterate this process until you are satisfied with this residual image.

Once you are happy with the full model, do the following for each galaxy in the cutout (with ID ```{GID}```). Open the file ```imfit/{MAIN_ID}/indiv/imfit_bestfit_full_{GID}.dat```, and copy the corresponding best-fit values from ```imfit/{MAIN_ID}/imfit_bestfit.dat``` (NB: the order of the galaxies in this file is the same as in ```imfit/{MAIN_ID}/imfit_config.dat```, make sure not to mix them up!). Then execute the scripts ```imfit/{MAIN_ID}/indiv/make_{GID}.sh``` and ```imfit/{MAIN_ID}/indiv/make_conv_{GID}.sh```. These scripts will produce the updated models for each galaxy.

**Advice and instructions for changing IMFIT parameters.** Briefly, ```X0``` and ```Y0``` are the centroid coordinates, ```PA``` is the position angle (in degrees), ```ell``` is the ellipticity (```0``` is perfectly circular, ```1``` is an infinitely thin galaxy and should be avoided), ```n``` is the Sersic index, ```I_e``` is the surface brightness at the half-light radius, and ```r_e``` is the half-light radius (a radius of zero should be avoided). Each parameter can be specified as just a single value ```X```, or a value and two bounds as ```X XLOW,XUP```. The bounds can be used to restrict the parameter to a particular range. If you do not specify the bounds, the parameter can take any value. If the values of these parameters look fine, just change them a bit so IMFIT has a chance to converge to a better solution, or specify some bounds to exclude the fit solution that you thought was incorrect (for example a too large ellipticity, or a too small size, etc.).


# De-blending other images

This is done with the ```extract.sh``` script. It has a large number of arguments, which *must* be provided in the order below (see an example in ```get_all.sh```):

1. The name of the image, without the extension. For example, ```wfc3-f160w``` if the image is called ```wfc3-f160w.fits```.

2. The path to the corresponding PSF. For example ```wfc3-f160w-psf.fits```.

3. The over-sampling factor of the PSF. Unless you know otherwise, there is probably no over-sampling, so you should leave this to ```1.0```.

4. The noise RMS in the image, computed in regions where there are no sources. This must be done manually, by opening the image in DS9, choosing a region without detectable sources, creating the largest possible circular region that does not touch any source, and reading the standard deviation of the pixel values given by DS9. It is not necessary to get this value exactly right, since the code will perform Monte Carlo simulations for computing the flux uncertainty for the source you are interested in. But it is good to give it a realistic value, as it will determine the uncertainties for all the other sources in the image (which can be used later for flux calibration).

5. The shift of the image in the X direction. You should leave this to zero at first, and we will see later if and how to adjust it.

6. The shift of the image in the Y direction, see above.

7. The right ascension coordinate of the "region of interest" of the image. This should be the position of the source you are interested in. The code will only create Monte Carlo realizations of the image noise around this position, so it is important that it is centered on the region you care about.

8. The declination coordinate of the "region of interest", see above.

9. A ```0``` or ```1``` flag to turn on large-scale background subtraction. In some images (particularly deep ground-based images like Subaru) there is a substantial variation of the background level at different locations in the image. If this is not taken care of, this can lead to systematic errors in the flux measurements. If you can see such kind of varying background on the final residual image, you should turn this option on. The code will then mask the sources specified in the ```mask_all.reg``` DS9 region file, and compute the background level using a moving boxcar average. This can take some time. It will then subtract the background from the image (leaving behind a backup copy with the name ```...-orig.fits``` in case you want to undo this).

10. A ```0``` or ```1``` flag to enable masking of more galaxies at the edges. By default (if you leave this flag to ```0```), when doing the extraction the code will mask all the sources listed in ```mask_border.reg```. This is needed because usually the detection images (e.g., from Hubble) are smaller than the images you want to measure the fluxes (e.g., from Spizter); the sources at the edge of these images must be masked otherwise they will bias the fit. If you set this flag to ```1```, the code will instead mask the sources listed in ```mask_border_more.reg```. This is useful, for example, with deep ground-based optical images which may contains a lot of very faint blue galaxies that do not show up in the detection images, but that you cannot mask on other images like Spizter because this will mask essentially all the image.

11. A ```0``` or ```1``` flag to enable computing the image astrometry shift. Leave it to ```0``` for now, we will discuss this later.

12. The radius (in arcseconds) of the area to use for computing the astrometry shift. A typical good value is about one or two arcseconds.

13. The source ID (from the segmentation map) to use as reference for the astrometry. You can leave a placeholder ID of ```1``` until you actually do the astrometry alignment later.

14. The astrometry shift step size, in arcseconds. The code will adjust the astrometry by trying a number of position shifts on a grid, and this determines the step size of the grid. A good value is a fraction of one HST pixel, or less than 0.06 arcseconds. A smaller step size will result in more precise astrometry alignment, but it will also take longer to compute.

15. The maximum allowed astrometry shift, in arcseconds. Astrometry should be good to 0.1-0.2 arcseconds, so setting a limit of 0.2 to 0.3 arcseconds is recommended. Increasing this value will make the code take longer to compute the astrometric shift, and it also increases the risks of the automatic alignment picking the wrong solution. On the other hand, too small values may prevent the code from reaching the true astrometric shift that is required to align the images.

Before you run the extraction, there are a couple of things you have to do:

* Create a DS9 region file called ```mask_border.reg``` containing circular regions that will be used to mask the sources at the image borders, which are not part of the detection image. The recommended way to do this is to open the first image, create the regions to mask the border galaxies, and save the regions into the file. Then, open the second image, load the region file, and add new regions to cover any new source that you want to mask (do not remove any region). Save the updated regions back into the file, and repeat this for all images.

* Inspect each image for artifacts and other image defects. This can include star spikes, saturated pixels, image borders, regions with abnormal noise, or any other "weirdness" that you believe is not in the real sky. If you find any, use DS9 to create as many circular regions as needed to mask out all the artifacts in the image. Then, save the regions into a file called ```{IMG}-bad.reg``` (where ```{IMG}``` is the image name provided in parameter 1).

* For each image, create the positions that will be used for the Monte Carlo noise realizations as circular regions in DS9, and save them in a file called ```{IMG}-sim.reg```. The way this Monte Carlo procedure works is the following. At each position listed in the ```{IMG}-sim.reg``` file, the code will extract the pixels within a circular aperture centered on the position, and with a radius defined by the radius of the circular region (this radius should be chosen to be larger than the image PSF, or larger than the galaxy, whichever is largest). It will then translate these pixels to the center coordinate of the main galaxy (the coordinates specified in parameters 7 and 8 of ```extract.sh```), and sum them with the exiting pixels there. Effectively, this creates another realization of the image noise at the galaxy's position. The fit is done again on this "perturbed" image, and the flux of the galaxy is saved. Repeating this for all positions in the region file, the flux uncertainty is finally computed from the standard deviation of all the measured flux values for the galaxy. It is therefore important that:

 * There should be enough regions to compute the standard deviation. With 10,20,50,100 regions, the accuracy of the flux uncertainty is 25%,17%,10%,7% respectively. So 10 regions should be a minimum, but there is not much to gain by going beyond 50 regions. I typically choose between 15 and 25 regions.

 * The circular regions defined in the file should be mostly empty. The pixels will be extracted from the residual image, so in principle it is OK if a circular region contains a bit of a galaxy, but you should try to avoid it if possible. Otherwise, you may include parts of the residual image with higher noise (because of bad fit residuals, for example), and this will lead to an over-estimate of the flux uncertainty. For this reason, it is recommended to look at the residual for each image after a first pass of the extraction, and use them to refine the list of positions if needed, and then re-run the extraction with the updated list.

 * The circular regions should not overlap if possible, and if they do, the overlap should be minimal. Up 20% overlap is perfectly fine, and 50% overlap should be the maximum. Otherwise, the noise realizations will not be independent, and this will lead to an under-estimate of the flux uncertainty.


# Adjusting astrometry

In most cases the astrometry of the images is correct within 0.1 arcseconds. This is good enough for computing fluxes in optical-NIR images, where the resolution is good enough and blending is not a big issue, however when de-blending Spizter IRAC images these shifts can cause significant biases in the fluxes (because the de-blending can go very wrong).

To correct the astrometry, the procedure that is used here is to pick one particular source in the image, typically a bright and compact galaxy (but not a star, since stars can have proper motions and can saturate). This is the parameter 13 in ```extract.sh``` (the ID of the galaxy). This galaxy is identified in the segmentation map, which is used as the reference for the "true" astrometry (it does not matter if this is not really the "true" astrometry, as long as all the images are aligned to the same astrometry). The code then looks for this object in the other images, and tries to find the best fit by slightly adjusting the astrometry until the residuals are optimal. The astrometry is varied on a grid of small position shifts, controlled by the parameters 14 (step size) and 15 (maximum shift). To focus only on the "reference" galaxy used for the alignment, the code will only consider the residuals within a small radius around the galaxy, controlled by the parameter 12.

The recommended way to proceed is the following:

* First run the extraction without any astrometric shift and look at the residuals. If the residuals are imperfect, where you see a clear trend that all galaxies seem to be slightly offset from where the model thought they were, then it may be necessary to correct the astrometry.

* Pick a galaxy that is reasonably isolated, bright, compact, and with a well-defined center. Avoid spiral galaxies, clumpy galaxies, mergers. Avoid stars, which can have proper motions, and avoid point-like sources (stars or quasars) since the PSF features (e.g., spikes) will change dramatically from one image to the next and can perturb the alignment. The best sources for this are bright ellipticals, because they have a steep core and limited color gradients. If you have a bright galaxy that would be a perfect candidate, except that it has a neighbor, you can include this neighbor in the fit; to do this, write the parameter 13 of ```extract.sh``` as ```[id1,id2]``` where ```id1``` is the ID of the bright galaxy and ```id2``` is the ID of its neighbor (you can include as many neighbors as you like, but the bright galaxy must be listed first in the list).

* Choose a fit radius (parameter 12) small enough that, when drawing a circle around the bright galaxy, no other galaxy is included (except possibly the neighbors that you have listed in parameter 13, see above).

* Set the parameter 11 to ```1```, and run ```extract.sh```. The code will adjust the astrometry and will print out the best shift it found, for example:
```bash
best shift: -0.05 0.03
```

* Before using these values, check the residual (```{IMG}-shift-sub.fits```, where ```{IMG}``` is the image name specified in parameter 1) to make sure it looks fine. If not, try to adjust the grid step size, the maximum allowed shift, or try using another galaxy.

* When you are satisfied, copy the shift values found by the code in the parameters 5 and 6.

* Set the parameter 11 to ```0```, and run ```extract.sh``` again to re-do the flux measurements for all sources. Check the final residual to see if it has improved, and if so you can keep the new fluxes.
