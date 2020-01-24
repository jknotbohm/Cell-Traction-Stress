# READ ME FOR TRACTION FORCE MICROSCOPY AND MONOLAYER STRESS MICROSCOPY

*Repository for traction force microscopy and monolayer stress microscopy.* https://notbohm.ep.wisc.edu

*Written by Notbohm Research Group, University of Wisconsin-Madison.*

This document explains the Notbohm Research Group's procedures for analyzing 
cell tractions and monolayer stresses. Relevant publications to read:

- Dembo, M. and Wang, Y. L. Stresses at the cell-to-substrate interface during 
locomotion of fibroblasts. Biophys J, 76(4):2307–2316, 1999.
- Butler, J. P., Tolic-Norrelykke, I. M., Fabry, B., and Fredberg, J. J. Traction 
fields, moments, and strain energy that cells exert on their surroundings. Am J 
Physiol Cell Ph, 282(3):C595–C605, 2002.
- del Alamo, J. C., Meili, R., Alonso-Latorre, B., Rodriguez-Rodriguez, J., Aliseda, 
A., Firtel, R. A., and Lasheras, J. C. Spatio-temporal analysis of eukaryotic cell 
motility by improved force cytometry. P Natl Acad Sci USA, 104(33):13343–13348, 2007.
- Trepat, X., Wasserman, M. R., Angelini, T. E., Millet, E., Weitz, D. A., Butler, 
J. P., and Fredberg, J. J. Physical forces during collective cell migration. Nat Phys, 
5(6):426–430, 2009.
- Tambe, D. T., Hardin, C. C., Angelini, T. E., Rajendran, K., Park, C. Y., Serra-Picamal, 
X., Zhou, E. H., Zaman, M. H., Butler, J. P., Weitz, D. A., et al. Collective cell 
guidance by cooperative intercellular forces. Nat Mater, 10(6):469–475, 2011.
- Tambe, D. T., Croutelle, U., Trepat, X., Park, C. Y., Kim, J. H., Millet, E., Butler, 
J. P., and Fredberg, J. J. Monolayer stress microscopy: limitations, artifacts, and 
accuracy of recovered intercellular stresses. Plos One, 8(2):e55172, 2013.

## SOFTWARE REQUIREMENTS

All scripts are written in Matlab. They have been tested on Matlab 2019b on a Windows 10
machine, but earlier versions of Matlab and other operating systems are likely to work.

It is recommended to have several GB of RAM available, as some steps, especially digital
image correlation, use a lot of RAM. We use machines having at least 16 GB of RAM, though
8 GB of RAM is typically sufficient.

## SUMMARY OF DATA REQUIRED

Images of two types:

1. Images of cells. Phase contrast. These images are required for the script **find_boundary.m**.
They are used to make a domain, an 8-bit image giving 0s at locations without cells and 
values of 255 at locations with cells.
2. Images of fluorescent particles in compliant substrate beneath cells. A reference image
is needed in the undeformed, stress-free state. The reference image is typically acquired
at the end of the experiment: after imaging the cells for a period of time, the cells 
can be removed from the substrate using trypsin or some other means, allowing for a stress-
free image to be acquired.

The two sets of images are used to compute displacements using digital image correlation 
(also referred to as particle image velocimetry).

A file called **ExperimentalSettings.txt**. This file contains the following lines:

X  			--- pixel size (m)

X			--- substrate Young's modulus (Pa)

X			--- substrate Poisson's ratio (-)

X			--- substrate thickness (um)

X	 		--- monolayer Poisson's ratio (matters only slightly; Young's modulus is irrelevant)

X			--- monolayer thickness (m)

1 or 2		--- 1-strip; 2-hole/island

An example is included.

## SUMMARY OF SCRIPTS TO RUN

- Compute cell-induced substrate displacements using **Digital Image Correlation** on the  
images of fluorescent particles in the substrate. We often use FIDIC, written by members of Christian 
Franck's research group. Citation: Bar-Kochba E, Toyjanova J, Andrews E, Kim K-S, Franck C.
A Fast Iterative Digital Volume Correlation Algorithm for Large Deformations. 
Experimental Mechanics 55:261–274, 2015. Available from https://github.com/FranckLab/
- Identify regions in images containing cells using **find_boundary.m**
(Calls function **smooth2a.m**, available on Matlab file exchange.)
- Compute cell-substrate tracitions. There exist numerous packages on Github and other 
sources to compute tractions. Our version of the software is available upon request.
- Compute monolayer stresses using **run_stress_calculation.m**
(Calls functions **force_moment_equilibrium.m**, **compute_stress.m**)

Each script has detailed comments.

## SUMMARY OF OUTPUTS

- **Digital Image Correlation** outputs displacement data and corresponding grid points. In addition
to these, the subset size (w0) and spacing (d0) are required. These data are stored in a 
Matlab mat file that is used in other scripts.
- **find_boundary.m** outputs an image called domain.tif that is used in other scripts. If you
only want to compute cell-substrate tractions and not monolayer stresses, you may skip running this file.
- **run_stress_calculation.m** outputs a mat file containing stresses in the x-y coordinate system (Sxx, Syy,
Sxy) and principle stresses and orientation (S1, S2, pangle).

## POSTPROCESSING

Sample scripts are included in the directory Cell-Traction-Stress-Velocity-Plots 
( https://github.com/jknotbohm/Cell-Traction-Stress-Velocity-Plots ):

- **plot_displ_tractions.m** generates color plot of substrate displacements and and cell-substrate tractions
- **plot_stress.m** generates color plots of monolayer stresses

Further analysis beyond color plots will likely be required. Such analysis has to be 
custom written by the user.

## FORMATTING REQUIREMENTS

#### Digital Image Correlation
Different versions of Digital Image Correlation software have different image formatting requirements. 
The DIC software we use reads multipage tif files, where the different pages correspond to different time
points.

#### Traction Calculation
Format data as required by the traction software you are using.

#### Stress Calculation
Computing stresses with **run_stress_calculation.m** requires the following: 
- A mat file containing the following variables: **x, y**: 2D arrays of size (M, N) containing the gridpoints on 
which the DIC was computed. This can be made using Matlab's meshgrid command. Units: pix. 
**tx, ty**: Tractions in horizontal and vertical directions. 2D or 3D arrays of size (M, N, P) where P is the
number of time points. **d0***: Scalar grid point (subset) spacing used in image correlation.
- The **domain** file described above. This file is **required** for computing stresses.

## EXAMPLE WORKFLOWS

### Example: Computing stresses for a cell island

A workflow to compute stress may be as follows:

1. Make **ExperimentalSettings.txt** file.
2. Organize your data. Save your images in tif format. Use a multipage tif file for a time lapse data set. Use a separate
single-page tif file for the reference image, acquired after removing the cells from the substrate.
Put the images and **ExperimentalSettings.txt** into a single directory.
3. Run DIC using a cumulative comparison that compares all images to the reference one. 
Store outputs in a mat file containing w0, d0, x, y, u, v as described above.
4. Run **find_boundary.m**. (Required for computing stresses.)
for drift in the next step.
5. Compute tractions using code available online or request our version of the code.
6. Format output of traction as follows: 
**w0**: scalar, subset/window size. Units: pix.
**d0**: scalar, subset spacing. Units: pix.
**x, y**: 2D arrays containing the gridpoints on which the DIC was computed. This can be made using Matlab's meshgrid
command. Units: pix.
**u, v**: Displacements computed by image correlation in horizontal and vertical directions. 
These are 2D or 3D arrays of size (M, N, P) where M and N are the number of rows and columns, which must match the
size of x and y. Variable P corresponds to different time points. If there is only one time point, then the array is 
2D (i.e., if P=1).
**tx, ty**: Tractions in horizontal and vertical directions. 2D or 3D arrays of size (M, N, P) that matches the 
size of u and v.
7. Run **plot_displ_tractions.m** to bring up a color map of the computed displacements and tractions. Verify that the displacements
and tractions are reasonable before going to the next step.
8. Run **run_stress_calculation.m**.
9. Run **plot_stress.m** to bring up color maps of the computed monolayer stresses.

### Example: Computing cell velocities

To compute cell velocities, you can apply DIC to phase contrast images of the cells as follows.

1. Organize your data. Save all images as a multipage tif. There is no separate reference image.
2. Run the DIC using an incremental comparison that compares each image to the one preceeding it in the multipage tif file.
3. Divide displacement by the time between imaging to get the velocity.

**Important note:** Digital image correlation only works to compute cell velocities for positions in the image that are confluent, 
i.e., fully covered by cells. If you have images with individual cells, you will have to use a different method, such as fluorescently
labeling cell nuclei and tracking them.

## TIPS

- Organize different data sets by putting them in different directories. Do not change the file names of the outputs
of the different scripts (DIC, find_boundary.m, run_traction_finite.m, run_stress_calculation.m). If you keep your
file names consistent, you will be able to run batch analyses.
- Most of the scripts in this repository have additional commands that were used in development or are useful for 
debugging. These commands are commented out and typically have notes indicating what they do. If you are having trouble
or get unexpected results, read through the script and identify the different commended lines to see if they help you
to address your problem.























