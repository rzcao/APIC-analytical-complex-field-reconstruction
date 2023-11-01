# Angular Ptychograhic Imaging with Closed-form method
This is a simple tutorial for Angular Ptychograhic Imaging with Closed-form method (APIC). We provide an example for performing reconstruction using the data included here. Interactive images of a few reconstruction results can be found [here](https://rzcao.github.io/APIC_Results/).

The preprint of APIC can be found on [arXiv](https://doi.org/10.48550/arXiv.2309.00755).

## Run the code with default settings
To get started, simply replace the folder's name in the main code `APIC_reconstruction.m` using the name where the data is. Assume we want to reconstruct the thyroid sample which was imaged using a highly aberrated imaging system, which is inside a folder named "reducedData". Then, we modify the code as
``` matlab
folderName = 'reducedData';
```
As there is only one file inside the `reducedData` folder whose name contains "Thyroid", we can ask the program to find data with name "Thyroid":
``` matlab
fileNameKeyword = 'Thyroid';
```
If we want to reconstruct the Siemens star target at a defocus distance of 16 µm, we can set `fileNameKeyword` to `'Siemens'`. However, we can find there are two files whose name contains this keyword. As a result, we can choose to specify the full name of the data or use an additional keyword. Here is an example of using an additional keyword (Note that the additional keyword must be selected from the words showing after the primary keyword. In the following example, the file has name '... Siemens ... defocus16 ...'):
``` matlab
fileNameKeyword   = 'Siemens';
additionalKeyword = 'defocus16';
```

Let's say we want to reconstruct a patch with size length of 256 pixels, we set `useROI` to be `true` and let `ROILength` equals 256:
``` matlab
useROI    = true;
ROILength = 256;
```
You can then run the code and see the reconstruction result using APIC.

## Parameters
We note that parameter tunning is unnecessary as APIC is an analytical method. However, there are several options that you can choose in performing the reconstruction.
### Basic parameters 
1. `matchingTol`: By default, the program automatically determines the NA-matching measurements according to the calibrated illumination NA (`na_calib`, a vector). `matchingTol` specifies the maximal absolute deviation allowed for a measurement to be considered as NA-matching measurement. The number of NA matching measurements found by the program will be displayed in MATLAB's command window. Consider to change this threshold if the number mismatches with your experiment.
2. `useROI`: When it is set to `false`, the program uses the entire set in the reconstruction. It is recommended to set to `true` as APIC scales badly with respect to the patch sizes. A good practice is conducting reconstruction using multiple patches and stiching them together to obtain a larger reconstruction coverage.
3. `ROILength`: This parameter is used only when `useROI` is `true`. It specifies the patch sizes used in the reconstruction. It is preferable to set this to be below 256.
4. `useAbeCorrection`: Whether to enable aberration correction. It is always recommended to set to `true`. We keep this parameter so that one can see the influence of the aberration if we do not take aberration into consideration.
5. `paddingHighRes`: To generate a high-resolution image, upsampling is typically requried due to the requirement of Nyquist sampling. `paddingHighRes` tells the program the upsampling ratio.

### Important options for APIC's sub-functions
#### Analytical complex field reconstruction using NA-matching measurements
Function `recFieldKK.m` is called to reconstruct the complex field reconstruction using NA-matching measurements. Required inputs for this function are the images and the illumination k-vectors (kx<sub>i</sub> ,ky<sub>i</sub>). When reconstructing multiple images, the first input should be a 3D matrix and the second input should be a 2D matrix that consists of vectors. The last dimension of the images and the first dimension of the k-vectors denote different measurements. Example with the minimal number of inputs:
``` matlab
im = measurements; % When the size of the acquired image is N-by-N and we have 2
                   % images, the size of im is N*N*2.
k  = [kx_1, ky_1;
      kx_2, ky_2]; % in pixels in the spatial frequency domain
reconstructedField = recFieldKK(im, k);
```
Two important optional arguments for `recFieldKK` are `'CTF'` and `'norm'`.
1. `'CTF'`: Incorporate the (absolute) coherent transfer function for noise reduction. This will force the regions that are not covered by the CTF in reconstructed spectrum to be zero.
2. `'norm'`: Whether to use normalization in the reconstruction. We assume that the amplitude of zero-frequency in the sample's spectrum does not change with repect to the tilted illumination. Thus, when `norm` is set to `true`, we force the zero-frequency component in reconstructed spectrums to have the same amplitude (This is essentially correcting for illumination intensity variation).

Example of using this two optional arguments is shown below
``` matlab
reconstructedFTField = recFieldKK(im, k,'CTF', CTF, 'norm', true);
```
#### Aberration extraction
Function `findAbeFromOverlap.m` is used to extract the aberration of an imaging system. Simply supply the `reconstructedFTField` generated by `recFieldKK`, the same illumination k-vector and `CTF` to this function to get the aberration (the second output is the Zernike coefficient of the extracted aberration).
``` matlab
[CTF_abe, zernikeCoeff] = findAbeFromOverlap(reconstructedFTField, k, CTF);
```
`'weighted'` can be used as an optional argument. When it is `true`, the program focuses more on places with larger weights and tends to ignore places with smaller weights, which is mainly designed to take signal-to-noise ratio (SNR) into consideration. When it is `false`, the function does not emphasize on any particular places and treats the small and large signals equally. It is set to `true` by default.

After we obtain the aberration and the complex spectrums sampled under NA-matching angle illumination, the aberration of the reconstructed spectrums is corrected and the corrected spectrums are then stitched. After stitching, we obtain `ftRecons` (Fourier transform of the reconstructed complex field) and `maskRecons` (mask where we put `true` in places that have been reconstructed).

#### Field reconstruction with darkfield measurements
Function `recFieldFromKnown.m` is used to do complex field reconstruction using darkfield measurements. Essentially, the sampled spectrum in one darkfield measurement is assumed to consist of two part: the known spectrum and the unkonwn spectrum. This function retrieves the unknown spectrum from the known spectrum. When the aberration of the system is extracted, we can use the following to call this function
``` matlab
[ftRecons,maskRecons] = recFieldFromKnown(im_DF, k_DF, ftRecons, maskRecons, CTF_abe);
```
Important optional arguments for this function are `'drift'`, `'threshold'`, and `'intensity correction'`.
1. `'drift'`: Whether to consider calibration error in illumination angle. When set to `true`, this function provides more robust results. However, this results in assuming a smaller known spectrum, which would in turn require a larger overlap ratio.
2. `'threshold'`: When the overlap ratio of the data is large, we can use this argument to improve the SNR of the reconstruction. The program calculates the same spectrum multiple times using different measurements and then averages over these independent reconstructions. It is recommended to choose a threshold number in between `[0,0.4]`. A larger threshold number leads to a longer reconstruction time.
3. `'intensity correction'`: When set to `true`, the function automatically compensates for illumination intensity differences.

One example of using these additional arguments is shown below
``` matlab
[ftRecons,maskRecons] = recFieldFromKnown(im_DF, k_DF, ftRecons, maskRecons, CTF_abe, ...
                                          'drift', true,'threshold', 0.3, 'intensity correction', true);
```

## FAQs
1. How is the illumination k-vector defined?

>In all of our functions, we assume the positive direction of `kx` points down on a screen and the positive direction of `ky` points right. This is defined based on how the coordinates are define in `imagesc`, which is used to show an image in MATLAB. The normal illumination should have [0, 0] as its illumination k-vector. 
>
>One important thing to note is that in our dataset a different definition is used for the convenience of most FPM algorithms.


2. What to do if the program tells me that 'the matrix is close to singular'?

>It is likely to happen when program finds the overlap ratio is insufficient for (at least one of) the darkfield measurements. Check if you have set `'drift'` to `true` when calling [`recFieldFromKnown.m`](#field-reconstruction-with-darkfield-measurements) when you aim to use a small dataset to do reconstruction. When working with a small dataset whose overlap ratio is small, you may want to set `'drift'` to `false`. If you want to enable the `'drift'` option, consider to increase the overlap ratio by using more tilted illuminations.

3. How should I order my measurements?

> Ordering is not necessary in APIC. APIC automatically sorts the measurements based on their illumination k-vector. However, one thing you need to check is that the number of NA-matching measurements found by the program matches your experiment. Check [basic parameter](#basic-parameters) for more information if the numbers do not match.
