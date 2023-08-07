# Angular Ptychograhic Imaging with Closed-form method
This is the data for Angular Ptychograhic Imaging with Closed-form method (APIC). We provide an example for performing reconstruction using the data included here. Interactive example of the reconstruction result can be found [here](https://rzcao.github.io/APIC_Results/).

To get started, simply replace the folder's name in the main code `APIC_reconstruction.m` using the name where the data is. Assume we want to reconstruct the thyroid sample which was imaged using a highly aberrated imaging system, which is inside a folder named "Data". Then, we modify the code as
``` matlab
folderName = 'Data';
```
As there is only one file inside the `Data` folder whose name contains "HEpath", we can ask the program to find data with name "HEpath":
``` matlab
fileNameKeyword = 'HEpath';
```
Let's say we want to reconstruct a patch with size length of 256 pixels, we set `useROI` to be `true` and let `ROILength` equals 256:
``` matlab
useROI    = true;
ROILength = 256;
```
You can then run the code and see the reconstruction result using APIC.

# Parameters
We note that parameter tunning is unnecessary as APIC is an analytical method. However, there are several options that you can choose in performing the reconstruction.
### Basic parameters 
1. `matchingTol`: By default, the program automatically determines the NA-matching measurements according to the calibrated illumination NA (`na_calib`, a vector). `matchingTol` specifies the maximal absolute deviation allowed for a measurement to be considered as NA-matching measurement. The number of NA matching measurements found by the program will be displayed in MATLAB's command window. Consider to change this threshold if the number mismatches with your experiment.
2. `useROI`: When it is set to `false`, the program uses the entire set in the reconstruction. It is recommended to set to `true` as APIC scales badly with respect to the patch sizes. A good practice is conducting reconstruction using multiple patches and stiching them together to obtain a larger reconstruction coverage.
3. `ROILength`: This parameter is used only when `useROI` is `true`. It specifies the patch sizes used in the reconstruction. It is preferrable to set this to be below 256.
4. `useAbeCorrection`: Whether to enable aberration correction. It is always recommended to set to `true`. We keep this parameter so that one can see the influence of the aberration if we do not take aberration into consideration.
5. `paddingHighRes`: To generate a high-resolution image, upsampling is typically requried due to the requirement of Nyquist sampling. `paddingHighRes` tells the program the upsampling ratio.

### Important options for APIC's sub-functions
#### Analytical complex field reconstruction using NA-matching measurements
Function `recFieldKK.m` is called to reconstruct the complex field reconstruction using NA-matching measurements. Required inputs for this function are the images and the illumination k-vectors (kx<sub>i</sub> ,ky<sub>i</sub>). When reconstructing multiple images, the first input should be a 3D matrix and the second input should be a 2D vector, and the last dimension of both demotes different measurements. Example with the minimal number of inputs:
``` matlab
im = measurement;
k = [kx, ky]; % in pixels in the spatial frequency domain
reconstructedField = recFieldKK(im, k);
```
Two important optional arguments for the complex field reconstruction with NA-matching measurements are `CTF` and `norm`.
1. `CTF`: Incorporate the (absolute) coherent transfer function for noise reduction. This will force the regions that are not covered by the CTF in reconstructed spectrum to be zero.
2. `norm`: Whether to use normalization in the reconstruction. We assume that the amplitude of zero-frequency in the sample's spectrum does not change with repect to the tilted illumination. Thus, when `norm` is set to `true`, we force the zero-frequency component in reconstructed spectrums to have the same amplitude (This is essentially correcting for illumination intensity variation).
Example of using this two optional arguments is shown below
``` matlab
reconstructedField = recFieldKK(im, k,'CTF', CTF, 'norm', true);
```
