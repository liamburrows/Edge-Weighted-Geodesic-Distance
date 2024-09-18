This is an implementation of the geodesic distance in both 2D and 3D first defined in [1] for selective segmentation. The code is implemented in c++ with usage in MATLAB. The geodesic distance involves solving the Eikonal equation. This repositry implements the fast sweeping method [2].

First compile the mex files by running 'compile_mex.m'

Then geodesic distance can be calculated on a 2D or 3D (greyscale) image by running geodesic2D(image) or geodesic3D(image) respectively. 

Edge weighted geodesic distance was first defined in [1] for selective segmentation, and has since been used in a number of works including [1,


[1] Roberts, Michael, Ke Chen, and Klaus L. Irion. "A convex geodesic selective model for image segmentation." Journal of Mathematical Imaging and Vision 61.4 (2019): 482-503.

[2] Zhao, Hongkai. "A fast sweeping method for eikonal equations." Mathematics of computation 74.250 (2005): 603-627.

[3] Roberts, Michael, and Jack Spencer. "Chan–vese reformulation for selective image segmentation." Journal of Mathematical Imaging and Vision 61.8 (2019): 1173-1196.

[4] Burrows, Liam, et al. "Reproducible kernel Hilbert space based global and local image segmentation." Inverse Problems & Imaging (2020).

[5] Burrows, Liam, Ke Chen, and Francesco Torella. "A deep image prior learning algorithm for joint selective segmentation and registration." International Conference on Scale Space and Variational Methods in Computer Vision. Cham: Springer International Publishing, 2021.

[6] Burrows, Liam, Ke Chen, and Francesco Torella. "Using deep image prior to assist variational selective segmentation deep learning algorithms." 17th International Symposium on Medical Information Processing and Analysis. Vol. 12088. SPIE, 2021.
