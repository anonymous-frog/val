# CNN-Accelerated Adaptive Video Block Compressive Sensing for Resource Constrained Scenarios

## Some results
![](images/Side-by-Side.png)

## Requirements

* The code was developed on a Linux machine running Linux version 18.04.2 LTS using Matlab 2018a. It was also tested on Matlab 2020a on Windows.
* An Nvidia GPU is required. We use a GTX 1080Ti card with the appropriate CUDA drivers.
* Conda 4.8.4 was used to create a Python 3.7.9 environment.
* Pytorch 1.1.0 was installed with GPU support.
* MatConvNet [2] is required for the DnCNN denoiser used in IDA. We installed matconvnet-1.0-beta25 from [3].
* The DnCNN [4] networks we use were trained by Metzler et al. and were originally downloaded with the D-AMP toolbox [5]. However *we include the trained nets and the software required to run DnCNN in the zip file.
* The DAIN [6] software is required and can be downloaded from [7].
* The VidSet6 videos can be downloaded from [8].

## Installing the code

* The following folders need to be populated
	* matconvnet-1.0-beta25 [2][3]
	* DAIN-master 
	* Data
* Download matconvnet-1.0-beta25 from [3].
* Follow the instructions in the package to install it in the sub-directory matconvnet-1.0-bet25.
* Compile it with GPU support.
* Download DAIN from [7].
* Follow the instructions to install DAIN into the subdirectory DAIN-master.
* Compile the code following the instructions in the package.
* Download the cif videos, in y4m format, from [8], and store them in the Data/Videos subdirectory.
* Move the file Dain_vfi.py to the DAIN-master subdirectory and the installation is complete.

## Running the code

* Running the default ```VAL_VFI.m``` in Matlab will generate the VAL-VFI result with GOP=4 and Key CR=0.5 in the table in [1]. This code will run DAIN as an external program and exchange files with it through the DAIN-master/VAL subdirectory structure.
* Setting CS_ref = 2 and CS_nonref = 2 generates the VAL-IDA-VFI results in [1].
* We also provide ```VAL_VFI_no_DAIN_2.m``` that can generate the VAL-VFI result with GOP=4 and Key CR=0.5 by accessing pre-interpolated non-key frames.
* ```VAL_VFI_no_DAIN_1.m``` is used to generate the non-Key frames used by ```VAL_VFI_no_DAIN_2.m``` and ```VAL_VFI_no_DAIN_no_DnCNN.m```. It requires the installation of the DAIN package.
* ```VAL_VFI_no_DAIN_no_DnCNN.m``` does not require DAIN nor DnCNN packages to be installed.

## References 

* [1] CNN-Accelerated Adaptive Video Block Compressive Sensing, Anonymous CVPR 2021 submission, Paper ID 7307
* [2] A. Vedaldi and K. Lenc.  Matconvnet – convolutional neural networks for Matlab. In Proceeding of the ACM Int. Conf. on Multimedia, 2015.
* [3] https://www.vlfeat.org/matconvnet/. 
* [4] Kai Zhang, Wangmeng Zuo, Yunjin Chen, Deyu Meng, and Lei Zhang. Beyond a gaussian denoiser: Residual learning of deep CNN for image denoising. IEEE Transactions on Image Processing, 26(7):3142–3155, 2017.
* [5] D-AMP Toolbox. https://github.com/ricedsp/d-amptoolbox, 08-Jul-2017.
* [6] Wenbo  Bao,  Wei-Sheng  Lai,  Chao  Ma,  Xiaoyun  Zhang,Zhiyong Gao,  and Ming-Hsuan Yang.   Depth-aware videoframe interpolation. In Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition, pages 3703–3712, 2019.
* [7] Wenbo, B., DAIN, https://github.com/baowenbo/DAIN.
* [8] https://media.xiph.org/video/derf/.
