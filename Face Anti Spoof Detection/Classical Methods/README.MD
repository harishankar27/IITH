# Classical Methods

This folder contains classical face anti-spoofing approaches based on image texture and quality assessment:
- `LBPdatageneration.ipynb, LBP.ipynb`,: Extracts multiscale Local Binary Pattern features, creating histograms of micro-textures to differentiate real and spoofed faces.
- `IQMcolor.ipynb`: Computes a suite of Image Quality Measures (IQMs), assembling statistical features from input images to train classifier models for spoof detection.

These methods are benchmarked using Linear Discriminant Analysis classifiers and compared across different datasets. Theoretical references include seminal works on LBP and IQA for biometric security. Classical methods excel on familiar data but experience significant accuracy drops in cross-dataset scenarios.

References:
1) LBP: J. Määttä, A. Hadid and M. Pietikäinen, "Face spoofing detection from single images using micro-texture analysis," 2011 International Joint Conference on Biometrics (IJCB), Washington, DC, USA, 2011, pp. 1-7, doi: 10.1109/IJCB.2011.6117510.
2) LBP: T. Ojala, M. Pietikainen and T. Maenpaa, "Multiresolution gray-scale and rotation invariant texture classification with local binary patterns," in IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 24, no. 7, pp. 971-987, July 2002, doi: 10.1109/TPAMI.2002.1017623.
3) IQA: J. Galbally, S. Marcel and J. Fierrez, "Image Quality Assessment for Fake Biometric Detection: Application to Iris, Fingerprint, and Face Recognition," in IEEE Transactions on Image Processing, vol. 23, no. 2, pp. 710-724, Feb. 2014, doi: 10.1109/TIP.2013.2292332.


