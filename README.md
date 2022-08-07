# Topological data analysis (TDA) and non-homogeneous Poisson process (NHPP) model for vascular structural analysis using high-resolution CUBIC 3D image data

---

This pipeline is a workflow for analyzing high-resolution CUBIC 3D image data and involves: (1) extraction of signals from CUBIC 3D image data, (2) generation of 3D geometic features from the signals, and (3) evaluation and visualization of the structual differences of vasculatures between samples.

With the following tutorials you can

- Extract the signals from CUBIC 3d image data (Python code)
- Perform persistent homology to generate the persistent diagrams from the signals ([R code](https://github.com/nagoya-sysbiol/cubic_analysis/blob/main/tutorials/ph.ipynb))
- Fit non-homogeneous Poission process (NHPP) model to extract the topological features from the signals (R code)
- Evaluate and visualize the structural differences of vasculatures between samples using sliced Wasserstein kernel and multi-dimensional scaling (MDS) ([R code](https://github.com/nagoya-sysbiol/cubic_analysis/blob/main/tutorials/swk.ipynb))

This analysis pipeline is described in:

Kei Takahashi, Ko Abe*, Shimpei I Kubota*, Noriaki Fukatsu*, Yasuyuki Morishita, Yasuhiro Yoshimatsu, Satoshi Hirakawa, Yoshiaki Kubota, Tetsuro Watabe, Shogo Ehata, Hiroki R Ueda, Teppei Shimamura†, Kohei Miyazono† An analysis modality for vascular structures combining tissue-clearing technology and topological data analysis, Nature Communications (2022).

*Authors contributed equally

†Corresponding author

