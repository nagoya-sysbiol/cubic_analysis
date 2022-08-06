# Topological data analysis (TDA) and non-homogeneous Poisson process (NHPP) model for vascular structural analysis using high-resolution CUBIC 3D image data

---

This pipeline is a workflow for analyzing high-resolution CUBIC 3D image data and involves: (1) extraction of signals from CUBIC 3D image data, (2) generation of 3D geometic features from the signals, and (3) evaluation of the structual differences of vasculatures between samples.

With this pipeline you can

- Extract the signals from CUBIC 3d image data 
- Perform persistent homology to generate the persistent diagrams from the signals
- Fit he non-homogeneous Poission process (NHPP) model to extract the topological features from the signals
- Evaluate the structural differences of vasculatures between samples using Wasserstein Kernel.
- Visualize and compare sample differences from a topological perspective using multi-dimensional scaling (MDS)

This analysis pipeline is described in:

Kei Takahashi, Ko Abe*, Shimpei I Kubota*, Noriaki Fukatsu*, Yasuyuki Morishita, Yasuhiro Yoshimatsu, Satoshi Hirakawa, Yoshiaki Kubota, Tetsuro Watabe, Shogo Ehata, Hiroki R Ueda, Teppei Shimamura†, Kohei Miyazono† An analysis modality for vascular structures combining tissue-clearing technology and topological data analysis, Nature Communications (2022).

*Authors contributed equally

†Corresponding author

