# Topological data analysis (TDA) and non-homogeneous Poisson process (NHPP) model for vascular structural analysis using high-resolution CUBIC 3D image data

---

This pipeline is a workflow for analyzing high-resolution CUBIC 3D image data and involves: (1) extraction of signals from CUBIC 3D image data, (2) generation of 3D geometic features from the signals, and (3) evaluation and visualization of the structual differences of vasculatures between samples.

With the following tutorials you can

- Extract the signals from CUBIC 3d image data ([Python code](https://github.com/nagoya-sysbiol/cubic_analysis/blob/main/tutorials/se.ipynb))
- Perform persistent homology to generate the persistent diagrams from the signals ([R code](https://github.com/nagoya-sysbiol/cubic_analysis/blob/main/tutorials/ph.ipynb))
- Fit non-homogeneous Poission process (NHPP) model to extract the topological features from the signals (R code)
- Evaluate and visualize the structural differences of vasculatures between samples using sliced Wasserstein kernel and multi-dimensional scaling (MDS) ([R code](https://github.com/nagoya-sysbiol/cubic_analysis/blob/main/tutorials/swk.ipynb))

The following R codes give reproducible results of Takahashi et al., Nat Commun (2022):

- [Figure 4b and Supplementary Figures 4a and 4B](https://github.com/nagoya-sysbiol/cubic_analysis/blob/main/scripts/fig4b_suppl_fig4a_suppl_fig4b.ipynb)
- [Figures 5b and 5c](https://github.com/nagoya-sysbiol/cubic_analysis/blob/main/scripts/fig5b_fig5c.ipynb)
- [Figure 6d and Supplementary Figure 9](https://github.com/nagoya-sysbiol/cubic_analysis/blob/main/scripts/fig6d_suppl_fig9.ipynb)
- [Figure 6e](https://github.com/nagoya-sysbiol/cubic_analysis/blob/main/scripts/fig6e.ipynb)
- [Figure 7g](https://github.com/nagoya-sysbiol/cubic_analysis/blob/main/scripts/fig7g.ipynb)
- [Figure 7e and Supplementary Figure 10a](https://github.com/nagoya-sysbiol/cubic_analysis/blob/main/scripts/fig7e_suppl_fig10a.ipynb)
- [Supplementary Figure 7](https://github.com/nagoya-sysbiol/cubic_analysis/blob/main/scripts/suppl_fig7.ipynb)
- [Supplementary Figure 10b](https://github.com/nagoya-sysbiol/cubic_analysis/blob/main/scripts/suppl_fig10.ipynb)

This analysis pipeline is described in:

Kei Takahashi, Ko Abe*, Shimpei I Kubota*, Noriaki Fukatsu*, Yasuyuki Morishita, Yasuhiro Yoshimatsu, Satoshi Hirakawa, Yoshiaki Kubota, Tetsuro Watabe, Shogo Ehata, Hiroki R Ueda, Teppei Shimamura†, Kohei Miyazono† An analysis modality for vascular structures combining tissue-clearing technology and topological data analysis, Nature Communications (2022).

*Authors contributed equally

†Corresponding author

