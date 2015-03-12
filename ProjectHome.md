**ssa-py is a singular spectrum analysis (SSA) suite for Python!**

This suite is written for use with Python 2.7, and the wizard-like routine in ssa.py has an input-file handling routine that may be peculiar to Windows.

The root or primary ssa routines (in ssa\_root.py) were originally developed in Matlab by Eric Breitenberger circa 1995, and later extended with his permission by Tongying Shun circa 1997 and Karsten Sedmera circa 2002 under the direction of Dr. Christopher Duffy. Karsten Sedmera finished translating both the core functions and the wizard-like routine into Python in 2013. With the exception of functions that were not ideal for Python, the credits for each translated function are provided in each function's doc-string.

The main sss.py file is designed to act like a wizard that queries the user for all necessary data and parameters. Visit the [source/browse](https://code.google.com/p/ssa-py/source/browse) tab to download and/or to request permission to contribute to this project, or the [Installation](http://code.google.com/p/ssa-py/wiki/Installation) tab to learn how to install, or the [tutorial](https://code.google.com/p/ssa-py/wiki/Tutorial) to start using it.

**General SSA Instructions (adapted from Tongying Shun):**

This routine performs Singular Spectrum Analysis for one time series. It is related to principal components analysis in that it identifies principal components within a time series, i.e. similar to what you might expect from Fourier analysis. There are many good textbooks that describe the pros and cons and limitations of this technique that you need to read in order to apply the results correctly.

Before performing SSA on a time series, you ideally should select and try a few values for the  embedding dimension, M. Usually, M should not be greater than 1/3 of the length of the time series. The ratios of the principal components are represented in a normalized eigenvalue spectrum to aid in quantifying how "dominant" each empirically orthogonal function (EOF) is.

The period or frequency of each EOF is estimated from the oscillatory pairs (i.e. two phase-shifted eigenvectors). The usual goal of SSA is to use only a
small subset of the principal components from this analysis to reconstruct the time series. Some excellent introductions to the SSA method include:

  1. Article by Vautard, Yiou, and Ghil, Physica D 58, 95-126, 1992.
  1. Book by J.B. Elsner and A.A. Tsonis, [Singular Spectrum Analysis - A New Tool in Time Series Analysis](http://books.google.com/books/about/Singular_Spectrum_Analysis.html?id=pHsGF9WIBxkC).

![https://ssa-py.googlecode.com/git/Test1/Step010b_ReconstructionOriginalResidualPlot&ExplainedVariance.png](https://ssa-py.googlecode.com/git/Test1/Step010b_ReconstructionOriginalResidualPlot&ExplainedVariance.png)