#Tutorial using Test1 Example.

# Introduction #

This tutorial uses the provided test1.txt dataset to illustrate how to use ssa.py, which is a wizard-like tool for performing singular spectrum analysis (SSA) written for Python 2.7. It assumes that you have already installed all of the prerequisite Python libraries.

# Tutorial via Example: \Test1 #

The "Test1" dataset is an example of a monthy-averaged streamflow record, where the streamflow measurements in cubic-feet per second (cfs) are given in a single column in the text file. The ssa.py script expects only **one column of regular-interval measurements, and does not provide a means to fill any missing time-steps in your data**.

Thus ssa.py only refers to the implicit time unit of each measurement in generic terms (e.g. incrementally, cycles per time-step). However, all of the figure windows allow you to pan, zoom, export, and otherwise customize their appearance for you specific needs. The script also permits you to save any interim results at the end of your runtime session to datasets that you can later customize for your reporting needs, e.g. via Excel.

## Step 1: Select a time-series, a txt file ##

When you first run ssa.py, it will give you some preliminary information about SSA, and then open a file-selection dialogue. Simply select an appropriately-formatted text file to begin.

![https://ssa-py.googlecode.com/git/Test1/Step01_SelectATimeSeriesFile.png](https://ssa-py.googlecode.com/git/Test1/Step01_SelectATimeSeriesFile.png)

Once you have selected a text file, ssa.py will report certain statistics about the time series, and then plot the raw time-series for you. Recall again that you can pan, zoom, export, and otherwise manipulate this figure for your session-documentation needs.

![https://ssa-py.googlecode.com/git/Test1/Step02_RawTimeSeriesPlot.png](https://ssa-py.googlecode.com/git/Test1/Step02_RawTimeSeriesPlot.png)

Once you close this figure, ssa.py will optionally permit you to resample or time-average the raw measurements before beginning the actual SSA routine.

![https://ssa-py.googlecode.com/git/Test1/Step03_OptionalTimeResampling-Averaging.png](https://ssa-py.googlecode.com/git/Test1/Step03_OptionalTimeResampling-Averaging.png)

If you choose to resample or time-average the raw data, ssa.py will prompt you for the necessary parameters, and then automatically use this modified data for the subsequent SSA routine.

## Step 2: Enter M,  Evaluate the resultant Eigen-spectrum ##

Once leave the previous data-manipulation step, ssa.py will then prompt you for the SSA embedding dimension, M. M is time-lag or window over which SSA computes the autocovariance within the dataset. Reading about principal components analysis and autocovariance will help you understand how to properly set the embedding dimension, but the commonly accepted upper limit for M is between 1/3 to 1/2 the total time series lenght, N. SSA will not reliably be able to resolve the frequency of components that have a wavelength significantly larger than M, and may not accurately reflect the significance of all the components when longer trends exist in the datset (which relates to the SSA assumption of stationarity).

![https://ssa-py.googlecode.com/git/Test1/Step04_EnterAnEmbeddingDimension_A-Time-Lag.png](https://ssa-py.googlecode.com/git/Test1/Step04_EnterAnEmbeddingDimension_A-Time-Lag.png)

Once you specify M, you will be presented with 2 ways to view the eigenvalues, which you will use to find significant components within the starting time series. You will be able to select a different M-value within the same runtime-session if you decide that your first attempt was not ideal.

The first and optional plot is the traditional ranked eigenvalue graph. This script makes this plot optional because we do not find it as valuable for identifying dominant higher-frequency components as the next graph. The following shows an external example using this plot to identify (usually) low-frequency components.

![http://www.nature.com/srep/2011/110914/srep00091/images_article/srep00091-f1.jpg](http://www.nature.com/srep/2011/110914/srep00091/images_article/srep00091-f1.jpg)

After the ranked eigenvalue plot, ssa.py will then use a Fast Fourier Transform method to calculate the most dominant frequency of each eigenvector. It will then plot the normalized eigenvalue (or spectral power) versus eigenvector frequency (in cycles per time step), which we will refer to as an eigenspectrum graph.

![https://ssa-py.googlecode.com/git/Test1/Step05_Eigenspectrum_NormalizeEigenvalueOrPowerVersusFrequency.png](https://ssa-py.googlecode.com/git/Test1/Step05_Eigenspectrum_NormalizeEigenvalueOrPowerVersusFrequency.png)

The first thing to note about this graph for the present example is that most of the eigenvalues in this plot conform to a "noise envelope" or "noise floor" that slopes to the right (i.e. toward higher-frequency or smaller wavelength). This left-skew in the noise floor is commonly referred to a "red-noise" pattern, indicating that the phenomenon or process that influenced this measurement prefers or allows-to-pass more low-frequency variation than high-frequency variation (i.e. behaving like a a low-pass filter or sieve).

The next thing to note about this graph are the few components (roughly 6 pairs of eigenvalues) that appear above the rest of the eigenvalues. These higher eigenvalues explain more of the total measurement variance than the rest. We will subsequently try to select these dominant components for further analysis. However, to our knowledge, there is no known significance to the few eigenvalues that are significantly below the rest.

## Step 3: Select eigenvalue pairs for evaluation ##

Once you close, ssa.py will then ask you if you want to select components from the eigenspectrum or not. If not, it then ask you if you want to consider a different M-value before terminating your analysis session. However, if you do opt to select components, it will then ask you if you want to fit a "red-noise" or power-law function to noise floor. Whichever fitting function you select, ssa.py will then give you another eigenspecturum graph. However, this eigenspectrum graph will be interactive in that it will be waiting for you to click two points within the graph/figure to define a line. Upon your second-click (a point within the graph), ssa.py will then draw a line between your two points. When you close this figure, ssa.py will then use that line (or if you specified one) to select the components that are above this line and fit your remaining components with the noise-floor-fitting method that you chose beforehand. If however, you didn't specify a line, ssa.py will simply select the highest (i.e. in spectral power) 20 eigenvalue components instead.

![https://ssa-py.googlecode.com/git/Test1/Step06_NoiseFloorFittingFunction&ComponentSelectionViaLine.png](https://ssa-py.googlecode.com/git/Test1/Step06_NoiseFloorFittingFunction&ComponentSelectionViaLine.png)

Once you have selected components this way, ssa.py will then open two new figures. The first will be an eigenspectrum graph with selected components circled and with the fitted noise-floor between dashed and area-shaded 95-percentile confidence intervals. The second graph will plot the eigenvectors for the selected components with labels indicating the dominant wavelength.

The second graph allows you to see how cleanly SSA was able to extract each presumably-dominant component. It also allows you to see if the selected components are really pairs of eigenvalues. Pairing is important because a pair of similarly-positioned eigenvalues indicates an authentic oscillation or periodicity in the raw dataset, as opposed to a trend (often extracted by one or an odd number of harmonic eigenvectors). Any non-stationary trend should be removed before trusting SSA to reliably extract other significant components, BTW.

![https://ssa-py.googlecode.com/git/Test1/Step07_FittedNoiseFloorWithSelectedComponents&SelectedEigenvectorsWithPeriod.png](https://ssa-py.googlecode.com/git/Test1/Step07_FittedNoiseFloorWithSelectedComponents&SelectedEigenvectorsWithPeriod.png)

Once you close these two figures, ssa.py will ask you a series of other optional questions for further analysis. First it will ask you if you wish to plot the eigenspectrum in semi-log scale or if you want to select other components.

![https://ssa-py.googlecode.com/git/Test1/Step08_VerifyComponentSelection&EigenspectrumScale.png](https://ssa-py.googlecode.com/git/Test1/Step08_VerifyComponentSelection&EigenspectrumScale.png)

Assuming you don't go back to select other components, ssa.py will then offer to plot any eigenvector pages of your choice.

After that, it will then allow you to plot any eigenvector-pairs against one another (i.e. as a parametric plot). This step requires you to select any eigenvalues that you want to plot by remembering their rank in the previous plots.

![https://ssa-py.googlecode.com/git/Test1/Step09_OptionalPairedEigenvectorPlots.png](https://ssa-py.googlecode.com/git/Test1/Step09_OptionalPairedEigenvectorPlots.png)

The last analysis that ssa.py will offer is to plot the reconstructed time series for any of your selected components against the original time series and with a residual plot (difference between the original and the reconstruction). The reconstruction steps offer several options for how to plot everything, and the ability to plot different sets of components (i.e. again by specifying the rank of the desired components). Each time it plots a reconstruction, ssa.py will report statistics for how much of the original variance is explained by the currently plotted reconstruction.

![https://ssa-py.googlecode.com/git/Test1/Step010b_ReconstructionOriginalResidualPlot&ExplainedVariance.png](https://ssa-py.googlecode.com/git/Test1/Step010b_ReconstructionOriginalResidualPlot&ExplainedVariance.png)

Once you exit the reconstructions, ssa.py will give you the option to save the latest version of any SSA variables to disk (i.e. as text files). That is, if you have selected components or reconstructions repeatedly, only the last-selected version of those variables can be saved to disk.

![https://ssa-py.googlecode.com/git/Test1/Step011_OptionalSave-A-Variable-eg-reconstructions.png](https://ssa-py.googlecode.com/git/Test1/Step011_OptionalSave-A-Variable-eg-reconstructions.png)