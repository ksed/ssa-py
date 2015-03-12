#Explains installation prerequisites, and execution.

# Installation #

Like most Python scripts, there are no formal installation steps beyond unzipping the project files to a folder of your choice. However, to run the main script, ssa.py, you must first have or install at least (I think) Python 2.6.5 or later as well as the [Numpy](http://www.numpy.org), [Scipy](http://www.scipy.org), and [Matplotlib](http://matplotlib.org) libraries (all of which are free and highly recommended for any science-related endeavor).

These scripts were furthermore designed in a Windows environment using Python 2.6.5 and 2.7.5, and thus all of the file-selection and output routines may have unintentional Windows-related sympathies. The code should be fairly easy to port into Python 3.x, when the time comes, as I have read several resources regarding ways to code in the Python 2.x line with the 3.x line in mind.

The README file includes the GNU GPL v3 license details for this code. The COPYING file includes citation details.

# Execution #

ssa.py is designed to interact with the user similar to a wizard. Once you have all of the prerequisites installed, all you have to do is double-click the ssa.py file in file explorer to get started. However, if you want to save a log of what the script reports, you will need to run the ssa.py script in an IDE (e.g. Python's default IDLE application). Otherwise, all of the things the script reports will disappear when the script finishes without letting you copy/save it to a log file.

Furthermore, every plot window that the script opens can be saved to a variety of image formats using the disk-like button at the bottom of each plot-window. You can also use other buttons on each plot window to zoom and pan into and around the plotted data, as well as see the coordinate under your mouse cursor. One of the eigen-spectrum plots instructs you to click two points on the graph to form a line, which it will use to select points (i.e. above the line). If you fail to select any points before closing the window, the script will select the first 20 ranked eigenvalues.  The last portion of the script allows you to save certain variables to csv-txt files, though it does not currently explain any of the variable names. I will post a legend that explains all of the available variables on this site soon.