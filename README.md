# Geneva-analysis
**This software is highly beta version so many improvements have to be made**

This software is a graphical interface built with pyqt4. All plots are done with the pyqt4 interface of matplotlib.

## Installation
### dependencies
1. python 2.7
2. PyQt4 (any version should work)
3. numpy 1.9.2 or higher (should work with lower versions but not tested)
4. scipy 0.15.1 or higher(should work with lower versions but not tested)
5. matplotlib 1.4.3 or higher

The development and tests have been done with the anaconda python distribution on windows and mac. The software should run on both.
It has not been tested on linux.

### Running

#### On windows:
you can either run it from the command prompt:

``` 
>> cd folder_containing_main.py_source_file 
```

``` 
>> python main.py 
``` 

Or you can open the main.py python file in the IDE spyder and press the run button. 
I sugget to set the run configuration as "run in a external terminal"

#### On MacOs :
I suggest to use the first method since it seems that the "run in an external terminal" is not working on spyder for macos.

## Basic things to know to use the software

The Gui interface is organized in 2 main tabs:

1. Calibration
2. Analysis

### Calibration

To begin, in menu, go to File->select folder. A window will ask you to choose the path to the folder where 
are stored the images of the scan.

**Be sure there are an even number of images in the folder (pumped + ref). Indeed it happened that when 
redoing a scan with same number, new images are added and then it can end up with an odd number of images. 
Erase then the first image. You will be able after to know where the real scan begins with the software**

Now, the calibration tab should be initialized with scan images. I will let you playing with the sotware 
since nothing is really complicated. Please just report any bug you see or anything which could be improved.

The folder contain a /XUVSpec subfolder which contain a set of images which I did the test with. It can be 
usefull to use the same set for discovering the software since all the functions should be running.

Finally, the minimum thing you have to do to go to the analysis tab is the select an area on the Up-left image
and press *crop Image*, doing the calibration (energy and time delay) if you want, and then click on *Proceed to analysis*

### Analysis

The analysis consists in doing some normalization, reference removing and moving average. I'll let you play with it again,
it should be intuitive the way it works.

When finished, fill the date and Scan name boxes and then press *Save Analysis*. The current displayed graph is saved as
3 text files (time axis, enegy axis and scan matrix) in:

```
folder_containing_main.py/Datas/date/scan name/current/
```

All the other processed and raw scans are saved in:
```
folder_containing_main.py/Datas/date/scan name/other/
```

I hope it will be usefull for the present and future data analysis. I will continue to improve the code so, the best way of 
being always up to date is to clone the reposito ry



