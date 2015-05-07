# Spinal Cord Recovery Analysis in Bitplane Imaris
Imaris Xtension to automate the analysis of zebrafish spinal cord recovery following injury to the nerve cord. In general this XTension requires live confocal or multiphoton images of a nuclei marker where Imaris spots analysis and tracking are using to generate tracks for the cells through time.

## General Requirements
* Imaris and XT Module
* Python 2.7.*

You will also need to install several packages (see below):

* numpy
* pandas
* matplotlib
* scipy
* seaborn

## Installation
You will need Python. On Windows I recommend installing the [Anaconda](http://continuum.io/downloads) distribution as it makes it easier to get all the dependencies. For *nix systems, you can generally use the distributed version of Python or install it using the your systems package manager. On OS X, I recommend [Homebrew](http://brew.sh).

If you are willing, I also highly recommend getting familiar with virtual environments for Python, whether that be [virtualenv](https://virtualenv.pypa.io/) or the Anaconda flavour. I personally prefer to set up a virtualenv with all the dependencies that I need for analysis and iPython and then use that virtual env python executable in Imaris (see below).

### Configure Imaris XT module
You will need to configure the XT module in Imaris to point to your Python executable:

1. Edit --> Preferences --> Custom Tools
2. Point to 'Python application' field to your Python installation. Important: If you created a virtualenv for Imaris, you need to point Imaris to that.

### Installing the SCR XTension
There are several ways one could do this, the aim is to get xtscr.py file into a place where Imaris can find it. Perhaps the easiest way would be to clone this repo using [git](http://git-scm.com):

```
cd ~/path/to/personal/imaris/xtensions
git clone https://github.com/keithschulze/xttgmmspotimport.git
```

or to download the zip from github and put in the place where you keep your Imaris XTensions (or create a place/folder if one doesn't exist). If you creat a new folder for extensions, then add that path into the 'XTension folders' box in the XT modules 'Custom Tools' config (see above).

Package requirements can be installed into your Python environment using Pip with the requirements file provided in this repository:

```
pip install -r requirements.txt
```

If you use a specific virtualenv for Imaris, make sure you activate this first before running the install.

## Usage
You will need to do 2 things before running the extension:

1. Run spot detection and tracking for cells in the nerve cord around the lesion site. Be sure to restrict the ROI of the spot detection to the nerve cord i.e., don't detect or track cells outside the nerve cord. You may need to manually remove cells outside the nerve cord.
2. Place 2 measurement points (mp; in the same measurement points object in Imaris). The first mp needs to be placed in approximately the centre of the lesion in the nerve cord. The second mp needs to be placed in the centre of the notochord, approximately perpendicular to the nerve cord/1st mp. This is used as a reference point for the of cell movement angle.

The extension can then be run by invoking the `Image Processing -> SCR` menuitem.
