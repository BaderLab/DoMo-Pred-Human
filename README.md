# DoMo-Pred-Human

The DoMo-Pred command line tool is implemented using Python 2.7 and C++. 

It is available for download under the GNU LGPL license from: http://www.baderlab.org/Software/DoMo-Pred

## Requirements 
* Python == 2.7
* Cython == 0.22
* setuptools == 17.1.1
* numpy == 1.9.2
* scipy == 0.16.0
* nwalign == 0.3.1 (provided with source code) (https://pypi.python.org/pypi/nwalign/?)
* networkx == 1.11
* sklearn == 0.0
* dill

## USAGE
```
python2.7 run_pwm.py 
```
* (set input or output paths in this file)


## COMPILE/INSTALL

1. Most of the code is written in Python and does not require any compilation or installation.
1. PWM scanning module (Peptide/PWMsearch/source/) is written in C++ to speed up the scanning process. This modules needs to be complied separately in order to be used by the pipeline.
1. Structural contact peptide feature uses nwalign python package for Needleman-Wunsch global sequence alignment.
This package needs to be installed and can be found here: Peptide/Structure/source/NW/ or https://pypi.python.org/pypi/nwalign/?
  * To install nwalign use the follwoing command: python2.7 setup.py 

### For compiling on Ubuntu:
1. g++ -fPIC -c PWMSeaech.cpp -I/usr/include/python2.7 -lpython2.7 -o PWMSearch.o
1. g++ -shared PWMSearch.o -o PWMSearch.so

### For compiling on Mac OS:
1.  g++ -fPIC -c PWMSearch.cpp -I/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7 -o PWMSearch.o
1.  g++ -shared -framework Python PWMSearch.o -o PWMSearch.so

### Using Macports version of Python make sure to use Macport Python library:
1.  g++-mp-4.8 -fPIC -c PWMSearch.cpp -I/opt/local/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7/ -o PWMSearch.o
1. Step2: g++-mp-4.8 -shared -L/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/config -ldl -framework CoreFoundation -lpython2.7 PWMSearch.o -o PWMSearch.so


