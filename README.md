# MachUp

MuPy is a Python (2.? - 3.?) library for the design and analysis of 
fixed-wing aircraft. This includes things like calculating lift, 
drag, pitching moments, and stability derivatives. 

At the heart of MuPy is a modern numerical lifting-line algorithm that
rapidly predicts flow over multiple lifting surfaces and can 
incorperate viscous effects. For a detailed explanation of the theory 
refer to...

The following code demonstrates how MuPy might be used in a 
Python script:

```python
import mupy

#Generate a new airplane object
new_plane = mupy.Plane(inputs...)
#Add main wing
new_plane.addWing(inputs...)
#Add vertical tail
new_plane.addWing(inputs...)
#Add horizontal tail
new_plane.addWing(inputs...)

#Generate lifting line model for airplane
myModel = mupy.createLLModel(new_plane)

#Generate solution and store in results
results = myModel.solve()

#Access results
print(results.Lift_Coeff)
#Save .stl file of airplane for viewing in an stl viewer
new_plane.saveSTL()
```

## Features

*Easy user interface
*Fast
*Incorperates viscous effects
*Handles multiple lifting surfaces that have sweep, dihedral, and twist
*Additional libraries available for...

## Documentation

Documentation can be found at [machup.readthedocs.io](machup.readthedocs.io)
or by using the built in Python help() function to consult the docstrings. 

## Installation

MuPy packages are available on PyPi and Conda and can be installed 
using the following commands respectively. 

'pip install MuPy'

'conda install MuPy'

### Prerequisites

* Python version (2.? - 3.?)
* Scipy/Numpy version (2.? - 3.?)

### Getting the Source Code

The source code can be found at [https://github.com/usuaero/MachUp](https://github.com/usuaero/MachUp)

You can either download the source as a ZIP file and extract the contents, or 
clone the MachUp repository using Git. If your system does not already have a 
version of Git installed, you will not be able to use this second option unless 
you first download and install Git. If you are unsure, you can check by typing 
`git --version` into a command prompt.

#### Downloading source as a ZIP file

1. Open a web browser and navigate to [https://github.com/usuaero/MachUp](https://github.com/usuaero/MachUp)
2. Make sure the Branch is set to `Master`
3. Click the `Clone or download` button
4. Select `Download ZIP`
5. Extract the downloaded ZIP file to a local directory on your machine

#### Cloning the Github repository

1. From the command prompt navigate to the directory where MachUp will be installed
2. `git clone https://github.com/usuaero/MachUp`

##Support
Contact ????@usu.edu with any questions.

##License
This project is licensed under the ??? license. See LICENSE file for more information. 
