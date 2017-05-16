# MachUp

A Numerical Lifting-Line Algorithm

## Installation

There are currently no installation packages available for MachUp. The source 
code can be downloaded from the USU Aero Lab's GitHub page and built manually. 
See instructions below.

### Prerequisites

* MinGW (Windows only)
* gcc version 4.9 or higher
* cmake version 3.5 or higher

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

### Building the executable

#### Linux / Mac OSX

1. From the command prompt navigate to the MachUp directory
2. `cmake [-Dndv=<ndv>]`
3. `make`
4. If successful, the executable (MachUp.out) will be created in the MachUp/bin/ directory

#### Windows

1. From the command prompt navigate to the MachUp directory
2. `cmake -G "MinGW Makefiles" [-Dndv=<ndv>]`
3. `mingw32-make`
4. If successful, the executable (MachUp.out.exe) will be created in the MachUp/bin/ directory

### Testing the executable

1. From the command prompt navigate to MachUp/examples/FlyingWing/
2. Copy the executable from bin/ to MachUp/examples/FlyingWing/
3. `./MachUp.out input.json` (Linux / Mac OSX) or `MachUp.out input.json` (Windows)


