# Installing the FPR Code
If you want, you can skip section 2 and just go to the GitHub page and download the zip, and unzip it

## 1. Preparing Your Computer
### 1a. For Mac
This section sets up the Homebrew package manager and installs git. If you have both of these, skip to the next section

#### Install Homebrew
First thing you should do is install Homebrew. This is a package manager for macOS, which is pretty useful.
  - Type the following commands into the terminal (I recommend using [iTerm2](https://iterm2.com) instead of macOS's preinstalled terminal)
    ```bash
    xcode-select --install # this installs Xcode, macOS’s development tools
    ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)” # this actually installs Homebrew
    brew doctor # checks to make sure everything installed properly
    brew install caskroom/cask/brew-cask # installs Cask, which allows you to install packages which are not part of the official Homebrew repository
    ```

#### Install Git
Now let’s install Git. This is what you use to download the code from GitHub and keep them updated. Also install make, which is needed to compile the code
```bash
brew install git
brew install cmake # only needed for the vector version, as of 07/25/18
brew install make
```

### 1b. For Linux
  - Your Linux system already has a package manager (`apt` for Ubuntu/Debian, `rpm` for Red Hat/Cent OS), and likely already has git and make. If it does not, type:
    ```
    sudo apt-get update # just to make sure the repository information is up-to-date
    sudo apt-get install git make cmake
    ```
  - If your system already has these installed, it will just skip them (though it may be useful to upgrade them using `sudo apt-get upgrade $PKG`)

## 2. Installing the Code
### 2a. Get the files from GitHub
  - Move to the Desktop (or whatever directory you want to use as your working directory) with `cd ~/Desktop`
  - Download the files from the repository:
    ```
    git clone https://github.com/mjohn218/fpr_code.git # this clones the master branch of the FPR code. This could take a while
    ```
    if you want to use a branch other than `master`, use the `--branch $BRANCH` flag

### 2b. Compile the code
#### Non-vector version
  - Move to the working directory and make the executables:
    ```
    cd ~/Desktop/fpr_code/FPR_CELL
    make # this compiles. Don’t worry about the warnings you get. They’re due to the large number of unused variables in the old code. 
    ```

#### Vector version
  - The vector version of the code uses CMake to automatically make the Makefile:
    ```
    cd ~/Desktop/fpr_code/vec_vers
    mkdir build # make a scratch directory to hold some temporary files
    cd build
    cmake .. # use cmake to make the Makefile
    make # make the executables
    ```
    

## 3. Running the code
### 3a. Non-vector version
- The executable is now in `~/Desktop/fpr_code/FPR_CELL/bin`
  - To run, it needs several required files, **parms.inp** and **rxn_rev.inp**, as well as **info** files for each specie
    - **parms.inp**: contains most of the information the code needs, like the box size, number of molecules, reactions, timestep size, etc.
    - **rxn_rev.inp**: contains all the reactions. Columns are, in order from left to right, reaction #, type (B = binding/association, U = unbinding/disassociation), reactant #1 (interface), reactant #2 (interface), product (interface), coupled reaction (0 = no, 1 = yes), binding radius, rate
    - It will also need an **.info** file for each specie in the starting configuration, which should be added at the bottom of the **parms.inp** file. If they are not, the code will crash. This file contains all the information specific to that specie, such as its interfaces and their coordinates
    - See below ([INCLUDE SECTION]) for notes about file format.
  - There are some examples, in `~/Desktop/fpr_code/FPR_CELL/examples/NEW_INPUT`. The ones you likely want to use/look at are `/AA` and `/AB`, which are homonuclear and heteronuclear diatomic "reactions" of hard spheres. 
- To run a simulation, type make sure the executable and the two inp files are in the same directory, and type `./rd_gen_reweight_NOPBC_complex parms.inp rxn_rev.inp > out.dat&"`
  - the > out dumps the output into a file called out.dat and the ampersand (&) gives you the terminal back while it’s running
  - if you want to kill the job while it’s running, type `jobs` to check what is running in that terminal, then `jobs -9 %$ID`, where `$ID` is the integer ID of the job you want to kill (it's probably 1)

### 3b. Vector-version
Note: This is still a work in progress, and quite subject to change.
- The executable is now in `~/Desktop/fpr_code/vec_vers/build`
- The code needs several required files to run, which are provided by using command line flags
  - **input file**: This file includes information about the system itself, such as the box dimensions, requested number of iterations, timestep, and reaction information. Included using the `--parmfile` flag [NOTE: going to change this to `-f`].
  - **mol file**: These files include all information specific to the initial species, such as their interfaces and coordinates. They are searched for automatically when the code reads the input file
