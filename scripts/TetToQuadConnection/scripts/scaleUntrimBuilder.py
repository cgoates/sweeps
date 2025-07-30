"""
This script is designed to successfully download and build the 
ScaleUntrim library, which can be found at this link:
https://github.com/colbyj427/edited-scale-untrim.git
"""

import platform
import shutil
import subprocess
from pathlib import Path
import sys

def findOS():
    try:
        OS = platform.system().lower()
        return OS
    except:
        print("Error finding OS.")
        sys.exit(-1)

def findPackageManager(os):
    os = os.lower()
    if os == "darwin" or os == "linux":
        if checkSystemToolIsInstalled("brew"):
            return "brew"
    elif os == "windows":
        if checkSystemToolIsInstalled("choco"):
            return "choco"
    return None

def checkSystemToolIsInstalled(dep):
    """
    This works for a program that creates output
    when called on terminal, such as git, cmake, or homebrew.
    dep: str - the name of the tool to check.
    """
    if shutil.which(dep):
        return True
    return False

def checkLibraryIsInstalled(manager, library):
    """
    Check if a library is installed.
    library: str - the name of the library to check.
    """
    if shutil.which(manager):
        command = [manager, "list", library]
        try:
            result = subprocess.run(command, check=True, capture_output=True, text=True)
            if result.returncode == 0:
                return True
            else:
                return False
        except subprocess.CalledProcessError:
            return False
    return "Given package manager is not installed."

def runCommand(command):
    """
    Run a command in the terminal. 
    Print any output.
    Command: List[str] - each item in the list is  an argument in the command
    """
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e.cmd}")
        print(f"Return code: {e.returncode}")
        print(f"OUTPUT: {e.output}")
        print(f"ERROR OUTPUT: {e.stderr}\n")
        sys.exit(-1)
    except:
        print(f"Error running command: {command}")
        sys.exit(-1)

def makeDirectory(path):
    """
    Make a directory at the given path.
    path: str - the path to the directory to create.
    """
    path = Path(path)
    path.mkdir()

def cloneRepository(repoUrl, targetDirectory=None):
    """
    Clone a github repository to the given target, or to the current directory if no target is specified.
    repoUrl: str - the url of the repo.
    targetDirectory: str - the target direcotry to clone the repo to.
    """
    command = ["git", "clone", repoUrl]
    if targetDirectory:
        command.append(targetDirectory)
    try:
        runCommand(command)
    except:
        print(f"Error running command: {command}")
        sys.exit(-1)

def installLibrary(manager, package):
    """
    Install a library using homebrew or another package manager. Error catching is implemented in the runCommand function.
    package: str - the name of the library that the package manager will recognize.
    """
    if checkSystemToolIsInstalled(manager):
        
        command = [manager, "install", package]
        try:
            runCommand(command)
        except:
            print(f"Error installing package {package} with {manager}.")
            sys.exit(-1)

def installLibraries(manager, libraries):
    """
    Install a list of packages using the given package manager.
    libraries: List[str] - a list of libraries to install.
    """
    for library in libraries:
        if not checkLibraryIsInstalled(manager, library):
            installLibrary(manager, library)
            print(f"{library} has been installed.")

def installSystemTools(manager, tools):
    """
    Check the list of tools input using the checkSystemToolIsInstalled function.
    Then install the tool if the user says yes.
    tools: List[str] - a list of dependencies to check
    """
    for tool in tools:
        if not checkSystemToolIsInstalled(tool):
            print(f"{tool} is not installed. It will be installed now.")
            usrIn = input("Enter y to continue or n to cancel program: ")
            if usrIn == "n":
                print("Exiting...")
                sys.exit(0)
            installLibrary(manager, tool)
            print(f"{tool} has been installed.")

def buildScaleUntrim():
    os = findOS()
    manager = findPackageManager(os)
    if not manager:
        print("No package manager found, please install homebrew or chocolatey to use this script.")
        return 
    # print("DEBUG: CHECKING DEPENDENCIES")
    installSystemTools(manager, ["cmake", "git", "make"])
    # print("DEBUG: FINISHED CHECKING DEPENDENCIES")
    installLibraries(manager, ["eigen", "boost", "OpenCascade"])
    # print("DEBUG: ABOUT TO BUILD")
    cloneRepository("https://github.com/colbyj427/edited-scale-untrim.git", "ScaleUntrim")
    makeDirectory("ScaleUntrim/build")
    runCommand(["cmake", "-B", "ScaleUntrim/build", "-S", "ScaleUntrim"])
    runCommand(["make", "-C", "ScaleUntrim/build"])
    makeDirectory("ScaleUntrim/build/tempDir")
    print("ScaleUntrim has been built successfully.")
    print("To run, update the setting.config file in the ScaleUntrim folder:\n1. Change the temp_dir to match your full working path to the tempDir directory in the build folder (This will be where the quad mesh is output).\n2. Change the magnitude value to your desired value, likely between 1 and 2, such as 1.4")
    print("3. Run the executable from the build folder with the command:\n  ./quadriflow <inputFilePath> <outputFilePath> ../setting.config")
    print("Currently the value entered for the output file doesn't change anything.")
    print("The resulting quadrilateral mesh is stored in build/tempDir/quad.vtk")

def runExample():
    """
    Run an example of the ScaleUntrim program.
    """
    print("Running example...")
    runCommand(["./ScaleUntrim/build/quadriflow", "./ScaleUntrim/hookBase.obj", "./ScaleUntrim/setting.config"])
    print("Example finished.")
    print("The resulting quadrilateral mesh is stored in ScaleUntrim/build/tempDir/quad.vtk")

buildScaleUntrim()