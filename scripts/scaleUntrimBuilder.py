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
import io

def findOS():
    try:
        OS = platform.system().lower()
        return OS
    except:
        print("Error finding OS.")
        sys.exit(-1)

def findPackageManager(os):
    os = os.lower()
    if os == "darwin":
        if checkSystemToolIsInstalled("brew"):
            return "brew"
    if os == "linux":
        if checkSystemToolIsInstalled("apt"):
            return "apt"
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
        if manager == "apt":
            command = [manager, "install", "-y", package]
        else:
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
            installLibrary(manager, tool)
            print(f"{tool} has been installed.")

def moveDirectoryFromSweeps(source, target, sourceFromSweeps=True):
    """
    Move a directory from a source to a target within the sweeps repository. 
    The path to be used as an argument should be written starting from the sweeps repository.
    For example, to move a folder in the sweeps/scripts directory to the sweeps/deps directory,
    the source would be "scripts/exampleFolder" and the target would be "deps/exampleFolder".
    source: str - the source folder that is being moved.
    target: str - the target folder where the source is being moved to.
    """
    # Create the paths by finding the path to sweeps and appending the arguments to it.
    sweepsPath = Path(__file__).resolve().parent.parent
    if sourceFromSweeps:
        source = Path(sweepsPath, source)
    target = Path(sweepsPath, target)
    oldPath = target / "ScaleUntrim"
    if oldPath.exists():
        try:
            shutil.rmtree(oldPath)
        except:
            print("Could not find the outdated file to remove.")
            raise(FileNotFoundError)
    runCommand(["mv", str(source), str(target)])
    return

def writeConfigFile():
    """
    Write the necessary config file so the built scaleUntrim library can be immediately run.
    """
    sweeps_path = Path(__file__).resolve().parent.parent
    config_path = sweeps_path / "deps" / "ScaleUntrim" / "setting.config"
    f = open(config_path, "r")
    lines = f.readlines()
    f.close()
    temp_dir_path = sweeps_path / "deps" / "ScaleUntrim" / "build" / "tempDir"
    f = open(config_path, "w")
    first_line = "temp_dir: " + str(temp_dir_path) + "/\n"
    lines[0] = first_line
    output = ""
    for line in lines:
        output += line
    f.write(output)
    f.close()

def buildScaleUntrim():
    os = findOS()
    manager = findPackageManager(os)
    if not manager:
        print("No package manager found, please install homebrew or chocolatey to use this script.")
        return 
    if manager == "apt":
        runCommand(["apt-get", "update"])
    installSystemTools(manager, ["cmake", "git", "make"])
    if os == "linux":
        openCascade = ["libocct-foundation-dev", "libocct-modeling-data-dev", "libocct-modeling-algorithms-dev", "libocct-ocaf-dev", "libocct-data-exchange-dev", "libocct-visualization-dev", "libocct-draw-dev"]
        installLibraries(manager, ["libeigen3-dev", "libboost-all-dev", "libtbb-dev"] + openCascade)
    else:
        installLibraries(manager, ["eigen", "boost", "OpenCascade"])
    cloneRepository("https://github.com/colbyj427/edited-scale-untrim.git", "ScaleUntrim")
    makeDirectory("ScaleUntrim/build")
    runCommand(["cmake", "-B", "ScaleUntrim/build", "-S", "ScaleUntrim"])
    print("DEBUG: RUNNING MAKE")
    runCommand(["make", "-C", "ScaleUntrim/build"])
    print("DEBUG: MAKING TEMP DIR")
    makeDirectory("ScaleUntrim/build/tempDir")
    print("DEBUG: MOVING TO DEPS")
    moveDirectoryFromSweeps("ScaleUntrim", "deps/", sourceFromSweeps=False)
    print("DEBUG: WRITING CONFIG")
    writeConfigFile()
    print("ScaleUntrim has been placed in the sweeps/deps folder and built successfully.")
    print("To run, update the setting.config file in the ScaleUntrim folder:\n1. Change the temp_dir to match your full working path to the tempDir directory in the build folder (This will be where the quad mesh is output).\n2. Change the magnitude value to your desired value, likely between 1 and 2, such as 1.4")
    print("3. Run the executable from the build folder with the command:\n  ./quadriflow <inputFilePath> ../setting.config")
    print("The resulting quadrilateral mesh is stored in build/tempDir/quad.vtk")

def runExample():
    """
    Run an example of the ScaleUntrim program.
    """
    sweepsPath = Path(__file__).resolve().parent.parent
    quadriflowPath = Path(sweepsPath, "deps", "ScaleUntrim", "build", "QuadriFlow")
    hookPath = Path(sweepsPath, "deps", "ScaleUntrim", "hookBase.obj")
    configPath = Path(sweepsPath, "deps", "ScaleUntrim", "setting.config")
    print("Running example...")
    runCommand([str(quadriflowPath), str(hookPath), str(configPath)])
    print("Example finished.")
    print("The resulting quadrilateral mesh is stored in ScaleUntrim/build/tempDir/quad.vtk")

buildScaleUntrim()
# runExample()