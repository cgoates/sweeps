This project requires CMake 3.24 or newer.
In the build directory of choice (e.g. `<git root>/build`) run
```
cmake <path/to/git/root> -DCMAKE_BUILD_TYPE=Release
make -j10
```
Dependencies will be downloaded and built automatically.
To run tests, run `ctest -j10` in the build directory after building.

To use the python API, after building the codebase, include the following in your python script:

```
import sys
sys.path.insert(0, '<path/to/build/directory>/src/api')
import splines
```

Documentation for the splines api can be found in python using `help(splines)`.
