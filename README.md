In the build directory of choice (e.g. `<git root>/build`) run `cmake path/to/git/root -DCMAKE_BUILD_TYPE=Release`, then `make -j10` or similar.
Dependencies will be downloaded and built automatically.
To run tests, run `ctest -j10` in the build directory after building.
