# Deborah

This fork includes many changes to the original code base, but the features are mostly unchanged.
The build system is now CMake, and the code is compliant to the C++20 standard.

Among the additions, the program can now output the output arrival curve related to the LUDB scenario, which is useful to integrate deborah into FIFO network studies.

## Reminder on CMake

On usual configurations, to build compile this program you would need to use the following commands

```
mkdir build
cd build
cmake ..
make
```

On Windows, in the last step run `cmake --build .` instead of `make`.

Alternatively, most IDEs will allow you to open the project via the `CMakeLists.txt` and take care of the build process themselves.