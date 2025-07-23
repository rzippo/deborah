# DEBORAH: A Tool for Worst-Case Analysis of FIFO Tandems

```
@InProceedings{10.1007/978-3-642-16558-0_15,
author="Bisti, Luca and Lenzini, Luciano and Mingozzi, Enzo and Stea, Giovanni",
title="DEBORAH: A Tool for Worst-Case Analysis of FIFO Tandems",
booktitle="Leveraging Applications of Formal Methods, Verification, and Validation",
year="2010",
doi={https://doi.org/10.1007/978-3-642-16558-0_15}
}
```

This is a fork of the original code base to make it compilable with modern tools.
The build system is now CMake, and the code is compliant to the C++20 standard.

The features are mostly unchanged, including the tool interface.

Among the additions, the program can now output the output arrival curve related to the LUDB scenario, which is useful to integrate `deborah` into FIFO network studies.

## How to use

> TODO

### Edge-cases

`deborah` makes a few assumptions on the input, which may lead to incorrect results when they are not verified.
See [here](/scenarios/edge-cases/) for more details and examples.

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