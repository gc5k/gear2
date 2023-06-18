# Genetic analysis repository (GEAR) for biobank-scale genetic data

## Getting Started

Check [GEAR wiki](https://github.com/gc5k/gear2/wiki) for kicking off.
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.


## Built With

* [Eigen](http://eigen.tuxfamily.org/) - The Linear algebra library for C++ (Eigen3)
* [Boost](http://boost.org) - Boost library for C++ （1.72.0 or above)
In addition, both eigen and boost should be installed first, and the cmake will automatically find the both libraries. Of note, some of Boost's library should be built first.

The following packages are required on a Linux machine to compile and use the software package.

```
g++
cmake
make
```

### Installing

Installing gear is fairly simple. Just issue the following commands on a Linux machine

```
git clone https://github.com/gc5k/gear2.git
cd gear2
mkdir build
cd build
cmake ..
make
```

## Acknowledgments
The mailman core is adapted from [ProPCA project](https://github.com/aman71197/ProPCA),

