# Legendre_functions

Code Information
================

For a maximum degree $N$, the following code computes the fully normalized Legendre functions $\bar{P}_{nm}(\cos(\theta))$ as a function of co-latitude $\theta$.

There are two source codes provided for the computation: one in C and one in Fortran.

The source codes have the following names:

- **geomalf.c**
- **geomalf.f90**

Code in C
---------

Before running the *geomalf.c* code, you must insert the values for the maximum degree $N$ and the co-latitude $\theta$ in lines **36** and **37**, respectively.
For example:

```c
N = 6; /* The maximum degree */
theta = 30.0; /* The co-latitude in degrees */
```

Enter the following command at the Linux command line to compile the code:

```bash
gcc -Wall geomalf.c -o geomalfc -lm
```

This command will build a single executable with the name *geomalfc*. Enter the following command in the Linux command line to execute the program:

```bash
./geomalfc
```

Below are the results:

```
N = 6
Co-latitude = 3.00000000e+01 deg
The number of Legendre Functions = 28

n         m           Pnm
0         0           1.0000000000000000e+00
1         0           1.5000000000000000e+00
1         1           8.6602540378443849e-01
2         0           1.3975424859373691e+00
2         1           1.6770509831248424e+00
2         2           4.8412291827592685e-01
3         0           8.5923294280422069e-01
3         1           2.2277546150777021e+00
3         2           1.1092649593311776e+00
3         3           2.6145625829189834e-01
4         0           7.0312500000000708e-02
4         1           2.3107045394749197e+00
4         2           1.7818666695701446e+00
4         3           6.7928328497763002e-01
4         4           1.3865811991639704e-01
5         0          -7.4051002865529225e-01
5         1           1.8565375211351955e+00
5         2           2.2993847894939750e+00
5         3           1.2465314425264324e+00
5         4           3.9826512815546256e-01
5         5           7.2712931519477580e-02
6         0          -1.3485606821315501e+00
6         1           9.5021287641141217e-01
6         2           2.4747031178290504e+00
6         3           1.8559287053259688e+00
6         4           8.1047568870385034e-01
6         5           2.2704605589840926e-01
6         6           3.7841009316401289e-02
```

---

Code in Fortran
---------------

Similarly, you need to enter the values for the co-latitude $\theta$ and maximum degree $N$ in lines **43** and **44**, respectively, in order to run the *geomalf.f90* code.

For instance:

```fortran
NX=6 ! The maximum degree, less than NMAX
theta=30.0d0 ! The co-latitude in degrees
```

In the Linux command line, enter the following command to compile this code:

```bash
gfortran -Wall geomalf.f90 -o geomalf90 -lm
```

One executable with the name *geomalf90* will be produced when you type this command.

In the Linux command line, enter the following command to execute the program:

```bash
./geomalf90
```

The outcomes are displayed below:

```
                               N=                   6
                      theta(deg)=      3.00000000E+01
The number of Legendre Functions=                  28
         n         m                        Pnm
         0         0     1.0000000000000000E+00
         1         0     1.5000000000000000E+00
         1         1     8.6602540378443849E-01
         2         0     1.3975424859373691E+00
         2         1     1.6770509831248424E+00
         2         2     4.8412291827592685E-01
         3         0     8.5923294280422069E-01
         3         1     2.2277546150777021E+00
         3         2     1.1092649593311774E+00
         3         3     2.6145625829189834E-01
         4         0     7.0312500000000708E-02
         4         1     2.3107045394749197E+00
         4         2     1.7818666695701446E+00
         4         3     6.7928328497763002E-01
         4         4     1.3865811991639704E-01
         5         0    -7.4051002865529225E-01
         5         1     1.8565375211351955E+00
         5         2     2.2993847894939745E+00
         5         3     1.2465314425264324E+00
         5         4     3.9826512815546267E-01
         5         5     7.2712931519477580E-02
         6         0    -1.3485606821315501E+00
         6         1     9.5021287641141217E-01
         6         2     2.4747031178290504E+00
         6         3     1.8559287053259688E+00
         6         4     8.1047568870385034E-01
         6         5     2.2704605589840926E-01
         6         6     3.7841009316401289E-02
```

*For individuals who are unfamiliar with the Linux command line, using an Integrated Development Environment (IDE) can facilitate the compilation and execution process.
