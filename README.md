# Ulam-spirals
Makes Ulam spirals of a sequence.

![sample](/sample.png)

# introduction
I knew I had to make Ulam spirals when I saw [Numberphile’s video](https://www.youtube.com/watch?v=iFuR97YcSLM) about them. I wrote Python code to generate these where colors alternate with each concentric square. A visit to [a page of Matthew Conroy](https://www.madandmoonly.com/doctormatt/mathematics/ulamSpirals/ulamSpirals.htm) gave me the notion to map sequences other than primes too. Factorials, the Fibonacci sequence, and Mersenne primes make for sparse Ulam spirals, but I like them for that fact that they show which sequences are spread far out.

I made a few I found at the [OEIS browse page](https://oeis.org/Sbrowse.html) when I watched for sequences that would fit this type, which, as it turned out, Ulam spirals transmuted from sequences to QR codes. :P

* [A030513](https://oeis.org/A030513) Numbers with 4 divisors.
* [A037003](https://oeis.org/A037003) Positions of the digit 4 in the decimal expansion of 𝜋.
* [A050704](https://oeis.org/A050704) Composite numbers k with the property that k minus the sum of the prime factors of k is prime.

# packages
This requires a few packages.
* mpmath
* numpy
* PIL

# run the script
_Colors_ is a color pallette.

```
colors = [[[130, 106, 237],
           [200, 121, 255],
           [255, 183, 255],
           [202, 255, 138],
           [115, 169, 66],
           [83, 141, 34]]]
```

_Sides_ is the number pixels wide/tall for the PNG square. An odd number is required since the Ulam spiral requires a center pixel, which means the same number of pixels to the left, to the right, above, and underneath the center pixel.

```
sides = [1001]
```

_Modes_ is a bit more complicated for Ulam spirals in that it requares a function written to return True or False for whether a number is in a sequence or not.

```
modes = ['a050704']
```

I wrote the following modes for Ulam spirals. A030513, A050704, A037003 test for OEIS integer sequences.
* prime
* triangular
* square
* pentagonal
* hexagonal
* heptagonal
* octogonal
* nonagonal
* decagonal
* hendecagonal
* dodecagonal
* fibonacci
* factorial
* mersenne_prime
* a030513
* a050704
* a037003


An image subdirectory is created in the same folder as the script for the images. To run the code.

```
colors = [[[130, 106, 237],
           [200, 121, 255],
           [255, 183, 255],
           [202, 255, 138],
           [115, 169, 66],
           [83, 141, 34]]]

sides = [1001]

modes = ['a050704']

us = UlamSpiral(sides, modes, colors)
```

For rendering multiple color palettes, image resolutions, and modes.

```
colors = [[[0, 127, 95],
           [43, 147, 72],
           [85, 166, 48],
           [128, 185, 24],
           [170, 204, 0],
           [191, 210, 0],
           [212, 215, 0],
           [221, 223, 0],
           [238, 239, 32],
           [255, 255, 63]],

          [[194, 0, 251],
           [236, 8, 104],
           [252, 47, 0],
           [236, 125, 16],
           [255, 188, 10]]]

sides = [101, 401, 601, 1001]

modes = ['a030513', 'a037003']

us = UlamSpiral(sides, modes, colors)
```
