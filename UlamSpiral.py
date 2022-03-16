import math
import mpmath
import numpy as np
from PIL import Image
import os


class EvenDimensionError(Exception):
    # The required resolution for the image is a square with a center pixel that
    # has the same number of pixels to the left, to the right, above, and
    # underneath, which precludes any even number of pixels.
    pass

class NumberError(Exception):
    # For a test of prime numbers, to make sure the number is a non-negative
    # integer.
    pass

class UlamSpiral:


    '''
    Makes Ulam spirals of any arbitrary sequence. A True/False function is
    required for each sequence.
    '''


    SIDE = None # The length of the square.
    image_array = [] # Pixel color data.
    # If the PNG directory is not made, make it.
    IMAGE_DIR = os.path.dirname(os.path.realpath(__file__))
    IMAGE_DIR = os.path.join(IMAGE_DIR, 'PNG')
    if not os.path.exists(IMAGE_DIR):
        os.makedirs(IMAGE_DIR)
    COLORS = None # Color palettes.
    is_prime_list = [] # A list for primes so recalculating is not necessary.
    is_not_prime_list = [] # To avoid recalculation.
    image_size = 1600 # For scaling the final image.
    image_size = (image_size, image_size) # The final image resolution.
    DIGITS_OF_PI = None # To store digits of pi.


    def __init__(self,
                 sides=[41],
                 modes=['prime'],
                 colors=[[255, 255, 255]],
                 debug_tests=False):

        '''
        Make a basic introduction of the image essentails and specific
        initializations for sequences that require it.
        '''

        if 'a037003' in modes or debug_tests:
            self.bake_pi(max(sides)**2+1)

        if 'a050704' in modes:
            for x in range(max(sides)**2+1):
                is_prime = self.is_prime(x)
                if x%1000 == 0:
                    print('primes ' + str(x) + ' of ' + str(max(sides)**2))

        if debug_tests:
            self.debug_tests()

        for i, mode in enumerate(modes):
            self.COLORS = colors[i%len(colors)]
            for side in sides:
                if side%2 == 0:
                    raise EvenDimensionError

                self.SIDE = side
                self.CENTER = int((self.SIDE-1)/2+1)
                self.CENTER = [self.CENTER, self.CENTER]

                self.IMAGE_F = str(mode) + ' {:,}'.format(self.SIDE**2) + \
                               '.png'
                self.IMAGE_PATH = os.path.join(self.IMAGE_DIR, self.IMAGE_F)

                self.image_array = np.zeros((self.SIDE, self.SIDE, 3), \
                                       dtype=np.uint8)

                self.calc_pixels(mode)
                self.write_image(mode)


    def calc_pixels(self, mode):

        '''
        Follow the path requirements of an Ulam spiral. With each step, test the
        pixel for the sequence.
        '''

        cursor = self.CENTER.copy() # The center pixel.
        # For the color palettes, keep count the concentric squares. That number
        # decides, for a pixel who passes the sequence test, which color from
        # the palettes to choose.
        num_square = 0

        for x in range(1, self.SIDE**2+1):
            if x%100 == 0:
                print(str(mode) + ' {:,}'.format(x) + ' of ' + \
                      '{:,}'.format(self.SIDE**2))
            if x > 1:
                if cursor == [self.CENTER[0]+num_square,
                              self.CENTER[1]+num_square]:
                    cursor[1] += 1
                    num_square += 1
                elif cursor == [self.CENTER[0]+num_square,
                                self.CENTER[1]-num_square]:
                    cursor[1] += 1
                elif cursor == [self.CENTER[0]-num_square,
                                self.CENTER[1]-num_square]:
                    cursor[0] += 1
                elif cursor == [self.CENTER[0]-num_square,
                                self.CENTER[1]+num_square]:
                    cursor[1] -= 1
                elif cursor[1] == self.CENTER[0]+num_square:
                    cursor[0] -= 1
                elif cursor[0] == self.CENTER[1]-num_square:
                    cursor[1] -= 1
                elif cursor[1] == self.CENTER[0]-num_square:
                    cursor[0] += 1
                elif cursor[0] == self.CENTER[1]+num_square:
                    cursor[1] += 1

            self.test_pixel(cursor, num_square, x, mode)


    def test_pixel(self, cursor, num_square, x, mode):

        '''
        The hub for tests of a pixel's presence in a sequence. t_f is True or
        False depending on if it is or is not a part of the sequence. If true,
        change the pixel color.
        '''

        t_f = None
        if mode == 'prime':
            t_f = self.is_prime(x)
        elif mode == 'triangular':
            t_f = self.is_triangular(x)
        elif mode == 'square':
            t_f = self.is_square(x)
        elif mode == 'pentagonal':
            t_f = self.is_pentagonal(x)
        elif mode == 'hexagonal':
            t_f = self.is_hexagonal(x)
        elif mode == 'heptagonal':
            t_f = self.is_heptagonal(x)
        elif mode == 'octogonal':
            t_f = self.is_octogonal(x)
        elif mode == 'nonagonal':
            t_f = self.is_nonagonal(x)
        elif mode == 'decagonal':
            t_f = self.is_decagonal(x)
        elif mode == 'hendecagonal':
            t_f = self.is_hendecagonal(x)
        elif mode == 'dodecagonal':
            t_f = self.is_dodecagonal(x)
        elif mode == 'fibonacci':
            t_f = self.is_fibonacci(x)
        elif mode == 'factorial':
            t_f = self.is_factorial(x)
        elif mode == 'mersenne_prime':
            t_f = self.is_mersenne_prime(x)
        elif mode == 'a030513':
            t_f = self.is_a030513(x)
        elif mode == 'a050704':
            t_f = self.is_a050704(x)
        elif mode == 'a037003':
            t_f = self.is_a037003(x)

        if t_f:
            color = self.COLORS[num_square % len(self.COLORS)]
            self.image_array[cursor[0]-1, cursor[1]-1][0] = color[0]
            self.image_array[cursor[0]-1, cursor[1]-1][1] = color[1]
            self.image_array[cursor[0]-1, cursor[1]-1][2] = color[2]


    def write_image(self, mode):

        '''
        Write the finalized pixel color values to a PNG.
        '''

        image = Image.fromarray(self.image_array)
        image = image.resize(self.image_size, Image.NEAREST)
        image.save(self.IMAGE_PATH)


    def is_prime(self, x):

        '''
        Return True if x is prime and False otherwise.
        '''

        if x in self.is_prime_list:
            return True
        elif x in self.is_not_prime_list:
            return False

        if not isinstance(x, int) or x < 0:
            raise NumberError

        if x==0 or x == 1:
            self.is_not_prime_list.append(x)
            return False

        if x == 2:
            self.is_prime_list.append(x)
            return True

        for y in range(2, math.floor(math.sqrt(x))+1):
            if x%y == 0:
                self.is_not_prime_list.append(x)
                return False

        self.is_prime_list.append(x)
        return True


    def is_triangular(self, x):

        '''
        Return True if x is triangular and False otherwise.
        '''

        for y in range(1, x+1):
            t = (y*(y+1)) / 2
            if t == x:
                return True
            elif t > x:
                return False


    def is_square(self, x):

        '''
        Return True if x is square and False otherwise.
        '''

        for y in range(1, x+1):
            s = y**2
            if s == x:
                return True
            elif s > x:
                return False


    def is_pentagonal(self, x):


        '''
        Return True if x is pentagonal and False otherwise.
        '''

        for y in range(1, x+1):
            p = (y*(3*y-1)) / 2
            if p == x:
                return True
            elif p > x:
                return False


    def is_hexagonal(self, x):

        '''
        Return True if x is hexagonal and False otherwise.
        '''

        for y in range(1, x+1):
            h = y*(2*y-1)
            if h == x:
                return True
            elif h > x:
                return False


    def is_heptagonal(self, x):

        '''
        Return True if x is heptagonal and False otherwise.
        '''

        for y in range(1, x+1):
            h = (y*(5*y-3)) / 2
            if h == x:
                return True
            elif h > x:
                return False


    def is_octogonal(self, x):

        '''
        Return True if x ia octogonal an False otherwise.
        '''

        for y in range(1, x+1):
            o = y*(3*y-2)
            if o == x:
                return True
            elif o > x:
                return False


    def is_nonagonal(self, x):

        '''
        Return True if x is nonagonal and False otherwise.
        '''

        for y in range(1, x+1):
            n = (y*(7*y-5)) / 2
            if n == x:
                return True
            elif n > x:
                return False


    def is_decagonal(self, x):

        '''
        Return True if x is decagonal and False otherwise.
        '''

        for y in range(1, x+1):
            d = 4*y**2 - 3*y
            if d == x:
                return True
            elif d > x:
                return False


    def is_hendecagonal(self, x):

        '''
        Return True if x is hendecagonal and False otherwise.
        '''

        for y in range(1, x+1):
            h = (9*y**2 - 7*y) / 2
            if h == x:
                return True
            elif h > x:
                return False


    def is_dodecagonal(self, x):

        '''
        Return True if x is dodecagonal and False otherwise.
        '''

        for y in range(1, x+1):
            d = 5*y**2 - 4*y
            if d == x:
                return True
            elif d > x:
                return False


    def is_fibonacci(self, x):

        '''
        Return True for numbers in the Fibonacci sequence and False otherwise.
        '''

        f1 = 0
        f2 = 1
        while True:
            f = f1 + f2
            if f == x:
                return True
            elif f > x:
                return False
            else:
                f1 = f2
                f2 = f


    def is_factorial(self, x):

        '''
        Return True for factorials and False otherwise.
        '''

        for y in range(1, x+1):
            f = 1
            for z in reversed(range(1, y+1)):
                f = f*z
            if f == x:
                return True
            elif f > x:
                return False

    def is_mersenne_prime(self, x):

        '''
        Return True for Mersenne primes and False otherwise.
        '''

        for y in range(1, x+1):
            m = 2**y-1
            if m == x and self.is_prime(x):
                return True
            elif m > x:
                return False

    def is_a030513(self, x):
        '''
        Return True for numbers in A030513 and False otherwise.
        https://oeis.org/A030513
        Numbers with 4 divisors.
        '''

        divisors = []

        for y in range(1, x+1):
            if x%y == 0:
                divisors.append(y)
            if len(divisors) > 4:
                return False

        if len(divisors) == 4:
            return True
        else:
            return False

    def is_a050704(self, x):
        '''
        Return True for numbers in A050704 and False otherwise.
        https://oeis.org/A050704
        Composite numbers k with the property that k minus the sum of the
        prime factors of k is prime. 
        '''

        if x == 1:
            return False

        primes = []
        prime_factors = []

        if self.is_prime(x):
            primes.append(x)
        else:
            primes = [p for p in self.is_prime_list if p <= x/2]

        d = x
        for prime in primes:
            while True:
                if d%prime == 0:
                    prime_factors.append(prime)
                    d //= prime
                else:
                    break

        sum_of_prime_factors = 0

        for prime_factor in prime_factors:
            sum_of_prime_factors += prime_factor

        k_minus_sum_of_prime_factors = x - sum_of_prime_factors

        if self.is_prime(k_minus_sum_of_prime_factors):
            return True
        else:
            return False

    def is_a037003(self, x):
        '''
        Return True for numbers in A037003 and False otherwise.
        https://oeis.org/A037003
         Positions of the digit '4' in the decimal expansion of Pi.
        '''

        if self.DIGITS_OF_PI[x-1] == '4':
            return True
        else:
            return False

    def bake_pi(self, num_digits):

        '''
        Set the value of DIGITS_OF_PI to the decimal expansion of pi for any
        arbitrary length.
        '''

        mpmath.mp.dps = num_digits
        pi = mpmath.mp.pi
        self.DIGITS_OF_PI = str(pi)[2:]

    def debug_tests(self):

        '''
        Outputs a list of numbers in the sequences for verification.
        '''

        list_={'prime': [],
               'triangular': [],
               'square': [],
               'pentagonal': [],
               'hexagonal': [],
               'heptagonal': [],
               'hexagonal': [],
               'octogonal': [],
               'nonagonal': [],
               'decagonal': [],
               'hendecagonal': [],
               'dodecagonal': [],
               'fibonacci': [],
               'factorial': [],
               'mersenne_prime': [],
               'a030513': [],
               'a050704': [],
               'a037003': []}
        for x in range(1, 100):
            if self.is_prime(x):
                list_['prime'].append(x)
            if self.is_triangular(x):
                list_['triangular'].append(x)
            if self.is_square(x):
                list_['square'].append(x)
            if self.is_pentagonal(x):
                list_['pentagonal'].append(x)
            if self.is_hexagonal(x):
                list_['hexagonal'].append(x)
            if self.is_heptagonal(x):
                list_['heptagonal'].append(x)
            if self.is_octogonal(x):
                list_['octogonal'].append(x)
            if self.is_nonagonal(x):
                list_['nonagonal'].append(x)
            if self.is_decagonal(x):
                list_['decagonal'].append(x)
            if self.is_hendecagonal(x):
                list_['hendecagonal'].append(x)
            if self.is_dodecagonal(x):
                list_['dodecagonal'].append(x)
            if self.is_fibonacci(x):
                list_['fibonacci'].append(x)
            if self.is_factorial(x):
                list_['factorial'].append(x)
            if self.is_mersenne_prime(x):
                list_['mersenne_prime'].append(x)
            if self.is_a030513(x):
                list_['a030513'].append(x)
            if self.is_a050704(x):
                list_['a050704'].append(x)
            if self.is_a037003(x):
                list_['a037003'].append(x)
        input(list_)

if __name__ == '__main__':
    '''
    modes = ['prime',
             'triangular',
             'square',
             'pentagonal',
             'hexagonal',
             'heptagonal',
             'octogonal',
             'nonagonal',
             'decagonal',
             'hendecagonal',
             'dodecagonal',
             'fibonacci',
             'factorial']

    colors = [[[0, 95, 115],
               [10, 147, 150],
               [148, 210, 189],
               [233, 216, 166],
               [238, 155, 0],
               [202, 103, 2],
               [187, 62, 3],
               [174, 32, 18],
               [155, 34, 38]],

               [[247, 37, 133],
                [181, 23, 158],
                [114, 9, 183],
                [86, 11, 173],
                [72, 12, 168],
                [58, 12, 163],
                [63, 55, 201],
                [67, 97, 238],
                [72, 149, 239],
                [76, 201, 240]],

               [[249, 65, 68],
                [243, 114, 44],
                [248, 150, 30],
                [249, 132, 74],
                [249, 199, 79],
                [144, 190, 109],
                [67, 170, 139],
                [77, 144, 142],
                [87, 117, 144],
                [39, 125, 161]],

               [[234, 105, 139],
                [213, 93, 146],
                [192, 82, 153],
                [172, 70, 161],
                [151, 58, 168],
                [130, 47, 175]],

               [[227, 242, 253],
                [187, 222, 251],
                [144, 202, 249],
                [100, 181, 246],
                [66, 165, 245],
                [33, 150, 243]],

               [[204, 88, 3],
                [226, 113, 29],
                [255, 149, 5],
                [255, 182, 39],
                [255, 201, 113]],

               [[140, 179, 105],
                [244, 226, 133],
                [244, 162, 89],
                [91, 142, 125],
                [188, 75, 81]],

               [[88, 107, 164],
                [50, 67, 118],
                [245, 221, 144],
                [246, 142, 95],
                [247, 108, 94]],

               [[204, 213, 174],
                [233, 237, 201],
                [254, 250, 224],
                [250, 237, 205],
                [212, 163, 115]],

               [[34, 124, 157],
                [23, 195, 178],
                [255, 203, 119],
                [254, 249, 239],
                [254, 109, 115]]]
    '''

    '''
    colors = [[[253, 255, 252],
               [46, 196, 182],
               [231, 29, 54],
               [255, 159, 28]]]

    modes = ['mersenne_prime']
    '''

    '''
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
    '''

    colors = [[[130, 106, 237],
               [200, 121, 255],
               [255, 183, 255],
               [202, 255, 138],
               [115, 169, 66],
               [83, 141, 34]]]

    sides = [601, 1001]

    modes = ['a050704']

    us = UlamSpiral(sides, modes, colors, debug_tests=False)
