#!/usr/bin/env python3

'''
MT19937_64 Python implementation of the 64-bit Mersenne Twister
Pseudorandom Number Generator with period 2^19937 - 1.

Author: Patrick Kelly
Email: patrickyunen@gmail.com
Last revised: October 5, 2020
'''

import sys

class MT19937_64:
    def __init__(self, seed):
        self.w = 64
        self.n = 340
        self.m = 156
        self.r = 59
        self.A = 0xb5026f5aa96619e9
        self.u = 29
        self.d = 0x5555555555555555
        self.s = 17
        self.b = 0x71d67fffeda60000
        self.t = 37
        self.c = 0xfff7eee000000000
        self.l = 43
        self.f = 6364136223846793005
        self.mask_64 = 0xffffffffffffffff
        self.upper_mask = 0xffffffff80000000
        self.lower_mask = 0x7fffffff

        # Create the state array
        self.MT = [0] * self.n

        # Initialize the state array
        self.index = self.n
        self.MT[0] = seed
        for i in range(1, self.n):
            temp = (self.f * (self.MT[i - 1] ^ (self.MT[i - 1] >> (self.w - 2))) + i)
            self.MT[i] = self.int_64(temp)

    def extract(self):
        # Get the next number in the state array.
        if self.index >= self.n:
            if self.index > self.n:
                raise InititionError('State array was not seeded')
                exit(0)
            self.twist() # Generate new array elements when needed
        return self.temper()

    def temper(self):
        # Tempering transformation for improved k-distribution
        y = self.MT[self.index]
        y ^= ((y >> self.u) & self.d)
        y ^= ((y << self.s) & self.b)
        y ^= ((y << self.t) & self.c)
        y ^= (y >> self.l)
        self.index = (self.index + 1) % self.n
        return self.int_64(y)

    def twist(self):
        for i in range(self.n):
            x = (self.MT[i] & self.upper_mask) + (self.MT[(i+1) % self.n] & self.lower_mask)
            xA = x >> 1
            if (x % 2) != 0:
                xA ^= self.A
            self.MT[i] = self.MT[(i + self.m) % self.n] ^ xA
        self.index = 0

    def int_64(self, number):
        # Truncates arbitrarily large number to 64 bits
        return number & self.mask_64


def main():
    length = int(sys.argv[1])
    seed = int(sys.argv[2])
    #mersenne = MT19937_64(143439545)
    mersenne = MT19937_64(seed)
    for i in range(length):
        number = mersenne.extract()
        normalized = '{:.16f}'.format(number/(2 ** 64))
        print(f"{'{:4d}'.format(i)}: {normalized}")

if __name__ == '__main__':
    main()
