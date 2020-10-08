#!/usr/bin/env python3

'''
Python implementation of a Twisted Generalized Feedback Shift
Register (TGFSR). The TGFSR has a maximum sequence length of
2 ^ (word_size * register_length) - 1. However, the maximum
length will only be attained with careful selection of the
trinomial_index and twist_vector parameters. The values
hard-coded here yield a sequence length of (2 ^ 400) - 1, or
approximately (10 ^ 121.2) - 1.

Author: Patrick Kelly
Email: patrickyunen@gmail.com
Last Updated: 10-8-2020
'''

register_length = 25
word_size = 16
trinomial_index = 11
twist_vector = 0xa875
initial_values = [0xaf6e, 0xe08e, 0x0caf, 0x6f13, 0xe5a0, 0xb79b,
                  0xbeb4, 0x9f5c, 0x5082, 0xeeda, 0x3a57, 0xb1b6,
                  0x8767, 0x98d5, 0x3a32, 0x398a, 0x5e57, 0x982b,
                  0xb960, 0x4d25, 0x3d6c, 0xd68b, 0x283a, 0x324e,
                  0xf51d]


class TGFSR:

    def __init__(self, initial_value, m, A):
        self.n = len(initial_value) #register size
        self.state = initial_value #register state
        self.m = m #index of trinomial coefficient Cm
        self.l = 0 #index of current word
        self.a = A #twist vector

    def increment(self):
        if (self.state[self.l] & 1) == 0:
            self.state[self.l] = self.state[(self.l+self.m)%self.n] ^ (self.state[self.l] >> 1)
        else:
            self.state[self.l] = self.state[(self.l+self.m)%self.n] ^ (self.state[self.l] >> 1)\
                                 ^ self.a

        self.l = (self.l + 1) % self.n

    def output(self):
        return self.state[self.l]


x = TGFSR(initial_values, trinomial_index, twist_vector)
out = ''
for k in range(10000):
    out += '{0:04x}'.format(x.output())
    x.increment()

print(out)
