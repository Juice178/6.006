#!/usr/bin/env python2.7

import unittest
from dnaseqlib import *

### Utility classes ###

# Maps integer keys to a set of arbitrary values.
class Multidict:
    # Initializes a new multi-value dictionary, and adds any key-value
    # 2-tuples in the iterable sequence pairs to the data structure.
    def __init__(self):
        self.holder = dict()

    # Associates the value v with the key k.
    def put(self, k, v):
        self.holder.setdefault(k, [])
        self.holder[k].append(v)


    # Gets any values that have been associated with the key k; or, if
    # none have been, returns an empty sequence.
    def get(self, k):
        try:
            return self.holder[k]
        except:
            return []



# Given a sequence of nucleotides, return all k-length subsequences
# and their hashes.  (What else do you need to know about each
# subsequence?)


# Outer loop : for each nucleotide in the nucleotides, append the nucleotide
# to the list called holder.
# Case1): the length of the holder is less than k and is_True is True,
# then append the val to the holder
# Case2): If the length of the holder equals k,
# convert it to the string, yield the hash value of it
# and turn is_True into False .
# Case3): if is_True is False , set next_val to the val
# in the next iteration from the outerloop)
# append the element(val) to the holder and
# call slide function to compute hash value

def subsequenceHashes(seq, k):
    holder = []
    is_True = True
    index = 0
    count = 0


    for val in seq:
        # Case1): not enough nucoltides (less than k)
        if len(holder) < k and is_True:
            holder.append(val)

        if len(holder) == k and is_True:
            sub_seq = ''.join(holder)

            rh = RollingHash(sub_seq)
            hash_val = rh.current_hash()

            yield sub_seq, hash_val
            is_True = False
            current_val = val

        if is_True is False:
            for val in seq:
                next_val = val
                holder.append(next_val)
                prev = holder.pop(0)
                hash_val = rh.slide(prev, next_val)
                sub_seq = ''.join(holder)
                yield sub_seq, hash_val




# Similar to subsequenceHashes(), but returns one k-length subsequence
# every m nucleotides.  (This will be useful when you try to use two
# whole data files.)

# case 1): if the length of the holder does not equal k, append k
# to the holder/
# case 2): if the lenght of the holder is k and count does not equal m ,
# increment count and remove the first elemetn in the holder
# case 3): if the length of the holder is k and count equals m,
# convert the holder to a string and yield hash value and decrement count
def intervalSubsequenceHashes(seq, k, m):
    holder = []
    count = 0
    index = 0
    for val in seq:
        # case 1) length of holder is less than k
        if len(holder) < k:
            holder.append(val)

        if len(holder) == k:
            count += 1

        # case 2) length of holder equals k and count equals m
        if len(holder) == k and count == m:
            sub_seq = ''.join(holder)
            rh = RollingHash(sub_seq)
            hash_val = rh.current_hash()
            yield sub_seq, hash_val
            count = 0
            holder.pop(0)

        # case 3) length of holder equals k and count is not m
        if len(holder) == k and count < m:
            count += 1
            holder.pop(0)


# Searches for commonalities between sequences a and b by comparing
# subsequences of length k.  The sequences a and b should be iterators
# that return nucleotides.  The table is built by computing one hash
# every m nucleotides (for m >= k).
def getExactSubmatches(a, b, k, m):
    matches = []
    count = 0
    index_a = 0
    index_b = 0
    a_values = dict()
    b_values = dict()
    # for a_values and b_values, keys are hash values, values are indexes
    # for each key can exist multiple-values(a key holds multple values)
    for a_val in intervalSubsequenceHashes(a, k, m):
        a_values.setdefault(a_val[1], [])
        a_values[a_val[1]].append(index_a)
        index_a += 1
    for b_val in subsequenceHashes(b, k):
        b_values.setdefault(b_val[1], [])
        b_values[b_val[1]].append(index_b)
        index_b += 1



    # use a dictionary for effienct runnning time (hash)
    # if val in a_values also exists in b_valeus
    # append the key(index) to matches(a list)
    for val in a_values:
        if val in b_values:
            # since a list is unhashable , append the element to matches(a list)
                for i in a_values[val]:
                    for j in b_values[val]:
                        matches.append((i, j))
    return matches


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print 'Usage: {0} [file_a.fa] [file_b.fa] [output.png]'.format(sys.argv[0])
        sys.exit(1)

    # The arguments are, in order: 1) Your getExactSubmatches
    # function, 2) the filename to which the image should be written,
    # 3) a tuple giving the width and height of the image, 4) the
    # filename of sequence A, 5) the filename of sequence B, 6) k, the
    # subsequence size, and 7) m, the sampling interval for sequence
    # A.
    compareSequences(getExactSubmatches, sys.argv[3], (500,500), sys.argv[1], sys.argv[2], 8, 100)
