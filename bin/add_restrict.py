#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pipes
import shutil
import bisect
import sys
import click


# Modified from pairtools_restrict, fixed the issue with first fragment of each chromosome
# and this implementation is two times faster

@click.command()
@click.argument(
    'pairs_path',
    type=str,
    required=False)
@click.option(
    '-f', '--frags',
    type=str,
    required=True,
    help='a tab-diliminated BED file with the positions of restriction fragments '
         '(chrom, start, end). Can be generated using cooler digest.'
         'The file must be sorted according to the location on genome')
@click.option(
    '-o', "--output",
    type=str,
    default="",
    help='output .pairs/.pairsam file.')
@click.option(
    '-t', "--threads",
    type=int,
    default=3,
    show_default=True,
    help='Number of threads used by input decompressing.'
)
def restrict(pairs_path, frags, output, **kwargs):
    '''Assign restriction fragment information to pairs.

    Identify the restriction fragments that got ligated into a Hi-C molecule.
    The index and the start and the end of the fragment will be appended in the pairs

    PAIRS_PATH : input .pairs/.pairsam file. If the file ends with .gz, pbgzip is required.
    '''
    add_restrict_frag(pairs_path, frags, output, **kwargs)


def read_pairs(pairs_path, nproc=1):
    if pairs_path.endswith('.gz'):
        if shutil.which('pbgzip') is None:
            raise ValueError('pbgzip is not found, cannot decompress input'
                             )
        else:
            t = pipes.Template()
            t.append('pbgzip -dc -n {}'.format(nproc), '--')
            f = t.open(pairs_path, 'r')
        return f
    else:
        return open(pairs_path, 'r')


def split_header_body(instream, comment_char='#'):
    header = []
    comment_byte = comment_char.encode()
    # get peekable buffer for the instream
    inbuffer = instream.buffer
    current_peek = inbuffer.peek()
    while current_peek.startswith(comment_byte):
        # consuming a line from buffer guarantees
        # that the remainder of the buffer starts
        # with the beginning of the line.
        line = inbuffer.readline()
        # append line to header, since it does start with header
        header.append(line.decode().strip())
        # peek into the remainder of the instream
        current_peek = inbuffer.peek()
    # apparently, next line does not start with the comment
    # return header and the instream, advanced to the beginning of the data
    return header, instream


def read_frags(frags):
    rfrags_dict = dict()
    # frag_np = numpy.loadtxt(frags, dtype=None, delimiter='\t',
    #                         comments='#', dtype=None, usecols=(0, 2),
    #                         names=['chrom', 'end'])

    with open(frags, "r") as frag_fh:
        for line in frag_fh:
            chrom, _start, end = line.rstrip().split("\t")
            # sizes file is 0-based, should be converted to 1-based
            end = int(end) + 1
            if chrom not in rfrags_dict:
                rfrags_dict[chrom] = [end]
            else:
                bisect.insort(rfrags_dict[chrom], end)
    # for chrom in rfrags_dict:
    #     rfrags_dict[chrom] = np.array(rfrags_dict[chrom])

    return rfrags_dict


def add_restrict_frag(pairs_path, frags, output, **kwargs):
    instream = (read_pairs(pairs_path,
                           nproc=kwargs.get('threads'))
                if pairs_path else sys.stdin)

    outstream = (open(output, "w")
                 if output else sys.stdout)

    rfrags = read_frags(frags)

    header, body = split_header_body(instream)
    header.append(
        "# The indice and locations of restriction fragments were appended in the end of pairs")
    outstream.writelines((line + '\n' for line in header))

    for line in body:
        cols = line.rstrip().split("\t")
        chrom1, pos1 = cols[1], int(
            cols[2])
        rfrag_idx1, rfrag_start1, rfrag_end1 = locate_frag(
            rfrags, chrom1, pos1)
        chrom2, pos2 = cols[3], int(
            cols[4])
        rfrag_idx2, rfrag_start2, rfrag_end2 = locate_frag(
            rfrags, chrom2, pos2)
        cols += [str(rfrag_idx1), str(rfrag_start1), str(rfrag_end1)]
        cols += [str(rfrag_idx2), str(rfrag_start2), str(rfrag_end2)]
        outstream.write("\t".join(cols) + "\n")
    if pairs_path:
        instream.close()
    if output:
        outstream.close()


# runtime 648.88s user 1227.64s system 714% cpu 4:22.76 total
def bsearch(left, right, target, nums):
    p = (left + right) // 2
    if target < nums[0]:
        return 0
    elif target >= nums[-1]:
        return len(nums)
    elif target == nums[p]:
        return p + 1
    elif (target < nums[p] and target > nums[p - 1]):
        return p
    elif (target > nums[p] and target < nums[p + 1]):
        return p + 1
    elif target > nums[p]:
        return bsearch(p, right, target, nums)
    else:
        return bsearch(left, p, target, nums)


# runtime 557.13s user 1059.06s system 704% cpu 3:49.51 total
def binsearch(sort_list, val):
    low = 0
    high = len(sort_list)
    while True:
        mid = (low + high) // 2
        if val >= sort_list[-1]:
            return high
        elif val < sort_list[0]:
            return 0
        elif (sort_list[mid] > val and sort_list[mid - 1] <= val):
            return mid
        elif sort_list[mid] > val and sort_list[mid - 1] > val:
            high = mid
        else:
            low = mid


def locate_frag(rfrags_dict, chrom, pos):
    rfrags_list = rfrags_dict[chrom]
    # runtime 1025.49s user 2132.81s system 687% cpu 7:39.58 total
    # idx = np.searchsorted(rfrags_list, pos, side='right')
    # runtime 562.33s user 1075.73s system 671% cpu 4:04.09 total
    # idx = bisect.bisect_right(rfrags_list, pos)
    #idx = bsearch(0, len(rfrags_list), pos, rfrags_list)
    idx = binsearch(rfrags_list, pos)
    if idx == 0:
        return 0, 0, rfrags_list[0]
    else:
        return idx, rfrags_list[idx - 1], rfrags_list[idx]


if __name__ == '__main__':
    restrict()
