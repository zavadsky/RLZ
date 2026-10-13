Unified RLZ paper benchmark
===========================

This package reproduces Table 1, Table 2, and the data underlying Figure 1
of the DCC paper from one executable and one file list.

The executable intentionally exposes only these three reproducibility modes;
legacy hashing, codec, parser-comparison, and diagnostic modes were removed.

Input
-----

files.csv contains one input path per row. Relative paths are resolved relative
to files.csv. The seed is stored in the second column:

  File;seed;uniform/LMC segment length;L
  hp1.txt;the ;300;24
  bible.txt;the ;300;24
  dickens.txt;the ;300;24
  ecoli.txt;t;300;24

Missing values default to segment length 300 and L=24. Table 1 always uses
these paper settings.

The FM-index dictionary constructor reads its seed from column 2. Its length
correction is selected automatically: b=0 for texts below 1,000,000 bytes and
b=1 otherwise. A caller may override b and the seed through the third and
fourth arguments of build_suffix_array_dictionary.

Modes
-----

  ./rlz_benchmark files.csv

Creates table1.csv, table2.csv, and figure1_table.csv. All rows of files.csv
are processed by all three experiments.

  ./rlz_benchmark files.csv table1

Creates table1.csv. For each input, it finds the smallest EulerTigs dictionary
with k>=3 whose size exceeds 1.4% of the text. That exact size is then used by
Uniform, LMC, and our FM-index-based character-by-character method. LMC tests
ell=0, 0.5, and 1 and reports the variant with the longest average RLZ phrase.
Dictionary construction time excludes RLZ parsing, as in the paper.

  ./rlz_benchmark files.csv table2

Creates table2.csv. It chooses D=2^j with the smallest j satisfying
D>0.014|T|. It reports our Bit-Fix, Bit-VB, EF-Fix, and EF-VB sizes, Uniform
Bit-Fix and EF-VB sizes, and the complete Liao/Zlib scheme. The statistical
and Uniform dictionaries use the latest leftmost-offset Liao parser augmented
with a succinct RMQ.

  ./rlz_benchmark files.csv figure1

Creates figure1_table.csv. Each row is one curve and columns 8,...,4096 are
microseconds per query, averaged over 10,000 random positions. It reports our
fast Bit-Fix and compact EF-VB configurations, level-22 independently framed
Zstd matched to each size within 0.2%, and the closest 4 KiB*2^j block choice used for the
SquashFS comparison. The selected frame/block size, compressed size, and
compression ratio are included in every row.

The SquashFS comparison tests the two legal block sizes adjacent to the matched
Zstd Seekable frame size. If a target lies outside the size range achievable by
Seekable Zstd, the program reports that the target cannot be reached.

Figure 1 mode loads libzstd dynamically (libzstd.so.1), so no zstd header is
needed at compile time. Seek-table overhead is counted explicitly. SquashFS
rows use independently Zstd-compressed blocks and its restricted block sizes;
the compact filesystem metadata accounting is included in the reported size.

Build
-----

Requirements: GCC with C++20, SDSL, divsufsort/divsufsort64, zlib, Microsoft
GSL (`gsl/span`), and the Zstd runtime library. PASTA is included with the
program in the `pasta` directory and does not need to be installed separately.

On Ubuntu, install the required packages with:

  sudo apt update
  sudo apt install build-essential libsdsl-dev libdivsufsort-dev \
      zlib1g-dev libmsgsl-dev libzstd1

Here `libsdsl-dev` provides the SDSL headers and library,
`libdivsufsort-dev` provides both `libdivsufsort` and `libdivsufsort64`, and
`zlib1g-dev` provides the zlib headers and library. `libmsgsl-dev` is the
Microsoft Guidelines Support Library required by PASTA; it is unrelated to
GNU GSL.

If Microsoft GSL is copied locally instead of installed as a package, place
it so that the following file exists:

  third_party/GSL/include/gsl/span

The build script recognizes system installations and these local layouts:

  pasta/bit_vector/bit_vector.hpp
  third_party/bit_vector/include/pasta/...
  third_party/GSL/include/gsl/...

Run:

  chmod +x build.sh
  ./build.sh

All times are wall-clock measurements. Run on an otherwise idle machine using
the same CPU-frequency policy for every method.
