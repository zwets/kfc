# kfc - k-mer frequency counter

Ultra-fast k-mer counter

Home: https://github.com/zwets/kfc

 
## Quick Start

Assuming you are on a GNU/Linux system, and you have just cloned `kfc` from
GitHub, here are the steps to run it:

* Requirements

  - To build `kfc` you need a C++ compiler and GNU `make`.  Run `c++ --version`
  and `make --version` to check that you have these.

  - Built-in support for gzipped files requires Boost Iostreams, available on
  Debian/Ubuntu as the `libboost-iostreams-dev` package.  Note though that
  until `kfc` is multi-threaded, using the shell's inherent parallelism is
  naturally more efficient: `gunzip -c file.fa.gz | kfc`.


* Build

      # Compile the kfc source code
      cd src
      make

      # Run the unit tests
      make test

* Install

  There is no need to install `kfc` in a specific place, but you may want to
  add it to the path:

      # For convenience in the examples below, put kfc on the path
      PATH="$PATH:$PWD"  # when in the ./src directory

  Alternatively, copy or symlink `kfc` to your `~/bin` or `~/.local/bin` directory.

      ln -sft ~/bin "$PWD/kfc"     # when in the ./src directory

* Run

      # kfc has self-contained usage instructions
      kfc --help

      # Example: codon count (canonical, ignoring strandedness)
      kfc -k 3 src/unit-test/data/ecoli.fa.gz

      # Example: codon count (single stranded)
      kfc -k 3 -s src/unit-test/data/ecoli.fa.gz


## FAQ

#### What are canonical (double-stranded) k-mers?

By default, `kfc` treats k-mers and their reverse complement as the same k-mer.
For instance, if 3-mers `tgc` and `gca` both occur in the input, then they are
reported as two occurrences of the 3-mer `gca`.

This makes sense, as seeing either implies the presence of the other (unless
we're explicitly dealing with single-stranded DNA):


    <--cgt---           <--acg---
       |||     is just     |||     from a funny angle
    ---gca-->           ---tgc-->


Clearly, `tgc` and `gca` refer to the same base pair sequence, and the choice
of name (`gca` or `tgc`) is arbitrary.  `kfc` simply picks whichever of the
two has *a* or *c* as its middle base, and calls the the 'canonical' kmer.
The requirement that k is odd guarantees that there ís a middle base.

If you do not want this 'destranded' behaviour, use option `-s` to treat the
input as single stranded DNA.  In this mode, only the k-mers which are
literally seen in the input are reported.  Note that this makes the k-mer
space twice as large as in the default mode.  In this mode k need not be odd.


#### What are k-mer numbers (s-code, c-code)?

`kfc` encodes each k-mer as a integral number.  For the sake of naming things,
kfc calls this number an s-code (when in single stranded mode), or c-code (in
default canonical mode).

The s-code is the unsigned integral number obtained by encoding each base of
the k-mer as two bits (a = 00, c = 01, g = 10, t = 11).  The c-code is the
same except that the middle base (which must be a=00 or c=01, see above) is
encoded as a single bit.


---

#### Licence

kfc - k-mer frequency counter  
Copyright (C) 2019  Marco van Zwetselaar <io@zwets.it>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

