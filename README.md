# kfc - k-mer frequency counter

Ultra-fast k-mer counter *UNDER CONSTRUCTION*

Home: https://github.com/zwets/kfc

 
## Quick Start

Assuming you are on a GNU/Linux system, and you have just cloned `kfc` from
GitHub, here are the steps to run it:

* Requirements

  To build `kcst` you need a C++ compiler and GNU `make`.  Run `c++ --version`
  and `make --version` to check that you have these.

  Support for gzipped files requires Boost Iostreams, available on Debian
  and Ubuntu as the `libboost-iostreams-dev` package.

* Build

      # Compile the khc source code
      cd src
      make

      # Run the unit tests
      make test

* Install

  There is no need to install `khc` in a specific place, but you may want to
  add it to the path:

      # For convenience in the examples below, put kfc on the path
      PATH="$PATH:$PWD"  # when in the ./src directory

  Alternatively, copy `kfc` to your `~/bin` or `~/.local/bin` directory.

* Run

      # kfc has self-contained usage instructions
      kfc --help

      # Example: codon count (canonical, ignoring strandedness)
      kfc -k 3 src/unit-test/data/test.fa.gz

      # Example: codon count (single stranded)
      kfc -k 3 -s src/unit-test/data/test.fa.gz


## FAQ

#### What are canonical k-mers?

By default, `kfc` treats k-mers and their reverse complement as the same k-mer.
For instance, if 3-mers `tgc` and `gca` both occur in the input, then they are
reported as two occurrences of the 3-mer `gca`.

This makes sense, as seeing either implies the presence of the other (unless
we're explicitly dealing with single-stranded DNA):


    <--cgt---           <--acg---
       |||     is just     |||     from a funny angle
    ---gca-->           ---tgc-->


Clearly, `tgc` and `gca` refer to the same _base pair sequence_, and the choice
to name this either `gca` or `tgc` is arbitrary.  `kfc` simply picks whichever
of the two has *a* or *c* as its middle base, and calls this the _'canonical'_
kmer.  (Note: the requirement that k be odd guarantees that there ís a middle
base!)


#### What are k-mer numbers?

`kfc` encodes each k-mer as an integral number.  @@TODO@@


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

