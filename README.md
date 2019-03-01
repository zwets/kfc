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

      # Example: codon count (with merged reverse complements)
      kfc -k 3 src/unit-test/data/test.fa.gz

      # Example: codon count (not merging reverse complements)
      kfc -k 3 -s src/unit-test/data/test.fa.gz

---

#### License

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

