#!/bin/sh
make maintainer-clean
./get_submodules.sh
./autogen.sh
./configure --prefix=/usr/local --with-osl=system --with-osl-prefix=/usr/local

#./configure --prefix=$HOME/usr --with-osl=bundled
make
