sudo apt-get update
sudo apt-get dist-upgrade
sudo apt-get upgrade
sudo apt-get build-dep dynare
sudo apt-get install git kdiff3 autoconf automake bison build-essential doxygen flex gettext gfortran git latex2html latex-beamer libboost-graph-dev libgsl0-dev liblapack3 liblapack-dev liblapacke libmatio-dev libslicot-dev libslicot-pic libsuitesparse-dev libtool make nano texi2html texinfo texlive texlive-extra-utils texlive-formats-extra texlive-latex-extra texlive-publishers unzip wget zip
sudo apt-get install binutils-mingw-w64 binutils-mingw-w64-i686 binutils-mingw-w64-x86-64 g++-mingw-w64 g++-mingw-w64-i686 g++-mingw-w64-x86-64 gcc-mingw-w64 gcc-mingw-w64-base gcc-mingw-w64-i686 gcc-mingw-w64-x86-64 gdb-mingw-w64 gdb-mingw-w64-target gfortran-mingw-w64 gfortran-mingw-w64-i686 gfortran-mingw-w64-x86-64 mingw-w64 mingw-w64-common mingw-w64-i686-dev mingw-w64-tools mingw-w64-x86-64-dev
sudo apt-get autoremove
git config --global merge.tool kdiff3
sudo mkdir -p /usr/local/lib/mingw64
cd /usr/local/lib/mingw64
sudo wget http://www.dynare.org/build/dynare-mingw64-libs.zip
sudo unzip dynare-mingw64-libs.zip
sudo rm dynare-mingw64-libs.zip
cd ~
git clone http://git.savannah.gnu.org/r/gsl.git
cd ~/gsl
git pull
autoreconf -si
./configure CFLAGS="-O3" CXXFLAGS="-O3" --host=x86_64-w64-mingw32 --prefix=/usr/local/lib/mingw64/gsl
make clean
make
sudo make install
cd ~/dynare/
git pull --recurse-submodules
git fetch upstream --recurse-submodules
git checkout master
git merge upstream/master
git mergetool
git commit -m "merge"
git push
autoreconf -si
./configure CFLAGS="-O3 -DCUDA=1" CXXFLAGS="-O3 -DCUDA=1" --host=x86_64-w64-mingw32 --enable-openmp --with-boost=/usr/local/lib/mingw64/boost --with-blas=/usr/local/lib/mingw64/blas/libopenblas.a --with-lapack=/usr/local/lib/mingw64/lapack/liblapack.a --with-gsl=/usr/local/lib/mingw64/gsl --with-matio=/usr/local/lib/mingw64/matio --with-slicot=/usr/local/lib/mingw64/slicot --with-matlab=/mnt/host/cygdrive/c/Progra~1/MATLAB/R2016a MATLAB_VERSION=R2016a --disable-octave
make clean
make all html
