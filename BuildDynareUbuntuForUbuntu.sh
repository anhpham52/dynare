sudo add-apt-repository -s "deb http://archive.ubuntu.com/ubuntu xenial main restricted universe multiverse"
sudo add-apt-repository -s "deb http://archive.ubuntu.com/ubuntu xenial-updates main restricted universe multiverse"
sudo add-apt-repository -s "deb http://security.ubuntu.com/ubuntu xenial-security main restricted universe multiverse"
sudo apt-get update
sudo apt-get dist-upgrade
sudo apt-get upgrade
sudo apt-get build-dep dynare
sudo apt-get install git kdiff3 autoconf automake bison build-essential doxygen flex gettext gfortran git latex2html latex-beamer libboost-graph-dev libgsl0-dev liblapack3 liblapack-dev liblapacke libmatio-dev libslicot-dev libslicot-pic libsuitesparse-dev libtool make nano texi2html texinfo texlive texlive-extra-utils texlive-formats-extra texlive-latex-extra texlive-publishers unzip wget zip binutils-mingw-w64 binutils-mingw-w64-i686 binutils-mingw-w64-x86-64 g++-mingw-w64 g++-mingw-w64-i686 g++-mingw-w64-x86-64 gcc-mingw-w64 gcc-mingw-w64-base gcc-mingw-w64-i686 gcc-mingw-w64-x86-64 gdb-mingw-w64 gdb-mingw-w64-target gfortran-mingw-w64 gfortran-mingw-w64-i686 gfortran-mingw-w64-x86-64 mingw-w64 mingw-w64-common mingw-w64-i686-dev mingw-w64-tools mingw-w64-x86-64-dev texlive-fonts-extra matlab-support
sudo apt-get autoremove
#sudo mkdir /usr/local/MATLAB/R2016a/sys/os/glnxa64.bak/
#sudo mv /usr/local/MATLAB/R2016a/sys/os/glnxa64/* /usr/local/MATLAB/R2016a/sys/os/glnxa64.bak/
#sudo mv /usr/local/MATLAB/R2016a/sys/os/glnxa64.bak/libiomp5.so /usr/local/MATLAB/R2016a/sys/os/glnxa64/

cd ~/dynare/
git pull --recurse-submodules
git submodule foreach git pull --recurse-submodules
git fetch upstream --recurse-submodules
git submodule foreach git fetch upstream --recurse-submodules
git checkout master
git merge upstream/master
git commit -m "merge"
git push
autoreconf -sfi
./configure CFLAGS="-O3" CXXFLAGS="-O3" MATLAB_MEX_CFLAGS="-O3" MATLAB_MEX_CXXFLAGS="-O3" --with-matlab=/usr/local/MATLAB/R2017a MATLAB_VERSION=R2017a
make clean
make
make pdf
make html
