sudo add-apt-repository -s "deb http://archive.ubuntu.com/ubuntu xenial main restricted universe multiverse"
sudo add-apt-repository -s "deb http://archive.ubuntu.com/ubuntu xenial-updates main restricted universe multiverse"
sudo add-apt-repository -s "deb http://security.ubuntu.com/ubuntu xenial-security main restricted universe multiverse"
sudo apt-get update
sudo apt-get dist-upgrade
sudo apt-get upgrade
sudo apt-get build-dep dynare
sudo apt-get install git kdiff3 autoconf automake bison build-essential doxygen flex gettext gfortran git latex2html latex-beamer libboost-graph-dev libgsl0-dev liblapack3 liblapack-dev liblapacke libmatio-dev libslicot-dev libslicot-pic libsuitesparse-dev libtool make nano texi2html texinfo texlive texlive-extra-utils texlive-formats-extra texlive-latex-extra texlive-publishers unzip wget zip
sudo apt-get install binutils-mingw-w64 binutils-mingw-w64-i686 binutils-mingw-w64-x86-64 g++-mingw-w64 g++-mingw-w64-i686 g++-mingw-w64-x86-64 gcc-mingw-w64 gcc-mingw-w64-base gcc-mingw-w64-i686 gcc-mingw-w64-x86-64 gdb-mingw-w64 gdb-mingw-w64-target gfortran-mingw-w64 gfortran-mingw-w64-i686 gfortran-mingw-w64-x86-64 mingw-w64 mingw-w64-common mingw-w64-i686-dev mingw-w64-tools mingw-w64-x86-64-dev texlive-fonts-extra
sudo apt-get install ca-certificates make mingw-w64 gcc xz-utils p7zip-full bzip2 patch git gfortran autoconf
sudo apt-get install parallel libtool pkg-config libssl-dev libcurl4-openssl-dev flex bison gfortran libsuitesparse-dev texlive texlive-publishers texlive-extra-utils texlive-formats-extra texlive-latex-extra texlive-math-extra texlive-fonts-extra latex-beamer texinfo texi2html latex2html doxygen gcc-multilib g++-multilib gfortran-multilib nsis zip libboost-dev
sudo apt-get autoremove

git pull --recurse-submodules
git submodule foreach git pull --recurse-submodules
git fetch upstream --recurse-submodules
git submodule foreach git fetch upstream --recurse-submodules
git checkout master
git merge upstream/master
git commit -m "merge"
git push

cd dynare-build/libs
make build
cd ..
./build.sh
cd ..
