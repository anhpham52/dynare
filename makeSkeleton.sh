#!/bin/bash
releaseFolder='4.5.6';
gitFolder=${PWD##*/};

cd ..

rm releaseListing.txt
rm gitListing.txt
rm diffListings.txt

cd "$releaseFolder"
find . -type f 1> ../releaseListing.txt
cd ..

cd "$gitFolder"
find . -type f 1> ../gitListing.txt
cd ..

diff --new-line-format="" --unchanged-line-format="" <(sort releaseListing.txt) <(sort gitListing.txt) 1> diffListings.txt

rm dynareReleaseSkeleton.tar

cd "$releaseFolder"
tar -cf ../dynareReleaseSkeleton.tar -T ../diffListings.txt
cd ..

rm -rf dynareReleaseSkeleton
mkdir dynareReleaseSkeleton
cd dynareReleaseSkeleton

tar -xf ../dynareReleaseSkeleton.tar

cd "../$gitFolder"
