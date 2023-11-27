# morpho-kdtree
 Fast KD-Tree data structure for Morpho
 
 (NOTE: This package can be imported using `import fastkdtree`, to avoid clash with the existing Morpho-based `kdtree` module.)


## Installation 

To install this, clone this repository onto your computer in any convenient place:

    git clone https://github.com/joshichaitanya3/morpho-kdtree.git

then add the location of this repository to your .morphopackages file.

    cd morpho-kdtree
    pwd >> ~/.morphopackages 

This extension needs to be compiled as follows:

    mkdir lib
    cd src
    make install

The package can be loaded into morpho using the `import` keyword.

    import fastkdtree
