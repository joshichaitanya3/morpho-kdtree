[comment]: # (Fastkdtree module help)
[version]: # (0.0.1)

# fastkdtree
[tagfastkdtree]: # (fastkdtree)

The `fastkdtree` module provides a C-based implementation of the KD-Tree data-structure. This is intended to replicate the Morpho-based `kdtree` module that currently ships with Morpho. It is restricted to 3D for now, but will be made general soon.

To use the module, first import it:

    import fastkdtree

The main object it provides is the `CKDTree` (look up the help for it).

# CKDTree
[tagckdtree]: # (ckdtree)

The `CKDTree` object is the equivalent of the `KDTree` object in the Morpho module. It can be initialized with list of 3D points. These points can themselves be either a `List` or a `Matrix`:

    var ctree = CKDtree(points)

The `CKDTree` object provides a number of methods, similar to `KDTree`.

[showsubtopics]: # (subtopics)

## Head
The `head` method returns a `CKDTreeNode` object (see `help` for `CKDTreeNode`) corresponding to the head of the tree.


## Search
The `search` method searches for nodes whose location lies within the bounding box provided by the user.
Input:
    `query`: List of 2-element lists containing the lower and upper bounds, eg `[[xl, xh], [yl, yh], [zl, zh]]`, where `l` and `h` correspond to lower and upper bounds in the respective dimensions. 

Returns
    `points`: a list of `CKDTreeNode` objects

Example:
    var points = ctree.search([[-5,5], [-4,4], [-1,1]])

## Ismember
The `ismember` method tries to find a node with a given location.
Input: 
* `point`: 3-element `List` or `Matrix` object containing the location

Output:
* `node`: The `CKDTreeNode` corresponding to the location if found, and `nil` otherwise.

## Insert
The `insert` method inserts a node in the tree.

Input:
* `location`: 3D `List` or `Matrix` corresponding to the new location
* `id` (optional): Index of the point. Will be assigned `0` if not probided.

Output:
* `node`: The `CKDTreeNode` object corresponding to the new node. If the location already exists, then the corresponding existing node is returned.

## Nearest
The `nearest` method returns the node in the tree nearest to the input point.

Input:
* `location`: 3D `List` or `Matrix` corresponding to the new location

Output:
* `node`: The `CKDTreeNode` object corresponding to the node nearest to the location.

## Maxdepth
The `maxdepth` method returns the maximum depth of a KD-Tree.

## CKDTreeNode
[tagckdtreenode]: # (ckdtreenode)

An object corresponding to a single node in a k-d tree. To get the location of the node, use the `location` method:

    var m = node.location()

To get the index, use the `id` method:

    var id = node.id()

You can access the left and right nodes using the corresponding methods:
    
    var l = node.left()
    var r = node.right()
