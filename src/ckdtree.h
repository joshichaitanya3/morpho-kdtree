#include <stdio.h>
#include <morpho/morpho.h>
#include <morpho/builtin.h>

extern objecttype objectkdtreenodetype;
#define OBJECT_KDTREE_NODE objectkdtreenodetype

extern objecttype objectkdtreetype;
#define OBJECT_KDTREE objectkdtreetype

/* This will probably remoed once `point` is changed to an objectmatrix */
#define kdtree_dimension 3  // Number of dimensions

typedef struct objectkdtreenode objectkdtreenode;
struct objectkdtreenode {
    object obj;
    double point[kdtree_dimension]; // Coordinates of the point in k-dimensional space
    int id; // ID of the vertex
    double tol;
    objectkdtreenode* left;
    objectkdtreenode* right;
};

/* --> The method below didn't work because of self-reference in the typedef struct*/
// typedeft struct {
//     object obj;
//     double point[kdtree_dimension]; // Coordinates of the point in k-dimensional space
//     int id; // ID of the vertex
//     double tol;
//     struct objectkdtreenode* left;
//     struct objectkdtreenode* right;
// };

typedef struct {
    object obj;
    int dimension;
    double tol;
    objectkdtreenode* head;
} objectkdtree;

// List structure to store nodes in range
DECLARE_VARRAY(kdtreenode, objectkdtreenode)

/** Tests whether an object is a objectkdtreenode */
#define MORPHO_ISKDTREENODE(val) object_istype(val, OBJECT_KDTREE_NODE)
/** Gets the object as a objectkdtreenode */
#define MORPHO_GETKDTREENODE(val) ((objectkdtreenode *) MORPHO_GETOBJECT(val))

/** Tests whether an object is a objectkdtree */
#define MORPHO_ISKDTREE(val) object_istype(val, OBJECT_KDTREE)
/** Gets the object as a objectkdtreenodelist */
#define MORPHO_GETKDTREE(val) ((objectkdtree *) MORPHO_GETOBJECT(val))

/** Gets the vertex id from a nodel */
#define MORPHO_GETKDTREENODEID(val) (((objectkdtreenode *) MORPHO_GETOBJECT(val))->id)

/* Prototypes for constructors*/
// objectkdtreenode *object_newkdtreenode(double point[kdtree_dimension], int id, struct objectkdtreenode* left, struct objectkdtreenode* right);
// objectkdtreenode *object_newkdtreenode(double point[kdtree_dimension], int id, objectkdtreenode* left, objectkdtreenode* right);
objectkdtreenode *object_newkdtreenode(double point[kdtree_dimension], int id);

objectkdtree *object_newkdtree(objectlist* points, int depth);

/* **********************************************************************
 * KD-Tree veneer class
 * ********************************************************************** */
    
#define KDTREE_CLASSNAME                  "CKDTree"
/** Constructor function for KDTree */
value kdtree_constructor(vm *v, int nargs, value *args);

/* -------------------------------------------------------
 * KDTree error messages
 * ------------------------------------------------------- */

#define KDTREE_CONSTRUCTOR                "KDTreeConstructorError"
#define KDTREE_CONSTRUCTOR_MSG            "KDTree(points, depth) constructor should be called with a list of 3d points and an integer depth"

#define KDTREE_SORTDIM                    "SortDimError"
#define KDTREE_SORTDIM_MSG                "Could not access the elements of the provided dimension."

/* **********************************************************************
 * objectkdtree utility functions
 * ********************************************************************** */

bool kdtree_valuetodoublearray(value val, double *pt);

int comparator(const void* p, const void* q);

// Function to build a kd-tree from a given list of points
objectkdtreenode* kdtree_build(objectlist* points, int depth);

void kdtree_printnode(vm *v, objectkdtreenode* node);
void kdtree_printtreefromnode(vm *v, objectkdtreenode* node, int depth);
void kdtree_print(vm *v, objectkdtree* tree);

void fastkdtree_initialize(void);
void fastkdtree_finalize(void);
