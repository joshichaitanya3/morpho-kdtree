/** @file kdree.c
 *  @author Chaitanya Joshi
 *
 *  @brief Implements the kd-tree data structure
 */

#include <stdio.h>
#include <morpho.h>
#include <common.h>
#include <morpho/list.h>
#include <math.h>
#include "ckdtree.h"

objecttype objectkdtreetype;
objecttype objectkdtreenodetype;


/* **********************************************************************
 * Utility functions
 * ********************************************************************** */
// DEFINE_VARRAY(kdtreenode, objectkdtreenode)


/* **********************************************************************
 * objectkdtreenode definitions
 * ********************************************************************** */

void objectkdtreenode_printfn(object *obj, void *v) {
    morpho_printf(v, "<KDTreeNode>");
}

// void objectkdtreenode_freefn(object *obj) {
//     objectkdtreenode *node = (objectkdtreenode *) obj;
//     kdtree_freetree(node);
// }

// void objectlist_markfn(object *obj, void *v) {
//     objectlist *c = (objectlist *) obj;
//     morpho_markvarrayvalue(v, &c->val);
// }

size_t objectkdtreenode_sizefn(object *obj) {
    return sizeof(objectkdtreenode);
}

objecttypedefn objectkdtreenodedefn = {
    .printfn=objectkdtreenode_printfn,
    .markfn=NULL,
    .freefn=NULL,
    .sizefn=objectkdtreenode_sizefn
};

/* todo: change the double[] to objectmatrix, thus allowing variable dimension*/
objectkdtreenode *object_newkdtreenode(double point[kdtree_dimension], int id) {
    objectkdtreenode *new = (objectkdtreenode *) object_new(sizeof(objectkdtreenode), OBJECT_KDTREE_NODE);
    if (new) {
        for (int i = 0; i < kdtree_dimension; i++) {
            new->point[i] = point[i];
        }
        new->id = id;
        new->tol = MORPHO_EPS;
        new->left = NULL;
        new->right = NULL;
    }

    return new;
}

/* **********************************************************************
 * objectkdtree definitions
 * ********************************************************************** */

void objectkdtree_printfn(object *obj, void *v) {
    morpho_printf(v, "<CKDTree>");
}

// void objectkdtreenode_freefn(object *obj) {
//     objectkdtreenode *node = (objectkdtreenode *) obj;
//     kdtree_freetree(node);
// }

// void objectlist_markfn(object *obj, void *v) {
//     objectlist *c = (objectlist *) obj;
//     morpho_markvarrayvalue(v, &c->val);
// }

size_t objectkdtree_sizefn(object *obj) {
    return sizeof(objectkdtree);
}

objecttypedefn objectkdtreedefn = {
    .printfn=objectkdtree_printfn,
    .markfn=NULL,
    .freefn=NULL,
    .sizefn=objectkdtree_sizefn
};

objectkdtree *object_newkdtree(objectlist* points, int depth) {

    int np = list_length(points);
    objectkdtree *new = (objectkdtree *) object_new(sizeof(objectkdtree)+np*sizeof(objectkdtreenode), OBJECT_KDTREE);
    if (new) {
        new->dimension = kdtree_dimension; // Keeping this hard-coded for now
        new->tol = MORPHO_EPS;
        new->head = kdtree_build(points, 0);
    }

    return new;
}

/* **********************************************************************
 * objectkdnode utility functions
 * ********************************************************************** */


/* **********************************************************************
 * objectkdtree utility functions
 * ********************************************************************** */

/* Comparator function that can sort the list of points along an axis (comparator_axis).
Allows for list of matrices or list of lists */
static vm * comparator_vm;
static int comparator_axis;
static bool comparator_err;

int comparator(const void* p, const void* q) {
    value al = *(value *) p;
    value ar = *(value *) q;

    value l, r;
    double outval;
    if (MORPHO_ISMATRIX(al)) {
        // printf("Left is a matrix.\n");
        objectmatrix *ml = MORPHO_GETMATRIX(al);
        comparator_err = !matrix_getelement(ml, comparator_axis, 0, &outval);
        // printf("Left outval is %f \n", outval);
        l = MORPHO_FLOAT(outval);
    }
    else if (MORPHO_ISLIST(al)) {
        // printf("Left is a list.\n");
        objectlist *ml = MORPHO_GETLIST(al);
        comparator_err = !list_getelement(ml, comparator_axis, &l);
        // printf("Left outval is %f \n", MORPHO_GETFLOATVALUE(l));
        // printf("Left outval is ");
        // morpho_printvalue(comparator_vm, l);
        // printf("\n");
        
    }
    else {
        comparator_err = true;
        return 0;
    }

    if (MORPHO_ISMATRIX(ar)) {
        // printf("Right is a matrix.\n");
        objectmatrix *mr = MORPHO_GETMATRIX(ar);
        comparator_err = !matrix_getelement(mr, comparator_axis, 0, &outval);
        // printf("Right outval is %f \n", outval);
        r = MORPHO_FLOAT(outval);
    }   
     else if (MORPHO_ISLIST(ar)) {
        // printf("Right is a list.\n");
        objectlist *mr = MORPHO_GETLIST(ar);
        comparator_err = !list_getelement(mr, comparator_axis, &r);
        // printf("Right outval is %f \n", MORPHO_GETFLOATVALUE(r));
        // printf("Right outval is ");
        // morpho_printvalue(comparator_vm, r);
        // printf("\n");   
    } 
    else {
        comparator_err = true;
        return 0;
    }

    return -morpho_comparevalue(l, r);

}
/* Probably going to be obsolete once the points in the node become of objectmatrix type*/
bool kdtree_valuetodoublearray(value val, double *pt) {

    value out;
    bool success = true;
    if (MORPHO_ISMATRIX(val)) {
        objectmatrix *m = MORPHO_GETMATRIX(val);
        // printf("Left is a matrix.\n");
        for (int i=0; i<kdtree_dimension; i++) {
            success = (success && matrix_getelement(m, i, 0, pt+i));
            // pt[i] = MORPHO_GETFLOATVALUE(out);
        }
    }
    else if (MORPHO_ISLIST(val)) {
        // printf("Left is a list.\n");
        objectlist *m = MORPHO_GETLIST(val);
        for (int i=0; i<kdtree_dimension; i++) {
            success = (success && list_getelement(m, i, &out));
            pt[i] = MORPHO_GETFLOATVALUE(out);
        }
    }

    return success;
}

// Function to build a kd-tree from a given list of points
objectkdtreenode* kdtree_build(objectlist* points, int depth) {
    int np = list_length(points);
    if (np==0) {
        objectkdtreenode* out = NULL;
        return out;
    }

    int axis = depth % kdtree_dimension;

    // Sort points to find median along current axis
    qsort(points->val.data, points->val.count, sizeof(value), comparator);
    int ipivot = ceil((np-1)/2);

    // Create lists for points on the left and right of this plane
    objectlist * left = object_newlist(ipivot, NULL);
    if (left && ipivot>=1) {
        for (unsigned int i=0; i<ipivot; i++) {
            left->val.data[i]=points->val.data[i];
        }
        left->val.count=ipivot;
    }

    objectlist * right = object_newlist(np-ipivot-1, NULL);
    if (right && ((np-ipivot-1)>=1)) {
        for (unsigned int i=0; i<(np-ipivot-1); i++) {
            right->val.data[i]=points->val.data[ipivot+1+i];
        }
        right->val.count=np-ipivot-1;
    }

    // Create node from the point at ipivot
    value ptval = points->val.data[ipivot];
    double pt[kdtree_dimension];
    bool success = kdtree_valuetodoublearray(ptval, &pt);
    
    objectkdtreenode* node = object_newkdtreenode(pt, 0);

    node->left = kdtree_build(left, depth+1);

    node->right = kdtree_build(right, depth+1);

    // Clear memory allocated for the temporary left and right lists
    varray_valueclear(&left->val);
    varray_valueclear(&right->val);

    return node;
}

// bool kdtree_ismember(pt) { // Tests if a point is in the tree.
//     var depth = 0

//     for (var node = self.head, next; !isnil(node); ) {
//       var pivot = mod(depth, self.dimension)
//       if (abs(pt[pivot] - node.location[pivot])<self.tol &&
//                 (pt-node.location).norm()<self.tol) {
//         return node
//       } else if (pt[pivot] < node.location[pivot]) {
//         node = node.left
//       } else {
//         node = node.right
//       }
//       depth+=1
//     }
//     return false
//   }

void kdtree_printnode(vm *v, objectkdtreenode* node) {
    morpho_printf(v, "[%g, %g, %g]\n", node->point[0], node->point[1], node->point[2]);
    
}

void kdtree_printtreefromnode(vm *v, objectkdtreenode* node, int depth) {
    morpho_printf(v, "Depth %d: ", depth);
    for (int i=0; i<depth; i++) {
        morpho_printf(v, " ");
    }
    kdtree_printnode(v, node);
    if (node->left) kdtree_printtreefromnode(v, node->left, depth+1);
    if (node->right) kdtree_printtreefromnode(v, node->right, depth+1);
}

void kdtree_print(vm *v, objectkdtree* tree) {
    objectkdtreenode* head = tree->head;
    kdtree_printtreefromnode(v, head, 0);
}

/* **********************************************************************
 * KD-Tree veneer class
 * ********************************************************************** */
   
// value kdtreenode_constructor(vm *v, int nargs, value *args) {
//     value out=MORPHO_NIL;
//     objectlist* points = NULL;
//     objectkdtreenode *new=NULL;
//     int id;

//     if (nargs==2 && MORPHO_ISLIST(MORPHO_GETARG(args, 0)) && MORPHO_ISINTEGER(MORPHO_GETARG(args, 1))) {
//         points = list_clone(MORPHO_GETLIST((MORPHO_GETARG(args, 0)))); 

//         depth = MORPHO_GETINTEGERVALUE(MORPHO_GETARG(args, 1));
//         new = object_newkdtree(points, depth);
//     }
//     else morpho_runtimeerror(v, KDTREE_CONSTRUCTOR);
//     if (comparator_err) morpho_runtimeerror(v, KDTREE_SORTDIM);
//     if (new) {
//         out=MORPHO_OBJECT(new);
//         morpho_bindobjects(v, 1, &out);


// }

/* **********************************************************************
 * KD-Tree veneer class
 * ********************************************************************** */
    
/** Constructor function for CKDTree */
value kdtree_constructor(vm *v, int nargs, value *args) {
    comparator_vm=v;
    value out=MORPHO_NIL;
    objectlist* points = NULL;
    objectkdtree *new=NULL;
    int depth;

    if (nargs==2 && MORPHO_ISLIST(MORPHO_GETARG(args, 0)) && MORPHO_ISINTEGER(MORPHO_GETARG(args, 1))) {
        points = list_clone(MORPHO_GETLIST((MORPHO_GETARG(args, 0)))); // Important to clone this, since we sort it while building.
        depth = MORPHO_GETINTEGERVALUE(MORPHO_GETARG(args, 1));
        new = object_newkdtree(points, depth);
    }
    else morpho_runtimeerror(v, KDTREE_CONSTRUCTOR);
    if (comparator_err) morpho_runtimeerror(v, KDTREE_SORTDIM);
    if (new) {
        out=MORPHO_OBJECT(new);
        morpho_bindobjects(v, 1, &out);

    }

    return out;
}

value KDTree_print(vm *v, int nargs, value *args) {
    value self = MORPHO_SELF(args);
    if (!MORPHO_ISKDTREE(self)) return Object_print(v, nargs, args);
    
    objectkdtree *tree=MORPHO_GETKDTREE(MORPHO_SELF(args));
    kdtree_print(v, tree);
    return MORPHO_NIL;
}

// // Function to insert a point into the kd-tree
// objectkdtreenode* kdtree_insert(objectkdtreenode* root, double point[kdtree_dimension], int id, int depth) {
//     if (root == NULL) {
//         return object_newkdtreenode(point, id);
//     }

//     int pivot = depth % kdtree_dimension;
//     double delta = 0;
//     for (int i = 0; i < kdtree_dimension; i++) {
//         delta += pow((point[i] - root->point[i]),2);
//     }
//     delta = pow(delta, 0.5);
//     if (fabs(point[pivot] - root->point[pivot])<root->tol &&
//                 delta<root->tol) {
//         printf("Warning: duplicate node\n");
//         return root;
//     } else if (point[pivot] < root->point[pivot]) {
//         root->left = kdtree_insert(root->left, point, id, depth + 1);
//     } else {
//         root->right = kdtree_insert(root->right, point, id, depth + 1);
//     }

//     return root;
// }


// // Function to concatenate two node lists
// kdtreenodelist* kdtree_concatenatelists(kdtreenodelist* list1, kdtreenodelist* list2) {
//     if (list1 == NULL) {
//         return list2;
//     }

//     kdtreenodelist* current = list1;
//     while (current->next != NULL) {
//         current = current->next;
//     }

//     current->next = list2;

//     return list1;
// }

// // Function to free the memory allocated for the kd-tree
// void kdtree_freetree(kdtreenode* root) {
//     if (root != NULL) {
//         kdtree_freetree(root->left);
//         kdtree_freetree(root->right);
//         MORPHO_FREE(root);
//     }
// }

// // Function to free the memory allocated for a node list
// void kdtree_freenodelist(kdtreenodelist* list) {
//     while (list != NULL) {
//         kdtreenodelist* temp = list;
//         list = list->next;
//         MORPHO_FREE(temp);
//     }
// }

// int main() {
//     kdtreenode* root = NULL;

//     // Insert points into the kd-tree
//     root = kdtree_insert(root, (double[]){2, 3}, 0, 0);
//     root = kdtree_insert(root, (double[]){5, 4}, 1, 0);
//     root = kdtree_insert(root, (double[]){9, 6}, 2, 0);
//     root = kdtree_insert(root, (double[]){4, 7}, 3, 0);
//     root = kdtree_insert(root, (double[]){8, 1}, 4, 0);
//     root = kdtree_insert(root, (double[]){7, 2}, 5, 0);

//     // Perform range search
//     printf("Points within range (4, 2) to (8, 6):\n");
//     kdtreenodelist* result = kdtree_rangesearch(root, (double[]){4, 2}, (double[]){8, 6}, 0);

//     // Print the result
//     kdtreenodelist* current = result;
//     while (current != NULL) {
//         printf("Point: (");
//         for (int i = 0; i < kdtree_dimension; i++) {
//             printf("%.2f", current->node->point[i]);
//             if (i < kdtree_dimension - 1) {
//                 printf(", ");
//             }
//         }
//         printf(")\n");

//         current = current->next;
//     }

//     // Free the memory allocated for the kd-tree and result list
//     kdtree_freetree(root);
//     kdtree_freenodelist(result);

//     return 0;
// }

// // Function to create a new node list
// objectkdtreenodelist *object_newkdtreenodelist(objectkdtreenode* node) {
// // NodeList* createNodeList(Node* node) {
//     // NodeList* newNode = (NodeList*)malloc(sizeof(NodeList));
//     objectkdtreenodelist *new = (objectkdtreenodelist *) object_new(sizeof(objectkdtreenodelist), OBJECT_KDTREE_NODE_LIST);
//     if (new != NULL) {
//         new->node = node;
//         new->next = NULL;
//     }
//     return new;
// }
MORPHO_BEGINCLASS(CKDTree)
MORPHO_METHOD(MORPHO_PRINT_METHOD, KDTree_print, BUILTIN_FLAGSEMPTY)
MORPHO_ENDCLASS


/* **********************************************************************
 * Initialization and finalization
 * ********************************************************************** */

objecttype objectkdtreenodetype;
objecttype objectkdtreetype;

void fastkdtree_initialize(void) {
    // Define list objecttype
    objectkdtreenodetype=object_addtype(&objectkdtreenodedefn);
    objectkdtreetype=object_addtype(&objectkdtreedefn);
    
    // Locate the Object class to use as the parent class of Range
    objectstring objname = MORPHO_STATICSTRING(OBJECT_CLASSNAME);
    value objclass = builtin_findclass(MORPHO_OBJECT(&objname));
    
    // List constructor function
    builtin_addfunction(KDTREE_CLASSNAME, kdtree_constructor, BUILTIN_FLAGSEMPTY);
    
    // List constructor function
    value kdtreeclass=builtin_addclass(KDTREE_CLASSNAME, MORPHO_GETCLASSDEFINITION(CKDTree), objclass);
    object_setveneerclass(OBJECT_KDTREE, kdtreeclass);
    
    // List error messages
    morpho_defineerror(KDTREE_CONSTRUCTOR, ERROR_HALT, KDTREE_CONSTRUCTOR_MSG);
    morpho_defineerror(KDTREE_SORTDIM, ERROR_HALT, KDTREE_SORTDIM_MSG);
}

void fastkdtree_finalize(void) {
}
