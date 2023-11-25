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
objectkdtreenode *object_newkdtreenode(value ptval, int id) {
    objectkdtreenode *new = (objectkdtreenode *) object_new(sizeof(objectkdtreenode), OBJECT_KDTREENODE);
    double point[kdtree_dimension];
    bool success = kdtree_valuetodoublearray(ptval, &point);
    if (new && success) {
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

bool kdtree_doublearraytomatrix(double pt[kdtree_dimension], objectmatrix *m) {

}
// Function to compute the l-norm distance between two points 
double kdtree_norm(double* pt1, double* pt2, double l) {
    double dist = 0.0;
    for (int i = 0; i < kdtree_dimension; i++) {
        dist += pow((pt1[i] - pt2[i]),l);
    }
    dist = pow(dist , 1.0/l);

    return dist;
}

// Function to build a kd-tree from a given list of points
objectkdtreenode* kdtree_build(objectlist* points, int depth) {
    int np = list_length(points);
    if (np==0) {
        objectkdtreenode* out = NULL;
        return out;
    }

    comparator_axis  = depth % kdtree_dimension;
    printf("np is now %d... \n", np);
    printf("Depth is now %d... \n", depth);
    printf("Comparator axis is now %d... \n", comparator_axis);
    // Sort points to find median along current axis
    qsort(points->val.data, points->val.count, sizeof(value), comparator);
    int ipivot = round((double)(np-1)/2);
    printf("ipivot is now %d... \n", ipivot);

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
    
    objectkdtreenode* node = object_newkdtreenode(ptval, 0);

    node->left = kdtree_build(left, depth+1);

    node->right = kdtree_build(right, depth+1);

    // Clear memory allocated for the temporary left and right lists
    varray_valueclear(&left->val);
    varray_valueclear(&right->val);

    return node;
}

// Tests if a point is in the tree and returns the point if true
objectkdtreenode* kdtree_ismember(objectkdtree* tree, value ptval) { 
    int depth = 0;
    double pt[kdtree_dimension];
    bool success = kdtree_valuetodoublearray(ptval, &pt);
    objectkdtreenode* node = tree->head;
    int pivot;
    double dist;
    while(node) {
        pivot = depth % kdtree_dimension;
        dist = kdtree_norm(pt, node->point, 2);
        if (fabs(pt[pivot] - node->point[pivot])<tree->tol && dist<tree->tol) {
            return node;
        }
        else if (pt[pivot] < node->point[pivot]) {
            node = node->left;
        } 
        else {
            node = node->right;
        }
        depth += 1;
    }
    return NULL;
  }

// Function to insert a point into the kd-tree
objectkdtreenode* kdtree_insert(objectkdtree* tree, value ptval, int id) { // Inserts a point into the tree
    int depth = 0;
    double pt[kdtree_dimension];
    bool success = kdtree_valuetodoublearray(ptval, &pt);
    objectkdtreenode* node = tree->head;
    objectkdtreenode* next;

    int pivot;
    double dist;
    // for (var node = self.head, next; node!=nil; node=next) {
    while (node) {
        pivot = depth % kdtree_dimension;
        dist = kdtree_norm(pt, node->point, 2);
        if (fabs(pt[pivot] - node->point[pivot])<tree->tol && dist<tree->tol) {
            printf("Warning: duplicate node!\n");
            return node;
        }
        else if (pt[pivot] < node->point[pivot]) {
            next = node->left;
            if (!next) {
            //   node.left = KDTreeNode(pt, nil, nil)
                node->left = object_newkdtreenode(ptval, id);
                return node->left;
            }
        } else {
            next = node->right;
            if (!next) {
                // node.right = KDTreeNode(pt, nil, nil)
                node->right = object_newkdtreenode(ptval, id);
                return node->right;
            }
        }
        depth += 1;
        node = next;
    }   
}

void kdtree_printnode(vm *v, objectkdtreenode* node) {
    morpho_printf(v, "[%g, %g, %g]", node->point[0], node->point[1], node->point[2]);
    
}

void kdtree_printtreefromnode(vm *v, objectkdtreenode* node, int depth) {
    morpho_printf(v, "Depth %d: ", depth);
    for (int i=0; i<depth; i++) {
        morpho_printf(v, " ");
    }
    kdtree_printnode(v, node);
    morpho_printf(v, "\n");
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
   
value kdtreenode_constructor(vm *v, int nargs, value *args) {
    value out=MORPHO_NIL;
    objectkdtreenode *new=NULL;
    int id;
    double point[kdtree_dimension];
    value ptval = MORPHO_NIL;
    if (nargs==1 && (MORPHO_ISLIST(MORPHO_GETARG(args, 0)) || MORPHO_ISMATRIX(MORPHO_GETARG(args, 0)))) {
        ptval = MORPHO_GETARG(args, 0);
        new = object_newkdtreenode(ptval, id);
    }
    else morpho_runtimeerror(v, KDTREENODE_CONSTRUCTOR);
    if (new) {
        out=MORPHO_OBJECT(new);
        morpho_bindobjects(v, 1, &out);
    }

    return out;
}

value KDTreeNode_print(vm *v, int nargs, value *args) {
    value self = MORPHO_SELF(args);
    if (!MORPHO_ISKDTREENODE(self)) return Object_print(v, nargs, args);
    
    objectkdtreenode *node=MORPHO_GETKDTREENODE(MORPHO_SELF(args));
    kdtree_printnode(v, node);
    return MORPHO_NIL;
}

value KDTreeNode_left(vm *v, int nargs, value *args) {
    value out=MORPHO_NIL;
    objectkdtreenode* new = NULL;
    objectkdtreenode *node=MORPHO_GETKDTREENODE(MORPHO_SELF(args));
    new = node->left;
    if (new) {
        out = MORPHO_OBJECT(new);
        morpho_bindobjects(v, 1, &out);
    }
    return out;
}

value KDTreeNode_right(vm *v, int nargs, value *args) {
    value out=MORPHO_NIL;
    objectkdtreenode* new = NULL;
    objectkdtreenode *node=MORPHO_GETKDTREENODE(MORPHO_SELF(args));
    new = node->right;
    if (new) {
        out = MORPHO_OBJECT(new);
        morpho_bindobjects(v, 1, &out);
    }
    return out;
}

MORPHO_BEGINCLASS(CKDTreeNode)
MORPHO_METHOD(MORPHO_PRINT_METHOD, KDTreeNode_print, BUILTIN_FLAGSEMPTY),
MORPHO_METHOD(KDTREENODE_LEFT_METHOD, KDTreeNode_left, BUILTIN_FLAGSEMPTY),
MORPHO_METHOD(KDTREENODE_RIGHT_METHOD, KDTreeNode_right, BUILTIN_FLAGSEMPTY)
MORPHO_ENDCLASS

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

value KDTree_ismember(vm *v, int nargs, value *args) {
    objectkdtree *tree=MORPHO_GETKDTREE(MORPHO_SELF(args));
    objectkdtreenode* node;
    value ptval = MORPHO_NIL;
    value out=MORPHO_NIL;
    if (nargs==1 && (MORPHO_ISLIST(MORPHO_GETARG(args, 0)) || MORPHO_ISMATRIX(MORPHO_GETARG(args, 0)))) {
        ptval = MORPHO_GETARG(args, 0);
        node = kdtree_ismember(tree, ptval);
    }
    if (node) {
        out=MORPHO_OBJECT(node);
        morpho_bindobjects(v, 1, &out);
    }
    return out;
}

value KDTree_insert(vm *v, int nargs, value *args) {
    objectkdtree *tree=MORPHO_GETKDTREE(MORPHO_SELF(args));
    objectkdtreenode* node;
    value ptval = MORPHO_NIL;
    value out=MORPHO_NIL;
    int id=0;
    if (nargs==1 && (MORPHO_ISLIST(MORPHO_GETARG(args, 0)) || MORPHO_ISMATRIX(MORPHO_GETARG(args, 0)))) {
        ptval = MORPHO_GETARG(args, 0);
        node = kdtree_insert(tree, ptval, id);
    }
    else if (nargs==2 && (MORPHO_ISLIST(MORPHO_GETARG(args, 0)) || MORPHO_ISMATRIX(MORPHO_GETARG(args, 0))) && MORPHO_ISINTEGER(MORPHO_GETARG(args, 1))) {
        ptval = MORPHO_GETARG(args, 0);
        id = MORPHO_GETINTEGERVALUE(MORPHO_GETARG(args, 1));
        node = kdtree_insert(tree, ptval, id);
    }
    if (node) {
        out=MORPHO_OBJECT(node);
        morpho_bindobjects(v, 1, &out);
    }
    return out;
}

value KDTree_head(vm *v, int nargs, value *args) {
    value out=MORPHO_NIL;
    objectkdtree *tree=MORPHO_GETKDTREE(MORPHO_SELF(args));
    objectkdtreenode* new = tree->head;
    if (new) {
        out = MORPHO_OBJECT(new);
        morpho_bindobjects(v, 1, &out);
    }
    return out;
}

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
MORPHO_METHOD(MORPHO_PRINT_METHOD, KDTree_print, BUILTIN_FLAGSEMPTY),
MORPHO_METHOD(KDTREE_HEAD_METHOD, KDTree_head, BUILTIN_FLAGSEMPTY),
MORPHO_METHOD(KDTREE_ISMEMBER_METHOD, KDTree_ismember, BUILTIN_FLAGSEMPTY),
MORPHO_METHOD(KDTREE_INSERT_METHOD, KDTree_insert, BUILTIN_FLAGSEMPTY)
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
    builtin_addfunction(KDTREENODE_CLASSNAME, kdtreenode_constructor, BUILTIN_FLAGSEMPTY);
    
    value kdtreeclass=builtin_addclass(KDTREE_CLASSNAME, MORPHO_GETCLASSDEFINITION(CKDTree), objclass);
    object_setveneerclass(OBJECT_KDTREE, kdtreeclass);
    
    value kdtreenodeclass=builtin_addclass(KDTREENODE_CLASSNAME, MORPHO_GETCLASSDEFINITION(CKDTreeNode), objclass);
    object_setveneerclass(OBJECT_KDTREENODE, kdtreenodeclass);
    
    // List error messages
    morpho_defineerror(KDTREE_CONSTRUCTOR, ERROR_HALT, KDTREE_CONSTRUCTOR_MSG);
    morpho_defineerror(KDTREENODE_CONSTRUCTOR, ERROR_HALT, KDTREENODE_CONSTRUCTOR_MSG);
    morpho_defineerror(KDTREE_SORTDIM, ERROR_HALT, KDTREE_SORTDIM_MSG);
}

void fastkdtree_finalize(void) {
}
