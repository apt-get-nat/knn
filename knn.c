/* Description: Library implementation of the KDtree and
                KNNsearch, which finds M nearest xpt neighbors
                for each ypt in M(log N) time.
	Author: Nathaniel Mathews,
	        based on a nearest-neighbor algorithm implemented
	        at rosettacode.org/wiki/K-d_tree
	Dates modified:
	    Last modified 2016-05-31
	          v.0.5 -            - external hooks completed
	        v.0.4.2 - 2016-06-01 - major (M>1) bug fix
	        v.0.4.1 - 2016-05-31 - various algorithm bug fixes
	          v.0.4 - 2016-03-23 - generalization to N
	                               nearest neighbors
	          v.0.3 - 2016-05-22 - external hook prototypes
	          v.0.2 - 2016-05-14 - ported to this project
	          v.0.1 - 2015-05-06 - rosetta implementation
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_DIM 3

struct kd_node_t{
	double x[MAX_DIM];
	struct kd_node_t *left, *right;
	int index;
};

inline double dist(struct kd_node_t *a, struct kd_node_t *b, int dim) {
	double t, d = 0;
	while (dim--){
		t = a->x[dim] - b->x[dim];
		d += t * t;
	}
	return sqrt(d);
}
inline void swap(struct kd_node_t *x, struct kd_node_t *y) {
	double tmp[MAX_DIM];
	int tmpi;
	memcpy(tmp, x->x, sizeof(tmp));
	memcpy(x->x, y->x, sizeof(tmp));
	memcpy(y->x, tmp, sizeof(tmp));
	tmpi = x->index;
	x->index = y->index;
	y->index = tmpi;
}

struct kd_node_t* find_median(struct kd_node_t *start, struct kd_node_t *end, int idx) {
	if(end <= start) return NULL;
	if(end == start+1) return start;

	struct kd_node_t *p, *store, *md = start + (end - start) / 2;
	double pivot;
	while(1) {
		pivot = md->x[idx];

		swap(md, end - 1);
		for(store = p = start; p < end; p++) {
			if(p->x[idx] < pivot) {
				if(p != store) swap(p, store);
				store++;
			}
		}
		swap(store, end - 1);

		if(store->x[idx] == md->x[idx]) return md;

		if(store > md) end = store;
		else start = store;
	}
}

struct kd_node_t* make_tree(struct kd_node_t *t, int len, int i, int dim) {
	struct kd_node_t *n;
	if(!len) return 0;
	if((n = find_median(t, t + len, i))) {
		i = (i+1) % dim;
		n->left = make_tree(t, n-t, i, dim);
		n->right = make_tree(n+1, t+len-(n+1), i, dim);
	}
	return n;
}

int visited;

void nearest(struct kd_node_t *root, struct kd_node_t *nd, int i, 
		int dim, int *best, double *best_dist, int M) {
	double d, dx, dx2; int j;
	if (!root) return;

	d = dist(root, nd, dim);
	dx = root->x[i] - nd->x[i];
	dx2 = dx*dx;

	visited++;

	for(j=M-1; j>=0; j--) {
		if( d < best_dist[j] && j<M-1) {
			best_dist[j+1] = best_dist[j];
			best[j+1] = best[j];
		}
		if( ( ( d>=best_dist[j-1] && best_dist[j-1]>0 ) || j==0) 
				&& (d < best_dist[j] || best_dist[j] == 0)) {
			best_dist[j] = d;
			best[j] = root->index;
		}
		if(d >= best_dist[M-1] && best_dist[M-1] != 0) {
			j = 0;
		}
	}

	if(!best_dist[0]) return;
	if(++i >= dim) i = 0;
	nearest(dx > 0 ? root->left : root->right, nd, i, dim, best, best_dist, M);
	if(dx2 >= best_dist[M-1] && best_dist[M-1] != 0) return;
	nearest(dx > 0 ? root->right : root->left, nd, i, dim, best, best_dist, M);
}

struct kd_node_t *PACKAGETREE;

int SETUP( int argc, void *argv[] ) {
	// Builds a kdtree and stores it in PACKAGETREE for later use.
	//
	// ARGUMENTS: 3 expected.
	// 0. dim is the dimension of the problem (usually 3)
	// 1. int N is the number of data locations.
	// 2. double[] X is a dim-wrapped array of data locations.
	//
	// Returns 0 if improper number of arguments, 1 otherwise.

	if(argc != 3) {
		printf("ERROR: SETUP expected 3 arguments, but got %i\n", argc);
		return 0;
	}
	int dim = argv[0];
	int N = argv[1];
	double *X = argv[2]; // size N

	int i, j;
	struct kd_node_t *tree, *root;
	tree = (struct kd_node_t*) calloc(N, sizeof(struct kd_node_t));

	for (i = 0; i < N; i++) {
		for(j = 0; j < dim; j++) {
			tree[i].x[j] = X[i*dim+j];
		}
		tree[i].index = i;
	}

	PACKAGETREE = make_tree(tree, N, 0, dim);

	return 1;
}

int ARRAYSEARCH( int argc, void *argv[] ) {
	// Searches PACKAGETREE for the nearest elements to those
	// in an inputted array.
	//
	// ARGUMENTS: 6 expected.
	// 0. dim is the dimension of the problem (usually 3).
	// 1. int N is the number of locations to find nearest neighbors of.
	// 2. double[] Y is the array of locations.
	// 3. int M is the number of nearest neighbors to each requested.
	// 4. int[] K is a pre-allocated array of length N*M.
	// 5. double[] D is a pre-allocated array of length N*M.
	//
	// K and d will be permuted by the function. Upon return of 1,
	// * K will contain the indices in the X vector (see SETUP) of the
	//   nearest neighbors, ordered as the array Y with length M looping.
	// * D will contain the distances to the elements of K.
	// A return of 0 indicates failure based on arguments inputted.

	if(argc != 6) {
		printf("ERROR: SETUP expected 6 arguments, but got %i\n", argc);
		return 0;
	}
	int dim = argv[0];
	int  N = argv[1];
	double *Y = argv[2]; // size dim*N
	int M = argv[3];
	int *K = argv[4];    // size N*M
	double *D = argv[5]; // size N*M
	int i, j;

	struct kd_node_t testNode;
	struct kd_node_t *root, *million;
	double best_dist[M]; int found[M];

	for(i = 0; i < N; i++) {
		visited = 0;
		memset(found, 0, M*sizeof(int));
		memset(best_dist, 0, M*sizeof(double));
		double Yloc[dim];
		for(j=0; j<dim; j++) testNode.x[j] = Y[dim*i+j];
		nearest(PACKAGETREE, &testNode, 0, dim, found, best_dist, M);
		for(j = 0; j < M; j++) {
			K[M*i + j] = found[j]; D[M*i + j] = best_dist[j];
		}
	}
	return 1;
}

int INSPECT_helper( struct kd_node_t *node, int dim ) {
	// Provides a printout of the PACKAGETREE.
	//
	// ARGUMENTS: 0 expected if called from outside,
	//            1 (next node) if called recursively.
	int i;
	if(!node) {
		printf("--^");
		return;
	}
	printf("%i: ", node->index);
	for(i=0; i<dim; i++) {
		printf("%f ", node->x[i]);
	}
	printf("\n/"); INSPECT_helper(node->left,  dim);
	printf("\n\\"); INSPECT_helper(node->right, dim);
	return;
}
int INSPECT(int argc, void *argv[]) {
	// Provides a printout of the PACKAGETREE.
	//
	// ARGUMENTS: 1 expected.
	// 0. dim is the number of dimensions.
	if(argc == 1) {
		int dim = argv[0];
		printf(" "); INSPECT_helper(PACKAGETREE, dim); printf("\n");
		return 1;
	} else {
		return 0;
	}
}

/*

// Main function is for c compilation. DLL compilation doesn't need it.

int main(void) {
	// For best results, comment out either the first or second case.

	// simple case: M = 1
	int i;
	double wp_array[7][2] = {{0,1},{0,2},{0,3},{0,4},{0,5},{0,6},{0,7}};
	struct kd_node_t wp[7];

	for(i=0;i<sizeof(wp_array)/sizeof(wp_array[0]);i++) {
		wp[i].x[0] = wp_array[i][0]; wp[i].x[1] = wp_array[i][1]; wp[i].index = i;
	}

	struct kd_node_t testNode = {{0,3.7}};
	struct kd_node_t *root, *million;
	double best_dist[1]; int found[1];

	root = make_tree(wp, sizeof(wp) / sizeof(wp[1]),0,2);
	visited = 0;
	memset(found, 0, sizeof(found));
	memset(best_dist, 0, m*sizeof(double));

	nearest(root, &testNode, 0, 2, found, best_dist, 1);

	printf("\n\nData set:\n");
	for(i=0; i<sizeof(wp)/sizeof(wp[1]); i++) {
		printf("(%g,%g), ", wp[i].x[0],wp[i].x[1]);
	}

	printf("\n>> WP tree\nsearching for (%g, %g)\n"
		   "found (%g, %g) dist %g\nseen %d nodes\n\n",
		   testNode.x[0], testNode.x[1],
		   wp_array[found[0]][0], wp_array[found[0]][1], sqrt(best_dist[0]), visited
		  );
	printf("----------------------------------\n");
	printf("root:  (%g,%g)\n",
		root[0].x[0], root[0].x[1]);
	printf("(%g,%g)        (%g,%g)\n",
		root[0].left->x[0], root[0].left->x[1],
		root[0].right->x[0], root[0].right->x[1]);
	printf("(%g,%g)(%g,%g)  (%g,%g)(%g,%g)\n",
		root[0].left->left->x[0], root[0].left->left->x[1],
		root[0].left->right->x[0], root[0].left->right->x[1],
		root[0].right->left->x[0], root[0].right->left->x[1],
		root[0].right->right->x[0], root[0].right->right->x[1]);

	// Case 2: M = 2

	int i; int m = 2;
	double wp_array[7][2] = {{0,1},{0,2},{0,3},{0,4},{0,5},{0,6},{0,7}};
	struct kd_node_t wp[7];

	for(i=0;i<sizeof(wp_array)/sizeof(wp_array[0]);i++) {
		wp[i].x[0] = wp_array[i][0]; wp[i].x[1] = wp_array[i][1]; wp[i].index = i;
	}

	struct kd_node_t testNode = {{0,3.7}};
	struct kd_node_t *root, *million;
	double best_dist[m]; int found[m];

	root = make_tree(wp, sizeof(wp) / sizeof(wp[1]),0,2);
	visited = 0;
	memset(found, 0, m*sizeof(int));
	memset(best_dist, 0, m*sizeof(double));

	nearest(root, &testNode, 0, 2, found, best_dist, m);

	printf("\n\nData set:\n");
	for(i=0; i<sizeof(wp)/sizeof(wp[1]); i++) {
		printf("(%g,%g), ", wp[i].x[0],wp[i].x[1]);
	}

	printf("\n>> WP tree\nsearching for (%g, %g)\n",
		   testNode.x[0], testNode.x[1]);
	for(i=0; i<m; i++) {
		printf(" %i. found (%g, %g) dist %g\n",
			   i, wp_array[found[i]][0],
			   wp_array[found[i]][1], sqrt(best_dist[i])
		  );
	}
	printf("----------------------------------\n");
	printf("root:  (%g,%g)\n",
		root[0].x[0], root[0].x[1]);
	printf("(%g,%g)        (%g,%g)\n",
		root[0].left->x[0], root[0].left->x[1],
		root[0].right->x[0], root[0].right->x[1]);
	printf("(%g,%g)(%g,%g)  (%g,%g)(%g,%g)\n",
		root[0].left->left->x[0], root[0].left->left->x[1],
		root[0].left->right->x[0], root[0].left->right->x[1],
		root[0].right->left->x[0], root[0].right->left->x[1],
		root[0].right->right->x[0], root[0].right->right->x[1]);

	return 0;
}
*/


























