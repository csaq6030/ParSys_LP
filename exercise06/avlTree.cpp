#include <cstdlib>
#include <iostream>
#include <cmath>
#include <climits>
#include <random>
#include <algorithm>
#include <omp.h>
#include <vector>

using namespace std;

struct node{
    unsigned int value;
    node * left;
    node * right;
    int height;
};

typedef node  * nodeptr;

class avlTree {
    
private:
    
    nodeptr root;
    
    int get_height(nodeptr & p) {
        if(p != NULL)
            return p->height;
        else
            return 0;
    }
    
//        T1, T2, T3 and T4 are subtrees.
//               z                                      y
//              / \                                   /   \
//             y   T4      Right Rotate (z)          x      z
//            / \          - - - - - - - - ->      /  \    /  \
//           x   T3                               T1  T2  T3  T4
//          / \
//        T1   T2
    nodeptr left_left_rotation(nodeptr &p1) {
        nodeptr p2 = p1->left;
        p1->left = p2->right;
        p2->right = p1;
        p1->height = max(get_height(p1->left), get_height(p1->right))+1;
        p2->height = max(get_height(p2->left), p1->height)+1;
        return p2;
    }

    
//         z                                y
//        /  \                            /   \
//       T1   y     Left Rotate(z)       z      x
//           /  \   - - - - - - - ->    / \    / \
//          T2   x                     T1  T2 T3  T4
//              / \
//             T3  T4
    nodeptr right_right_rotation(nodeptr &p1) {
        nodeptr p2 = p1->right;
        p1->right = p2->left;
        p2->left = p1;
        p1->height = max(get_height(p1->left), get_height(p1->right))+1;
        p2->height = max(get_height(p2->right), p1->height)+1;
        return p2;
    }

    
//         z                               z                           x
//        / \                            /   \                        /  \
//       y   T4  Left Rotate (y)        x    T4  Right Rotate(z)    y      z
//      / \      - - - - - - - - ->    /  \      - - - - - - - ->  / \    / \
//    T1   x                          y    T3                    T1  T2 T3  T4
//        / \                        / \
//      T2   T3                    T1   T2
    nodeptr left_right_rotation(nodeptr &p1) {
        p1->left = right_right_rotation(p1->left);
        return left_left_rotation(p1);
    }

//       z                            z                            x
//      / \                          / \                          /  \
//    T1   y   Right Rotate (y)    T1   x      Left Rotate(z)   z      y
//        / \  - - - - - - - - ->     /  \   - - - - - - - ->  / \    / \
//       x   T4                      T2   y                  T1  T2  T3  T4
//      / \                              /  \
//    T2   T3                           T3   T4
    nodeptr right_left_rotation(nodeptr &p1){
        p1->right = left_left_rotation(p1->right);
        return right_right_rotation(p1);
    }
    
    void insert(unsigned int value, nodeptr & p){
        if(p == NULL){
            p = new node;
            p->value = value;
            p->left = NULL;
            p->right = NULL;
            p->height = 0;
            return;
        }
        
        if(value == p->value)
            return;
        else if (value < p->value) {
            insert(value, p->left);
            if(get_height(p->left) - get_height(p->right) == 2) {
                if(value < p->left->value)
                    p = left_left_rotation(p);
                else
                    p = left_right_rotation(p);
            }
        }
        else {
            insert(value, p->right);
            if ((get_height(p->right) - get_height(p->left)) == 2) {
                if(value > p->right->value)
                    p = right_right_rotation(p);
                else
                    p = right_left_rotation(p);
            }
        }
    }
    
    int get_value_left (nodeptr & p) {
        if (p->left != NULL) {
            return p->left->value;
        }
        else
            return -1;
    }
    
    int get_value_right (nodeptr & p) {
        if (p->right != NULL)
            return p->right->value;
        else
            return numeric_limits<int>::max();
    }
    
    int check_tree (nodeptr & p) {
        int out = 0;
        if (get_height(p->left) - get_height(p->right) < 2 && get_height(p->right) - get_height(p->left) < 2) {
            if (get_value_left(p) < (int) p->value  && (int) p->value < get_value_right(p)) {
                if (p->left != NULL)
                    out += check_tree (p->left);
                if (p->right != NULL)
                    out += check_tree (p->right);
            } else
                out++;
        } else
            out++;
        
        
        return out;
    }
    
    void storeInorder(nodeptr node, vector<unsigned int> *inorder) {
        if (node == NULL)
            return;
        
        //left child
        storeInorder(node->left, inorder);
        
        inorder->push_back(node->value);
        
        //right child
        storeInorder(node->right, inorder);
    }
    
    void insert_unsafe(unsigned int value, nodeptr & p, int depth){
        if(p == NULL){
            p = new node;
            p->value = value;
            p->left = NULL;
            p->right = NULL;
            p->height = depth;
            return;
        }
        
        if(value == p->value)
            return;
        else if (value < p->value)
            insert_unsafe(value, p->left , depth + 1);
        else
            insert_unsafe(value, p->right, depth + 1);
        
    }
    
    void printInorder(nodeptr node) {
        if (node == NULL)
            return;
        
        /* first recur on left child */
        printInorder(node->left);
        
        cout << node->value << " ";
        //printf("%d ", node->value);
        
        /* now recur on right child */
        printInorder(node->right);
    }
    
    
    
public:
    
    avlTree() : root(NULL){}
    
    ~avlTree() { delete root; }
    
    void insert(vector <unsigned int> values) {
        ;
        if(root == NULL) {
            root = new node;
            root->value = values[0];
            root->right = NULL;
            root->left = NULL;
            root->height = 0;
            values.erase(values.begin());
            
        }
        
        for(size_t j = 0; j<values.size(); j++)
            insert(values[j], root);
    }
    
    bool check () {
        if (root != NULL) {
            return 0 == check_tree(root);
        } else
            return true;
    }
    
    void storeInorder (vector<unsigned int> *inorder) {
        if (root != NULL) {
            storeInorder(root, inorder);
        }
    }
    
    
    nodeptr sortedArrayToAVL(vector<unsigned int> arr, int start, int end) {
        if (start > end)
            return NULL;
        
        /* Get the middle element and make it root */
        int mid = (start + end) / 2;
        nodeptr tmp = new node;
        tmp->value = arr[mid];
        tmp->right = NULL;
        tmp->left = NULL;
        tmp->height = 0;
        
        
        #pragma omp parallel
        {
            
            #pragma omp single
            {
                /* Recursively construct the left subtree and make it
                left child of root */
                #pragma omp task
                tmp->left =  sortedArrayToAVL(arr, start, mid-1);
                
                /* Recursively construct the right subtree and make it
                right child of root */
                #pragma omp task
                tmp->right = sortedArrayToAVL(arr, mid+1, end);
            
                #pragma omp taskwait
                tmp->height = max(get_height(tmp->left),get_height(tmp->right)) + 1;
            }
        }
        
        return tmp;
    }
    
    avlTree(vector<unsigned int> arr) : root(sortedArrayToAVL(arr, 0, arr.size() - 1)){}
    
    void printInorder() {
        printInorder(root);
        cout << endl;
    }
    
    
};


vector<unsigned int> thread_helper(vector<unsigned int> arr, vector<unsigned int>::const_iterator first, size_t block_size) {
    avlTree tmp;
    vector<unsigned int>::const_iterator last = first + block_size;
    vector<unsigned int> newVec(first, last);
    
    tmp.insert(newVec);
    
    vector<unsigned int> testMe;
    vector<unsigned int> *testMePtr = &testMe;
    tmp.storeInorder(testMePtr);
    //cout << "Helper: " << testMe.size() << endl;
    return testMe;
    
}

std::vector<unsigned int> flatten(const std::vector<std::vector<unsigned int>>& v) {
    std::size_t total_size = 0;
    for (const auto& sub : v)
        total_size += sub.size();
    
    std::vector<unsigned int> result;
    result.reserve(total_size);
    for (const auto& sub : v)
        result.insert(result.end(), sub.begin(), sub.end());
    
    return result;
}

int main(int argc, char* argv[]) {
    if ( argc != 3 ) {
        cout << "Usage: size number-of-threads" << endl;
        return EXIT_FAILURE;
    }
    
    const size_t size = atoi(argv[1]); // array size
    const int threads = atoi(argv[2]); // number of threads
    
    
    vector <unsigned int> value;
    value.reserve(size);
    
    #pragma omp parallel
    {
        random_device rd;
        mt19937 gen(rd() + omp_get_thread_num());
        uniform_int_distribution<unsigned int> dis(0, size / 8);
        
        #pragma omp for
        for (size_t i = 0; i < size; i++) {
            value[i] = dis(gen);
        }
    }
    
    if (threads > 1) {
        //omp_set_nested(1);
        if (omp_get_nested() != 1) {
            cout << "Nested omp not supported" << endl;
        }
        
        const size_t block_size = size / threads;
        //cout << block_size << endl;
        
        vector<vector<unsigned int>> results;
        results.reserve(threads);
        
        #pragma omp parallel num_threads(threads)
        {
            int id = omp_get_thread_num();
            //cout << id << endl;
            if (id != threads - 1) {
                //cout << thread_helper(value, value.begin() + block_size * id, block_size).size() << endl;
                results.push_back(thread_helper(value, value.begin() + block_size * id, block_size));
            }
            else {
                results.push_back(thread_helper(value, value.begin() + block_size * id + (size % threads), block_size));
            }
        }
        
        vector<unsigned int> flat = flatten(results);
        
        //cout << flat.size() << endl;
        sort(flat.begin(), flat.end());
        flat.erase(unique(flat.begin(), flat.end()), flat.end());
        cout << flat.size() << endl;
        
        
        for (int i = 0; i < flat.size() - 1 ; i++) {
            if (flat[i] >= flat[i + 1])
                cout << "error with sorting: " << flat[i] << ">" << flat[i + 1]<< endl;
        }
        
        avlTree out(flat);
        
        cout << out.check() << endl;
        
        if(out.check())
            return EXIT_SUCCESS;
        else
            return EXIT_FAILURE;
    }
    else {
        avlTree out(value);
        
        cout << out.check() << endl;
        
        if(out.check())
            return EXIT_SUCCESS;
        else
            return EXIT_FAILURE;
    }
    
    
    
    
    /*
    avlTree tree;
    vector <unsigned int> test = {1,564,232,576, 576,23232,45646};
    tree.insert(test);
    
    cout << tree.check() << endl;
    
    tree.printInorder();
    
    vector<unsigned int> testMe;
    vector<unsigned int> *testMePtr = &testMe;
    tree.storeInorder(testMePtr);
    
    cout << testMe.size() << "foo" << endl;
    
    for (size_t i = 0; i < testMe.size(); i++)
        cout << testMe[i] << " ";
    cout << endl;
    
    cout << "mark 1" << endl;
    avlTree treeTwo(testMe);
    cout << "mark 2" << endl;
    
    treeTwo.printInorder();
    cout << treeTwo.check() << endl;
    */
    
    
}
    
    
