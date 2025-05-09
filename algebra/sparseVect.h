#ifndef SPARSEVECT_H
#define SPARSEVECT_H

/** \file sparseVect.h
 * \brief sparse vector
a write sparse vector is a collection of v_coeff, which is a couple composed of an index and a double value
to populate with coefficients the sparse vector, use insert method
If several v_coeff have the same index, then they are automatically summed up.
a read sparse vector is a std::vector of v_coeff, built by its constructor
 **/

#include <vector>
#include <iostream>
#include <algorithm>

namespace algebra
{
/** elementary node */
template<typename T>
struct node_t
    {
    T *data;
    node_t *left;
    node_t *right;
    };

/** tree to store sparse vector coefficient */
template<typename T>
class Tree
{
public:
    Tree(): root(NULL) {}
    ~Tree() { destroy_tree(root); }

    void insert(T data);
    node_t<T> *search(int _i) const;
    void inorder_print();
    void postorder_print();
    void preorder_print();
    void inorder_insert(std::vector<T> &v);

private:
    void destroy_tree(node_t<T> *leaf);
    void insert(T data, node_t<T> *leaf);
    node_t<T> *search(int _i, node_t<T> *leaf) const;
    void inorder_print(node_t<T> *leaf);
    void postorder_print(node_t<T> *leaf);
    void preorder_print(node_t<T> *leaf);
    void inorder_insert(node_t<T> *leaf, std::vector<T> &v);
    node_t<T> *root;
};

/** recursive template function to delete the tree */
template<typename T>
void Tree<T>::destroy_tree(node_t<T> *leaf)
    {
    if(leaf != NULL)
        {
        destroy_tree(leaf->left);
        destroy_tree(leaf->right);
        delete leaf->data;
        delete leaf;
        }
    }

/** leaf inserter */
template<typename T>
void Tree<T>::insert(T data, node_t<T> *leaf)
    {
    if (data._i == leaf->data->_i)
        {
        leaf->data->add(data.getVal());
        return;
        }

    if (data._i < leaf->data->_i)
        {
        if (leaf->left != NULL) { insert(data, leaf->left); }
        else
            {
            leaf->left = new node_t<T>;
            if (!leaf->left) exit(1);
            T* data_ptr = new T(data._i, data.getVal());
            if (!data_ptr) exit(1);
            leaf->left->data  = data_ptr;
            leaf->left->left  = NULL;
            leaf->left->right = NULL;
            }
        return;
        }

    if (data._i > leaf->data->_i)
        {
        if (leaf->right != NULL) { insert(data, leaf->right); }
        else
            {
            leaf->right = new node_t<T>;
            if (!leaf->right) exit(1);
            T* data_ptr = new T(data._i, data.getVal());
            if (!data_ptr) exit(1);
            leaf->right->data  = data_ptr;
            leaf->right->right = NULL;
            leaf->right->left  = NULL;
            }
        return;
        }
    }

/** coeff inserter */
template<typename T>
void Tree<T>::insert(T data)
    {
    if(root != NULL)
        { insert(data, root); }
    else
        {
        root = new node_t<T>;
        if (!root) exit(1);
        root->data = new T(data._i, data.getVal());
        root->left = NULL;
        root->right = NULL;
        }
    }

/** search template function for the first occurence of index _i from position leaf */
template<typename T>
node_t<T> *Tree<T>::search(int _i, node_t<T> *leaf) const
    {
    if(leaf != NULL)
        {
        if(_i == leaf->data->_i)
            { return leaf; }

        if(_i < leaf->data->_i)
            { return search(_i, leaf->left); }
        else
            { return search(_i, leaf->right); }
        }
    else
        { return NULL; }
    }

/** search template function for the first occurence of index _i in the whole tree */
template<typename T>
node_t<T> *Tree<T>::search(int _i) const { return search(_i, root); }

/** printing function */
template<typename T>
void Tree<T>::inorder_print()
    {
    inorder_print(root);
    std::cout << std::endl;
    }

/** printing function */
template<typename T>
void Tree<T>::inorder_print(node_t<T> *leaf)
    {
    if(leaf != NULL)
        {
        inorder_print(leaf->left);
        std::cout << "{" << leaf->data->_i<< ":"<<leaf->data->getVal() << "},";
        inorder_print(leaf->right);
        }
    }

/** printing function */
template<typename T>
void Tree<T>::postorder_print()
    {
    postorder_print(root);
    std::cout << std::endl;
    }

/** printing function */
template<typename T>
void Tree<T>::postorder_print(node_t<T> *leaf)
    {
    if(leaf != NULL)
        {
        inorder_print(leaf->left);
        inorder_print(leaf->right);
        std::cout << "{" << leaf->data->_i<< ":"<<leaf->data->getVal() << "},";
        }
    }

/** printing function */
template<typename T>
void Tree<T>::preorder_print()
    {
    preorder_print(root);
    std::cout << std::endl;
    }

/** printing function */
template<typename T>
void Tree<T>::preorder_print(node_t<T> *leaf)
    {
    if(leaf != NULL)
        {
        std::cout << "{" << leaf->data->_i<< ":"<< leaf->data->getVal() << "},";
        inorder_print(leaf->left);
        inorder_print(leaf->right);
        }
    }

/** inserter for a bunch of coefficient from position root */
template<typename T>
void Tree<T>::inorder_insert(std::vector<T> &v) { inorder_insert(root, v); }

/** inserter for a bunch of coefficient from position leaf */
template<typename T>
void Tree<T>::inorder_insert(node_t<T> *leaf, std::vector<T> &v)
    {
    if(leaf != NULL)
        {
        inorder_insert(leaf->left, v);
        v_coeff data{(int)(leaf->data->_i), leaf->data->getVal()}; //CT : cast to int to silent narrowing
        v.push_back(data);
        inorder_insert(leaf->right, v);
        }
    }

/** overloading of << for printing functions */
template<typename T>
std::ostream& operator << (std::ostream& os, const std::vector<T>& v)
    {
    os << "[";
    for (auto it = v.begin(); it != v.end(); ++it)
        { os << "{" << it->_i << ":" << it->getVal() << "},"; }
    os << "]";
    return os;
    }

/**
\class w_sparseVect
it is a container for v_coeff, in writing mode, using template<v_coeff> tree
*/
class w_sparseVect
{
    friend class r_sparseVect;
public:
    /** inserter with a coefficient */
    inline void insert(v_coeff coeff) {tree.insert(coeff);}

    /** return true if the coefficient exists */
    inline bool exist(const int &idx) const { return (tree.search(idx) != NULL); }

    /** getter for the value of a coefficient of index idx, if several coeffs have the same index then it returns the value of the first occurence */
    inline double getVal(const int &idx)
        {
        node_t< v_coeff > *node=tree.search(idx);
        if (node==NULL) return 0.0;
        return node->data->getVal();
        }

private:
    /** coeffs container */
    Tree< v_coeff > tree;
}; // end class w_sparseVect


/**
\class r_sparseVect
read sparse vector : it is a std::vector container for v_coeff, in reading mode
*/
class r_sparseVect: public std::vector<v_coeff>
{
public:
    /** default constructor */
    r_sparseVect(): std::vector<v_coeff>() {}

    /** constructor from a write sparse vector */
    r_sparseVect(w_sparseVect &v): std::vector<v_coeff>() { collect(v); }

    /** return true if the coefficient exists */
    inline bool exist(const int idx) const
        {
        return ( std::find_if(begin(),end(),[this,&idx](v_coeff coeff){return (coeff._i == idx);}) != end() );
        }

    /** collect method is sorting all v_coeffs, eventually with redundant indices, and is summing coeffs with same indices. It removes the coeffs that have been summed. */
    inline void collect(w_sparseVect &v)
        {
        clear();
        v.tree.inorder_insert(*this);
        }

    /** getter for the value of a coefficient of index idx
    if several coeffs have the same index then it returns the value of the first occurence
    return zero if coeffciet of index idx does not exist
     */
    inline double getVal(const int idx) const
        {
        double val(0);
        auto it = std::find_if(begin(),end(),[this,&idx](v_coeff coeff){return (coeff._i == idx); } );
        if (it != end()) val = it->getVal();
        return val;
        }

    /** scalar product */
    inline double dot(const std::vector<double> & X) const
        {
        double val(0);
        for(auto it=begin();it!=end();++it)
            { if(it->_i < (int)(X.size()) ) { val += it->getVal()*X[it->_i]; } }
        return val;
        }

    /** printing function */
    inline void print(std::ostream & flux) const
        {
        flux<<'{';
        std::for_each(begin(),end(), [&flux](const v_coeff &c)
            { flux << '{' << c._i << ':' << c.getVal() <<'}';});
        flux<<"}\n";
        }
}; // end class r_sparseVect

} // end namespace algebra

#endif
