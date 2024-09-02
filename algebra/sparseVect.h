#ifndef SPARSEVECT_H
#define SPARSEVECT_H

/** \file sparseVect.h 
 * \brief sparse vector
a sparse vector is a collection of v_coeff, which is a couple composed of an index and a value
to populate with coefficients the sparse vector, use insert method
If several v_coeff have the same index, then they are automatically summed up.
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
sparse vector : it is a container for v_coeff, in reading mode
*/
class r_sparseVect
{
public:
	/** default constructor */
	inline r_sparseVect() {collected = false;}

	/** constructor */
	inline r_sparseVect(w_sparseVect &v) { collect(v); }

    /** return true if the coefficient exists */
	inline bool exist(const int &idx) const
		{ 
		auto it = std::find_if(x.begin(),x.end(),[this,&idx](v_coeff coeff){return (coeff._i == idx); } ); 
		return (it != x.end());
		}

	/** collect method is sorting all v_coeffs, eventually with redundant indices, and is summing coeffs with same indices. It removes the coeffs that have been summed. */
	inline void collect(w_sparseVect &v)
		{
        x.clear();
        v.tree.inorder_insert(x);	
		collected = true;		
		}

	/** getter for emptyness of the container of the coeffs */
	inline bool isEmpty(void) const {return x.empty();} 

	/** getter for collected */
	inline bool isCollected(void) const {return collected;}

	/** getter for the value of a coefficient of index idx, if several coeffs have the same index then it returns the value of the first occurence */
	inline double getVal(int idx) const
		{
		double val(0);
		auto it = std::find_if(x.begin(),x.end(),[this,&idx](v_coeff coeff){return (coeff._i == idx); } ); 
		if (it != x.end()) val = it->getVal();		
		return val;		
		}

    /** getter for a reference to the value of a coefficient */
    inline double & getValRef(int idx)
		{
		auto it = std::find_if(x.begin(),x.end(),[this,&idx](v_coeff coeff){return (coeff._i == idx); } ); 
		return it->valRef();// carefull might be out of bounds when it == x.end()	
		}

	/** setter for the value of a coefficient of index idx, all coeffs must have a unique idx, call collect() method before if needed */
	inline void setVal(const int idx,const double val)
		{
		if (collected)
			{
			auto it = std::find_if(x.begin(),x.end(),[this,&idx](v_coeff coeff){return (coeff._i == idx); } );
			if (it != x.end()) it->setVal(val);
			}
		}

	/** scalar product */
    inline double dot(const std::vector<double> & X) const
    	{
    	double val(0);
    	if (!isCollected()) {std::cout << "warning : cannot dot on an uncollected sparseVect" << std::endl;exit(1);}
    	else
    		{
    		for(auto it=x.begin();it!=x.end();++it)
    			{ if(it->_i < (int)(X.size()) ) { val += it->getVal()*X[it->_i]; } }
    		}	
    	return val;
    	}

	/** printing function */
	inline void print(std::ostream & flux) const
	{ flux<<'{'; std::for_each(x.begin(),x.end(), [&flux](const v_coeff &c){ flux << '{' << c._i << ':' << c.getVal() <<'}';}); flux<<"}\n"; }

        /** iterators */
        std::vector< v_coeff >::iterator begin() { return x.begin(); }
        std::vector< v_coeff >::iterator end()   { return x.end();   }
        std::vector< v_coeff >::const_iterator begin() const { return x.begin(); }
        std::vector< v_coeff >::const_iterator end()   const { return x.end();   }
        std::vector< v_coeff >::const_iterator cbegin() const { return x.cbegin(); }
        std::vector< v_coeff >::const_iterator cend()   const { return x.cend();   }

private:
	/** coeffs container */
    std::vector< v_coeff > x;

    /** if true the coeffs have been collected */
    bool collected;
}; // end class r_sparseVect


/** operator<< for r_sparseVect */
inline std::ostream & operator<<(std::ostream & flux, r_sparseVect const& v) {v.print(flux); return flux;}

} // end namespace algebra

#endif
