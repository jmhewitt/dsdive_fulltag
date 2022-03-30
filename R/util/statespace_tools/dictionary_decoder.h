#ifndef STATESPACE_TOOLS_DICTIONARY_DECODER_H
#define STATESPACE_TOOLS_DICTIONARY_DECODER_H

// [[Rcpp::depends(RcppEigen)]]

#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;

/**
 * The templating for the state space methods is useful because it only requires 
 * an iterator over the matrices and likelihood functions required.  The things
 * could be explicitly placed in a std::vector, or they can be placed in more
 * sophisticated classes.  All that is required is the ability to work with an 
 * iterator and return objects compatible with Eigen::MatrixXd addition and 
 * multiplication.
 */

/**
 * Sequentially decompress EntryType obj's into corresponding dictionary entries
 * of pairs (EntryType current, EntryType next)
 *
 * @tparam Dict
 * @tparam EntryType
 * @tparam EntryIdType
 */
template<typename Dict, typename EntryType, typename EntryIdType>
struct MarkovDictIter {

    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = EntryType;
    using pointer = value_type*;
    using reference = value_type&;

    using entry_pair = std::pair<EntryIdType, EntryIdType>;

    MarkovDictIter<Dict, EntryType, EntryIdType>(
        EntryIdType* it, Dict *dict
    ) : m_it(it+1), m_it_last(it), m_dict(dict) { }

    // return dictionary entry associated w/current pos. of a pointer into x
    reference operator*() { return (*m_dict)[entry_pair(*m_it_last, *m_it)]; }
    pointer operator->() { return &((*m_dict)[entry_pair(*m_it_last, *m_it)]); }

    // prefix increment
    MarkovDictIter<Dict, EntryType, EntryIdType>& operator++() {
        m_it_last = m_it;
        m_it++;
        return *this;
    }
    // postfix increment
    MarkovDictIter<Dict, EntryType, EntryIdType> operator++(int) {
        MarkovDictIter<Dict, EntryType, EntryIdType> tmp = *this;
        ++(*this);
        return tmp;
    }

    friend bool operator== (
        const MarkovDictIter<Dict, EntryType, EntryIdType>& a,
        const MarkovDictIter<Dict, EntryType, EntryIdType>& b
    ) {
        return a.m_it_last == b.m_it;
    };
    friend bool operator!= (
        const MarkovDictIter<Dict, EntryType, EntryIdType>& a,
        const MarkovDictIter<Dict, EntryType, EntryIdType>& b
    ) {
        return a.m_it_last != b.m_it;
    };

private:

    EntryIdType *m_it, *m_it_last;
    Dict *m_dict;

};

/**
 * Sequentially decompress EntryType obj's into corresponding dictionary entries
 * @tparam Dict
 * @tparam EntryType
 * @tparam EntryIdType
 */
template<typename Dict, typename EntryType, typename EntryIdType>
struct DictionaryIterator {

    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = EntryType;
    using pointer = value_type*;
    using reference = value_type&;

    DictionaryIterator<Dict, EntryType, EntryIdType>(
        EntryIdType* it, Dict *dict
    ) : m_it(it), m_dict(dict) { }

    // return dictionary entry associated w/current pos. of a pointer into x
    reference operator*() { return (*m_dict)[*m_it]; }
    pointer operator->() { return &((*m_dict)[*m_it]); }

    // prefix increment
    DictionaryIterator<Dict, EntryType, EntryIdType>& operator++() {
        m_it++;
        return *this;
    }
    // postfix increment
    DictionaryIterator<Dict, EntryType, EntryIdType> operator++(int) {
        DictionaryIterator<Dict, EntryType, EntryIdType> tmp = *this;
        ++(*this);
        return tmp;
    }

    friend bool operator== (
        const DictionaryIterator<Dict, EntryType, EntryIdType>& a,
        const DictionaryIterator<Dict, EntryType, EntryIdType>& b
    ) {
        return a.m_it == b.m_it;
    };
    friend bool operator!= (
        const DictionaryIterator<Dict, EntryType, EntryIdType>& a,
        const DictionaryIterator<Dict, EntryType, EntryIdType>& b
    ) {
        return a.m_it != b.m_it;
    };

private:

    EntryIdType* m_it;
    Dict *m_dict;

};

template<typename Container, typename ContainerValue>
struct SubscriptIterator {

        using iterator_category = std::forward_iterator_tag;
        using difference_value = std::ptrdiff_t;
        using value_type = ContainerValue;
        using pointer = value_type*;
        using reference = value_type&;

        SubscriptIterator(Container & c, std::size_t ind) :
            m_container(c), m_ind(ind) { }

        // return mapped element associated w/current pos. of a pointer into x
        reference operator*() { return m_container[m_ind]; }
        pointer operator->() { return &(m_container[m_ind]); }

        // prefix increment
        SubscriptIterator& operator++() {
            m_ind++;
            return *this;
        }

        // postfix increment
        SubscriptIterator operator++(int) {
            SubscriptIterator tmp = *this;
            ++(*this);
            return tmp;
        }

        friend bool operator== (
            const SubscriptIterator& a,
            const SubscriptIterator& b
        ) {
            return a.m_ind == b.m_ind;
        }

        friend bool operator!= (
            const SubscriptIterator& a,
            const SubscriptIterator& b
        ) {
            return a.m_ind != b.m_ind;
        }

    private:

        Container & m_container;
        std::size_t m_ind;

    };

/**
 * Decompress a sequence of objects stored via dictionary compression
 * 
 * DictionaryDecoder provides a custom iterator that generates a sequence of 
 * EntryType objects corresponding to elements of an array of EntryIdType
 * objects.  In context, DictionaryDecoder decompresses a sequence of EntryType
 * objects that were previously compressed into an array of EntryIdType
 * objects using dictionary compression.
 * 
 * Dict is a "dictionary" class that behaves like a key-value store. Dict
 * objects associate an EntryIdType object (i.e., unsigned int) with an 
 * EntryType object (i.e., Eigen::VectorXd) using the function
 * operator[](EntryIdType), which may be overloaded.  Examples of simple 
 * Dict objects include std::map and std::vector objects, which do not require
 * operator overloading.  More complicated Dict objects can be designed that
 * return an Eigen::Map<Eigen::VectorXd> object for each column of an 
 * Eigen::MatrixXd object.
 */
template<typename Dict, typename EntryType, typename EntryIdType = std::size_t,
         typename Iterator = DictionaryIterator<Dict, EntryType, EntryIdType> >
class DictionaryDecoder {
  
  Dict *dictionary;

  // addresses of the start and (one past) end of a compressed data array, x
  EntryIdType *xbegin, *xend;
  
public:

  /**
   * DictionaryDecoder as wrapper for externally-compressed data
   * @param d dictionary for sequence to decompress
   * @param x std::vector containing EntryIdType array
   */
  DictionaryDecoder(Dict &d, std::vector<EntryIdType> &x) :
    dictionary(&d), xbegin(&(*x.begin())), xend(&(*x.end()))  { }

  /**
   * DictionaryDecoder as wrapper for externally-compressed data
   * @param d dictionary for sequence to decompress
   * @param x address of the start of EntryIdType array
   * @param n number of elements in array
   */
  DictionaryDecoder(Dict &d, EntryIdType* x, unsigned int n) :
    dictionary(&d), xbegin(x), xend(x+n)  { }

  typedef Iterator iterator;
  
  Iterator begin() { return Iterator(xbegin, dictionary); };
  Iterator end() { return Iterator(xend, dictionary); };
  
};

/**
 *  Wrapper to access columns of an Eigen::MatrixXd object via 
 *  Eigen::Map<Eigen::VectorXd> objects.  Facilitates dictionary decompression 
 *  without creating temporary copies of matrix columns.
 *  
 *  Assumes matrices are stored in column-major format.
 */
struct ColumnMapper {

  double *data;
  unsigned int rows;
  
  Eigen::Map<Eigen::VectorXd> v;
  
  ColumnMapper(Eigen::MatrixXd & mat) : data(mat.data()), rows(mat.rows()),
    v(data, rows) {}

  ColumnMapper(double* mat, unsigned int nrows) : data(mat), rows(nrows),
    v(data, rows) {}
  
  Eigen::Map<Eigen::VectorXd>& operator[](unsigned int colind) {
    // change the mapped array using "placement new" syntax, see:
    // https://eigen.tuxfamily.org/dox/group__TutorialMapClass.html#TutorialMapPlacementNew
    new (&v) Eigen::Map<Eigen::VectorXd>(data + rows * colind, rows);
    return v;
  }

  // support decompression via arbitrary inputs that can cast to int
  template<typename Scalar>
  Eigen::Map<Eigen::VectorXd>& operator[](Scalar colind) {
      return (*this)[(unsigned int)(colind)];
  }

};

/**
 * Wrapper to access matrices stored as a 3-dimensional column major array
 * via Eigen::Map<Eigen::MatrixXd> objects.  Facilitates dictionary
 * decompression without creating temporary copies of matrices stored in tensor
 * format.
 *
 * Assumes data is stored in column-major format, which implies the associated
 * indices access matrix coefficients via [row,col,matrix].
 */
struct MatrixMapper {

    double *data;
    unsigned int rows, cols, nelem, nmats;

    typedef Eigen::Map<Eigen::MatrixXd> MappedMat;
    MappedMat m;

    MatrixMapper(double* mat_array, int nrows, int ncols) : data(mat_array),
        rows(nrows), cols(ncols), nelem(rows*cols), nmats(0),
        m(data, rows, cols) { }

    MatrixMapper(double* mat_array, int nrows, int ncols, int nmatrices) :
        data(mat_array), rows(nrows), cols(ncols), nelem(rows*cols),
        nmats(nmatrices), m(data, rows, cols) { }

    MappedMat& operator[](unsigned int matind) {
        new (&m) MappedMat(data + nelem * matind, rows, cols);
        return m;
    }

    // support decompression via arbitrary inputs that can cast to int
    template<typename Scalar>
    MappedMat& operator[](Scalar matind) {
        return (*this)[(unsigned int)(matind)];
    }

    typedef SubscriptIterator<MatrixMapper, MappedMat> iterator;

    iterator begin() { return iterator(*this, 0); }
    iterator end() { return iterator(*this, nmats); }
};

/**
 * Wrapper to access columns of a 3d array via Eigen::Map<Eigen::VectorXd>
 * objects.  Facilitates dictionary decompression without creating temporary
 * copies of matrix elements.
 *
 * Assumes the source array is stored in column-major format with indices
 * (x,y,z) with dimensions (m,n,p) which implies there are m values for
 * x = 0,...,m-1, n values for y = 0,...,n-1, and p values for z = 0,...,p-1.
 *
 * The Eigen::Map<Eigen::VectorXd> that will be returned will iterate over all
 * x values (0,y,z),...,(m-1,y,z) in a "column" for a specific value of y and z.
 * The values of y and z are wrapped in an std::pair object, passed as an
 * argument to operator[]
 */
struct ColumnMapper3 {

    double *data;
    unsigned int m, n, nelem;

    Eigen::Map<Eigen::VectorXd> v;

    ColumnMapper3(
        double* mat_array, unsigned int rows, unsigned int cols
    ) : data(mat_array), m(rows), n(cols), nelem(m*n), v(data, m) { }

    Eigen::Map<Eigen::VectorXd>& operator[](
        const std::pair<unsigned int, unsigned int>& ind
    ) {
        new (&v) Eigen::Map<Eigen::VectorXd>(
            data + m * ind.first + nelem * ind.second, m
        );
        return v;
    }
};

#endif