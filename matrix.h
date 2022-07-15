#pragma  once
#ifndef MATRIX_H
#define MATRIX_H

#include <algorithm>
#include <time.h>

template <typename MatrixRow>
class MatrixRowIterator
{
    public:
    using value_type = typename MatrixRow::value_type;
    using pointer = value_type*;
    using reference = value_type&;

    MatrixRowIterator(pointer ptr):m_ptr(ptr){}

    MatrixRowIterator& operator++()
    {
        m_ptr++;
        return *this;
    }

MatrixRowIterator operator++(int)
{
    MatrixRowIterator it = *this;
    ++(this);
    return it;
}

MatrixRowIterator& operator--()
{
    m_ptr--;
    return *this;
}

MatrixRowIterator operator--(int)
{
    MatrixRowIterator it = *this;
    --(this);
    return it;
}

pointer operator->()
{
    return m_ptr;
}

reference operator*()
{
    return *m_ptr;
}

bool operator==(MatrixRowIterator other)
{
    return this->m_ptr == other.m_ptr;
}

bool operator!=(MatrixRowIterator other)
{
    return this->m_ptr != other.m_ptr;
}

private:
pointer m_ptr;

};


#endif 