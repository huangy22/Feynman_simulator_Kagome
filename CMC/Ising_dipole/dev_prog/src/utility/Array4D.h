
//
//  vector.h
//  Feynman_Simulator
//
//  Created by yuan on 10/14/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include <sstream>
#include <initializer_list>

using namespace std;

template <class T, size_t dim1, size_t dim2, size_t dim3, size_t dim4>
class Array4D
{
public:
    Array4D()
        : buffer(dim1*dim2*dim3*dim4)
    {
    }

    inline T& at(unsigned int x1, unsigned int x2, unsigned int x3, unsigned int x4)
    {
        return buffer[x4*dim3*dim2*dim1 + x3*dim2*dim1 + x2*dim1 + x1];
    }

    inline const T& at(unsigned int x1, unsigned int x2, unsigned int x3, unsigned int x4) const
    {
        return buffer[x4*dim3*dim2*dim1 + x3*dim2*dim1 + x2*dim1 + x1];
    }

private:
    std::vector<T> buffer;
};
