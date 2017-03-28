
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

template <class T, size_t width, size_t height>
class Array2D
{
public:
    //const int width = W;
    //const int height = H;

    Array2D()
        : buffer(width*height)
    {
    }

    inline T& at(unsigned int x, unsigned int y)
    {
        return buffer[y*width + x];
    }

    inline const T& at(unsigned int x, unsigned int y) const
    {
        return buffer[y*width + x];
    }

private:
    std::vector<T> buffer;
};
