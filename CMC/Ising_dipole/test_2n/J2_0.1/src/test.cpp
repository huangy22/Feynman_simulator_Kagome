//
//  test.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/4/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "test.h"
#include "utility/dictionary.h"

using namespace std;

#define TEST(func)                  \
    {                               \
        if (EXIT_SUCCESS != func()) \
            exit(0);                \
    }
int RunTest()
{
    LOGGER_CONF("test.log", "test", Logger::file_on | Logger::screen_on, INFO, INFO);
    //    TEST(TestDictionary);

    return 0;
}
