#pragma once

#include <complex>

#include "MCC_GalileoCodes.h"

template <class T>
class E1_SearchReference_Pilot
{
public:

    // n в диапазоне [0, 8183]
    T operator() (int code, int n)
    {
        return (MCC_C1_CODE[code][n / 64] & (1 << (31 - ((n / 2) % 32))) ? -1 : 1) * ((n % 2) ? 1 : -1);
    }
};

template <class T>
class E1_SearchReference_Data
{
public:

    // n в диапазоне [0, 8183]
    T operator() (int code, int n)
    {
        return (MCC_B1_CODE[code][n / 64] & (1 << (31 - ((n / 2) % 32))) ? -1 : 1) * ((n % 2) ? 1 : -1);
    }
};