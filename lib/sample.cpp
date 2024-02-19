#include "sample.h"

int rfd_decode(uint8_t val) {
    int res = 0;
    val &= 3;
    res = (val == 0) ? 1 : (val == 1) ? 3 : (val == 2) ? -3 : -1;
    return res;
//    return res > 0 ? 1 : -1; // 1 бит
}
