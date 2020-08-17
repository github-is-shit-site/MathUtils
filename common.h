#pragma once

template <typename T>
inline T abs(const T& val)
{
    T neg = -val;
    return neg > val ? neg : val;
}
