#include <string>
#include "utils.h"

std::string trans(std::string a)
{
    std::string b = "";
    int len = a.length();
    for(int i = len - 1; i >= 0; i--)
    {
        if(a[i] == 'A')
            b = b + 'T';
        if(a[i] == 'G')
            b = b + 'C';
        if(a[i] == 'C')
            b = b + 'G';
        if(a[i] == 'T')
            b = b + 'A';
        if(a[i] == 'N')
            b = b + 'N';
    }

    return b;
}