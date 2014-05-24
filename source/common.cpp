#include <cstdlib>
#include <iostream>
#include <stdint.h>

using namespace std;

void* reallocAndZero(void* ptr, size_t oldSize, size_t size){
    ptr = realloc(ptr, size);
    if(ptr == NULL)
        cerr << "ERROR: cant realloc" << endl;
    if(size > oldSize)
        for(uint64_t i = oldSize; i < size; i++)
            ((uint8_t*)ptr)[i] = 0;
    return ptr;
}
