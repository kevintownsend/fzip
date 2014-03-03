#if !defined(ARITHMETIC_FLOATING_POINT_COMPRESSION_H)
#define ARITHMETIC_FLOATING_POINT_COMPRESSION_H
#include <stdint.h>
typedef struct{
    uint64_t count;
    uint64_t codes;
    uint64_t codeLengths;
    uint64_t codeStream;
    uint64_t leftOverStream;
    uint64_t size;
}compressedArithmeticFloatingPointPackage_t;

typedef struct{
    uint64_t count;
    uint64_t* codes;
    uint8_t* codeLengths;
    uint16_t* codeStream;
    void* leftOverStream;
    uint64_t leftOverSize;
}compressedArithmeticFloatingPoint_t;

compressedArithmeticFloatingPointPackage_t* compressArithmeticFloatingPoint(double* data, uint64_t count, compressedArithmeticFloatingPoint_t* cfp);
compressedArithmeticFloatingPointPackage_t* compressArithmeticFloatingPoint(double* data, uint64_t count);
void decompressArithmeticFloatingPoint(double* dest, compressedArithmeticFloatingPoint_t* source);
void decompressArithmeticFloatingPoint(double* dest, compressedArithmeticFloatingPointPackage_t* source);

struct Investor{
    uint64_t value;
    uint32_t suitNumber;
    uint16_t depth;
    uint16_t level;
};
struct JilesPackage_t{
    uint64_t count;
    uint64_t lastLevel;
    uint64_t suits;
    uint64_t finalRepeats;
    uint64_t codeStream;
    uint64_t leftOverStream;
    uint64_t size;
};
struct Jiles_t{
    uint64_t count;
    Investor* suits;
    uint64_t suitCount;
    void* codeStream;
    void* leftOverStream;
};
JilesPackage_t* compressJilesFloatingPoint(double* data, uint64_t count);
double* decompressJilesFloatingPoint(JilesPackage_t* package);
#endif
