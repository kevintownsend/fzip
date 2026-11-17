/**
 * Copyright 2014 Kevin Townsend
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef FZIP_H
#define FZIP_H

#include <algorithm>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

typedef unsigned long long ull;

static_assert(sizeof(ull) == 8, "ull must be 64 bits");
static_assert(sizeof(double) == 8, "IEEE-754 64-bit double required");

namespace fzip_detail {

static const std::uint8_t kMagic[4] = {'F', 'Z', 'I', 'P'};
static const std::uint16_t kFormatVersion = 2;
static const std::size_t kHeaderSize = 44;
static const int kMaxAllowedCodeLength = 15;
static const int kMaxCommonBits = 13;

inline ull doubleToBits(double value) {
    ull bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    return bits;
}

inline double bitsToDouble(ull bits) {
    double value = 0.0;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
}

inline bool startsWith(const std::string& s, const char* prefix) {
    std::size_t len = std::strlen(prefix);
    return s.size() >= len && s.compare(0, len, prefix) == 0;
}

inline bool parseIntStrict(const std::string& s, int& out) {
    if (s.empty()) {
        return false;
    }
    char* endPtr = NULL;
    errno = 0;
    long val = std::strtol(s.c_str(), &endPtr, 10);
    if (errno != 0 || endPtr == s.c_str() || *endPtr != '\0') {
        return false;
    }
    if (val < std::numeric_limits<int>::min() || val > std::numeric_limits<int>::max()) {
        return false;
    }
    out = static_cast<int>(val);
    return true;
}

inline int ceilLog2(std::size_t n) {
    int bits = 0;
    while ((std::size_t{1} << bits) < n) {
        ++bits;
    }
    return bits;
}

inline void putU16LE(std::uint8_t* dst, std::uint16_t v) {
    dst[0] = static_cast<std::uint8_t>(v & 0xFFU);
    dst[1] = static_cast<std::uint8_t>((v >> 8) & 0xFFU);
}

inline void putU32LE(std::uint8_t* dst, std::uint32_t v) {
    dst[0] = static_cast<std::uint8_t>(v & 0xFFU);
    dst[1] = static_cast<std::uint8_t>((v >> 8) & 0xFFU);
    dst[2] = static_cast<std::uint8_t>((v >> 16) & 0xFFU);
    dst[3] = static_cast<std::uint8_t>((v >> 24) & 0xFFU);
}

inline void putU64LE(std::uint8_t* dst, std::uint64_t v) {
    for (int i = 0; i < 8; ++i) {
        dst[i] = static_cast<std::uint8_t>((v >> (8 * i)) & 0xFFU);
    }
}

inline std::uint16_t getU16LE(const std::uint8_t* src) {
    return static_cast<std::uint16_t>(src[0]) |
           (static_cast<std::uint16_t>(src[1]) << 8);
}

inline std::uint32_t getU32LE(const std::uint8_t* src) {
    return static_cast<std::uint32_t>(src[0]) |
           (static_cast<std::uint32_t>(src[1]) << 8) |
           (static_cast<std::uint32_t>(src[2]) << 16) |
           (static_cast<std::uint32_t>(src[3]) << 24);
}

inline std::uint64_t getU64LE(const std::uint8_t* src) {
    std::uint64_t v = 0;
    for (int i = 0; i < 8; ++i) {
        v |= static_cast<std::uint64_t>(src[i]) << (8 * i);
    }
    return v;
}

inline std::uint32_t crc32Update(std::uint32_t crc, const std::uint8_t* data, std::size_t len) {
    crc = ~crc;
    for (std::size_t i = 0; i < len; ++i) {
        crc ^= data[i];
        for (int b = 0; b < 8; ++b) {
            crc = (crc & 1U) ? ((crc >> 1) ^ 0xEDB88320U) : (crc >> 1);
        }
    }
    return ~crc;
}

inline void appendBits(std::vector<ull>& stream, ull& buffer, int& endBit, ull bits, int width) {
    if (width <= 0) {
        return;
    }
    if (width < 64) {
        bits &= (1ULL << width) - 1ULL;
    }
    buffer |= (bits << endBit);
    endBit += width;
    if (endBit >= 64) {
        stream.push_back(buffer);
        endBit -= 64;
        int consumed = width - endBit;
        buffer = (consumed < 64) ? (bits >> consumed) : 0ULL;
    }
}

inline ull readBits(const std::vector<ull>& stream, ull bitPos, int width) {
    if (width <= 0) {
        return 0ULL;
    }
    std::size_t wordIdx = static_cast<std::size_t>(bitPos / 64ULL);
    int bitOffset = static_cast<int>(bitPos % 64ULL);
    ull val = (wordIdx < stream.size()) ? (stream[wordIdx] >> bitOffset) : 0ULL;
    if (bitOffset > 0 && bitOffset + width > 64 && (wordIdx + 1) < stream.size()) {
        val |= stream[wordIdx + 1] << (64 - bitOffset);
    }
    if (width < 64) {
        val &= (1ULL << width) - 1ULL;
    }
    return val;
}

inline void packWordsToBytes(const std::vector<ull>& words,
                             ull bitLength,
                             std::vector<std::uint8_t>& outBytes) {
    std::size_t byteCount = static_cast<std::size_t>((bitLength + 7ULL) / 8ULL);
    outBytes.assign(byteCount, 0U);
    for (std::size_t i = 0; i < byteCount; ++i) {
        std::size_t wordIdx = i / 8;
        int byteShift = static_cast<int>((i % 8) * 8);
        if (wordIdx < words.size()) {
            outBytes[i] = static_cast<std::uint8_t>((words[wordIdx] >> byteShift) & 0xFFU);
        }
    }
}

inline void unpackBytesToWords(const std::uint8_t* bytes,
                               std::size_t byteCount,
                               std::vector<ull>& outWords) {
    std::size_t wordCount = (byteCount + 7) / 8;
    outWords.assign(wordCount, 0ULL);
    for (std::size_t i = 0; i < byteCount; ++i) {
        std::size_t wordIdx = i / 8;
        int byteShift = static_cast<int>((i % 8) * 8);
        outWords[wordIdx] |= static_cast<ull>(bytes[i]) << byteShift;
    }
}

inline bool checkTrailingZeroBits(const std::vector<std::uint8_t>& bytes, ull bitPos) {
    std::size_t byteLen = bytes.size();
    ull expectedBytes = (bitPos + 7ULL) / 8ULL;
    if (expectedBytes != static_cast<ull>(byteLen)) {
        return false;
    }
    int remBits = static_cast<int>(static_cast<ull>(byteLen) * 8ULL - bitPos);
    if (remBits > 0 && byteLen > 0) {
        std::uint8_t lastByte = bytes[byteLen - 1];
        std::uint8_t usedMask = static_cast<std::uint8_t>((1U << (8 - remBits)) - 1U);
        if ((lastByte & static_cast<std::uint8_t>(~usedMask)) != 0U) {
            return false;
        }
    }
    return true;
}

inline bool writeAll(const void* data, std::size_t elemSize, std::size_t count, FILE* out) {
    if (count == 0) {
        return true;
    }
    return std::fwrite(data, elemSize, count, out) == count;
}

inline std::uint16_t bitReverseLow(ull code, int len) {
    std::uint16_t lsbCode = 0;
    for (int b = 0; b < len; ++b) {
        if ((code >> (len - 1 - b)) & 1ULL) {
            lsbCode |= static_cast<std::uint16_t>(1U << b);
        }
    }
    return lsbCode;
}

struct ByteSource {
    FILE* file;
    const std::uint8_t* mem;
    std::size_t memSize;
    std::size_t memPos;
    std::uint32_t crc;

    explicit ByteSource(FILE* f) : file(f), mem(NULL), memSize(0), memPos(0), crc(0U) {}
    ByteSource(const std::uint8_t* data, std::size_t size)
        : file(NULL), mem(data), memSize(size), memPos(0), crc(0U) {}

    void resetCrc() { crc = 0U; }

    bool readExact(std::uint8_t* dst, std::size_t count) {
        if (count == 0) {
            return true;
        }
        if (file != NULL) {
            if (std::fread(dst, 1, count, file) != count) {
                return false;
            }
        } else {
            if (mem == NULL || count > memSize - memPos) {
                return false;
            }
            std::memcpy(dst, mem + memPos, count);
            memPos += count;
        }
        crc = crc32Update(crc, dst, count);
        return true;
    }

    bool readExactVec(std::vector<std::uint8_t>& out, std::size_t byteCount) {
        out.clear();
        if (byteCount == 0) {
            return true;
        }
        if (file == NULL) {
            if (mem == NULL || byteCount > memSize - memPos) {
                return false;
            }
            out.assign(mem + memPos, mem + memPos + byteCount);
            memPos += byteCount;
            crc = crc32Update(crc, &out[0], byteCount);
            return true;
        }
        const std::size_t kChunk = 65536;
        std::vector<std::uint8_t> chunk(std::min(byteCount, kChunk));
        std::size_t remaining = byteCount;
        while (remaining > 0) {
            std::size_t want = std::min(remaining, kChunk);
            if (std::fread(&chunk[0], 1, want, file) != want) {
                return false;
            }
            crc = crc32Update(crc, &chunk[0], want);
            out.insert(out.end(), chunk.begin(), chunk.begin() + static_cast<std::ptrdiff_t>(want));
            remaining -= want;
        }
        return true;
    }

    bool atEndIfMemory() const {
        return file != NULL || memPos == memSize;
    }
};

struct NoopFileCloser {
    void operator()(FILE*) const {}
};

struct FileCloser {
    void operator()(FILE* f) const {
        if (f != NULL) {
            std::fclose(f);
        }
    }
};

}  // namespace fzip_detail

struct FzipOptions {
    FILE* inputFile;
    FILE* outputFile;
    std::shared_ptr<FILE> inputHandle;
    std::shared_ptr<FILE> outputHandle;
    bool compress;
    bool valid;
    bool showHelp;
    int maxCodeLength;
    int targetCodeLength;
    int commonCodes;
    int log2CommonCodes;

    FzipOptions(int argc, char* argv[])
        : inputFile(stdin),
          outputFile(stdout),
          inputHandle(stdin, fzip_detail::NoopFileCloser()),
          outputHandle(stdout, fzip_detail::NoopFileCloser()),
          compress(true),
          valid(true),
          showHelp(false),
          maxCodeLength(9),
          targetCodeLength(8),
          commonCodes(8192),
          log2CommonCodes(13) {
        bool seenCompress = false;
        bool seenDecompress = false;
        bool hasTargetCodeLength = false;
        int nonOptionArgs = 0;

        for (int currentArg = 1; currentArg < argc; ++currentArg) {
            std::string arg(argv[currentArg]);
            if (arg == "-h" || arg == "--help") {
                showHelp = true;
                return;
            }
            if (fzip_detail::startsWith(arg, "--")) {
                if (fzip_detail::startsWith(arg, "--MaxCodeLength=")) {
                    std::string valStr = arg.substr(std::strlen("--MaxCodeLength="));
                    if (!fzip_detail::parseIntStrict(valStr, maxCodeLength) ||
                        maxCodeLength < 1 || maxCodeLength > fzip_detail::kMaxAllowedCodeLength) {
                        std::cerr << "Error: --MaxCodeLength must be an integer in [1, "
                                  << fzip_detail::kMaxAllowedCodeLength << "]\n";
                        valid = false;
                        return;
                    }
                } else if (fzip_detail::startsWith(arg, "--TargetCodeLength=")) {
                    std::string valStr = arg.substr(std::strlen("--TargetCodeLength="));
                    if (!fzip_detail::parseIntStrict(valStr, targetCodeLength) ||
                        targetCodeLength < 0 ||
                        targetCodeLength > fzip_detail::kMaxAllowedCodeLength) {
                        std::cerr << "Error: --TargetCodeLength must be an integer in [0, "
                                  << fzip_detail::kMaxAllowedCodeLength << "]\n";
                        valid = false;
                        return;
                    }
                    hasTargetCodeLength = true;
                } else {
                    std::cerr << "Error: unknown option '" << arg << "'\n";
                    valid = false;
                    return;
                }
            } else if (arg.size() > 1 && arg[0] == '-') {
                for (std::size_t i = 1; i < arg.size(); ++i) {
                    if (arg[i] == 'c') {
                        compress = true;
                        seenCompress = true;
                    } else if (arg[i] == 'd') {
                        compress = false;
                        seenDecompress = true;
                    } else if (arg[i] == 'h') {
                        showHelp = true;
                        return;
                    } else {
                        std::cerr << "Error: unknown option '-" << arg[i] << "'\n";
                        valid = false;
                        return;
                    }
                }
                if (seenCompress && seenDecompress) {
                    std::cerr << "Error: cannot specify both -c and -d\n";
                    valid = false;
                    return;
                }
            } else {
                if (nonOptionArgs == 0) {
                    FILE* f = std::fopen(arg.c_str(), "rb");
                    if (f == NULL) {
                        std::cerr << "Error: failed to open input file '" << arg
                                  << "': " << std::strerror(errno) << "\n";
                        valid = false;
                        return;
                    }
                    inputHandle.reset(f, fzip_detail::FileCloser());
                    inputFile = f;
                } else if (nonOptionArgs == 1) {
                    FILE* f = std::fopen(arg.c_str(), "wb");
                    if (f == NULL) {
                        std::cerr << "Error: failed to open output file '" << arg
                                  << "': " << std::strerror(errno) << "\n";
                        valid = false;
                        return;
                    }
                    outputHandle.reset(f, fzip_detail::FileCloser());
                    outputFile = f;
                } else {
                    std::cerr << "Error: too many positional arguments\n";
                    valid = false;
                    return;
                }
                nonOptionArgs++;
            }
        }

        if (!hasTargetCodeLength) {
            targetCodeLength = std::max(0, maxCodeLength - 1);
        } else if (targetCodeLength > maxCodeLength) {
            std::cerr << "Error: --TargetCodeLength (" << targetCodeLength
                      << ") cannot exceed --MaxCodeLength (" << maxCodeLength << ")\n";
            valid = false;
            return;
        }
    }
};

struct FzipCode {
    std::uint8_t codeLength;
    std::uint8_t isCommon;
    std::uint8_t prefixLength;
    std::uint16_t code;
    ull prefix;

    FzipCode() : codeLength(0), isCommon(0), prefixLength(0), code(0), prefix(0) {}

    void print() const {
        std::cerr << "codeLength: " << static_cast<int>(codeLength) << "\n"
                  << "isCommon: " << static_cast<int>(isCommon) << "\n"
                  << "prefixLength: " << static_cast<int>(prefixLength) << "\n"
                  << "code: " << code << "\n"
                  << "prefix: " << std::hex << prefix << std::dec << "\n";
    }
};

struct FzipEncoded {
    std::uint8_t maxCodeLength;
    std::uint8_t commonIndexBits;
    ull elementCount;
    std::vector<FzipCode> symbolTable;
    std::vector<ull> commons;
    std::vector<ull> codeStream;
    std::vector<ull> argumentStream;
    ull codeStreamBitLength;
    ull argumentStreamBitLength;

    FzipEncoded()
        : maxCodeLength(9),
          commonIndexBits(0),
          elementCount(0),
          codeStreamBitLength(0),
          argumentStreamBitLength(0) {}
};

inline void printBits(ull val) {
    for (int i = 0; i < 64; ++i) {
        std::cerr << ((val & (1ULL << i)) ? '1' : '0');
    }
    std::cerr << "\n";
}

inline bool readFloats(FILE* inputFile, std::vector<double>& out) {
    out.clear();
    if (inputFile == NULL) {
        return false;
    }
    std::uint8_t buf[8];
    while (true) {
        std::size_t n = std::fread(buf, 1, sizeof(buf), inputFile);
        if (n == sizeof(buf)) {
            ull bits = fzip_detail::getU64LE(buf);
            out.push_back(fzip_detail::bitsToDouble(bits));
        } else if (n == 0) {
            if (std::ferror(inputFile)) {
                std::cerr << "Error: failed reading input stream\n";
                return false;
            }
            return true;
        } else {
            std::cerr << "Error: input byte count is not a multiple of 8 (" << n
                      << " trailing bytes)\n";
            return false;
        }
    }
}

inline std::vector<double> readFloats(FILE* inputFile) {
    std::vector<double> ret;
    readFloats(inputFile, ret);
    return ret;
}

inline std::vector<FzipCode> createFzipCodes(const std::vector<std::pair<ull, ull> >& values,
                                             std::size_t targetCodeCount);

inline std::vector<FzipCode> huffmanCoding(
    const std::vector<std::pair<ull, FzipCode> >& codeFrequencies,
    int maxCodeLength = 9);

inline bool assignCanonicalCodes(std::vector<FzipCode>& symbols, int maxCodeLength) {
    if (symbols.empty()) {
        return true;
    }
    if (maxCodeLength < 1 || maxCodeLength > fzip_detail::kMaxAllowedCodeLength ||
        symbols.size() > (std::size_t{1} << maxCodeLength)) {
        return false;
    }
    if (symbols.size() == 1) {
        if (symbols[0].codeLength != 1) {
            return false;
        }
        symbols[0].code = 0;
        return true;
    }

    ull kraftSum = 0;
    const ull kraftTarget = 1ULL << maxCodeLength;
    ull nextCode = 0;
    int prevLen = 0;
    for (std::size_t i = 0; i < symbols.size(); ++i) {
        int len = static_cast<int>(symbols[i].codeLength);
        if (len < 1 || len > maxCodeLength || len < prevLen) {
            return false;
        }
        kraftSum += (1ULL << (maxCodeLength - len));
        if (kraftSum > kraftTarget) {
            return false;
        }
        nextCode <<= (len - prevLen);
        if (nextCode >= (1ULL << len)) {
            return false;
        }
        symbols[i].code = fzip_detail::bitReverseLow(nextCode, len);
        ++nextCode;
        prevLen = len;
    }
    return kraftSum == kraftTarget;
}

inline bool buildDecodeLut(const std::vector<FzipCode>& symbols,
                           int maxCodeLength,
                           std::vector<FzipCode>& lut) {
    if (maxCodeLength < 1 || maxCodeLength > fzip_detail::kMaxAllowedCodeLength) {
        return false;
    }
    const std::size_t lutSize = std::size_t{1} << maxCodeLength;
    lut.assign(lutSize, FzipCode());
    for (std::size_t i = 0; i < symbols.size(); ++i) {
        if (symbols[i].codeLength == 0 || symbols[i].codeLength > maxCodeLength) {
            return false;
        }
        std::size_t stride = std::size_t{1} << symbols[i].codeLength;
        for (std::size_t j = symbols[i].code; j < lutSize; j += stride) {
            if (lut[j].codeLength != 0) {
                return false;
            }
            lut[j] = symbols[i];
        }
    }
    return true;
}

inline bool fzipEncode(const std::vector<double>& rawStream,
                       FzipEncoded& encoded,
                       bool staticCommon = false,
                       int maxCodeLength = 9,
                       int targetCodeLength = 8,
                       int maxCommonBits = 13) {
    if (maxCodeLength < 1 || maxCodeLength > fzip_detail::kMaxAllowedCodeLength ||
        targetCodeLength < 0 || targetCodeLength > maxCodeLength ||
        maxCommonBits < 0 || maxCommonBits > fzip_detail::kMaxCommonBits) {
        std::cerr << "Error: invalid compression parameters\n";
        return false;
    }

    encoded.maxCodeLength = static_cast<std::uint8_t>(maxCodeLength);
    encoded.elementCount = static_cast<ull>(rawStream.size());
    if (!staticCommon) {
        encoded.commons.clear();
    }
    encoded.symbolTable.clear();
    encoded.codeStream.clear();
    encoded.argumentStream.clear();
    encoded.codeStreamBitLength = 0;
    encoded.argumentStreamBitLength = 0;
    encoded.commonIndexBits = 0;

    if (rawStream.empty()) {
        if (staticCommon && !encoded.commons.empty()) {
            const std::size_t maxCap =
                (maxCommonBits > 0) ? (std::size_t{1} << maxCommonBits) : 0;
            if (encoded.commons.size() > maxCap) {
                encoded.commons.resize(maxCap);
            }
            encoded.commonIndexBits =
                static_cast<std::uint8_t>(fzip_detail::ceilLog2(encoded.commons.size()));
        }
        return true;
    }

    std::map<ull, ull> frequencies;
    for (std::size_t i = 0; i < rawStream.size(); ++i) {
        frequencies[fzip_detail::doubleToBits(rawStream[i])]++;
    }

    std::vector<std::pair<ull, ull> > sortedByFrequency;
    sortedByFrequency.reserve(frequencies.size());
    for (std::map<ull, ull>::const_iterator it = frequencies.begin(); it != frequencies.end(); ++it) {
        sortedByFrequency.push_back(std::make_pair(it->second, it->first));
    }
    std::sort(sortedByFrequency.begin(), sortedByFrequency.end());

    std::size_t maxAllowedCodes = std::size_t{1} << maxCodeLength;
    std::size_t targetPrefixCodes = std::size_t{1} << targetCodeLength;
    if (targetPrefixCodes >= maxAllowedCodes) {
        targetPrefixCodes = maxAllowedCodes - 1;
    }

    std::vector<FzipCode> prefixCodes = createFzipCodes(sortedByFrequency, targetPrefixCodes);
    std::map<ull, std::size_t> mapToOrigPrefix;
    for (std::size_t i = 0; i < prefixCodes.size(); ++i) {
        mapToOrigPrefix[prefixCodes[i].prefix] = i;
    }

    const std::size_t commonCapacity =
        (maxCommonBits > 0) ? (std::size_t{1} << maxCommonBits) : 0;

    if (!staticCommon) {
        std::vector<std::pair<ull, ull> > candidates;
        for (std::vector<std::pair<ull, ull> >::const_reverse_iterator it =
                 sortedByFrequency.rbegin();
             it != sortedByFrequency.rend() && candidates.size() < commonCapacity; ++it) {
            ull freq = it->first;
            ull val = it->second;
            if (freq < 2ULL) {
                break;
            }
            std::map<ull, std::size_t>::const_iterator ub = mapToOrigPrefix.upper_bound(val);
            if (ub == mapToOrigPrefix.begin()) {
                std::cerr << "Error: value missing from prefix trie partition\n";
                return false;
            }
            std::size_t origIdx = std::prev(ub)->second;
            int prefLen = static_cast<int>(prefixCodes[origIdx].prefixLength);
            if (prefLen < 64 && freq * static_cast<ull>(64 - prefLen) > 64ULL) {
                candidates.push_back(std::make_pair(freq, val));
            }
        }

        int tentativeBits = fzip_detail::ceilLog2(candidates.size());
        for (std::size_t i = 0; i < candidates.size(); ++i) {
            ull freq = candidates[i].first;
            ull val = candidates[i].second;
            std::size_t origIdx = std::prev(mapToOrigPrefix.upper_bound(val))->second;
            int prefLen = static_cast<int>(prefixCodes[origIdx].prefixLength);
            int argSavedPerOcc = (64 - prefLen) - tentativeBits;
            if (argSavedPerOcc > 0 &&
                freq * static_cast<ull>(argSavedPerOcc) > 64ULL) {
                encoded.commons.push_back(val);
            }
        }
    } else if (encoded.commons.size() > commonCapacity) {
        encoded.commons.resize(commonCapacity);
    }

    int commonBits = fzip_detail::ceilLog2(encoded.commons.size());
    encoded.commonIndexBits = static_cast<std::uint8_t>(commonBits);

    std::unordered_set<ull> commonSet;
    std::unordered_map<ull, std::size_t> mapToCommon;
    for (std::size_t i = 0; i < encoded.commons.size(); ++i) {
        commonSet.insert(encoded.commons[i]);
        if (mapToCommon.find(encoded.commons[i]) == mapToCommon.end()) {
            mapToCommon[encoded.commons[i]] = i;
        }
    }

    FzipCode commonCode;
    commonCode.isCommon = 1;
    commonCode.prefixLength = static_cast<std::uint8_t>(64 - commonBits);

    const std::uint8_t commonThresholdPrefixLen = static_cast<std::uint8_t>(64 - commonBits);
    std::vector<ull> origFrequencies(prefixCodes.size() + 1, 0ULL);

    for (std::map<ull, ull>::const_iterator it = frequencies.begin(); it != frequencies.end(); ++it) {
        std::map<ull, std::size_t>::const_iterator ub = mapToOrigPrefix.upper_bound(it->first);
        if (ub == mapToOrigPrefix.begin()) {
            std::cerr << "Error: value missing from prefix trie partition\n";
            return false;
        }
        std::size_t origIdx = std::prev(ub)->second;
        if (!encoded.commons.empty() &&
            prefixCodes[origIdx].prefixLength <= commonThresholdPrefixLen &&
            commonSet.count(it->first)) {
            origFrequencies[prefixCodes.size()] += it->second;
        } else {
            origFrequencies[origIdx] += it->second;
        }
    }

    if (!staticCommon && origFrequencies[prefixCodes.size()] == 0ULL && !encoded.commons.empty()) {
        encoded.commons.clear();
        commonSet.clear();
        mapToCommon.clear();
        commonBits = 0;
        encoded.commonIndexBits = 0;
    }

    std::vector<std::pair<ull, FzipCode> > codeFrequencies;
    std::vector<int> origToActive(prefixCodes.size(), -1);
    int activeCommonIdx = -1;
    for (std::size_t i = 0; i < prefixCodes.size(); ++i) {
        if (origFrequencies[i] > 0) {
            origToActive[i] = static_cast<int>(codeFrequencies.size());
            codeFrequencies.push_back(std::make_pair(origFrequencies[i], prefixCodes[i]));
        }
    }
    if (origFrequencies[prefixCodes.size()] > 0) {
        activeCommonIdx = static_cast<int>(codeFrequencies.size());
        codeFrequencies.push_back(std::make_pair(origFrequencies[prefixCodes.size()], commonCode));
    }

    std::vector<FzipCode> rawSymbols = huffmanCoding(codeFrequencies, maxCodeLength);
    if (rawSymbols.empty()) {
        std::cerr << "Error: failed to construct bounded-length Huffman codes\n";
        return false;
    }

    // Reorder symbols into canonical (codeLength, original index) order so the serialized
    // symbol table directly yields identical canonical codes without storing codewords.
    std::vector<std::size_t> canonicalOrder(rawSymbols.size());
    for (std::size_t i = 0; i < rawSymbols.size(); ++i) {
        canonicalOrder[i] = i;
    }
    std::stable_sort(canonicalOrder.begin(), canonicalOrder.end(),
                     [&](std::size_t a, std::size_t b) {
                         return rawSymbols[a].codeLength < rawSymbols[b].codeLength;
                     });

    std::vector<std::size_t> activeToCanonical(rawSymbols.size());
    encoded.symbolTable.reserve(rawSymbols.size());
    for (std::size_t pos = 0; pos < canonicalOrder.size(); ++pos) {
        std::size_t origActiveIdx = canonicalOrder[pos];
        activeToCanonical[origActiveIdx] = pos;
        encoded.symbolTable.push_back(rawSymbols[origActiveIdx]);
    }
    if (!assignCanonicalCodes(encoded.symbolTable, maxCodeLength)) {
        std::cerr << "Error: canonical code assignment failed\n";
        return false;
    }

    int codeBufferEndBit = 0;
    int argumentBufferEndBit = 0;
    ull codeBuffer = 0;
    ull argumentBuffer = 0;

    for (std::size_t i = 0; i < rawStream.size(); ++i) {
        ull rawBits = fzip_detail::doubleToBits(rawStream[i]);
        std::map<ull, std::size_t>::const_iterator ub = mapToOrigPrefix.upper_bound(rawBits);
        if (ub == mapToOrigPrefix.begin()) {
            std::cerr << "Error: value missing from prefix map during encoding\n";
            return false;
        }
        std::size_t origIdx = std::prev(ub)->second;
        int activeIdx = origToActive[origIdx];
        if (!encoded.commons.empty() &&
            prefixCodes[origIdx].prefixLength <= commonThresholdPrefixLen &&
            commonSet.count(rawBits)) {
            activeIdx = activeCommonIdx;
        }
        if (activeIdx < 0 || static_cast<std::size_t>(activeIdx) >= activeToCanonical.size()) {
            std::cerr << "Error: unmapped symbol during encoding\n";
            return false;
        }

        const FzipCode& sym = encoded.symbolTable[activeToCanonical[static_cast<std::size_t>(activeIdx)]];
        fzip_detail::appendBits(encoded.codeStream, codeBuffer, codeBufferEndBit,
                                static_cast<ull>(sym.code),
                                static_cast<int>(sym.codeLength));

        ull argument = 0;
        int argumentLength = 0;
        if (sym.isCommon) {
            std::unordered_map<ull, std::size_t>::const_iterator cit = mapToCommon.find(rawBits);
            if (cit == mapToCommon.end()) {
                std::cerr << "Error: common value missing from index map\n";
                return false;
            }
            argument = static_cast<ull>(cit->second);
            argumentLength = commonBits;
        } else if (sym.prefixLength < 64) {
            argumentLength = 64 - static_cast<int>(sym.prefixLength);
            argument = (argumentLength < 64) ? (rawBits & ((1ULL << argumentLength) - 1ULL))
                                             : rawBits;
        }
        fzip_detail::appendBits(encoded.argumentStream, argumentBuffer, argumentBufferEndBit,
                                argument, argumentLength);
    }

    encoded.codeStreamBitLength = static_cast<ull>(encoded.codeStream.size()) * 64ULL +
                                  static_cast<ull>(codeBufferEndBit);
    encoded.argumentStreamBitLength = static_cast<ull>(encoded.argumentStream.size()) * 64ULL +
                                      static_cast<ull>(argumentBufferEndBit);

    if (codeBufferEndBit > 0) {
        encoded.codeStream.push_back(codeBuffer);
    }
    if (argumentBufferEndBit > 0) {
        encoded.argumentStream.push_back(argumentBuffer);
    }
    return true;
}

inline bool fzipDecode(const std::vector<FzipCode>& lut,
                       int maxCodeLength,
                       const std::vector<ull>& commons,
                       const std::vector<ull>& codeStream,
                       const std::vector<ull>& argumentStream,
                       ull codeStreamMaxBits,
                       ull argumentStreamMaxBits,
                       bool hasExactCount,
                       ull expectedElements,
                       std::vector<double>& rawStream,
                       ull* outCodeBitsUsed = NULL,
                       ull* outArgBitsUsed = NULL) {
    rawStream.clear();
    if (maxCodeLength < 1 || maxCodeLength > fzip_detail::kMaxAllowedCodeLength ||
        lut.size() != (std::size_t{1} << maxCodeLength)) {
        std::cerr << "Error: invalid decode lookup table\n";
        return false;
    }

    if ((hasExactCount && expectedElements == 0) || (!hasExactCount && codeStreamMaxBits == 0)) {
        if (argumentStreamMaxBits != 0) {
            if (!hasExactCount) {
                std::cerr << "Error: argument stream non-empty for empty code stream\n";
                return false;
            }
        }
        if (outCodeBitsUsed) *outCodeBitsUsed = 0;
        if (outArgBitsUsed) *outArgBitsUsed = 0;
        return true;
    }

    if (codeStreamMaxBits > static_cast<ull>(codeStream.size()) * 64ULL ||
        argumentStreamMaxBits > static_cast<ull>(argumentStream.size()) * 64ULL) {
        std::cerr << "Error: bitstream bit length exceeds buffer bounds\n";
        return false;
    }

    if (hasExactCount) {
        rawStream.reserve(static_cast<std::size_t>(expectedElements));
    }

    ull currCodeBit = 0;
    ull currArgBit = 0;
    ull decodedCount = 0;

    while (hasExactCount ? (decodedCount < expectedElements) : (currCodeBit < codeStreamMaxBits)) {
        if (currCodeBit >= codeStreamMaxBits) {
            std::cerr << "Error: code stream truncated before elementCount reached\n";
            return false;
        }
        ull peek = fzip_detail::readBits(codeStream, currCodeBit, maxCodeLength);
        const FzipCode& code = lut[static_cast<std::size_t>(peek)];
        if (code.codeLength == 0 || code.codeLength > maxCodeLength) {
            std::cerr << "Error: invalid Huffman code encountered\n";
            return false;
        }
        if (currCodeBit + code.codeLength > codeStreamMaxBits) {
            std::cerr << "Error: code stream truncated mid-codeword\n";
            return false;
        }
        currCodeBit += code.codeLength;

        if (code.prefixLength > 64) {
            std::cerr << "Error: invalid prefix length\n";
            return false;
        }
        int argBits = 64 - static_cast<int>(code.prefixLength);
        if (currArgBit + static_cast<ull>(argBits) > argumentStreamMaxBits) {
            std::cerr << "Error: argument stream truncated\n";
            return false;
        }
        ull arg = fzip_detail::readBits(argumentStream, currArgBit, argBits);
        currArgBit += static_cast<ull>(argBits);

        ull value = code.prefix;
        if (code.isCommon) {
            if (arg >= commons.size()) {
                std::cerr << "Error: common table index out of bounds\n";
                return false;
            }
            value = commons[static_cast<std::size_t>(arg)];
        } else if (argBits > 0) {
            value |= arg;
        }
        rawStream.push_back(fzip_detail::bitsToDouble(value));
        ++decodedCount;
    }

    if (!hasExactCount) {
        if (currCodeBit != codeStreamMaxBits) {
            std::cerr << "Error: code stream mismatch: " << currCodeBit << "/"
                      << codeStreamMaxBits << "\n";
            return false;
        }
        if (currArgBit != argumentStreamMaxBits) {
            std::cerr << "Error: argument stream mismatch: " << currArgBit << "/"
                      << argumentStreamMaxBits << "\n";
            return false;
        }
    }
    if (outCodeBitsUsed) *outCodeBitsUsed = currCodeBit;
    if (outArgBitsUsed) *outArgBitsUsed = currArgBit;
    return true;
}

inline bool fzipSerializeV2(const FzipEncoded& encoded, std::vector<std::uint8_t>& outBytes) {
    if (encoded.symbolTable.size() > std::numeric_limits<std::uint16_t>::max() ||
        encoded.commons.size() > std::numeric_limits<std::uint16_t>::max()) {
        std::cerr << "Error: table count exceeds v2 header limit\n";
        return false;
    }

    std::vector<std::uint8_t> payload;
    // Symbol table entries
    for (std::size_t i = 0; i < encoded.symbolTable.size(); ++i) {
        const FzipCode& sym = encoded.symbolTable[i];
        std::uint8_t flagsLen = static_cast<std::uint8_t>(
            (sym.codeLength & 0x0FU) | (sym.isCommon ? 0x80U : 0x00U));
        payload.push_back(flagsLen);
        payload.push_back(sym.prefixLength);
        if (!sym.isCommon && sym.prefixLength > 0) {
            int prefixBytes = (static_cast<int>(sym.prefixLength) + 7) / 8;
            for (int b = 0; b < prefixBytes; ++b) {
                payload.push_back(static_cast<std::uint8_t>((sym.prefix >> (56 - 8 * b)) & 0xFFU));
            }
        }
    }

    // Commons table
    std::uint8_t tmp8[8];
    for (std::size_t i = 0; i < encoded.commons.size(); ++i) {
        fzip_detail::putU64LE(tmp8, encoded.commons[i]);
        payload.insert(payload.end(), tmp8, tmp8 + 8);
    }

    // Packed bitstreams
    std::vector<std::uint8_t> codeBytes;
    std::vector<std::uint8_t> argBytes;
    fzip_detail::packWordsToBytes(encoded.codeStream, encoded.codeStreamBitLength, codeBytes);
    fzip_detail::packWordsToBytes(encoded.argumentStream, encoded.argumentStreamBitLength, argBytes);
    payload.insert(payload.end(), codeBytes.begin(), codeBytes.end());
    payload.insert(payload.end(), argBytes.begin(), argBytes.end());

    std::uint32_t payloadCrc =
        payload.empty() ? 0U : fzip_detail::crc32Update(0U, &payload[0], payload.size());

    outBytes.resize(fzip_detail::kHeaderSize + payload.size());
    std::uint8_t* hdr = &outBytes[0];
    std::memcpy(hdr, fzip_detail::kMagic, 4);
    fzip_detail::putU16LE(hdr + 4, fzip_detail::kFormatVersion);
    hdr[6] = encoded.maxCodeLength;
    hdr[7] = encoded.commonIndexBits;
    fzip_detail::putU64LE(hdr + 8, encoded.elementCount);
    fzip_detail::putU16LE(hdr + 16, static_cast<std::uint16_t>(encoded.symbolTable.size()));
    fzip_detail::putU16LE(hdr + 18, static_cast<std::uint16_t>(encoded.commons.size()));
    fzip_detail::putU64LE(hdr + 20, static_cast<ull>(codeBytes.size()));
    fzip_detail::putU64LE(hdr + 28, static_cast<ull>(argBytes.size()));
    fzip_detail::putU32LE(hdr + 36, payloadCrc);
    std::uint32_t headerCrc = fzip_detail::crc32Update(0U, hdr, 40);
    fzip_detail::putU32LE(hdr + 40, headerCrc);

    if (!payload.empty()) {
        std::memcpy(&outBytes[fzip_detail::kHeaderSize], &payload[0], payload.size());
    }
    return true;
}

inline bool fzipCompressBytes(const std::vector<double>& input,
                              std::vector<std::uint8_t>& outBytes,
                              int maxCodeLength = 9,
                              int targetCodeLength = 8,
                              int maxCommonBits = 13) {
    FzipEncoded encoded;
    if (!fzipEncode(input, encoded, false, maxCodeLength, targetCodeLength, maxCommonBits)) {
        return false;
    }
    return fzipSerializeV2(encoded, outBytes);
}

inline bool fzipCompress(const FzipOptions& options) {
    if (!options.valid || options.inputFile == NULL || options.outputFile == NULL) {
        return false;
    }
    std::vector<double> rawStream;
    if (!readFloats(options.inputFile, rawStream)) {
        return false;
    }
    std::vector<std::uint8_t> outBytes;
    if (!fzipCompressBytes(rawStream, outBytes, options.maxCodeLength, options.targetCodeLength,
                           options.log2CommonCodes)) {
        return false;
    }
    if (!fzip_detail::writeAll(outBytes.empty() ? NULL : &outBytes[0], 1, outBytes.size(),
                               options.outputFile) ||
        std::fflush(options.outputFile) != 0) {
        std::cerr << "Error: failed writing compressed output\n";
        return false;
    }
    return true;
}

inline bool fzipDeserializeV2(fzip_detail::ByteSource& src, std::vector<double>& rawStream) {
    rawStream.clear();

    std::uint8_t hdr[fzip_detail::kHeaderSize];
    if (!src.readExact(hdr, sizeof(hdr))) {
        std::cerr << "Error: failed reading fzip v2 header\n";
        return false;
    }
    if (std::memcmp(hdr, fzip_detail::kMagic, 4) != 0) {
        std::cerr << "Error: invalid fzip magic bytes\n";
        return false;
    }
    std::uint16_t version = fzip_detail::getU16LE(hdr + 4);
    if (version != fzip_detail::kFormatVersion) {
        std::cerr << "Error: unsupported fzip version (" << version << ")\n";
        return false;
    }
    std::uint32_t expectedHeaderCrc = fzip_detail::getU32LE(hdr + 40);
    if (fzip_detail::crc32Update(0U, hdr, 40) != expectedHeaderCrc) {
        std::cerr << "Error: fzip header checksum mismatch\n";
        return false;
    }

    int maxCodeLength = static_cast<int>(hdr[6]);
    int commonIndexBits = static_cast<int>(hdr[7]);
    ull elementCount = fzip_detail::getU64LE(hdr + 8);
    std::uint16_t symbolCount = fzip_detail::getU16LE(hdr + 16);
    std::uint16_t commonCount = fzip_detail::getU16LE(hdr + 18);
    ull codeStreamBytes = fzip_detail::getU64LE(hdr + 20);
    ull argumentStreamBytes = fzip_detail::getU64LE(hdr + 28);
    std::uint32_t expectedPayloadCrc = fzip_detail::getU32LE(hdr + 36);

    if (maxCodeLength < 1 || maxCodeLength > fzip_detail::kMaxAllowedCodeLength ||
        commonIndexBits < 0 || commonIndexBits > fzip_detail::kMaxCommonBits) {
        std::cerr << "Error: invalid fzip header parameters\n";
        return false;
    }
    if ((commonCount == 0 && commonIndexBits != 0) ||
        (commonCount > 0 && commonIndexBits != fzip_detail::ceilLog2(commonCount))) {
        std::cerr << "Error: invalid common table dimensions\n";
        return false;
    }
    if (static_cast<std::size_t>(symbolCount) > (std::size_t{1} << maxCodeLength)) {
        std::cerr << "Error: symbolCount exceeds 2^maxCodeLength\n";
        return false;
    }
    if (elementCount == 0) {
        if (symbolCount != 0 || commonCount != 0 ||
            codeStreamBytes != 0 || argumentStreamBytes != 0) {
            std::cerr << "Error: non-empty payload metadata for 0 elements\n";
            return false;
        }
        if (expectedPayloadCrc != 0U || !src.atEndIfMemory()) {
            std::cerr << "Error: invalid payload for empty stream\n";
            return false;
        }
        return true;
    }
    if (symbolCount == 0) {
        std::cerr << "Error: missing symbol table for non-empty stream\n";
        return false;
    }
    ull minCodeBytes = (elementCount / 8ULL) + ((elementCount % 8ULL) != 0ULL ? 1ULL : 0ULL);
    ull minElementsForArgs =
        (argumentStreamBytes / 8ULL) + ((argumentStreamBytes % 8ULL) != 0ULL ? 1ULL : 0ULL);
    if (codeStreamBytes < minCodeBytes || minElementsForArgs > elementCount) {
        std::cerr << "Error: inconsistent elementCount and bitstream sizes\n";
        return false;
    }
    if (codeStreamBytes > std::numeric_limits<std::size_t>::max() ||
        argumentStreamBytes > std::numeric_limits<std::size_t>::max() ||
        elementCount > std::numeric_limits<std::size_t>::max() / sizeof(double)) {
        std::cerr << "Error: stream size exceeds platform limits\n";
        return false;
    }

    src.resetCrc();
    std::vector<FzipCode> symbols(symbolCount);
    ull kraftSum = 0;
    const ull kraftTarget = 1ULL << maxCodeLength;
    int prevCodeLen = 0;
    for (std::size_t i = 0; i < symbolCount; ++i) {
        std::uint8_t meta[2];
        if (!src.readExact(meta, 2)) {
            std::cerr << "Error: truncated symbol table\n";
            return false;
        }
        if ((meta[0] & 0x70U) != 0U) {
            std::cerr << "Error: reserved symbol flags set\n";
            return false;
        }
        std::uint8_t codeLen = static_cast<std::uint8_t>(meta[0] & 0x0FU);
        bool isCommon = (meta[0] & 0x80U) != 0U;
        std::uint8_t prefixLen = meta[1];
        if (codeLen < 1 || codeLen > maxCodeLength || codeLen < prevCodeLen || prefixLen > 64) {
            std::cerr << "Error: invalid symbol metadata\n";
            return false;
        }
        if (isCommon && (commonCount == 0 ||
                         prefixLen != static_cast<std::uint8_t>(64 - commonIndexBits))) {
            std::cerr << "Error: invalid common symbol prefix length\n";
            return false;
        }
        prevCodeLen = codeLen;
        kraftSum += (1ULL << (maxCodeLength - codeLen));
        if (kraftSum > kraftTarget) {
            std::cerr << "Error: oversubscribed Huffman code lengths\n";
            return false;
        }

        ull prefix = 0;
        if (!isCommon && prefixLen > 0) {
            int prefixBytes = (static_cast<int>(prefixLen) + 7) / 8;
            std::uint8_t prefixBuf[8];
            if (!src.readExact(prefixBuf, static_cast<std::size_t>(prefixBytes))) {
                std::cerr << "Error: truncated symbol prefix\n";
                return false;
            }
            for (int b = 0; b < prefixBytes; ++b) {
                prefix |= static_cast<ull>(prefixBuf[b]) << (56 - 8 * b);
            }
            if (prefixLen < 64) {
                ull lowMask = (1ULL << (64 - prefixLen)) - 1ULL;
                if ((prefix & lowMask) != 0ULL) {
                    std::cerr << "Error: non-zero padding bits in symbol prefix\n";
                    return false;
                }
            }
        }
        symbols[i].codeLength = codeLen;
        symbols[i].isCommon = isCommon ? 1U : 0U;
        symbols[i].prefixLength = prefixLen;
        symbols[i].prefix = prefix;
    }

    if (!assignCanonicalCodes(symbols, maxCodeLength)) {
        std::cerr << "Error: invalid canonical Huffman symbol table\n";
        return false;
    }
    std::vector<FzipCode> lut;
    if (!buildDecodeLut(symbols, maxCodeLength, lut)) {
        std::cerr << "Error: failed to construct decode lookup table\n";
        return false;
    }

    std::vector<ull> commons(commonCount, 0ULL);
    if (commonCount > 0) {
        std::vector<std::uint8_t> commonBytes;
        if (!src.readExactVec(commonBytes, static_cast<std::size_t>(commonCount) * 8)) {
            std::cerr << "Error: truncated common values table\n";
            return false;
        }
        for (std::size_t i = 0; i < commonCount; ++i) {
            commons[i] = fzip_detail::getU64LE(&commonBytes[i * 8]);
        }
    }

    std::vector<std::uint8_t> codeBytes;
    std::vector<std::uint8_t> argBytes;
    if (!src.readExactVec(codeBytes, static_cast<std::size_t>(codeStreamBytes)) ||
        !src.readExactVec(argBytes, static_cast<std::size_t>(argumentStreamBytes))) {
        std::cerr << "Error: truncated bitstream payload\n";
        return false;
    }
    if (!src.atEndIfMemory()) {
        std::cerr << "Error: trailing bytes after fzip payload\n";
        return false;
    }
    if (src.crc != expectedPayloadCrc) {
        std::cerr << "Error: payload checksum mismatch\n";
        return false;
    }

    std::vector<ull> codeWords;
    std::vector<ull> argWords;
    fzip_detail::unpackBytesToWords(codeBytes.empty() ? NULL : &codeBytes[0],
                                    codeBytes.size(), codeWords);
    fzip_detail::unpackBytesToWords(argBytes.empty() ? NULL : &argBytes[0],
                                    argBytes.size(), argWords);

    ull usedCodeBits = 0;
    ull usedArgBits = 0;
    if (!fzipDecode(lut, maxCodeLength, commons, codeWords, argWords,
                    codeStreamBytes * 8ULL, argumentStreamBytes * 8ULL,
                    true, elementCount, rawStream, &usedCodeBits, &usedArgBits)) {
        return false;
    }

    if (!fzip_detail::checkTrailingZeroBits(codeBytes, usedCodeBits) ||
        !fzip_detail::checkTrailingZeroBits(argBytes, usedArgBits)) {
        std::cerr << "Error: invalid trailing bitstream padding or unconsumed bytes\n";
        return false;
    }
    return true;
}

inline bool fzipDecompressFromFile(FILE* inFile, std::vector<double>& rawStream) {
    rawStream.clear();
    if (inFile == NULL) {
        return false;
    }
    fzip_detail::ByteSource src(inFile);
    return fzipDeserializeV2(src, rawStream);
}

inline bool fzipDecompressBytes(const std::uint8_t* data,
                                std::size_t size,
                                std::vector<double>& rawStream) {
    rawStream.clear();
    if (data == NULL && size > 0) {
        return false;
    }
    fzip_detail::ByteSource src(data, size);
    return fzipDeserializeV2(src, rawStream);
}

inline bool fzipDecompress(const FzipOptions& options) {
    if (!options.valid || options.inputFile == NULL || options.outputFile == NULL) {
        return false;
    }
    std::vector<double> rawStream;
    if (!fzipDecompressFromFile(options.inputFile, rawStream)) {
        return false;
    }
    std::uint8_t buf[8];
    for (std::size_t i = 0; i < rawStream.size(); ++i) {
        fzip_detail::putU64LE(buf, fzip_detail::doubleToBits(rawStream[i]));
        if (std::fwrite(buf, 1, sizeof(buf), options.outputFile) != sizeof(buf)) {
            std::cerr << "Error: failed writing decompressed output\n";
            return false;
        }
    }
    if (std::fflush(options.outputFile) != 0) {
        std::cerr << "Error: failed flushing decompressed output\n";
        return false;
    }
    return true;
}

inline bool fzipCompress(std::vector<double>& rawStream,
                         std::vector<ull>& commons,
                         std::vector<FzipCode>& codes,
                         std::vector<ull>& codeStream,
                         std::vector<ull>& argumentStream,
                         ull& codeStreamBitLength,
                         ull& argumentStreamBitLength,
                         bool staticCommon = false,
                         int maxCodeLength = 9,
                         int targetCodeLength = 8,
                         int commonBits = 13) {
    FzipEncoded encoded;
    if (staticCommon) {
        encoded.commons = commons;
    }
    if (!fzipEncode(rawStream, encoded, staticCommon, maxCodeLength, targetCodeLength, commonBits)) {
        return false;
    }
    if (!buildDecodeLut(encoded.symbolTable, maxCodeLength, codes)) {
        return false;
    }
    commons.swap(encoded.commons);
    codeStream.swap(encoded.codeStream);
    argumentStream.swap(encoded.argumentStream);
    codeStreamBitLength = encoded.codeStreamBitLength;
    argumentStreamBitLength = encoded.argumentStreamBitLength;
    return true;
}

inline bool fzipDecompress(std::vector<double>& rawStream,
                           std::vector<ull>& commons,
                           std::vector<FzipCode>& codes,
                           std::vector<ull>& codeStream,
                           std::vector<ull>& argumentStream,
                           ull& codeStreamBitLength,
                           ull& argumentStreamBitLength) {
    if (codes.empty() && codeStreamBitLength == 0 && argumentStreamBitLength == 0) {
        rawStream.clear();
        return true;
    }
    std::size_t sz = codes.size();
    if (sz < 2 || (sz & (sz - 1)) != 0) {
        std::cerr << "Error: lookup table size must be a power of 2\n";
        return false;
    }
    int maxCodeLength = fzip_detail::ceilLog2(sz);
    return fzipDecode(codes, maxCodeLength, commons, codeStream, argumentStream,
                      codeStreamBitLength, argumentStreamBitLength,
                      false, 0ULL, rawStream);
}

inline std::vector<FzipCode> huffmanCoding(
    const std::vector<std::pair<ull, FzipCode> >& codeFrequencies,
    int maxCodeLength) {
    std::vector<FzipCode> ret;
    std::size_t n = codeFrequencies.size();
    if (n == 0 || maxCodeLength < 1 || maxCodeLength > fzip_detail::kMaxAllowedCodeLength ||
        n > (std::size_t{1} << maxCodeLength)) {
        return ret;
    }
    ret.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        ret.push_back(codeFrequencies[i].second);
    }
    if (n == 1) {
        ret[0].code = 0;
        ret[0].codeLength = 1;
        return ret;
    }

    struct PMNode {
        ull weight;
        int left;
        int right;
        int sym;
    };

    std::vector<std::size_t> sortedSyms(n);
    for (std::size_t i = 0; i < n; ++i) {
        sortedSyms[i] = i;
    }
    std::stable_sort(sortedSyms.begin(), sortedSyms.end(),
                     [&](std::size_t a, std::size_t b) {
                         return codeFrequencies[a].first < codeFrequencies[b].first;
                     });

    std::vector<PMNode> arena;
    arena.reserve(n * static_cast<std::size_t>(maxCodeLength) * 2);
    std::vector<int> baseCoins;
    baseCoins.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t sym = sortedSyms[i];
        ull w = std::max<ull>(1ULL, codeFrequencies[sym].first);
        PMNode node = {w, -1, -1, static_cast<int>(sym)};
        baseCoins.push_back(static_cast<int>(arena.size()));
        arena.push_back(node);
    }

    std::vector<int> curr = baseCoins;
    for (int level = 1; level < maxCodeLength; ++level) {
        std::vector<int> pkgs;
        pkgs.reserve(curr.size() / 2);
        for (std::size_t i = 0; i + 1 < curr.size(); i += 2) {
            PMNode pkg = {
                arena[static_cast<std::size_t>(curr[i])].weight +
                    arena[static_cast<std::size_t>(curr[i + 1])].weight,
                curr[i],
                curr[i + 1],
                -1
            };
            pkgs.push_back(static_cast<int>(arena.size()));
            arena.push_back(pkg);
        }
        std::vector<int> nextCurr;
        nextCurr.reserve(baseCoins.size() + pkgs.size());
        std::merge(baseCoins.begin(), baseCoins.end(),
                   pkgs.begin(), pkgs.end(),
                   std::back_inserter(nextCurr),
                   [&](int a, int b) {
                       return arena[static_cast<std::size_t>(a)].weight <
                              arena[static_cast<std::size_t>(b)].weight;
                   });
        curr.swap(nextCurr);
    }

    std::size_t takeCount = 2 * n - 2;
    if (curr.size() < takeCount) {
        return std::vector<FzipCode>();
    }

    std::vector<int> lengths(n, 0);
    std::vector<int> stack;
    stack.reserve(takeCount * 2);
    for (std::size_t i = 0; i < takeCount; ++i) {
        stack.push_back(curr[i]);
    }
    while (!stack.empty()) {
        int nodeIdx = stack.back();
        stack.pop_back();
        const PMNode& node = arena[static_cast<std::size_t>(nodeIdx)];
        if (node.sym >= 0) {
            lengths[static_cast<std::size_t>(node.sym)]++;
        } else {
            stack.push_back(node.left);
            stack.push_back(node.right);
        }
    }

    std::vector<std::size_t> canonicalOrder(n);
    for (std::size_t i = 0; i < n; ++i) {
        canonicalOrder[i] = i;
    }
    std::stable_sort(canonicalOrder.begin(), canonicalOrder.end(),
                     [&](std::size_t a, std::size_t b) {
                         return lengths[a] < lengths[b];
                     });

    ull nextCode = 0;
    int prevLen = 0;
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t idx = canonicalOrder[i];
        int len = lengths[idx];
        if (len < 1 || len > maxCodeLength) {
            return std::vector<FzipCode>();
        }
        nextCode <<= (len - prevLen);
        ret[idx].code = fzip_detail::bitReverseLow(nextCode, len);
        ret[idx].codeLength = static_cast<std::uint8_t>(len);
        ++nextCode;
        prevLen = len;
    }
    return ret;
}

struct PrefixNode {
    std::size_t leftChild;
    std::size_t rightChild;
    ull frequency;
    PrefixNode()
        : leftChild(std::numeric_limits<std::size_t>::max()),
          rightChild(std::numeric_limits<std::size_t>::max()),
          frequency(0) {}
};

inline void dfsPrefixTree(const std::vector<PrefixNode>& tree,
                          std::size_t curr,
                          const std::unordered_set<std::size_t>& partition,
                          std::vector<FzipCode>& codes,
                          int depth,
                          ull prefix) {
    const std::size_t kNull = std::numeric_limits<std::size_t>::max();
    if (curr == kNull) {
        return;
    }
    if (partition.count(curr)) {
        FzipCode code;
        code.prefixLength = static_cast<std::uint8_t>(depth);
        code.prefix = prefix;
        codes.push_back(code);
        return;
    }
    dfsPrefixTree(tree, tree[curr].leftChild, partition, codes, depth + 1, prefix);
    dfsPrefixTree(tree, tree[curr].rightChild, partition, codes, depth + 1,
                  prefix | (1ULL << (63 - depth)));
}

inline std::vector<FzipCode> createFzipCodes(const std::vector<std::pair<ull, ull> >& values,
                                             std::size_t targetCodeCount) {
    std::vector<FzipCode> ret;
    if (values.empty() || targetCodeCount == 0) {
        return ret;
    }

    const std::size_t kNull = std::numeric_limits<std::size_t>::max();
    std::vector<PrefixNode> tree(1);

    for (std::size_t i = 0; i < values.size(); ++i) {
        ull value = values[i].second;
        ull frequency = values[i].first;
        std::size_t currNode = 0;
        tree[currNode].frequency += frequency;
        for (int j = 0; j < 64; ++j) {
            if (value & (1ULL << (63 - j))) {
                if (tree[currNode].rightChild == kNull) {
                    std::size_t nextIdx = tree.size();
                    tree[currNode].rightChild = nextIdx;
                    tree.push_back(PrefixNode());
                }
                currNode = tree[currNode].rightChild;
            } else {
                if (tree[currNode].leftChild == kNull) {
                    std::size_t nextIdx = tree.size();
                    tree[currNode].leftChild = nextIdx;
                    tree.push_back(PrefixNode());
                }
                currNode = tree[currNode].leftChild;
            }
            tree[currNode].frequency += frequency;
        }
    }

    auto collapseChain = [&](std::size_t nodeIdx) {
        while (true) {
            std::size_t l = tree[nodeIdx].leftChild;
            std::size_t r = tree[nodeIdx].rightChild;
            if (l != kNull && r == kNull) {
                nodeIdx = l;
            } else if (l == kNull && r != kNull) {
                nodeIdx = r;
            } else {
                break;
            }
        }
        return nodeIdx;
    };

    struct HeapEntry {
        ull frequency;
        std::size_t nodeIdx;
        bool operator<(const HeapEntry& other) const {
            return frequency < other.frequency;
        }
    };

    std::priority_queue<HeapEntry> pq;
    std::vector<std::size_t> veryCommonValues;
    std::size_t root = collapseChain(0);
    pq.push(HeapEntry{tree[root].frequency, root});

    while (!pq.empty() && (pq.size() + veryCommonValues.size()) < targetCodeCount) {
        HeapEntry top = pq.top();
        pq.pop();
        std::size_t u = top.nodeIdx;
        std::size_t l = tree[u].leftChild;
        std::size_t r = tree[u].rightChild;
        if (l == kNull && r == kNull) {
            veryCommonValues.push_back(u);
        } else {
            std::size_t cl = collapseChain(l);
            std::size_t cr = collapseChain(r);
            pq.push(HeapEntry{tree[cl].frequency, cl});
            pq.push(HeapEntry{tree[cr].frequency, cr});
        }
    }

    std::unordered_set<std::size_t> partitionSet;
    while (!pq.empty()) {
        partitionSet.insert(pq.top().nodeIdx);
        pq.pop();
    }
    for (std::size_t i = 0; i < veryCommonValues.size(); ++i) {
        partitionSet.insert(veryCommonValues[i]);
    }

    dfsPrefixTree(tree, 0, partitionSet, ret, 0, 0ULL);
    return ret;
}

#endif
