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
#include <cmath>
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

inline bool writeAll(const void* data, std::size_t elemSize, std::size_t count, FILE* out) {
    if (count == 0) {
        return true;
    }
    return std::fwrite(data, elemSize, count, out) == count;
}

template <typename T>
inline bool readExactVec(FILE* in, std::vector<T>& out, std::size_t count) {
    out.clear();
    if (count == 0) {
        return true;
    }
    const std::size_t kChunk = 65536;
    std::vector<T> chunk(std::min(count, kChunk));
    std::size_t remaining = count;
    while (remaining > 0) {
        std::size_t want = std::min(remaining, kChunk);
        std::size_t got = std::fread(&chunk[0], sizeof(T), want, in);
        if (got != want) {
            return false;
        }
        out.insert(out.end(), chunk.begin(), chunk.begin() + static_cast<std::ptrdiff_t>(want));
        remaining -= want;
    }
    return true;
}

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
                        maxCodeLength < 1 || maxCodeLength > 9) {
                        std::cerr << "Error: --MaxCodeLength must be an integer in [1, 9]\n";
                        valid = false;
                        return;
                    }
                } else if (fzip_detail::startsWith(arg, "--TargetCodeLength=")) {
                    std::string valStr = arg.substr(std::strlen("--TargetCodeLength="));
                    if (!fzip_detail::parseIntStrict(valStr, targetCodeLength) ||
                        targetCodeLength < 0 || targetCodeLength > 9) {
                        std::cerr << "Error: --TargetCodeLength must be an integer in [0, 9]\n";
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

struct FzipHeader {
    ull r0;
    ull codeStreamBitLength;
    ull argumentStreamBitLength;
    ull r1[13];
    ull codesPointer;
    ull commonPointer;
    ull codeStreamPointer;
    ull argumentStreamPointer;
    ull size;
    ull r2[11];

    FzipHeader()
        : r0(0),
          codeStreamBitLength(0),
          argumentStreamBitLength(0),
          r1(),
          codesPointer(0),
          commonPointer(0),
          codeStreamPointer(0),
          argumentStreamPointer(0),
          size(0),
          r2() {}
};

static_assert(sizeof(FzipHeader) == 256, "FzipHeader must remain 256 bytes");

struct FzipCode {
    ull codeLength : 4, isCommon : 1, prefixLength : 7, code : 9;
    ull prefix;
    FzipCode()
        : codeLength(0), isCommon(0), prefixLength(0), code(0), prefix(0) {}
    void print() const {
        std::cerr << "codeLength: " << codeLength << "\n"
                  << "isCommon: " << isCommon << "\n"
                  << "prefixLength: " << prefixLength << "\n"
                  << "code: " << code << "\n"
                  << "prefix: " << std::hex << prefix << std::dec << "\n";
    }
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
    unsigned char buf[8];
    while (true) {
        std::size_t n = std::fread(buf, 1, sizeof(buf), inputFile);
        if (n == sizeof(buf)) {
            double val = 0.0;
            std::memcpy(&val, buf, sizeof(val));
            out.push_back(val);
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
                         int commonBits = 13);

inline bool fzipDecompress(std::vector<double>& rawStream,
                           std::vector<ull>& commons,
                           std::vector<FzipCode>& codes,
                           std::vector<ull>& codeStream,
                           std::vector<ull>& argumentStream,
                           ull& codeStreamBitLength,
                           ull& argumentStreamBitLength);

inline bool fzipCompress(const FzipOptions& options) {
    if (!options.valid || options.inputFile == NULL || options.outputFile == NULL) {
        return false;
    }
    std::vector<double> rawStream;
    if (!readFloats(options.inputFile, rawStream)) {
        return false;
    }
    std::vector<ull> commons;
    std::vector<FzipCode> codes;
    std::vector<ull> codeStream;
    std::vector<ull> argumentStream;
    ull codeStreamBitLength = 0;
    ull argumentStreamBitLength = 0;
    if (!fzipCompress(rawStream, commons, codes, codeStream, argumentStream,
                      codeStreamBitLength, argumentStreamBitLength, false,
                      options.maxCodeLength, options.targetCodeLength,
                      options.log2CommonCodes)) {
        return false;
    }

    FzipHeader header;
    header.codeStreamBitLength = codeStreamBitLength;
    header.argumentStreamBitLength = argumentStreamBitLength;
    std::size_t tmp = sizeof(FzipHeader);
    header.codesPointer = tmp;
    tmp += codes.size() * sizeof(FzipCode);
    header.commonPointer = tmp;
    tmp += commons.size() * sizeof(ull);
    header.codeStreamPointer = tmp;
    tmp += codeStream.size() * sizeof(ull);
    header.argumentStreamPointer = tmp;
    tmp += argumentStream.size() * sizeof(ull);
    header.size = tmp;

    if (!fzip_detail::writeAll(&header, sizeof(FzipHeader), 1, options.outputFile) ||
        !fzip_detail::writeAll(codes.empty() ? NULL : &codes[0], sizeof(FzipCode), codes.size(),
                               options.outputFile) ||
        !fzip_detail::writeAll(commons.empty() ? NULL : &commons[0], sizeof(ull), commons.size(),
                               options.outputFile) ||
        !fzip_detail::writeAll(codeStream.empty() ? NULL : &codeStream[0], sizeof(ull),
                               codeStream.size(), options.outputFile) ||
        !fzip_detail::writeAll(argumentStream.empty() ? NULL : &argumentStream[0], sizeof(ull),
                               argumentStream.size(), options.outputFile) ||
        std::fflush(options.outputFile) != 0) {
        std::cerr << "Error: failed writing compressed output\n";
        return false;
    }
    return true;
}

inline bool fzipDecompress(const FzipOptions& options) {
    if (!options.valid || options.inputFile == NULL || options.outputFile == NULL) {
        return false;
    }
    FzipHeader header;
    if (std::fread(&header, sizeof(FzipHeader), 1, options.inputFile) != 1) {
        std::cerr << "Error: failed reading fzip header\n";
        return false;
    }
    if (header.codesPointer != sizeof(FzipHeader) ||
        header.commonPointer < header.codesPointer ||
        header.codeStreamPointer < header.commonPointer ||
        header.argumentStreamPointer < header.codeStreamPointer ||
        header.size < header.argumentStreamPointer) {
        std::cerr << "Error: invalid fzip header layout\n";
        return false;
    }
    ull codesBytes = header.commonPointer - header.codesPointer;
    ull commonsBytes = header.codeStreamPointer - header.commonPointer;
    ull codeStreamBytes = header.argumentStreamPointer - header.codeStreamPointer;
    ull argumentStreamBytes = header.size - header.argumentStreamPointer;
    if (codesBytes % sizeof(FzipCode) != 0 ||
        commonsBytes % sizeof(ull) != 0 ||
        codeStreamBytes % sizeof(ull) != 0 ||
        argumentStreamBytes % sizeof(ull) != 0) {
        std::cerr << "Error: misaligned fzip section sizes\n";
        return false;
    }
    ull codesCount = codesBytes / sizeof(FzipCode);
    ull commonsCount = commonsBytes / sizeof(ull);
    ull codeWords = codeStreamBytes / sizeof(ull);
    ull argumentWords = argumentStreamBytes / sizeof(ull);
    if (codesCount != 512ULL || commonsCount > 65536ULL ||
        codeWords > std::numeric_limits<std::size_t>::max() / sizeof(ull) ||
        argumentWords > std::numeric_limits<std::size_t>::max() / sizeof(ull)) {
        std::cerr << "Error: out-of-range fzip table sizes\n";
        return false;
    }
    if (codeWords > std::numeric_limits<ull>::max() / 64ULL ||
        argumentWords > std::numeric_limits<ull>::max() / 64ULL ||
        header.codeStreamBitLength > codeWords * 64ULL ||
        header.argumentStreamBitLength > argumentWords * 64ULL) {
        std::cerr << "Error: bitstream bit length exceeds section size\n";
        return false;
    }

    std::vector<FzipCode> codes;
    std::vector<ull> commons;
    std::vector<ull> codeStream;
    std::vector<ull> argumentStream;
    if (!fzip_detail::readExactVec(options.inputFile, codes, static_cast<std::size_t>(codesCount)) ||
        !fzip_detail::readExactVec(options.inputFile, commons,
                                   static_cast<std::size_t>(commonsCount)) ||
        !fzip_detail::readExactVec(options.inputFile, codeStream,
                                   static_cast<std::size_t>(codeWords)) ||
        !fzip_detail::readExactVec(options.inputFile, argumentStream,
                                   static_cast<std::size_t>(argumentWords))) {
        std::cerr << "Error: truncated fzip payload\n";
        return false;
    }

    std::vector<double> rawStream;
    if (!fzipDecompress(rawStream, commons, codes, codeStream, argumentStream,
                        header.codeStreamBitLength, header.argumentStreamBitLength)) {
        return false;
    }
    if (!fzip_detail::writeAll(rawStream.empty() ? NULL : &rawStream[0], sizeof(double),
                               rawStream.size(), options.outputFile) ||
        std::fflush(options.outputFile) != 0) {
        std::cerr << "Error: failed writing decompressed output\n";
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
                         bool staticCommon,
                         int maxCodeLength,
                         int targetCodeLength,
                         int commonBits) {
    if (maxCodeLength < 1 || maxCodeLength > 9 ||
        targetCodeLength < 0 || targetCodeLength > maxCodeLength ||
        commonBits < 0 || commonBits > 13) {
        std::cerr << "Error: invalid compression parameters\n";
        return false;
    }
    const std::size_t lutSize = std::size_t{1} << 9;
    const std::size_t commonCapacity = (commonBits > 0) ? (std::size_t{1} << commonBits) : 0;

    if (!staticCommon) {
        commons.assign(commonCapacity, 0ULL);
    }
    codes.clear();
    codeStream.clear();
    argumentStream.clear();
    codeStreamBitLength = 0;
    argumentStreamBitLength = 0;

    if (rawStream.empty()) {
        codes.assign(lutSize, FzipCode());
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

    std::size_t activeCommonCount = std::min(commons.size(), commonCapacity);
    if (!staticCommon) {
        activeCommonCount = std::min(sortedByFrequency.size(), commonCapacity);
        std::size_t idx = 0;
        for (std::vector<std::pair<ull, ull> >::const_reverse_iterator it =
                 sortedByFrequency.rbegin();
             idx < activeCommonCount && it != sortedByFrequency.rend(); ++idx, ++it) {
            commons[idx] = it->second;
        }
    }

    std::unordered_set<ull> commonSet;
    std::unordered_map<ull, std::size_t> mapToCommon;
    for (std::size_t i = 0; i < activeCommonCount; ++i) {
        commonSet.insert(commons[i]);
        mapToCommon[commons[i]] = i;
    }

    FzipCode commonCode;
    commonCode.isCommon = 1;
    commonCode.prefixLength = static_cast<ull>(64 - commonBits);

    std::map<ull, std::size_t> mapToOrigPrefix;
    std::vector<ull> origFrequencies(prefixCodes.size() + 1, 0ULL);
    for (std::size_t i = 0; i < prefixCodes.size(); ++i) {
        mapToOrigPrefix[prefixCodes[i].prefix] = i;
    }

    const ull commonThresholdPrefixLen = static_cast<ull>(64 - commonBits);
    for (std::map<ull, ull>::const_iterator it = frequencies.begin(); it != frequencies.end(); ++it) {
        std::map<ull, std::size_t>::const_iterator ub = mapToOrigPrefix.upper_bound(it->first);
        if (ub == mapToOrigPrefix.begin()) {
            std::cerr << "Error: value missing from prefix trie partition\n";
            return false;
        }
        std::size_t origIdx = std::prev(ub)->second;
        if (activeCommonCount > 0 &&
            prefixCodes[origIdx].prefixLength <= commonThresholdPrefixLen &&
            commonSet.count(it->first)) {
            origFrequencies[prefixCodes.size()] += it->second;
        } else {
            origFrequencies[origIdx] += it->second;
        }
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

    codes = huffmanCoding(codeFrequencies, maxCodeLength);
    if (codes.empty()) {
        std::cerr << "Error: failed to construct bounded-length Huffman codes\n";
        return false;
    }
    for (std::size_t i = 0; i < codes.size(); ++i) {
        if (codes[i].codeLength == 0 || codes[i].codeLength > static_cast<ull>(maxCodeLength)) {
            std::cerr << "Error: invalid Huffman code length produced\n";
            return false;
        }
    }

    std::vector<FzipCode> lut(lutSize, FzipCode());
    for (std::size_t i = 0; i < codes.size(); ++i) {
        std::size_t stride = std::size_t{1} << codes[i].codeLength;
        for (std::size_t j = static_cast<std::size_t>(codes[i].code); j < lutSize; j += stride) {
            lut[j] = codes[i];
        }
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
        if (activeCommonCount > 0 &&
            prefixCodes[origIdx].prefixLength <= commonThresholdPrefixLen &&
            commonSet.count(rawBits)) {
            activeIdx = activeCommonIdx;
        }
        if (activeIdx < 0 || static_cast<std::size_t>(activeIdx) >= codes.size()) {
            std::cerr << "Error: unmapped symbol during encoding\n";
            return false;
        }

        const FzipCode& sym = codes[static_cast<std::size_t>(activeIdx)];
        fzip_detail::appendBits(codeStream, codeBuffer, codeBufferEndBit,
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
        fzip_detail::appendBits(argumentStream, argumentBuffer, argumentBufferEndBit,
                                argument, argumentLength);
    }

    codeStreamBitLength = static_cast<ull>(codeStream.size()) * 64ULL +
                          static_cast<ull>(codeBufferEndBit);
    argumentStreamBitLength = static_cast<ull>(argumentStream.size()) * 64ULL +
                              static_cast<ull>(argumentBufferEndBit);

    if (codeBufferEndBit > 0) {
        codeStream.push_back(codeBuffer);
    }
    if (argumentBufferEndBit > 0) {
        argumentStream.push_back(argumentBuffer);
    }

    codes.swap(lut);
    return true;
}

inline std::vector<FzipCode> huffmanCoding(
    const std::vector<std::pair<ull, FzipCode> >& codeFrequencies,
    int maxCodeLength) {
    std::vector<FzipCode> ret;
    std::size_t n = codeFrequencies.size();
    if (n == 0 || maxCodeLength < 1 || maxCodeLength > 15 ||
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
        ull lsbCode = 0;
        for (int b = 0; b < len; ++b) {
            if ((nextCode >> (len - 1 - b)) & 1ULL) {
                lsbCode |= (1ULL << b);
            }
        }
        ret[idx].code = lsbCode;
        ret[idx].codeLength = static_cast<ull>(len);
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
        code.prefixLength = static_cast<ull>(depth);
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

inline bool fzipDecompress(std::vector<double>& rawStream,
                           std::vector<ull>& commons,
                           std::vector<FzipCode>& codes,
                           std::vector<ull>& codeStream,
                           std::vector<ull>& argumentStream,
                           ull& codeStreamBitLength,
                           ull& argumentStreamBitLength) {
    rawStream.clear();
    if (codeStreamBitLength == 0) {
        if (argumentStreamBitLength != 0) {
            std::cerr << "Error: argument stream non-empty for empty code stream\n";
            return false;
        }
        return true;
    }
    if (codes.size() != 512) {
        std::cerr << "Error: invalid Huffman lookup table size (" << codes.size() << ")\n";
        return false;
    }
    if (codeStreamBitLength > static_cast<ull>(codeStream.size()) * 64ULL ||
        argumentStreamBitLength > static_cast<ull>(argumentStream.size()) * 64ULL) {
        std::cerr << "Error: bitstream bit length exceeds buffer bounds\n";
        return false;
    }

    ull currCodeStreamBit = 0;
    ull currArgumentStreamBit = 0;
    while (currCodeStreamBit < codeStreamBitLength) {
        ull peek = fzip_detail::readBits(codeStream, currCodeStreamBit, 9);
        const FzipCode& code = codes[static_cast<std::size_t>(peek)];
        if (code.codeLength == 0 || code.codeLength > 9) {
            std::cerr << "Error: invalid Huffman code length (" << code.codeLength << ")\n";
            return false;
        }
        if (currCodeStreamBit + code.codeLength > codeStreamBitLength) {
            std::cerr << "Error: code stream truncated\n";
            return false;
        }
        currCodeStreamBit += code.codeLength;

        if (code.prefixLength > 64) {
            std::cerr << "Error: invalid prefix length (" << code.prefixLength << ")\n";
            return false;
        }
        int argBits = 64 - static_cast<int>(code.prefixLength);
        if (currArgumentStreamBit + static_cast<ull>(argBits) > argumentStreamBitLength) {
            std::cerr << "Error: argument stream truncated\n";
            return false;
        }
        ull arg = fzip_detail::readBits(argumentStream, currArgumentStreamBit, argBits);
        currArgumentStreamBit += static_cast<ull>(argBits);

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
    }

    if (currCodeStreamBit != codeStreamBitLength) {
        std::cerr << "Error: code stream mismatch: " << currCodeStreamBit << "/"
                  << codeStreamBitLength << "\n";
        return false;
    }
    if (currArgumentStreamBit != argumentStreamBitLength) {
        std::cerr << "Error: argument stream mismatch: " << currArgumentStreamBit << "/"
                  << argumentStreamBitLength << "\n";
        return false;
    }
    return true;
}

#endif
