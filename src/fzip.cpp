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
#include "fzip/fzip.hpp"

#include <algorithm>
#include <cstring>
#include <iterator>
#include <limits>
#include <map>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace fzip {

namespace {

static_assert(sizeof(double) == 8, "IEEE-754 64-bit double required");
static_assert(std::numeric_limits<double>::is_iec559, "IEEE-754 double required");

constexpr std::uint8_t kMagic[4] = {'F', 'Z', 'I', 'P'};

struct FzipCode {
    std::uint8_t codeLength = 0;
    std::uint8_t isCommon = 0;
    std::uint8_t prefixLength = 0;
    std::uint16_t code = 0;
    std::uint64_t prefix = 0;
};

struct FzipEncoded {
    std::uint8_t maxCodeLength = 9;
    std::uint8_t commonIndexBits = 0;
    std::uint64_t elementCount = 0;
    std::vector<FzipCode> symbolTable;
    std::vector<std::uint64_t> commons;
    std::vector<std::uint64_t> codeStream;
    std::vector<std::uint64_t> argumentStream;
    std::uint64_t codeStreamBitLength = 0;
    std::uint64_t argumentStreamBitLength = 0;
};

inline std::uint64_t doubleToBits(double value) {
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    return bits;
}

inline double bitsToDouble(std::uint64_t bits) {
    double value = 0.0;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
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
           static_cast<std::uint16_t>(static_cast<std::uint16_t>(src[1]) << 8);
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

inline std::uint16_t bitReverseLow(std::uint64_t code, int len) {
    std::uint16_t lsbCode = 0;
    for (int b = 0; b < len; ++b) {
        if ((code >> (len - 1 - b)) & 1ULL) {
            lsbCode |= static_cast<std::uint16_t>(1U << b);
        }
    }
    return lsbCode;
}

inline void appendBits(std::vector<std::uint64_t>& stream,
                       std::uint64_t& buffer,
                       int& endBit,
                       std::uint64_t bits,
                       int width) {
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

inline std::uint64_t readBits(const std::vector<std::uint64_t>& stream,
                              std::uint64_t bitPos,
                              int width) {
    if (width <= 0) {
        return 0ULL;
    }
    std::size_t wordIdx = static_cast<std::size_t>(bitPos / 64ULL);
    int bitOffset = static_cast<int>(bitPos % 64ULL);
    std::uint64_t val = (wordIdx < stream.size()) ? (stream[wordIdx] >> bitOffset) : 0ULL;
    if (bitOffset > 0 && bitOffset + width > 64 && (wordIdx + 1) < stream.size()) {
        val |= stream[wordIdx + 1] << (64 - bitOffset);
    }
    if (width < 64) {
        val &= (1ULL << width) - 1ULL;
    }
    return val;
}

inline void packWordsToBytes(const std::vector<std::uint64_t>& words,
                             std::uint64_t bitLength,
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
                               std::vector<std::uint64_t>& outWords) {
    std::size_t wordCount = (byteCount + 7) / 8;
    outWords.assign(wordCount, 0ULL);
    for (std::size_t i = 0; i < byteCount; ++i) {
        std::size_t wordIdx = i / 8;
        int byteShift = static_cast<int>((i % 8) * 8);
        outWords[wordIdx] |= static_cast<std::uint64_t>(bytes[i]) << byteShift;
    }
}

inline bool checkTrailingZeroBits(const std::vector<std::uint8_t>& bytes, std::uint64_t bitPos) {
    std::size_t byteLen = bytes.size();
    std::uint64_t expectedBytes = (bitPos + 7ULL) / 8ULL;
    if (expectedBytes != static_cast<std::uint64_t>(byteLen)) {
        return false;
    }
    int remBits = static_cast<int>(static_cast<std::uint64_t>(byteLen) * 8ULL - bitPos);
    if (remBits > 0 && byteLen > 0) {
        std::uint8_t lastByte = bytes[byteLen - 1];
        std::uint8_t usedMask = static_cast<std::uint8_t>((1U << (8 - remBits)) - 1U);
        if ((lastByte & static_cast<std::uint8_t>(~usedMask)) != 0U) {
            return false;
        }
    }
    return true;
}

class ByteSource {
public:
    explicit ByteSource(std::FILE* f) : file_(f) {}
    ByteSource(const std::uint8_t* data, std::size_t size) : mem_(data), memSize_(size) {}

    void resetCrc() { crc_ = 0U; }
    std::uint32_t crc() const { return crc_; }

    bool readExact(std::uint8_t* dst, std::size_t count) {
        if (count == 0) {
            return true;
        }
        if (file_ != nullptr) {
            if (std::fread(dst, 1, count, file_) != count) {
                return false;
            }
        } else {
            if (mem_ == nullptr || count > memSize_ - memPos_) {
                return false;
            }
            std::memcpy(dst, mem_ + memPos_, count);
            memPos_ += count;
        }
        crc_ = crc32Update(crc_, dst, count);
        return true;
    }

    bool readExactVec(std::vector<std::uint8_t>& out, std::size_t byteCount) {
        out.clear();
        if (byteCount == 0) {
            return true;
        }
        if (file_ == nullptr) {
            if (mem_ == nullptr || byteCount > memSize_ - memPos_) {
                return false;
            }
            out.assign(mem_ + memPos_, mem_ + memPos_ + byteCount);
            memPos_ += byteCount;
            crc_ = crc32Update(crc_, out.data(), byteCount);
            return true;
        }
        constexpr std::size_t kChunk = 65536;
        std::vector<std::uint8_t> chunk(std::min(byteCount, kChunk));
        std::size_t remaining = byteCount;
        while (remaining > 0) {
            std::size_t want = std::min(remaining, kChunk);
            if (std::fread(chunk.data(), 1, want, file_) != want) {
                return false;
            }
            crc_ = crc32Update(crc_, chunk.data(), want);
            out.insert(out.end(), chunk.begin(), chunk.begin() + static_cast<std::ptrdiff_t>(want));
            remaining -= want;
        }
        return true;
    }

    bool atEndIfMemory() const {
        return file_ != nullptr || memPos_ == memSize_;
    }

private:
    std::FILE* file_ = nullptr;
    const std::uint8_t* mem_ = nullptr;
    std::size_t memSize_ = 0;
    std::size_t memPos_ = 0;
    std::uint32_t crc_ = 0U;
};

struct PrefixNode {
    std::size_t leftChild = std::numeric_limits<std::size_t>::max();
    std::size_t rightChild = std::numeric_limits<std::size_t>::max();
    std::uint64_t frequency = 0;
};

void dfsPrefixTree(const std::vector<PrefixNode>& tree,
                   std::size_t curr,
                   const std::unordered_set<std::size_t>& partition,
                   std::vector<FzipCode>& codes,
                   int depth,
                   std::uint64_t prefix) {
    constexpr std::size_t kNull = std::numeric_limits<std::size_t>::max();
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

std::vector<FzipCode> createFzipCodes(
    const std::vector<std::pair<std::uint64_t, std::uint64_t>>& values,
    std::size_t targetCodeCount) {
    std::vector<FzipCode> ret;
    if (values.empty() || targetCodeCount == 0) {
        return ret;
    }

    constexpr std::size_t kNull = std::numeric_limits<std::size_t>::max();
    std::vector<PrefixNode> tree(1);

    for (const auto& [frequency, value] : values) {
        std::size_t currNode = 0;
        tree[currNode].frequency += frequency;
        for (int j = 0; j < 64; ++j) {
            if (value & (1ULL << (63 - j))) {
                if (tree[currNode].rightChild == kNull) {
                    std::size_t nextIdx = tree.size();
                    tree[currNode].rightChild = nextIdx;
                    tree.emplace_back();
                }
                currNode = tree[currNode].rightChild;
            } else {
                if (tree[currNode].leftChild == kNull) {
                    std::size_t nextIdx = tree.size();
                    tree[currNode].leftChild = nextIdx;
                    tree.emplace_back();
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
        std::uint64_t frequency;
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
    for (std::size_t idx : veryCommonValues) {
        partitionSet.insert(idx);
    }

    dfsPrefixTree(tree, 0, partitionSet, ret, 0, 0ULL);
    return ret;
}

bool assignCanonicalCodes(std::vector<FzipCode>& symbols, int maxCodeLength) {
    if (symbols.empty()) {
        return true;
    }
    if (maxCodeLength < 1 || maxCodeLength > kMaxAllowedCodeLength ||
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

    std::uint64_t kraftSum = 0;
    const std::uint64_t kraftTarget = 1ULL << maxCodeLength;
    std::uint64_t nextCode = 0;
    int prevLen = 0;
    for (auto& sym : symbols) {
        int len = static_cast<int>(sym.codeLength);
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
        sym.code = bitReverseLow(nextCode, len);
        ++nextCode;
        prevLen = len;
    }
    return kraftSum == kraftTarget;
}

bool buildDecodeLut(const std::vector<FzipCode>& symbols,
                    int maxCodeLength,
                    std::vector<FzipCode>& lut) {
    if (maxCodeLength < 1 || maxCodeLength > kMaxAllowedCodeLength) {
        return false;
    }
    const std::size_t lutSize = std::size_t{1} << maxCodeLength;
    lut.assign(lutSize, FzipCode{});
    for (const auto& sym : symbols) {
        if (sym.codeLength == 0 || sym.codeLength > maxCodeLength) {
            return false;
        }
        std::size_t stride = std::size_t{1} << sym.codeLength;
        for (std::size_t j = sym.code; j < lutSize; j += stride) {
            if (lut[j].codeLength != 0) {
                return false;
            }
            lut[j] = sym;
        }
    }
    return true;
}

std::vector<FzipCode> huffmanCoding(
    const std::vector<std::pair<std::uint64_t, FzipCode>>& codeFrequencies,
    int maxCodeLength) {
    std::vector<FzipCode> ret;
    std::size_t n = codeFrequencies.size();
    if (n == 0 || maxCodeLength < 1 || maxCodeLength > kMaxAllowedCodeLength ||
        n > (std::size_t{1} << maxCodeLength)) {
        return ret;
    }
    ret.reserve(n);
    for (const auto& item : codeFrequencies) {
        ret.push_back(item.second);
    }
    if (n == 1) {
        ret[0].code = 0;
        ret[0].codeLength = 1;
        return ret;
    }

    struct PMNode {
        std::uint64_t weight;
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
        std::uint64_t w = std::max<std::uint64_t>(1ULL, codeFrequencies[sym].first);
        baseCoins.push_back(static_cast<int>(arena.size()));
        arena.push_back(PMNode{w, -1, -1, static_cast<int>(sym)});
    }

    std::vector<int> curr = baseCoins;
    for (int level = 1; level < maxCodeLength; ++level) {
        std::vector<int> pkgs;
        pkgs.reserve(curr.size() / 2);
        for (std::size_t i = 0; i + 1 < curr.size(); i += 2) {
            std::uint64_t w = arena[static_cast<std::size_t>(curr[i])].weight +
                              arena[static_cast<std::size_t>(curr[i + 1])].weight;
            pkgs.push_back(static_cast<int>(arena.size()));
            arena.push_back(PMNode{w, curr[i], curr[i + 1], -1});
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
        return {};
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

    std::uint64_t nextCode = 0;
    int prevLen = 0;
    for (std::size_t idx : canonicalOrder) {
        int len = lengths[idx];
        if (len < 1 || len > maxCodeLength) {
            return {};
        }
        nextCode <<= (len - prevLen);
        ret[idx].code = bitReverseLow(nextCode, len);
        ret[idx].codeLength = static_cast<std::uint8_t>(len);
        ++nextCode;
        prevLen = len;
    }
    return ret;
}

bool fzipEncode(const double* rawData,
                std::size_t count,
                const CompressOptions& options,
                FzipEncoded& encoded,
                std::string& err) {
    if (count > 0 && rawData == nullptr) {
        err = "input data pointer is null for non-zero element count";
        return false;
    }
    if (options.maxCodeLength < 1 || options.maxCodeLength > kMaxAllowedCodeLength ||
        options.targetCodeLength < 0 || options.targetCodeLength > options.maxCodeLength ||
        options.maxCommonBits < 0 || options.maxCommonBits > kMaxCommonBits) {
        err = "invalid compression parameters";
        return false;
    }

    encoded = FzipEncoded{};
    encoded.maxCodeLength = static_cast<std::uint8_t>(options.maxCodeLength);
    encoded.elementCount = static_cast<std::uint64_t>(count);

    if (count == 0) {
        return true;
    }

    std::map<std::uint64_t, std::uint64_t> frequencies;
    for (std::size_t i = 0; i < count; ++i) {
        frequencies[doubleToBits(rawData[i])]++;
    }

    std::vector<std::pair<std::uint64_t, std::uint64_t>> sortedByFrequency;
    sortedByFrequency.reserve(frequencies.size());
    for (const auto& [val, freq] : frequencies) {
        sortedByFrequency.emplace_back(freq, val);
    }
    std::sort(sortedByFrequency.begin(), sortedByFrequency.end());

    std::size_t maxAllowedCodes = std::size_t{1} << options.maxCodeLength;
    std::size_t targetPrefixCodes = std::size_t{1} << options.targetCodeLength;
    if (targetPrefixCodes >= maxAllowedCodes) {
        targetPrefixCodes = maxAllowedCodes - 1;
    }

    std::vector<FzipCode> prefixCodes = createFzipCodes(sortedByFrequency, targetPrefixCodes);
    std::map<std::uint64_t, std::size_t> mapToOrigPrefix;
    for (std::size_t i = 0; i < prefixCodes.size(); ++i) {
        mapToOrigPrefix[prefixCodes[i].prefix] = i;
    }

    const std::size_t commonCapacity =
        (options.maxCommonBits > 0) ? (std::size_t{1} << options.maxCommonBits) : 0;

    std::vector<std::pair<std::uint64_t, std::uint64_t>> candidates;
    for (auto it = sortedByFrequency.rbegin();
         it != sortedByFrequency.rend() && candidates.size() < commonCapacity; ++it) {
        std::uint64_t freq = it->first;
        std::uint64_t val = it->second;
        if (freq < 2ULL) {
            break;
        }
        auto ub = mapToOrigPrefix.upper_bound(val);
        if (ub == mapToOrigPrefix.begin()) {
            err = "internal error: value missing from prefix trie partition";
            return false;
        }
        std::size_t origIdx = std::prev(ub)->second;
        int prefLen = static_cast<int>(prefixCodes[origIdx].prefixLength);
        if (prefLen < 64 && freq * static_cast<std::uint64_t>(64 - prefLen) > 64ULL) {
            candidates.emplace_back(freq, val);
        }
    }

    int tentativeBits = ceilLog2(candidates.size());
    for (const auto& [freq, val] : candidates) {
        std::size_t origIdx = std::prev(mapToOrigPrefix.upper_bound(val))->second;
        int prefLen = static_cast<int>(prefixCodes[origIdx].prefixLength);
        int argSavedPerOcc = (64 - prefLen) - tentativeBits;
        if (argSavedPerOcc > 0 && freq * static_cast<std::uint64_t>(argSavedPerOcc) > 64ULL) {
            encoded.commons.push_back(val);
        }
    }

    int commonBits = ceilLog2(encoded.commons.size());
    encoded.commonIndexBits = static_cast<std::uint8_t>(commonBits);

    std::unordered_set<std::uint64_t> commonSet;
    std::unordered_map<std::uint64_t, std::size_t> mapToCommon;
    for (std::size_t i = 0; i < encoded.commons.size(); ++i) {
        commonSet.insert(encoded.commons[i]);
        mapToCommon.emplace(encoded.commons[i], i);
    }

    FzipCode commonCode;
    commonCode.isCommon = 1;
    commonCode.prefixLength = static_cast<std::uint8_t>(64 - commonBits);

    const std::uint8_t commonThresholdPrefixLen = static_cast<std::uint8_t>(64 - commonBits);
    std::vector<std::uint64_t> origFrequencies(prefixCodes.size() + 1, 0ULL);

    for (const auto& [val, freq] : frequencies) {
        auto ub = mapToOrigPrefix.upper_bound(val);
        if (ub == mapToOrigPrefix.begin()) {
            err = "internal error: value missing from prefix trie partition";
            return false;
        }
        std::size_t origIdx = std::prev(ub)->second;
        if (!encoded.commons.empty() &&
            prefixCodes[origIdx].prefixLength <= commonThresholdPrefixLen &&
            commonSet.count(val)) {
            origFrequencies[prefixCodes.size()] += freq;
        } else {
            origFrequencies[origIdx] += freq;
        }
    }

    if (origFrequencies[prefixCodes.size()] == 0ULL && !encoded.commons.empty()) {
        encoded.commons.clear();
        commonSet.clear();
        mapToCommon.clear();
        commonBits = 0;
        encoded.commonIndexBits = 0;
    }

    std::vector<std::pair<std::uint64_t, FzipCode>> codeFrequencies;
    std::vector<int> origToActive(prefixCodes.size(), -1);
    int activeCommonIdx = -1;
    for (std::size_t i = 0; i < prefixCodes.size(); ++i) {
        if (origFrequencies[i] > 0) {
            origToActive[i] = static_cast<int>(codeFrequencies.size());
            codeFrequencies.emplace_back(origFrequencies[i], prefixCodes[i]);
        }
    }
    if (origFrequencies[prefixCodes.size()] > 0) {
        activeCommonIdx = static_cast<int>(codeFrequencies.size());
        codeFrequencies.emplace_back(origFrequencies[prefixCodes.size()], commonCode);
    }

    std::vector<FzipCode> rawSymbols = huffmanCoding(codeFrequencies, options.maxCodeLength);
    if (rawSymbols.empty()) {
        err = "failed to construct bounded-length Huffman codes";
        return false;
    }

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
    if (!assignCanonicalCodes(encoded.symbolTable, options.maxCodeLength)) {
        err = "canonical code assignment failed";
        return false;
    }

    int codeBufferEndBit = 0;
    int argumentBufferEndBit = 0;
    std::uint64_t codeBuffer = 0;
    std::uint64_t argumentBuffer = 0;

    for (std::size_t i = 0; i < count; ++i) {
        std::uint64_t rawBits = doubleToBits(rawData[i]);
        auto ub = mapToOrigPrefix.upper_bound(rawBits);
        if (ub == mapToOrigPrefix.begin()) {
            err = "internal error: value missing from prefix map during encoding";
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
            err = "internal error: unmapped symbol during encoding";
            return false;
        }

        const FzipCode& sym =
            encoded.symbolTable[activeToCanonical[static_cast<std::size_t>(activeIdx)]];
        appendBits(encoded.codeStream, codeBuffer, codeBufferEndBit,
                   static_cast<std::uint64_t>(sym.code),
                   static_cast<int>(sym.codeLength));

        std::uint64_t argument = 0;
        int argumentLength = 0;
        if (sym.isCommon) {
            auto cit = mapToCommon.find(rawBits);
            if (cit == mapToCommon.end()) {
                err = "internal error: common value missing from index map";
                return false;
            }
            argument = static_cast<std::uint64_t>(cit->second);
            argumentLength = commonBits;
        } else if (sym.prefixLength < 64) {
            argumentLength = 64 - static_cast<int>(sym.prefixLength);
            argument = (argumentLength < 64) ? (rawBits & ((1ULL << argumentLength) - 1ULL))
                                             : rawBits;
        }
        appendBits(encoded.argumentStream, argumentBuffer, argumentBufferEndBit,
                   argument, argumentLength);
    }

    encoded.codeStreamBitLength =
        static_cast<std::uint64_t>(encoded.codeStream.size()) * 64ULL +
        static_cast<std::uint64_t>(codeBufferEndBit);
    encoded.argumentStreamBitLength =
        static_cast<std::uint64_t>(encoded.argumentStream.size()) * 64ULL +
        static_cast<std::uint64_t>(argumentBufferEndBit);

    if (codeBufferEndBit > 0) {
        encoded.codeStream.push_back(codeBuffer);
    }
    if (argumentBufferEndBit > 0) {
        encoded.argumentStream.push_back(argumentBuffer);
    }
    return true;
}

bool fzipDecode(const std::vector<FzipCode>& lut,
                int maxCodeLength,
                const std::vector<std::uint64_t>& commons,
                const std::vector<std::uint64_t>& codeStream,
                const std::vector<std::uint64_t>& argumentStream,
                std::uint64_t codeStreamMaxBits,
                std::uint64_t argumentStreamMaxBits,
                std::uint64_t expectedElements,
                std::vector<double>& rawStream,
                std::uint64_t& outCodeBitsUsed,
                std::uint64_t& outArgBitsUsed,
                std::string& err) {
    rawStream.clear();
    outCodeBitsUsed = 0;
    outArgBitsUsed = 0;

    if (maxCodeLength < 1 || maxCodeLength > kMaxAllowedCodeLength ||
        lut.size() != (std::size_t{1} << maxCodeLength)) {
        err = "invalid decode lookup table";
        return false;
    }
    if (expectedElements == 0) {
        return true;
    }
    if (codeStreamMaxBits > static_cast<std::uint64_t>(codeStream.size()) * 64ULL ||
        argumentStreamMaxBits > static_cast<std::uint64_t>(argumentStream.size()) * 64ULL) {
        err = "bitstream bit length exceeds buffer bounds";
        return false;
    }

    rawStream.reserve(static_cast<std::size_t>(expectedElements));

    std::uint64_t currCodeBit = 0;
    std::uint64_t currArgBit = 0;
    for (std::uint64_t decoded = 0; decoded < expectedElements; ++decoded) {
        if (currCodeBit >= codeStreamMaxBits) {
            err = "code stream truncated before elementCount reached";
            return false;
        }
        std::uint64_t peek = readBits(codeStream, currCodeBit, maxCodeLength);
        const FzipCode& code = lut[static_cast<std::size_t>(peek)];
        if (code.codeLength == 0 || code.codeLength > maxCodeLength) {
            err = "invalid Huffman code encountered";
            return false;
        }
        if (currCodeBit + code.codeLength > codeStreamMaxBits) {
            err = "code stream truncated mid-codeword";
            return false;
        }
        currCodeBit += code.codeLength;

        if (code.prefixLength > 64) {
            err = "invalid prefix length";
            return false;
        }
        int argBits = 64 - static_cast<int>(code.prefixLength);
        if (currArgBit + static_cast<std::uint64_t>(argBits) > argumentStreamMaxBits) {
            err = "argument stream truncated";
            return false;
        }
        std::uint64_t arg = readBits(argumentStream, currArgBit, argBits);
        currArgBit += static_cast<std::uint64_t>(argBits);

        std::uint64_t value = code.prefix;
        if (code.isCommon) {
            if (arg >= commons.size()) {
                err = "common table index out of bounds";
                return false;
            }
            value = commons[static_cast<std::size_t>(arg)];
        } else if (argBits > 0) {
            value |= arg;
        }
        rawStream.push_back(bitsToDouble(value));
    }

    outCodeBitsUsed = currCodeBit;
    outArgBitsUsed = currArgBit;
    return true;
}

bool fzipSerializeV2(const FzipEncoded& encoded,
                     std::vector<std::uint8_t>& outBytes,
                     std::string& err) {
    if (encoded.symbolTable.size() > std::numeric_limits<std::uint16_t>::max() ||
        encoded.commons.size() > std::numeric_limits<std::uint16_t>::max()) {
        err = "table count exceeds v2 header limit";
        return false;
    }

    std::vector<std::uint8_t> payload;
    for (const auto& sym : encoded.symbolTable) {
        auto flagsLen = static_cast<std::uint8_t>(
            (sym.codeLength & 0x0FU) | (sym.isCommon ? 0x80U : 0x00U));
        payload.push_back(flagsLen);
        payload.push_back(sym.prefixLength);
        if (!sym.isCommon && sym.prefixLength > 0) {
            int prefixBytes = (static_cast<int>(sym.prefixLength) + 7) / 8;
            for (int b = 0; b < prefixBytes; ++b) {
                payload.push_back(
                    static_cast<std::uint8_t>((sym.prefix >> (56 - 8 * b)) & 0xFFU));
            }
        }
    }

    std::uint8_t tmp8[8];
    for (std::uint64_t c : encoded.commons) {
        putU64LE(tmp8, c);
        payload.insert(payload.end(), tmp8, tmp8 + 8);
    }

    std::vector<std::uint8_t> codeBytes;
    std::vector<std::uint8_t> argBytes;
    packWordsToBytes(encoded.codeStream, encoded.codeStreamBitLength, codeBytes);
    packWordsToBytes(encoded.argumentStream, encoded.argumentStreamBitLength, argBytes);
    payload.insert(payload.end(), codeBytes.begin(), codeBytes.end());
    payload.insert(payload.end(), argBytes.begin(), argBytes.end());

    std::uint32_t payloadCrc =
        payload.empty() ? 0U : crc32Update(0U, payload.data(), payload.size());

    outBytes.resize(kHeaderSize + payload.size());
    std::uint8_t* hdr = outBytes.data();
    std::memcpy(hdr, kMagic, 4);
    putU16LE(hdr + 4, kFormatVersion);
    hdr[6] = encoded.maxCodeLength;
    hdr[7] = encoded.commonIndexBits;
    putU64LE(hdr + 8, encoded.elementCount);
    putU16LE(hdr + 16, static_cast<std::uint16_t>(encoded.symbolTable.size()));
    putU16LE(hdr + 18, static_cast<std::uint16_t>(encoded.commons.size()));
    putU64LE(hdr + 20, static_cast<std::uint64_t>(codeBytes.size()));
    putU64LE(hdr + 28, static_cast<std::uint64_t>(argBytes.size()));
    putU32LE(hdr + 36, payloadCrc);
    putU32LE(hdr + 40, crc32Update(0U, hdr, 40));

    if (!payload.empty()) {
        std::memcpy(outBytes.data() + kHeaderSize, payload.data(), payload.size());
    }
    return true;
}

bool fzipDeserializeV2(ByteSource& src, std::vector<double>& rawStream, std::string& err) {
    rawStream.clear();

    std::uint8_t hdr[kHeaderSize];
    if (!src.readExact(hdr, sizeof(hdr))) {
        err = "failed reading fzip v2 header";
        return false;
    }
    if (std::memcmp(hdr, kMagic, 4) != 0) {
        err = "invalid fzip magic bytes";
        return false;
    }
    std::uint16_t version = getU16LE(hdr + 4);
    if (version != kFormatVersion) {
        err = "unsupported fzip format version";
        return false;
    }
    if (crc32Update(0U, hdr, 40) != getU32LE(hdr + 40)) {
        err = "fzip header checksum mismatch";
        return false;
    }

    int maxCodeLength = static_cast<int>(hdr[6]);
    int commonIndexBits = static_cast<int>(hdr[7]);
    std::uint64_t elementCount = getU64LE(hdr + 8);
    std::uint16_t symbolCount = getU16LE(hdr + 16);
    std::uint16_t commonCount = getU16LE(hdr + 18);
    std::uint64_t codeStreamBytes = getU64LE(hdr + 20);
    std::uint64_t argumentStreamBytes = getU64LE(hdr + 28);
    std::uint32_t expectedPayloadCrc = getU32LE(hdr + 36);

    if (maxCodeLength < 1 || maxCodeLength > kMaxAllowedCodeLength ||
        commonIndexBits < 0 || commonIndexBits > kMaxCommonBits) {
        err = "invalid fzip header parameters";
        return false;
    }
    if ((commonCount == 0 && commonIndexBits != 0) ||
        (commonCount > 0 && commonIndexBits != ceilLog2(commonCount))) {
        err = "invalid common table dimensions";
        return false;
    }
    if (static_cast<std::size_t>(symbolCount) > (std::size_t{1} << maxCodeLength)) {
        err = "symbolCount exceeds 2^maxCodeLength";
        return false;
    }
    if (elementCount == 0) {
        if (symbolCount != 0 || commonCount != 0 ||
            codeStreamBytes != 0 || argumentStreamBytes != 0) {
            err = "non-empty payload metadata for 0 elements";
            return false;
        }
        if (expectedPayloadCrc != 0U || !src.atEndIfMemory()) {
            err = "invalid payload for empty stream";
            return false;
        }
        return true;
    }
    if (symbolCount == 0) {
        err = "missing symbol table for non-empty stream";
        return false;
    }

    std::uint64_t minCodeBytes =
        (elementCount / 8ULL) + ((elementCount % 8ULL) != 0ULL ? 1ULL : 0ULL);
    std::uint64_t minElementsForArgs =
        (argumentStreamBytes / 8ULL) + ((argumentStreamBytes % 8ULL) != 0ULL ? 1ULL : 0ULL);
    if (codeStreamBytes < minCodeBytes || minElementsForArgs > elementCount) {
        err = "inconsistent elementCount and bitstream sizes";
        return false;
    }
    if (codeStreamBytes > std::numeric_limits<std::size_t>::max() ||
        argumentStreamBytes > std::numeric_limits<std::size_t>::max() ||
        elementCount > std::numeric_limits<std::size_t>::max() / sizeof(double)) {
        err = "stream size exceeds platform limits";
        return false;
    }

    src.resetCrc();
    std::vector<FzipCode> symbols(symbolCount);
    std::uint64_t kraftSum = 0;
    const std::uint64_t kraftTarget = 1ULL << maxCodeLength;
    int prevCodeLen = 0;

    for (std::size_t i = 0; i < symbolCount; ++i) {
        std::uint8_t meta[2];
        if (!src.readExact(meta, 2)) {
            err = "truncated symbol table";
            return false;
        }
        if ((meta[0] & 0x70U) != 0U) {
            err = "reserved symbol flags set";
            return false;
        }
        auto codeLen = static_cast<std::uint8_t>(meta[0] & 0x0FU);
        bool isCommon = (meta[0] & 0x80U) != 0U;
        std::uint8_t prefixLen = meta[1];
        if (codeLen < 1 || codeLen > maxCodeLength || codeLen < prevCodeLen || prefixLen > 64) {
            err = "invalid symbol metadata";
            return false;
        }
        if (isCommon &&
            (commonCount == 0 || prefixLen != static_cast<std::uint8_t>(64 - commonIndexBits))) {
            err = "invalid common symbol prefix length";
            return false;
        }
        prevCodeLen = codeLen;
        kraftSum += (1ULL << (maxCodeLength - codeLen));
        if (kraftSum > kraftTarget) {
            err = "oversubscribed Huffman code lengths";
            return false;
        }

        std::uint64_t prefix = 0;
        if (!isCommon && prefixLen > 0) {
            int prefixBytes = (static_cast<int>(prefixLen) + 7) / 8;
            std::uint8_t prefixBuf[8];
            if (!src.readExact(prefixBuf, static_cast<std::size_t>(prefixBytes))) {
                err = "truncated symbol prefix";
                return false;
            }
            for (int b = 0; b < prefixBytes; ++b) {
                prefix |= static_cast<std::uint64_t>(prefixBuf[b]) << (56 - 8 * b);
            }
            if (prefixLen < 64) {
                std::uint64_t lowMask = (1ULL << (64 - prefixLen)) - 1ULL;
                if ((prefix & lowMask) != 0ULL) {
                    err = "non-zero padding bits in symbol prefix";
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
        err = "invalid canonical Huffman symbol table";
        return false;
    }
    std::vector<FzipCode> lut;
    if (!buildDecodeLut(symbols, maxCodeLength, lut)) {
        err = "failed to construct decode lookup table";
        return false;
    }

    std::vector<std::uint64_t> commons(commonCount, 0ULL);
    if (commonCount > 0) {
        std::vector<std::uint8_t> commonBytes;
        if (!src.readExactVec(commonBytes, static_cast<std::size_t>(commonCount) * 8)) {
            err = "truncated common values table";
            return false;
        }
        for (std::size_t i = 0; i < commonCount; ++i) {
            commons[i] = getU64LE(&commonBytes[i * 8]);
        }
    }

    std::vector<std::uint8_t> codeBytes;
    std::vector<std::uint8_t> argBytes;
    if (!src.readExactVec(codeBytes, static_cast<std::size_t>(codeStreamBytes)) ||
        !src.readExactVec(argBytes, static_cast<std::size_t>(argumentStreamBytes))) {
        err = "truncated bitstream payload";
        return false;
    }
    if (!src.atEndIfMemory()) {
        err = "trailing bytes after fzip payload";
        return false;
    }
    if (src.crc() != expectedPayloadCrc) {
        err = "payload checksum mismatch";
        return false;
    }

    std::vector<std::uint64_t> codeWords;
    std::vector<std::uint64_t> argWords;
    unpackBytesToWords(codeBytes.data(), codeBytes.size(), codeWords);
    unpackBytesToWords(argBytes.data(), argBytes.size(), argWords);

    std::uint64_t usedCodeBits = 0;
    std::uint64_t usedArgBits = 0;
    if (!fzipDecode(lut, maxCodeLength, commons, codeWords, argWords,
                    codeStreamBytes * 8ULL, argumentStreamBytes * 8ULL,
                    elementCount, rawStream, usedCodeBits, usedArgBits, err)) {
        return false;
    }

    if (!checkTrailingZeroBits(codeBytes, usedCodeBits) ||
        !checkTrailingZeroBits(argBytes, usedArgBits)) {
        err = "invalid trailing bitstream padding or unconsumed bytes";
        return false;
    }
    return true;
}

bool readDoubles(std::FILE* in, std::vector<double>& out, std::string& err) {
    out.clear();
    if (in == nullptr) {
        err = "input stream is null";
        return false;
    }
    std::uint8_t buf[8];
    while (true) {
        std::size_t n = std::fread(buf, 1, sizeof(buf), in);
        if (n == sizeof(buf)) {
            out.push_back(bitsToDouble(getU64LE(buf)));
        } else if (n == 0) {
            if (std::ferror(in)) {
                err = "failed reading input stream";
                return false;
            }
            return true;
        } else {
            err = "input byte count is not a multiple of 8";
            return false;
        }
    }
}

bool writeAll(const void* data, std::size_t elemSize, std::size_t count, std::FILE* out) {
    if (count == 0) {
        return true;
    }
    return std::fwrite(data, elemSize, count, out) == count;
}

}  // namespace

Result<std::vector<std::uint8_t>> compress(
    const double* data,
    std::size_t count,
    const CompressOptions& options) {
    Result<std::vector<std::uint8_t>> result;
    FzipEncoded encoded;
    if (!fzipEncode(data, count, options, encoded, result.error)) {
        return result;
    }
    if (!fzipSerializeV2(encoded, result.value, result.error)) {
        result.value.clear();
        return result;
    }
    result.ok = true;
    return result;
}

Result<std::vector<std::uint8_t>> compress(
    const std::vector<double>& input,
    const CompressOptions& options) {
    return compress(input.data(), input.size(), options);
}

Result<std::vector<double>> decompress(const std::uint8_t* data, std::size_t size) {
    Result<std::vector<double>> result;
    if (data == nullptr && size > 0) {
        result.error = "input buffer pointer is null";
        return result;
    }
    ByteSource src(data, size);
    if (!fzipDeserializeV2(src, result.value, result.error)) {
        result.value.clear();
        return result;
    }
    result.ok = true;
    return result;
}

Result<std::vector<double>> decompress(const std::vector<std::uint8_t>& input) {
    return decompress(input.data(), input.size());
}

Status compressStream(std::FILE* in, std::FILE* out, const CompressOptions& options) {
    Status status;
    if (in == nullptr || out == nullptr) {
        status.error = "input or output stream is null";
        return status;
    }
    std::vector<double> rawStream;
    if (!readDoubles(in, rawStream, status.error)) {
        return status;
    }
    auto compressed = compress(rawStream, options);
    if (!compressed) {
        status.error = std::move(compressed.error);
        return status;
    }
    if (!writeAll(compressed.value.data(), 1, compressed.value.size(), out) ||
        std::fflush(out) != 0) {
        status.error = "failed writing compressed output";
        return status;
    }
    status.ok = true;
    return status;
}

Status decompressStream(std::FILE* in, std::FILE* out) {
    Status status;
    if (in == nullptr || out == nullptr) {
        status.error = "input or output stream is null";
        return status;
    }
    std::vector<double> rawStream;
    ByteSource src(in);
    if (!fzipDeserializeV2(src, rawStream, status.error)) {
        return status;
    }
    std::uint8_t buf[8];
    for (double val : rawStream) {
        putU64LE(buf, doubleToBits(val));
        if (std::fwrite(buf, 1, sizeof(buf), out) != sizeof(buf)) {
            status.error = "failed writing decompressed output";
            return status;
        }
    }
    if (std::fflush(out) != 0) {
        status.error = "failed flushing decompressed output";
        return status;
    }
    status.ok = true;
    return status;
}

}  // namespace fzip
