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

#include <sys/wait.h>
#include <unistd.h>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

namespace {

int gFailures = 0;
int gChecks = 0;

void checkImpl(bool cond, const char* expr, int line) {
    ++gChecks;
    if (!cond) {
        ++gFailures;
        std::cerr << "FAIL [line " << line << "]: " << expr << "\n";
    }
}

#define CHECK(expr) checkImpl(static_cast<bool>(expr), #expr, __LINE__)

std::uint32_t crc32(const std::uint8_t* data, std::size_t len) {
    std::uint32_t crc = ~0U;
    for (std::size_t i = 0; i < len; ++i) {
        crc ^= data[i];
        for (int b = 0; b < 8; ++b) {
            crc = (crc & 1U) ? ((crc >> 1) ^ 0xEDB88320U) : (crc >> 1);
        }
    }
    return ~crc;
}

void putU16LE(std::uint8_t* dst, std::uint16_t v) {
    dst[0] = static_cast<std::uint8_t>(v & 0xFFU);
    dst[1] = static_cast<std::uint8_t>((v >> 8) & 0xFFU);
}

void putU32LE(std::uint8_t* dst, std::uint32_t v) {
    dst[0] = static_cast<std::uint8_t>(v & 0xFFU);
    dst[1] = static_cast<std::uint8_t>((v >> 8) & 0xFFU);
    dst[2] = static_cast<std::uint8_t>((v >> 16) & 0xFFU);
    dst[3] = static_cast<std::uint8_t>((v >> 24) & 0xFFU);
}

void putU64LE(std::uint8_t* dst, std::uint64_t v) {
    for (int i = 0; i < 8; ++i) {
        dst[i] = static_cast<std::uint8_t>((v >> (8 * i)) & 0xFFU);
    }
}

double bitsToDouble(std::uint64_t bits) {
    double val = 0.0;
    std::memcpy(&val, &bits, sizeof(val));
    return val;
}

std::vector<std::uint8_t> buildV2Packet(std::uint8_t maxCl,
                                        std::uint8_t commonBits,
                                        std::uint64_t elementCount,
                                        std::uint16_t symbolCount,
                                        std::uint16_t commonCount,
                                        std::uint64_t codeStreamBytes,
                                        std::uint64_t argStreamBytes,
                                        const std::vector<std::uint8_t>& payload,
                                        bool overridePayloadCrc = false,
                                        std::uint32_t forcedPayloadCrc = 0) {
    std::uint32_t payloadCrc =
        overridePayloadCrc
            ? forcedPayloadCrc
            : (payload.empty() ? 0U : crc32(payload.data(), payload.size()));

    std::vector<std::uint8_t> out(fzip::kHeaderSize + payload.size(), 0U);
    std::uint8_t* hdr = out.data();
    hdr[0] = 'F';
    hdr[1] = 'Z';
    hdr[2] = 'I';
    hdr[3] = 'P';
    putU16LE(hdr + 4, fzip::kFormatVersion);
    hdr[6] = maxCl;
    hdr[7] = commonBits;
    putU64LE(hdr + 8, elementCount);
    putU16LE(hdr + 16, symbolCount);
    putU16LE(hdr + 18, commonCount);
    putU64LE(hdr + 20, codeStreamBytes);
    putU64LE(hdr + 28, argStreamBytes);
    putU32LE(hdr + 36, payloadCrc);
    putU32LE(hdr + 40, crc32(hdr, 40));

    if (!payload.empty()) {
        std::memcpy(out.data() + fzip::kHeaderSize, payload.data(), payload.size());
    }
    return out;
}

bool bitEqual(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        return false;
    }
    if (a.empty()) {
        return true;
    }
    return std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

void verifyRoundtrip(const std::vector<double>& input,
                     const fzip::CompressOptions& opts = fzip::CompressOptions{}) {
    auto comp = fzip::compress(input, opts);
    CHECK(comp.ok);
    if (!comp) {
        std::cerr << "  compress error: " << comp.error << "\n";
        return;
    }
    auto dec = fzip::decompress(comp.value);
    CHECK(dec.ok);
    if (!dec) {
        std::cerr << "  decompress error: " << dec.error << "\n";
        return;
    }
    CHECK(bitEqual(input, dec.value));
}

struct TempDir {
    std::string path;
    bool valid = false;

    TempDir() {
        char tmpl[] = "/tmp/fzip_test_XXXXXX";
        char* created = ::mkdtemp(tmpl);
        if (created != nullptr) {
            path = created;
            valid = true;
        }
    }

    ~TempDir() {
        if (valid) {
            std::string cmd = "rm -rf '" + path + "'";
            if (std::system(cmd.c_str()) != 0) {}
        }
    }
};

int runCmdExitCode(const std::string& cmd) {
    int status = std::system(cmd.c_str());
    if (status == -1 || !WIFEXITED(status)) {
        return -1;
    }
    return WEXITSTATUS(status);
}

bool writeRawBytes(const std::string& path, const void* data, std::size_t size) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (f == nullptr) {
        return false;
    }
    bool ok = (size == 0) || (std::fwrite(data, 1, size, f) == size);
    if (std::fclose(f) != 0) {
        ok = false;
    }
    return ok;
}

bool readRawBytes(const std::string& path, std::vector<std::uint8_t>& out) {
    out.clear();
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (f == nullptr) {
        return false;
    }
    std::uint8_t buf[4096];
    while (true) {
        std::size_t n = std::fread(buf, 1, sizeof(buf), f);
        if (n > 0) {
            out.insert(out.end(), buf, buf + n);
        }
        if (n < sizeof(buf)) {
            bool err = std::ferror(f) != 0;
            std::fclose(f);
            return !err;
        }
    }
}

void testRoundtripMatrix() {
    // Empty inputs
    {
        auto emptyVec = fzip::compress(std::vector<double>{});
        CHECK(emptyVec.ok);
        CHECK(emptyVec.value.size() == fzip::kHeaderSize);
        auto dec = fzip::decompress(emptyVec.value);
        CHECK(dec.ok && dec.value.empty());

        auto emptyPtr = fzip::compress(nullptr, 0);
        CHECK(emptyPtr.ok);
        CHECK(emptyPtr.value.size() == fzip::kHeaderSize);

        auto nullNonZero = fzip::compress(nullptr, 1);
        CHECK(!nullNonZero.ok && !nullNonZero.error.empty());
    }

    // Small & identical runs
    verifyRoundtrip({3.14159});
    verifyRoundtrip({1.0, -2.0});
    {
        std::vector<double> same100(100, 42.5);
        auto comp = fzip::compress(same100);
        CHECK(comp.ok);
        CHECK(comp.value.size() < 100);
        auto dec = fzip::decompress(comp.value);
        CHECK(dec.ok && bitEqual(same100, dec.value));
    }
    verifyRoundtrip(std::vector<double>(1000, -123.75));

    // Special IEEE-754 bit patterns
    {
        const std::uint64_t patterns[] = {
            0x0000000000000000ULL,  // +0.0
            0x8000000000000000ULL,  // -0.0
            0x0000000000000001ULL,  // smallest subnormal
            0x000fffffffffffffULL,  // largest subnormal
            0x7ff0000000000000ULL,  // +Inf
            0xfff0000000000000ULL,  // -Inf
            0x7ff8000000000001ULL,  // quiet NaN
            0x7ff0000000000001ULL,  // signaling NaN
            0xfff8000000000002ULL,  // negative NaN payload
        };
        std::vector<double> specials;
        for (int rep = 0; rep < 50; ++rep) {
            for (std::uint64_t p : patterns) {
                specials.push_back(bitsToDouble(p));
            }
        }
        verifyRoundtrip(specials);
    }

    // Cardinality boundaries across prefix-table and common-table limits
    const std::size_t counts[] = {255, 256, 257, 511, 512, 513, 8191, 8192, 8193, 100000};
    for (std::size_t n : counts) {
        std::vector<double> seq(n);
        for (std::size_t i = 0; i < n; ++i) {
            seq[i] = static_cast<double>(i) + 0.123;
        }
        verifyRoundtrip(seq);
    }

    // Skewed distribution (previously triggered negative bitLengthTarget loop)
    {
        std::vector<double> skewed(99000, 1.25);
        for (int i = 0; i < 300; ++i) {
            skewed.push_back(static_cast<double>(i) + 0.5);
        }
        verifyRoundtrip(skewed);
    }

    // Smallest prefix values (previously triggered prev(begin()) UB)
    {
        std::vector<double> smallPref;
        for (int i = 0; i < 50; ++i) {
            smallPref.push_back(bitsToDouble(1ULL));
            smallPref.push_back(bitsToDouble(2ULL));
        }
        for (int i = 0; i < 500; ++i) {
            smallPref.push_back(static_cast<double>(i) + 1000.0);
        }
        verifyRoundtrip(smallPref);
    }

    // Periodic telemetry & random doubles
    {
        std::vector<double> telem(5000);
        for (std::size_t i = 0; i < telem.size(); ++i) {
            telem[i] = static_cast<double>(i % 50) * 0.125 + ((i % 7 == 0) ? 1e-6 : 0.0);
        }
        auto comp = fzip::compress(telem);
        CHECK(comp.ok && comp.value.size() < telem.size() * sizeof(double) / 4);
        auto dec = fzip::decompress(comp.value);
        CHECK(dec.ok && bitEqual(telem, dec.value));
    }
    {
        std::mt19937_64 rng(42);
        std::uniform_real_distribution<double> dist(-1e12, 1e12);
        std::vector<double> rnd(20000);
        for (double& x : rnd) {
            x = dist(rng);
        }
        verifyRoundtrip(rnd);
    }

    // Option sweep & validation
    {
        std::vector<double> data(600);
        for (std::size_t i = 0; i < data.size(); ++i) {
            data[i] = static_cast<double>(i) + 0.25;
        }
        const int maxLens[] = {1, 2, 5, 9, 12, 15};
        for (int mcl : maxLens) {
            verifyRoundtrip(data, fzip::CompressOptions{mcl, 0, 13});
            verifyRoundtrip(data, fzip::CompressOptions{mcl, std::max(0, mcl - 1), 13});
            verifyRoundtrip(data, fzip::CompressOptions{mcl, mcl, 13});
            verifyRoundtrip(data, fzip::CompressOptions{mcl, std::max(0, mcl - 1), 0});
            verifyRoundtrip(data, fzip::CompressOptions{mcl, std::max(0, mcl - 1), 4});
        }

        CHECK(!fzip::compress(data, fzip::CompressOptions{0, 0, 13}).ok);
        CHECK(!fzip::compress(data, fzip::CompressOptions{16, 8, 13}).ok);
        CHECK(!fzip::compress(data, fzip::CompressOptions{9, -1, 13}).ok);
        CHECK(!fzip::compress(data, fzip::CompressOptions{9, 10, 13}).ok);
        CHECK(!fzip::compress(data, fzip::CompressOptions{9, 8, -1}).ok);
        CHECK(!fzip::compress(data, fzip::CompressOptions{9, 8, 14}).ok);
    }
}

void testStreamsAndAdversarial() {
    // Null stream checks
    CHECK(!fzip::compressStream(nullptr, stdout).ok);
    CHECK(!fzip::compressStream(stdin, nullptr).ok);
    CHECK(!fzip::decompressStream(nullptr, stdout).ok);
    CHECK(!fzip::decompressStream(stdin, nullptr).ok);

    // Stream roundtrip via tmpfile
    {
        std::vector<double> input{1.0, -0.0, 3.5, 1.0, 99.25};
        std::FILE* raw = std::tmpfile();
        std::FILE* comp = std::tmpfile();
        std::FILE* dec = std::tmpfile();
        CHECK(raw != nullptr && comp != nullptr && dec != nullptr);
        if (raw && comp && dec) {
            for (double d : input) {
                std::uint64_t bits = 0;
                std::memcpy(&bits, &d, sizeof(bits));
                std::uint8_t b[8];
                putU64LE(b, bits);
                CHECK(std::fwrite(b, 1, 8, raw) == 8);
            }
            std::rewind(raw);
            CHECK(fzip::compressStream(raw, comp).ok);
            std::rewind(comp);
            CHECK(fzip::decompressStream(comp, dec).ok);
            std::rewind(dec);
            std::vector<double> out(input.size());
            std::uint8_t b[8];
            for (std::size_t i = 0; i < input.size(); ++i) {
                CHECK(std::fread(b, 1, 8, dec) == 8);
                std::uint64_t bits = 0;
                for (int k = 0; k < 8; ++k) {
                    bits |= static_cast<std::uint64_t>(b[k]) << (8 * k);
                }
                out[i] = bitsToDouble(bits);
            }
            CHECK(std::fgetc(dec) == EOF);
            CHECK(bitEqual(input, out));
        }
        if (raw) std::fclose(raw);
        if (comp) std::fclose(comp);
        if (dec) std::fclose(dec);
    }

    // Non-8-byte aligned stream input
    {
        std::FILE* bad = std::tmpfile();
        std::FILE* out = std::tmpfile();
        CHECK(bad != nullptr && out != nullptr);
        if (bad && out) {
            const char junk[5] = {'a', 'b', 'c', 'd', 'e'};
            CHECK(std::fwrite(junk, 1, sizeof(junk), bad) == sizeof(junk));
            std::rewind(bad);
            CHECK(!fzip::compressStream(bad, out).ok);
        }
        if (bad) std::fclose(bad);
        if (out) std::fclose(out);
    }

    // Empty / NULL / truncated buffers
    CHECK(!fzip::decompress(nullptr, 0).ok);
    CHECK(!fzip::decompress(nullptr, 10).ok);
    CHECK(!fzip::decompress(std::vector<std::uint8_t>{}).ok);

    std::vector<double> sample;
    for (int i = 0; i < 300; ++i) {
        sample.push_back(static_cast<double>(i % 25) * 0.5);
    }
    auto validComp = fzip::compress(sample);
    CHECK(validComp.ok);
    const std::vector<std::uint8_t>& validBytes = validComp.value;

    for (std::size_t trunc = 0; trunc < fzip::kHeaderSize; ++trunc) {
        CHECK(!fzip::decompress(validBytes.data(), trunc).ok);
    }
    // Truncated payload
    if (validBytes.size() > fzip::kHeaderSize + 1) {
        CHECK(!fzip::decompress(validBytes.data(), validBytes.size() - 1).ok);
    }
    // Trailing garbage byte
    {
        auto extra = validBytes;
        extra.push_back(0x00);
        CHECK(!fzip::decompress(extra).ok);
    }
    // Single-bit flips across header & payload
    for (std::size_t pos : std::vector<std::size_t>{
             0, 3, 4, 6, 7, 8, 16, 18, 20, 28, 36, 40, validBytes.size() - 1}) {
        auto mutated = validBytes;
        mutated[pos] ^= 0x01U;
        CHECK(!fzip::decompress(mutated).ok);
    }

    // Unsupported version & bad magic
    {
        auto badMagic = validBytes;
        badMagic[0] = 'G';
        putU32LE(badMagic.data() + 40, crc32(badMagic.data(), 40));
        CHECK(!fzip::decompress(badMagic).ok);

        auto badVer = validBytes;
        putU16LE(badVer.data() + 4, 1);
        putU32LE(badVer.data() + 40, crc32(badVer.data(), 40));
        CHECK(!fzip::decompress(badVer).ok);
    }

    // Crafted valid-CRC headers & payloads
    // 1. elementCount == 0 with non-zero symbolCount / payload
    CHECK(!fzip::decompress(buildV2Packet(9, 0, 0, 1, 0, 0, 0, {})).ok);
    CHECK(!fzip::decompress(buildV2Packet(9, 0, 0, 0, 0, 0, 0, {}, true, 1U)).ok);

    // 2. Huge elementCount on 1-byte codeStream (minCodeBytes check)
    CHECK(!fzip::decompress(buildV2Packet(9, 0, 1000000000ULL, 1, 0, 1, 0,
                                          std::vector<std::uint8_t>{0x01, 64, 0, 0, 0, 0, 0, 0, 0, 0, 0}))
               .ok);

    // 3. argumentStreamBytes > 8 * elementCount
    CHECK(!fzip::decompress(buildV2Packet(9, 0, 1, 1, 0, 1, 16, {})).ok);

    // 4. Non-minimal commonIndexBits for commonCount == 1
    CHECK(!fzip::decompress(buildV2Packet(9, 1, 1, 1, 1, 1, 0, {})).ok);

    // 5. Huge forged codeStreamBytes on a 44-byte FILE* (must fail without OOM)
    {
        auto dosPacket = buildV2Packet(9, 0, 1, 1, 0, 100ULL * 1024 * 1024 * 1024, 0, {});
        std::FILE* f = std::tmpfile();
        std::FILE* out = std::tmpfile();
        CHECK(f != nullptr && out != nullptr);
        if (f && out) {
            CHECK(std::fwrite(dosPacket.data(), 1, dosPacket.size(), f) == dosPacket.size());
            std::rewind(f);
            CHECK(!fzip::decompressStream(f, out).ok);
        }
        if (f) std::fclose(f);
        if (out) std::fclose(out);
    }

    // 6. Invalid symbol parameters with valid CRCs
    // codeLength == 0
    CHECK(!fzip::decompress(buildV2Packet(9, 0, 1, 1, 0, 1, 0,
                                          {0x00, 64, 0, 0, 0, 0, 0, 0, 0, 0, 0}))
               .ok);
    // reserved flag bits 0x70 set
    CHECK(!fzip::decompress(buildV2Packet(9, 0, 1, 1, 0, 1, 0,
                                          {0x11, 64, 0, 0, 0, 0, 0, 0, 0, 0, 0}))
               .ok);
    // single symbol with codeLength == 2
    CHECK(!fzip::decompress(buildV2Packet(9, 0, 1, 1, 0, 1, 0,
                                          {0x02, 64, 0, 0, 0, 0, 0, 0, 0, 0, 0}))
               .ok);
    // prefixLength == 65 (> 64)
    CHECK(!fzip::decompress(buildV2Packet(9, 0, 1, 1, 0, 1, 0,
                                          {0x01, 65, 0}))
               .ok);
    // non-monotonic codeLength [2, 1]
    CHECK(!fzip::decompress(buildV2Packet(9, 0, 1, 2, 0, 1, 0,
                                          {0x02, 64, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0x01, 64, 0, 0, 0, 0, 0, 0, 0, 0, 0}))
               .ok);
    // oversubscribed Kraft sum [1, 1] for 2 symbols? wait: [1, 1] sums to 1.0!
    // Use 3 symbols with [1, 1, 2] -> oversubscribed!
    CHECK(!fzip::decompress(buildV2Packet(9, 0, 1, 3, 0, 1, 0,
                                          {0x01, 64, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0x01, 64, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0x02, 64, 0, 0, 0, 0, 0, 0, 0, 0, 0}))
               .ok);
    // incomplete Kraft sum [2, 2, 2] for 3 symbols
    CHECK(!fzip::decompress(buildV2Packet(9, 0, 1, 3, 0, 1, 0,
                                          {0x02, 64, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0x02, 64, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0x02, 64, 0, 0, 0, 0, 0, 0, 0, 0, 0}))
               .ok);
    // dirty non-zero bits in prefixLength == 12 (2 prefix bytes, low 4 bits must be 0)
    CHECK(!fzip::decompress(buildV2Packet(9, 0, 1, 1, 0, 1, 7,
                                          {0x01, 12, 0x3F, 0xF1, 0, 0, 0, 0, 0, 0, 0, 0}))
               .ok);
    // dirty trailing padding bit in codeStream (elemCount = 1 uses bit 0; bit 7 is dirty)
    CHECK(!fzip::decompress(buildV2Packet(9, 0, 1, 1, 0, 1, 0,
                                          {0x01, 64, 0, 0, 0, 0, 0, 0, 0, 0, 0x80}))
               .ok);
    // dirty trailing padding bit in argumentStream (prefixLength = 60 uses 4 bits; bit 7 is dirty)
    CHECK(!fzip::decompress(buildV2Packet(9, 0, 1, 1, 0, 1, 1,
                                          {0x01, 60, 0, 0, 0, 0, 0, 0, 0, 0, 0x00, 0x80}))
               .ok);
    // isCommon with commonCount == 0
    CHECK(!fzip::decompress(buildV2Packet(9, 0, 1, 1, 0, 1, 0, {0x81, 64, 0x00})).ok);
}

void testCli(const std::string& fzipBin) {
    TempDir dir;
    CHECK(dir.valid);
    if (!dir.valid) {
        return;
    }

    CHECK(runCmdExitCode(fzipBin + " --help >/dev/null") == 0);
    CHECK(runCmdExitCode(fzipBin + " -h >/dev/null") == 0);

    // Invalid CLI invocations -> exit 2
    CHECK(runCmdExitCode(fzipBin + " --Unknown=1 2>/dev/null") == 2);
    CHECK(runCmdExitCode(fzipBin + " --MaxCodeLength=0 2>/dev/null") == 2);
    CHECK(runCmdExitCode(fzipBin + " --MaxCodeLength=16 2>/dev/null") == 2);
    CHECK(runCmdExitCode(fzipBin + " --TargetCodeLength=10 2>/dev/null") == 2);
    CHECK(runCmdExitCode(fzipBin + " -c -d 2>/dev/null") == 2);
    CHECK(runCmdExitCode(fzipBin + " -cd 2>/dev/null") == 2);
    CHECK(runCmdExitCode(fzipBin + " -x 2>/dev/null") == 2);
    CHECK(runCmdExitCode(fzipBin + " a b c 2>/dev/null") == 2);

    // File / runtime errors -> exit 1
    std::string missingIn = dir.path + "/does_not_exist.bin";
    std::string outFz = dir.path + "/out.fz";
    CHECK(runCmdExitCode(fzipBin + " -c '" + missingIn + "' '" + outFz + "' 2>/dev/null") == 1);
    CHECK(runCmdExitCode(fzipBin + " -c '" + dir.path + "/no_dir/out.fz' '" + outFz + "' 2>/dev/null") == 1);

    std::string oddBin = dir.path + "/odd.bin";
    CHECK(writeRawBytes(oddBin, "12345", 5));
    CHECK(runCmdExitCode(fzipBin + " -c '" + oddBin + "' '" + outFz + "' 2>/dev/null") == 1);

    // File-to-file roundtrip with custom code lengths
    std::string inBin = dir.path + "/in.bin";
    std::string decBin = dir.path + "/dec.bin";
    std::string pipeBin = dir.path + "/pipe.bin";
    std::vector<std::uint8_t> rawBytes;
    for (int i = 0; i < 400; ++i) {
        double v = static_cast<double>(i % 20) * 0.25;
        std::uint64_t bits = 0;
        std::memcpy(&bits, &v, sizeof(bits));
        std::uint8_t tmp[8];
        putU64LE(tmp, bits);
        rawBytes.insert(rawBytes.end(), tmp, tmp + 8);
    }
    CHECK(writeRawBytes(inBin, rawBytes.data(), rawBytes.size()));
    CHECK(runCmdExitCode(fzipBin + " -c --MaxCodeLength=12 --TargetCodeLength=10 '" +
                         inBin + "' '" + outFz + "'") == 0);
    CHECK(runCmdExitCode(fzipBin + " -d '" + outFz + "' '" + decBin + "'") == 0);
    std::vector<std::uint8_t> decBytes;
    CHECK(readRawBytes(decBin, decBytes));
    CHECK(decBytes == rawBytes);

    // Stdin/stdout pipe roundtrip
    CHECK(runCmdExitCode(fzipBin + " -c < '" + inBin + "' | " + fzipBin + " -d > '" +
                         pipeBin + "'") == 0);
    std::vector<std::uint8_t> pipeBytes;
    CHECK(readRawBytes(pipeBin, pipeBytes));
    CHECK(pipeBytes == rawBytes);
}

}  // namespace

int main(int argc, char* argv[]) {
    std::string fzipBin = (argc > 1) ? argv[1] : "./fzip";
    testRoundtripMatrix();
    testStreamsAndAdversarial();
    testCli(fzipBin);

    if (gFailures != 0) {
        std::cerr << gFailures << " of " << gChecks << " checks failed.\n";
        return 1;
    }
    std::cout << "All " << gChecks << " checks passed.\n";
    return 0;
}
