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
#pragma once
#ifndef FZIP_FZIP_HPP_
#define FZIP_FZIP_HPP_

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

namespace fzip {

inline constexpr std::uint16_t kFormatVersion = 2;
inline constexpr std::size_t kHeaderSize = 44;
inline constexpr int kMaxAllowedCodeLength = 15;
inline constexpr int kMaxCommonBits = 13;

struct CompressOptions {
    int maxCodeLength = 9;
    int targetCodeLength = 8;
    int maxCommonBits = 13;
};

struct Status {
    bool ok = false;
    std::string error;

    explicit operator bool() const noexcept { return ok; }
};

template <typename T>
struct Result {
    bool ok = false;
    T value{};
    std::string error;

    explicit operator bool() const noexcept { return ok; }
};

Result<std::vector<std::uint8_t>> compress(
    const double* data,
    std::size_t count,
    const CompressOptions& options = CompressOptions{});

Result<std::vector<std::uint8_t>> compress(
    const std::vector<double>& input,
    const CompressOptions& options = CompressOptions{});

Result<std::vector<double>> decompress(const std::uint8_t* data, std::size_t size);

Result<std::vector<double>> decompress(const std::vector<std::uint8_t>& input);

Status compressStream(
    std::FILE* in,
    std::FILE* out,
    const CompressOptions& options = CompressOptions{});

Status decompressStream(std::FILE* in, std::FILE* out);

}  // namespace fzip

#endif  // FZIP_FZIP_HPP_
