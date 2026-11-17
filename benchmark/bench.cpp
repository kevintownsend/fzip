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

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace {

struct Dataset {
    std::string name;
    std::vector<double> data;
};

std::vector<Dataset> generateDatasets(std::size_t count) {
    std::vector<Dataset> datasets;

    datasets.push_back({"constant", std::vector<double>(count, 1.0)});

    {
        std::vector<double> v(count);
        for (std::size_t i = 0; i < count; ++i) {
            v[i] = static_cast<double>(i);
        }
        datasets.push_back({"monotonic_int", std::move(v)});
    }

    {
        std::vector<double> v(count);
        for (std::size_t i = 0; i < count; ++i) {
            v[i] = std::sin(static_cast<double>(i) * 0.001);
        }
        datasets.push_back({"sine_wave", std::move(v)});
    }

    {
        std::vector<double> v(count);
        for (std::size_t i = 0; i < count; ++i) {
            v[i] = static_cast<double>(i % 64) * 0.125 + ((i % 97 == 0) ? 1e-6 : 0.0);
        }
        datasets.push_back({"quantized_telemetry", std::move(v)});
    }

    {
        std::mt19937_64 rng(2026);
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        std::vector<double> v(count);
        for (double& x : v) {
            x = dist(rng);
        }
        datasets.push_back({"uniform_random", std::move(v)});
    }

    return datasets;
}

}  // namespace

int main() {
    constexpr std::size_t kElements = 250000;  // 2 MB raw per dataset
    constexpr int kIters = 3;
    const auto datasets = generateDatasets(kElements);

    std::cout << std::left << std::setw(22) << "Dataset"
              << std::right << std::setw(12) << "Raw (B)"
              << std::setw(15) << "Compressed (B)"
              << std::setw(10) << "Ratio"
              << std::setw(16) << "Comp (MB/s)"
              << std::setw(16) << "Decomp (MB/s)"
              << std::setw(11) << "Verified\n";
    std::cout << std::string(102, '-') << "\n";

    for (const auto& ds : datasets) {
        const std::size_t rawBytes = ds.data.size() * sizeof(double);
        const double rawMB = static_cast<double>(rawBytes) / (1024.0 * 1024.0);

        auto warmupComp = fzip::compress(ds.data);
        if (!warmupComp) {
            std::cerr << "Benchmark compress failed on " << ds.name << ": "
                      << warmupComp.error << "\n";
            return 1;
        }
        auto warmupDec = fzip::decompress(warmupComp.value);
        if (!warmupDec || warmupDec.value.size() != ds.data.size() ||
            std::memcmp(warmupDec.value.data(), ds.data.data(), rawBytes) != 0) {
            std::cerr << "Benchmark roundtrip verification failed on " << ds.name << "\n";
            return 1;
        }

        double compSec = 0.0;
        double decompSec = 0.0;
        std::size_t compressedBytes = warmupComp.value.size();

        for (int it = 0; it < kIters; ++it) {
            auto t0 = std::chrono::steady_clock::now();
            auto c = fzip::compress(ds.data);
            auto t1 = std::chrono::steady_clock::now();
            auto d = fzip::decompress(c.value);
            auto t2 = std::chrono::steady_clock::now();

            if (!c || !d || d.value.size() != ds.data.size() ||
                std::memcmp(d.value.data(), ds.data.data(), rawBytes) != 0) {
                std::cerr << "Benchmark verification failed on iteration " << it << "\n";
                return 1;
            }
            compressedBytes = c.value.size();
            compSec += std::chrono::duration<double>(t1 - t0).count();
            decompSec += std::chrono::duration<double>(t2 - t1).count();
        }

        double compMBs = (rawMB * kIters) / compSec;
        double decompMBs = (rawMB * kIters) / decompSec;
        double ratio = static_cast<double>(rawBytes) / static_cast<double>(compressedBytes);

        std::cout << std::left << std::setw(22) << ds.name
                  << std::right << std::setw(12) << rawBytes
                  << std::setw(15) << compressedBytes
                  << std::setw(9) << std::fixed << std::setprecision(2) << ratio << "x"
                  << std::setw(16) << std::fixed << std::setprecision(1) << compMBs
                  << std::setw(16) << std::fixed << std::setprecision(1) << decompMBs
                  << std::setw(10) << "OK\n";
    }
    return 0;
}
