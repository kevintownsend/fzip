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

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size) {
    auto dec = fzip::decompress(data, size);
    if (dec && dec.value.size() <= 100000) {
        auto comp = fzip::compress(dec.value);
        assert(comp.ok);
        auto roundtrip = fzip::decompress(comp.value);
        assert(roundtrip.ok);
        assert(roundtrip.value.size() == dec.value.size());
        if (!dec.value.empty()) {
            assert(std::memcmp(dec.value.data(), roundtrip.value.data(),
                               dec.value.size() * sizeof(double)) == 0);
        }
    }
    return 0;
}
