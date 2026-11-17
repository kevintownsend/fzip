# fzip

`fzip` is a lossless, bit-exact C++17 compression library and CLI tool for IEEE-754 64-bit (`double`) floating-point streams.

It partitions unique bit patterns in a prefix trie, extracts high-frequency repeated values into a compact dictionary (`commons`), assigns length-limited canonical Huffman codes via Package-Merge, and emits an explicit little-endian `.fz` v2 stream protected by header and payload CRC-32 checksums.

## Features

- **Bit-exact roundtrip:** Preserves every 64-bit IEEE-754 bit pattern identically (`-0.0`, subnormals, `±Inf`, and quiet/signaling `NaN` payloads).
- **Compact, portable `.fz` v2 wire format:** 44-byte explicit little-endian header with `"FZIP"` magic, format version, exact element count, minimal symbol & dictionary tables, byte-aligned bitstreams, and dual CRC-32 checksums.
- **Hardened decoder:** Rejects corrupted checksums, truncated inputs, non-canonical/oversubscribed Huffman trees, dirty padding bits, and out-of-range section counts without unbounded allocations.
- **Clean C++17 API + CLI:** Static library (`libfzip.a`) with explicit `fzip::Result<T>` / `fzip::Status` error handling and a thin `fzip` command-line tool.

## Project Layout

- `include/fzip/fzip.hpp` — Public C++17 API (`fzip::compress`, `fzip::decompress`, `fzip::compressStream`, `fzip::decompressStream`).
- `src/fzip.cpp` — Core codec implementation and v2 serializer/deserializer.
- `src/main.cpp` — CLI entry point.
- `tests/fzip_test.cpp` — Unit, stream, CLI, and adversarial format tests.
- `tests/fuzz_decompress.cpp` — LibFuzzer entry point (`LLVMFuzzerTestOneInput`).
- `benchmark/bench.cpp` — Self-contained compression ratio and throughput benchmark.

## Build & Test

Requirements: a C++17 compiler (`g++` or `clang++`) and GNU `make`.

```bash
# Build CLI binary (`./fzip`) and static library (`build/libfzip.a`)
make all lib

# Run full test suite (Release -O3)
make test

# Run full test suite under AddressSanitizer + UndefinedBehaviorSanitizer
make clean && make SANITIZE=1 test

# Build and run with clang++
make clean && make CXX=clang++ SANITIZE=1 test

# Run synthetic benchmark suite
make clean && make bench

# Stage/install binary, static library, and header
make install PREFIX=/usr/local
```

## CLI Usage

```text
Usage: fzip [-c | -d] [--MaxCodeLength=N] [--TargetCodeLength=N] [inputFile [outputFile]]
  -c                    Compress (default)
  -d                    Decompress
  --MaxCodeLength=N     Max Huffman code length in [1, 15] (default: 9)
  --TargetCodeLength=N  Target prefix tree log2 count in [0, MaxCodeLength]
  -h, --help            Show this help message
```

Compress and decompress files or streams:

```bash
./fzip -c input.bin output.fz
./fzip -d output.fz restored.bin

# Pipe via stdin / stdout
./fzip -c < input.bin | ./fzip -d > restored.bin
```

Input byte sizes must be multiples of 8 bytes (64-bit IEEE-754 doubles). Invalid CLI arguments exit with code `2`; I/O, format, and checksum errors exit with code `1`.

## C++ Library Example

```cpp
#include "fzip/fzip.hpp"
#include <cassert>
#include <vector>

int main() {
    std::vector<double> samples{1.0, 2.5, 1.0, 2.5, 42.0};
    auto compressed = fzip::compress(samples);
    assert(compressed.ok);

    auto restored = fzip::decompress(compressed.value);
    assert(restored.ok && restored.value == samples);
    return 0;
}
```

## Synthetic Benchmarks

Observed on Linux x86-64 (`g++ -O3`, 250,000 `double`s / 2,000,000 bytes per dataset via `make bench`):

| Dataset | Raw (B) | Compressed (B) | Ratio | Comp (MB/s) | Decomp (MB/s) | Verified |
|---|---:|---:|---:|---:|---:|---|
| `constant` | 2,000,000 | 31,304 | 63.89x | 1223.1 | 1129.3 | OK |
| `quantized_telemetry` | 2,000,000 | 197,522 | 10.13x | 389.8 | 892.5 | OK |
| `monotonic_int` | 2,000,000 | 1,688,944 | 1.18x | 4.2 | 179.1 | OK |
| `sine_wave` | 2,000,000 | 1,699,942 | 1.18x | 3.9 | 176.3 | OK |
| `uniform_random` | 2,000,000 | 1,689,005 | 1.18x | 3.4 | 179.5 | OK |

## License

Licensed under the [Apache License, Version 2.0](LICENSE).
