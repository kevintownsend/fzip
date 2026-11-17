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
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <string>

namespace {

struct FileCloser {
    void operator()(std::FILE* f) const noexcept {
        if (f != nullptr) {
            std::fclose(f);
        }
    }
};

using UniqueFile = std::unique_ptr<std::FILE, FileCloser>;

void printUsage(const char* progName) {
    std::cout << "Usage: " << progName
              << " [-c | -d] [--MaxCodeLength=N] [--TargetCodeLength=N]"
                 " [inputFile [outputFile]]\n"
                 "  -c                    Compress (default)\n"
                 "  -d                    Decompress\n"
                 "  --MaxCodeLength=N     Max Huffman code length in [1, "
              << fzip::kMaxAllowedCodeLength
              << "] (default: 9)\n"
                 "  --TargetCodeLength=N  Target prefix tree log2 count in [0, MaxCodeLength]\n"
                 "  -h, --help            Show this help message\n";
}

bool startsWith(const std::string& s, const char* prefix) {
    std::size_t len = std::strlen(prefix);
    return s.size() >= len && s.compare(0, len, prefix) == 0;
}

bool parseIntStrict(const std::string& s, int& out) {
    if (s.empty()) {
        return false;
    }
    char* endPtr = nullptr;
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

struct CliArgs {
    bool compress = true;
    bool showHelp = false;
    fzip::CompressOptions options;
    std::string inputPath;
    std::string outputPath;
};

bool parseCli(int argc, char* argv[], CliArgs& out) {
    bool seenCompress = false;
    bool seenDecompress = false;
    bool hasTargetCodeLength = false;
    int positionalCount = 0;

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "-h" || arg == "--help") {
            out.showHelp = true;
            return true;
        }
        if (startsWith(arg, "--")) {
            if (startsWith(arg, "--MaxCodeLength=")) {
                std::string valStr = arg.substr(std::strlen("--MaxCodeLength="));
                if (!parseIntStrict(valStr, out.options.maxCodeLength) ||
                    out.options.maxCodeLength < 1 ||
                    out.options.maxCodeLength > fzip::kMaxAllowedCodeLength) {
                    std::cerr << "Error: --MaxCodeLength must be an integer in [1, "
                              << fzip::kMaxAllowedCodeLength << "]\n";
                    return false;
                }
            } else if (startsWith(arg, "--TargetCodeLength=")) {
                std::string valStr = arg.substr(std::strlen("--TargetCodeLength="));
                if (!parseIntStrict(valStr, out.options.targetCodeLength) ||
                    out.options.targetCodeLength < 0 ||
                    out.options.targetCodeLength > fzip::kMaxAllowedCodeLength) {
                    std::cerr << "Error: --TargetCodeLength must be an integer in [0, "
                              << fzip::kMaxAllowedCodeLength << "]\n";
                    return false;
                }
                hasTargetCodeLength = true;
            } else {
                std::cerr << "Error: unknown option '" << arg << "'\n";
                return false;
            }
        } else if (arg.size() > 1 && arg[0] == '-') {
            for (std::size_t j = 1; j < arg.size(); ++j) {
                if (arg[j] == 'c') {
                    out.compress = true;
                    seenCompress = true;
                } else if (arg[j] == 'd') {
                    out.compress = false;
                    seenDecompress = true;
                } else if (arg[j] == 'h') {
                    out.showHelp = true;
                    return true;
                } else {
                    std::cerr << "Error: unknown option '-" << arg[j] << "'\n";
                    return false;
                }
            }
            if (seenCompress && seenDecompress) {
                std::cerr << "Error: cannot specify both -c and -d\n";
                return false;
            }
        } else {
            if (positionalCount == 0) {
                out.inputPath = arg;
            } else if (positionalCount == 1) {
                out.outputPath = arg;
            } else {
                std::cerr << "Error: too many positional arguments\n";
                return false;
            }
            positionalCount++;
        }
    }

    if (!hasTargetCodeLength) {
        out.options.targetCodeLength = std::max(0, out.options.maxCodeLength - 1);
    } else if (out.options.targetCodeLength > out.options.maxCodeLength) {
        std::cerr << "Error: --TargetCodeLength (" << out.options.targetCodeLength
                  << ") cannot exceed --MaxCodeLength (" << out.options.maxCodeLength << ")\n";
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char* argv[]) {
    CliArgs args;
    if (!parseCli(argc, argv, args)) {
        return 2;
    }
    if (args.showHelp) {
        printUsage(argc > 0 ? argv[0] : "fzip");
        return 0;
    }

    UniqueFile ownedInput;
    std::FILE* in = stdin;
    if (!args.inputPath.empty()) {
        ownedInput.reset(std::fopen(args.inputPath.c_str(), "rb"));
        if (!ownedInput) {
            std::cerr << "Error: failed to open input file '" << args.inputPath
                      << "': " << std::strerror(errno) << "\n";
            return 1;
        }
        in = ownedInput.get();
    }

    UniqueFile ownedOutput;
    std::FILE* out = stdout;
    if (!args.outputPath.empty()) {
        ownedOutput.reset(std::fopen(args.outputPath.c_str(), "wb"));
        if (!ownedOutput) {
            std::cerr << "Error: failed to open output file '" << args.outputPath
                      << "': " << std::strerror(errno) << "\n";
            return 1;
        }
        out = ownedOutput.get();
    }

    fzip::Status status = args.compress
                              ? fzip::compressStream(in, out, args.options)
                              : fzip::decompressStream(in, out);
    if (!status) {
        std::cerr << "Error: " << status.error << "\n";
        return 1;
    }

    if (ownedOutput && std::fclose(ownedOutput.release()) != 0) {
        std::cerr << "Error: failed closing output file '" << args.outputPath
                  << "': " << std::strerror(errno) << "\n";
        return 1;
    }
    return 0;
}
