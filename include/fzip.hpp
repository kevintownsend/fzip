
#include <stdio.h>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <vector>
#include <map>
#include <algorithm>
#include <set>
#include <iostream>

using namespace std;

#ifndef FZIP_H
#define FZIP_H
#define nullptr NULL

struct FzipOptions{
    FILE* inputFile;
    FILE* outputFile;
    bool compress;
    int maxCodeLength;
    int targetCodeLength;
    int commonCodes;
    int log2CommonCodes;
    FzipOptions(int argc, char* argv[]){
        compress = true;
        maxCodeLength = 9;
        targetCodeLength = this->maxCodeLength-1;
        commonCodes = 8192;
        log2CommonCodes = log2(this->commonCodes);
        int nonOptionArgs = 0;
        int currentArg = 1;
        while(currentArg != argc){
            string arg(argv[currentArg]);
            if(arg.substr(0,2) == "--"){
                if(arg.find("MaxCodeLength")){
                    this->maxCodeLength = atoi(arg.substr(arg.find("=")+1,arg.size()).c_str());
                } else if(arg.find("TargetCodeLength")){
                    this->targetCodeLength = atoi(arg.substr(arg.find("=")+1,arg.size()).c_str());
                }
            }else if(arg[0] == '-'){
                cerr << "found -" << endl;
                if(arg.find("c") != string::npos){
                    cerr << "found -c" << endl;
                    this->compress = true;
                }
                if(arg.find("d") != string::npos){
                    cerr << "found -d" << endl;
                    this->compress = false;
                }
            }else{
                if(nonOptionArgs == 0)
                    inputFile = fopen(arg.c_str(),"r");
                else if(nonOptionArgs == 1)
                    outputFile = fopen(arg.c_str(),"w");
                else
                    cerr << "this should not be reached" << endl;
                nonOptionArgs++;
            }
            currentArg++;
        }
        if(nonOptionArgs < 2)
            this->outputFile = stdout;
        if(nonOptionArgs < 1){
            cerr << "assigning input to stdin" << endl;
            this->inputFile = stdin;
        }
    }
};

typedef unsigned long long ull;

struct FzipHeader{
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
};

struct FzipCode{
    ull codeLength : 4, isCommon : 1, prefixLength : 7, code : 9;
    ull prefix;
    FzipCode(){
        codeLength = 0;
        isCommon = 0;
        prefixLength = 0;
        code = 0;
        prefix = 0;
    }
    void print(){
        cerr << "codeLength: " << codeLength << endl;
        cerr << "isCommon: " << isCommon << endl;
        cerr << "prefixLength: " << prefixLength << endl;
        cerr << "code: " << code << endl;
        cerr << "prefix: " << prefix << endl;
    }
};

void printBits(ull val){
    for(int i = 0; i < 64; ++i){
        if(val & 1LL << i)
            cerr << "1";
        else
            cerr << "0";
    }
    cerr << endl;
}


vector<double> readFloats(FILE* inputFile){
    vector<double> ret;
    double tmp;
    fread(&tmp, 8, 1, inputFile);
    while(!feof(inputFile)){
        ret.push_back(tmp);
        fread(&tmp, 8, 1, inputFile);
    }
    return ret;
}

bool fzipCompress(vector<double> &rawStream, vector<ull> &commons, vector<FzipCode> &codes, vector<ull> &codeStream, vector<ull> &argumentStream, ull &codeStreamBitLength, ull &argumentStreamBitLength, bool staticCommon = false);

bool fzipCompress(FzipOptions options){
    vector<ull> commons;
    vector<FzipCode> codes;
    vector<ull> codeStream;
    vector<ull> argumentStream;
    vector<double> rawStream = readFloats(options.inputFile);
    cerr << "size of rawStream: " << rawStream.size() << endl;
    cerr << "compressing values" << endl;
    ull codeStreamBitLength, argumentStreamBitLength;
    fzipCompress(rawStream, commons, codes, codeStream, argumentStream, codeStreamBitLength, argumentStreamBitLength); //TODO: options
    cerr << "what muther: " << codeStreamBitLength << endl;
    FzipHeader header;
    header.codeStreamBitLength = codeStreamBitLength;
    header.argumentStreamBitLength = argumentStreamBitLength;
    size_t tmp = sizeof(FzipHeader);
    header.codesPointer = tmp;
    tmp += codes.size() * sizeof(FzipCode);
    header.commonPointer = tmp;
    tmp += commons.size() * sizeof(ull);
    header.codeStreamPointer = tmp;
    tmp += codeStream.size() * sizeof(ull);
    header.argumentStreamPointer = tmp;
    tmp += argumentStream.size() * sizeof(ull);
    header.size = tmp;
    fwrite(&header, sizeof(FzipHeader), 1, options.outputFile);
    fwrite(&codes[0], sizeof(FzipCode), codes.size(), options.outputFile);
    fwrite(&commons[0], sizeof(ull), commons.size(), options.outputFile);
    fwrite(&codeStream[0], sizeof(ull), codeStream.size(), options.outputFile);
    fwrite(&argumentStream[0], sizeof(ull), argumentStream.size(), options.outputFile);
    return true;
};

bool fzipDecompress(vector<double> &rawStream, vector<ull> &commons, vector<FzipCode> &codes, vector<ull> &codeStream, vector<ull> &argumentStream, ull &codeStreamBitLength, ull &argumentStreamBitLength);

bool fzipDecompress(FzipOptions options){
    //TODO: 
    vector<ull> commons;
    vector<FzipCode> codes;
    vector<ull> codeStream;
    vector<ull> argumentStream;
    FzipHeader header;
    fread(&header, sizeof(FzipHeader), 1, options.inputFile);
    codes.resize((header.commonPointer - header.codesPointer) / sizeof(FzipCode));
    commons.resize((header.codeStreamPointer - header.commonPointer) / sizeof(ull));
    codeStream.resize((header.argumentStreamPointer - header.codeStreamPointer) / sizeof(ull));
    argumentStream.resize((header.size - header.argumentStreamPointer) / sizeof(ull));
    fread(&codes[0], sizeof(FzipCode), codes.size(), options.inputFile);
    fread(&commons[0], sizeof(ull), commons.size(), options.inputFile);
    fread(&codeStream[0], sizeof(ull), codeStream.size(), options.inputFile);
    fread(&argumentStream[0], sizeof(ull), argumentStream.size(), options.inputFile);

    cerr << "codes size: " << codes.size() << endl;
    cerr << "commons size: " << commons.size() << endl;
    cerr << "codeStream size: " << codeStream.size() << endl;
    cerr << "argumentStream size: " << argumentStream.size() << endl;

    vector<double> rawStream;
    fzipDecompress(rawStream, commons, codes, codeStream, argumentStream, header.codeStreamBitLength, header.argumentStreamBitLength);
    fwrite(&rawStream[0], sizeof(double), rawStream.size(), options.outputFile);
    //TODO: write rawStream
};


vector<FzipCode> createFzipCodes(vector<pair<ull, ull> > &values, int targetCodeCount);

vector<FzipCode> huffmanCoding(vector<pair<ull, FzipCode> > codeFrequencies);

bool fzipCompress(vector<double> &rawStream, vector<ull> &commons, vector<FzipCode> &codes, vector<ull> &codeStream, vector<ull> &argumentStream, ull &codeStreamBitLength, ull &argumentStreamBitLength, bool staticCommon){
    if(!staticCommon)
        commons.clear();
    codes.clear();
    codeStream.clear();
    argumentStream.clear();
    ull* rawStreamUll = (ull*)&rawStream[0];
    //int targetCodeCount = pow(2,8);
    cerr << "Sorting by repeats" << endl;
    //sort by repeats
    map<ull, ull> frequencies;
    for(auto it = rawStream.begin(); it != rawStream.end(); ++it){
        frequencies[*(ull*)&*it]++;
    }

    vector<pair<ull,ull> > sortedByFrequency;
    for(auto it = frequencies.begin(); it != frequencies.end(); ++it){
        sortedByFrequency.push_back(make_pair(it->second, it->first));
    }
    sort(sortedByFrequency.begin(), sortedByFrequency.end());

    //create first fzip codes
    cerr << "creating initial codes" << endl;
    codes = createFzipCodes(sortedByFrequency, pow(2,8));
    //codes = createFzipCodes(sortedByFrequency, 4);
    //TODO: look at initial codes
    cerr << "codes size: " << codes.size() << endl;

    //create commons
    cerr << "creating the commons" << endl;
    int i = 0;
    if(!staticCommon){
        commons.resize((int)pow(2,13));
        for(auto it = frequencies.begin(); i < commons.size() && it != frequencies.end(); ++i, ++it){
            commons[i] = it->second;
        }
    }
    FzipCode commonCode;
    commonCode.isCommon = 1;
    commonCode.prefixLength = 64 - 13;
    codes.push_back(commonCode);

    //TODO: update fzipcode frequencies
    cerr << "updating fzipcode frequencies" << endl;
    vector<pair<ull, FzipCode> > codeFrequencies;
    map<ull, int> mapToPrefixCode;
    i = 0;
    for(auto it = codes.begin(); it != codes.end(); ++it){
        if(!it->isCommon)
            mapToPrefixCode[it->prefix] = codeFrequencies.size();
        codeFrequencies.push_back(make_pair(0, *it));
        i++;
    }
    cerr << "creating commonSet" << endl;
    set<ull> commonSet(commons.begin(), commons.end());
    cerr << "actually updating code frequencies" << endl;
    for(auto it = frequencies.begin(); it != frequencies.end(); ++it){
        //find code
        int index = prev(mapToPrefixCode.upper_bound(it->first))->second;
        if(codes[index].prefixLength != 64 && commonSet.count(it->first))
            index = codeFrequencies.size()-1;
        codeFrequencies[index].first += it->second;
    }
    //use huffman codes
    cerr << "creating huffman codes" << endl;
    int total = 0;
    for(auto it = codeFrequencies.begin(); it != codeFrequencies.end(); ++it){
        total += it->first;
        //cerr << "freq before: " << it->first << endl;
        if(it->first == 0){
            codeFrequencies.erase(it);
            it--;
            continue;
        }
        if(it->first < rawStream.size() / pow(2,8)){
            cerr << "frequency too small" << endl;
            it->first = rawStream.size() / pow(2,8) + 1;
        }
        //cerr << "freq after: " << it->first << endl;
    }
    cerr << "total: " << total << endl;
    cerr << "check total: " << rawStream.size() << endl;
    //TODO: force code to obey my command
    cerr << "number of codes: " << codeFrequencies.size() << endl;
    cerr << "min code length: " << log2(codeFrequencies.size()) << endl;
    codes = huffmanCoding(codeFrequencies);
    for(auto it = codes.begin(); it != codes.end(); ++it){
        if(it->codeLength > 9){
            cerr << "code greater than 9" << endl;
            exit(1);
        }else if(it->codeLength == 0){
            cerr << "code length equals 0" << endl;
            exit(1);
        }
    }
    mapToPrefixCode.clear();
    i = 0;
    for(auto it = codes.begin(); it != codes.end(); ++it){
        if(!it->isCommon)
            mapToPrefixCode[it->prefix] = i;
        i++;
    }

    //TODO: create final fzip codes
    vector<FzipCode> newCodes;
    newCodes.resize(pow(2,9));
    for(int i = 0; i < codes.size(); i++){
        for(int j = codes[i].code; j < pow(2,9); j += pow(2,codes[i].codeLength))
            newCodes[j] = codes[i];
    }
    //TODO: create codeStream and argumentStream
    map<ull, int> mapToCommon;
    for(int i = 0; i < commons.size(); ++i){
        mapToCommon[commons[i]] = i;
    }
    int codeBufferEndBit = 0;
    int argumentBufferEndBit = 0;

    ull codeBuffer = 0;
    ull argumentBuffer = 0;
    for(int i = 0; i < rawStream.size(); ++i){
        //cerr << "compressing i: " << i << endl;
        int index = prev(mapToPrefixCode.upper_bound(rawStreamUll[i]))->second;
        //cerr << "index: " << index << "/" << codes.size() << endl;
        if(mapToPrefixCode.upper_bound(rawStreamUll[i]) == mapToPrefixCode.begin()){
            //cerr << "that's the problem" << endl;
            index = codes.size() - 1;
        }
        if(codes[index].prefixLength == 64 && codes[index].prefix == rawStreamUll[i]){
        } else if(commonSet.count(rawStreamUll[i]))
            index = codes.size() - 1;
//        cerr << "info: " << endl;
//        cerr << "codes size: " << codes.size() << endl;
//        cerr << "index: " << index << endl;
//        codes[index].print();
        if(codes[index].codeLength == 0){
            cerr << "fail" << endl;
            exit(1);
        }
        codeBuffer |= (ull)codes[index].code << codeBufferEndBit;
        codeBufferEndBit += codes[index].codeLength;
//        cerr << "putting in code: " << i << endl;
//        printBits(codes[index].code);
//        printBits(codeBuffer);
//        cerr << "wtf: " << codeBufferEndBit << endl;
//        cerr << "code length: " << codes[index].codeLength << endl;
        if(codeBufferEndBit >= 64){
            codeStream.push_back(codeBuffer);
            codeBufferEndBit -= 64;
            codeBuffer = codes[index].code >> (codes[index].codeLength - codeBufferEndBit);
        }
        ull argument;
        int argumentLength;
        if(codes[index].isCommon){
            argument = mapToCommon[rawStreamUll[i]];
            argumentLength = 13;
        }else if(codes[index].prefixLength == 64){
            argument = 0;
            argumentLength = 0;
        }else{
            ull mask = -1LL;
            argument = rawStreamUll[i] & ((ull)-1LL >> codes[index].prefixLength);
            argumentLength = 64 - codes[index].prefixLength;
        }
        argumentBuffer |= argument << argumentBufferEndBit;
        argumentBufferEndBit += argumentLength;
        if(argumentBufferEndBit >= 64){
            argumentStream.push_back(argumentBuffer);
            argumentBufferEndBit -= 64;
            argumentBuffer = argument >> (argumentLength - argumentBufferEndBit);
        }
    }
    codeStreamBitLength = codeStream.size() * 64 + codeBufferEndBit;
    argumentStreamBitLength = argumentStream.size() * 64 + argumentBufferEndBit;
    codeStream.push_back(codeBuffer);
    codeStream.push_back(0);
    argumentStream.push_back(argumentBuffer);
    for(int i = 0; i < 5; ++i)
        argumentStream.push_back(0);

    codes = newCodes;
    return true;
};

struct HuffmanNode{
    HuffmanNode* leftChild;
    HuffmanNode* rightChild;
    int frequency;
    FzipCode* codePtr;
    bool traversed;
    HuffmanNode(){
        traversed = 0;
    }
};

bool cmpHuffmanNode(HuffmanNode a, HuffmanNode b){
    return a.frequency < b.frequency;
}

void setHuffmanCodesDFS(HuffmanNode* curr, ull code, ull depth){
    if(curr->traversed){
        cerr << "loop encoutered" << endl;
        return;
    }
    curr->traversed = 1;
    if(depth > 9){
        cerr << "something is wrong" << endl;
        cerr << "depth: " << depth << endl;
        //return;
    }
    if(curr->codePtr != nullptr){
        //cerr << "child" << endl;
        //cerr << "codeptr: " << curr->codePtr << endl;
        curr->codePtr->code = code;
        curr->codePtr->codeLength = depth;
    }else{
        //cerr << "not child" << endl;
        setHuffmanCodesDFS(curr->leftChild, code, depth + 1);
        setHuffmanCodesDFS(curr->rightChild, code | (1 << depth), depth + 1);
    }
}

vector<FzipCode> huffmanCoding(vector<pair<ull, FzipCode> > codeFrequencies){
    cerr << "code frequencies size: " << codeFrequencies.size() << endl;
    if(codeFrequencies.size() == 0)
        codeFrequencies.push_back(make_pair(0, FzipCode()));
    if(codeFrequencies.size() == 1)
        codeFrequencies.push_back(make_pair(0, FzipCode()));
    vector<FzipCode> ret;
    for(int i = 0; i < codeFrequencies.size(); ++i){
        ret.push_back(codeFrequencies[i].second);
    }
    vector<HuffmanNode> tree;
    tree.resize(2 * ret.size());
    HuffmanNode* beginning = &tree[0];
    HuffmanNode* ending = beginning + ret.size();
    cerr << "creating initial nodes" << endl;
    for(int i = 0; i < ret.size(); ++i){
        (beginning + i)->frequency = codeFrequencies[i].first;
        (beginning + i)->codePtr = &ret[0] + i;
        (beginning + i)->leftChild = nullptr;
        (beginning + i)->rightChild = nullptr;
    }
    cerr << "creating intermediate nodes" << endl;
    while(beginning + 1 < ending){
        sort(beginning, ending, cmpHuffmanNode);
        ending->leftChild = beginning;
        ending->rightChild = beginning + 1;
        ending->frequency = beginning->frequency + (beginning+1)->frequency;
        ending->codePtr = nullptr;
        ending++;
        beginning += 2;
    }
    cerr << "size of tree: " << tree.size() << endl;
    cerr << "setting huffman codes" << endl;
    setHuffmanCodesDFS(ending - 1, 0, 0);
    return ret;
}

struct PrefixNode{
    PrefixNode* leftChild;
    PrefixNode* rightChild;
    int frequency;
    PrefixNode(){
        leftChild = nullptr;
        rightChild = nullptr;
        frequency = 0;
    }
    void clear(){
        leftChild = nullptr;
        rightChild = nullptr;
        frequency = 0;
    }
};

bool checkBounds(PrefixNode* begin, PrefixNode* end){
    for(PrefixNode* ptr = begin; ptr != end; ++ptr){
        if((ptr->leftChild < begin || ptr->leftChild >= end) && ptr->leftChild != nullptr)
            return false;
        if((ptr->rightChild < begin || ptr->rightChild >= end) && ptr->rightChild != nullptr)
            return false;
    }
    return true;
}

void dfsPrefixTree(PrefixNode* curr, set<PrefixNode*> &partition, vector<FzipCode> &codes, int depth, ull prefix);


bool sortByFrequency(PrefixNode* left, PrefixNode* right){
    return left->frequency < right->frequency;
}

vector<FzipCode> createFzipCodes(vector<pair<ull, ull> > &values, int targetCodeCount){
    //TODO: fix this function
    vector<FzipCode> ret;
    //Create tree
    cerr << "creating prefix tree" << endl;
    ull allocationSize = sizeof(PrefixNode) * (values.size() * 64 + 1);
    PrefixNode* tree = (PrefixNode*)malloc(allocationSize);
    while(tree == nullptr){
        cerr << "malloc failed reallocating" << endl;
        allocationSize /= 2;
        tree = (PrefixNode*)malloc(allocationSize);
    }
    PrefixNode* treeEnd = tree + (allocationSize / sizeof(PrefixNode));
    PrefixNode* freePtr = tree + 1;
    PrefixNode* root = tree;
    root->clear();
    for(int i = 0; i < values.size(); ++i){
        ull value = values[i].second;
        ull frequency = values[i].first;
        PrefixNode* currNode = root;
        currNode->frequency += frequency;
        for(int j = 0; j < 64; ++j){
            if(value & (1LL << (63 - j))){
                if(currNode->rightChild == nullptr){
                    freePtr->clear();
                    currNode->rightChild = freePtr++;
                }
                currNode = currNode->rightChild;
            }else{
                if(currNode->leftChild == nullptr){
                    freePtr->clear();
                    currNode->leftChild = freePtr++;
                }
                currNode = currNode->leftChild;
            }
            currNode->frequency += frequency;
        }
    }
    if(!checkBounds(root, freePtr))
        cerr << "bounds failure" << endl;
    else
        cerr << "bounds ok" << endl;
    //cerr << "root frequency: " << tree[0].frequ
    int total = tree[0].frequency;

    //Create partition
    cerr << "creating partition" << endl;
    vector<PrefixNode*> partition;
    vector<PrefixNode*> veryCommonValues;
    int currentCodeCount = 1;
    partition.push_back(root);
    cerr << "targetCodeCount: " << targetCodeCount << endl;
    int i = 0;
    while(currentCodeCount < targetCodeCount && partition.size() != 0){
        if(currentCodeCount < 0){
            cerr << "wtf mate" << endl;
            exit(1);
        }
        if(partition.back()->leftChild == nullptr && partition.back()->rightChild == nullptr){
            veryCommonValues.push_back(partition.back());
            partition.pop_back();
            //add to currentCodeCount
            cerr << veryCommonValues.back()->frequency << endl; // - 1;
            cerr << total << endl;
            cerr << (targetCodeCount * 1.0 * veryCommonValues.back()->frequency / total) << endl; // - 1;
            currentCodeCount += targetCodeCount * 1.0 * veryCommonValues.back()->frequency / total; // - 1;
        }else if(partition.back()->leftChild == nullptr){
            PrefixNode* tmp = partition.back()->rightChild;
            partition.pop_back();
            partition.push_back(tmp);
        }else if(partition.back()->rightChild == nullptr){
            PrefixNode* tmp = partition.back()->leftChild;
            partition.pop_back();
            partition.push_back(tmp);
        }else{
            PrefixNode* tmp0 = partition.back()->leftChild;
            PrefixNode* tmp1 = partition.back()->rightChild;
            partition.pop_back();
            partition.push_back(tmp0);
            partition.push_back(tmp1);
            currentCodeCount++;
        }
        sort(partition.begin(), partition.end(), sortByFrequency);
    }
    set<PrefixNode*> partitionSet;
    for(int i = 0; i < partition.size(); ++i){
        partitionSet.insert(partition[i]);
    }
    for(int i = 0; i < veryCommonValues.size(); ++i){
        partitionSet.insert(veryCommonValues[i]);
    }
    dfsPrefixTree(root, partitionSet, ret, 0, 0LL);
    free(tree);
    return ret;
}

void dfsPrefixTree(PrefixNode* curr, set<PrefixNode*> &partition, vector<FzipCode> &codes, int depth, ull prefix){
    if(curr == nullptr)
        return;
    if(partition.count(curr)){
        FzipCode code;
        code.prefixLength = depth;
        code.prefix = prefix;
        codes.push_back(code);
    }else{
        dfsPrefixTree(curr->leftChild, partition, codes, depth + 1, prefix);
        dfsPrefixTree(curr->rightChild, partition, codes, depth + 1, prefix | (1LL << (63 - depth)));
    }
}

bool fzipDecompress(vector<double> &rawStream, vector<ull> &commons, vector<FzipCode> &codes, vector<ull> &codeStream, vector<ull> &argumentStream, ull &codeStreamBitLength, ull &argumentStreamBitLength){
    //TODO: write
    ull currCodeStreamBit = 0;
    ull currArgumentStreamBit = 0;
    int i = 0;
    while(currCodeStreamBit < codeStreamBitLength){
        //cerr << "decoding index: " << i << endl;
        //cerr << "at bit: " << currCodeStreamBit << "/" << codeStreamBitLength << endl;
        ull buffer = codeStream[currCodeStreamBit / 64] >> (currCodeStreamBit % 64);
        if(currCodeStreamBit % 64 > 64 - 9)
            buffer |= codeStream[currCodeStreamBit / 64 + 1] << (64 - currCodeStreamBit % 64);
        FzipCode code = codes[buffer & (ull)-1LL >> 64 - 9];
        ull value = code.prefix;
        if(i == 0){
            cerr << "prefix: " << hex << value << endl;
            cerr << "prefix length: " << dec << code.prefixLength << endl;
        }
        currCodeStreamBit += code.codeLength;

        buffer = argumentStream[currArgumentStreamBit / 64] >> (currArgumentStreamBit % 64);
        if(currArgumentStreamBit % 64)
            buffer |= argumentStream[currArgumentStreamBit / 64 + 1] << (64 - currArgumentStreamBit % 64);
        if(code.prefixLength != 64)
            buffer &= (ull)-1LL >> (code.prefixLength);
        if(i == 0){
            cerr << "argument: " << hex << buffer << endl;
        }
        if(code.isCommon){
            value = commons[buffer];
        }else if(code.prefixLength != 64)
            value |= buffer;
        currArgumentStreamBit += 64 - code.prefixLength;
        rawStream.push_back(*(double*)&value);
        //cerr << "value: " << (*(double*)&value) << endl;
        if(i == 0){
            cerr << "value: " << hex << value << endl;
        }

        i++;
    }
    if(currCodeStreamBit != codeStreamBitLength)
        cerr << "overflow: " << currCodeStreamBit << "/" << codeStreamBitLength << endl;
}

#endif
