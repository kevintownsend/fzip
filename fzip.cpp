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
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <stdint.h>
#include "sortByRepeats.h"
#include "arithmeticFloatingPointCompression.h"
#include <bwt.hpp>
#include <map>

using namespace std;

int bwt(vector<uint64_t> &data);

struct fzip0aFmt_t{
    uint8_t options;
    uint8_t reserved1;
    uint8_t reserved2;
    uint8_t reserved3;
    uint8_t reserved4;
    uint8_t reserved5;
    uint8_t reserved6;
    uint8_t reserved7;
    uint64_t count;
    uint64_t bwtEndPtr;
    uint64_t jilesEncoding;
    uint64_t isRepeatArray;
    uint64_t size;
};
int fzipCompress();

int main(int argc, char* argv[]){
    if(argc == 1){
        fzipCompress();
    }else if(argc == 2){
        if("-d" == string(argv[1])){
            cerr << "decompress" << endl;
        }
    }else{
        cerr << "Usage: " << argv[0] << " < input.bin > output.fz" << endl;
    }
}
int fzipCompress(){
    ofstream logFile("fzip0aData", ofstream::app);
    int c;
    uint64_t value;
    vector<uint64_t> rawInput;
    while(true){
        value = 0;
        for(int i =0; i < 8; i++){
            c = getchar();
            if(c == -1)
                break;
            value |= ((uint64_t)c)<<(i*8);
        }
        if(c == -1)
            break;
        rawInput.push_back(value);
    }
    int count = rawInput.size();
    bwt(rawInput);
    int zeroDeltas = 0;
    for(int i = 0; i < rawInput.size()-1; i++){
        if(rawInput[i] == rawInput[i+1])
            zeroDeltas++;
    }
    cerr << "bwt Compress est: " << ((double)rawInput.size()/(rawInput.size()-zeroDeltas)) << endl;
    //TODO: bit array
    vector<bool> isRepeatArray;
    vector<uint64_t> bwtCompressed;
    bwtCompressed.push_back(rawInput[0]);
    isRepeatArray.push_back(false);
    for(int i = 1; i < rawInput.size(); i++){
        if(rawInput[i-1] == rawInput[i]){
            isRepeatArray.push_back(true);
        }else{
            isRepeatArray.push_back(false);
            bwtCompressed.push_back(rawInput[i]);
        }
    }
    //TODO: new array
    JilesPackage_t* theJilesPart = compressJilesFloatingPoint((double*)&rawInput[0], rawInput.size());
    cerr << "compression of just jiles: " << theJilesPart->size << endl;
    free(theJilesPart);

    theJilesPart = compressJilesFloatingPoint((double*)&bwtCompressed[0], bwtCompressed.size());
    cerr << "compression of bwt jiles: " << theJilesPart->size << endl;
    int finalSize = sizeof(fzip0aFmt_t) + theJilesPart->size + ((isRepeatArray.size()-1)/64)*8+8;
    cerr << "total size: " << finalSize << endl;
    cerr << "compression ratio: " << ((double)count*8/finalSize) << endl;
    logFile << "compression ratio: " << ((double)count*8/finalSize) << endl;
    fzip0aFmt_t* finalData = (fzip0aFmt_t*)malloc(finalSize);
    finalData->size = finalSize;
    finalData->count = count;
    finalData->bwtEndPtr = 0;//TODO: finish
    logFile.close();
    return 0;
}
/*
    //TODO: BWT transform
    int count = rawInput.size();
    uint64_t* commons = (uint64_t*)malloc(sizeof(uint64_t)*rawInput.size());
    int repeatCount = sortByRepeats((double*)&rawInput[0], rawInput.size(), (double*)commons,
    rawInput.size(), 2);
    map<uint64_t, int> commonsReverseMap;
    for(int i = 0; i < repeatCount; i++){
        commonsReverseMap[commons[i]] = i;
    }
    cerr << "repeat count: " << repeatCount << endl;
    //TODO: theory
    int indexWidth = 0;
    int tmp = repeatCount;
    while(tmp){
        tmp >>= 1;
        indexWidth++;
    }
    cerr << "index width: " << dec << indexWidth << endl;
    map<uint64_t, int> repeatMap;
    for(int i = 0; i < repeatCount; i++){
        repeatMap[commons[i]] = 0;
    }
    int repeatsInValues = 0;
    int uniques = 0;
    int indexBits = 0;
    for(int i = 0; i < count; i++){
        if(repeatMap.count(rawInput[i])){
            repeatMap[rawInput[i]]++;
            repeatsInValues++;
        }else{
            uniques++;
        }
    }
    uint64_t sizeInBits = count + indexWidth*repeatsInValues + 64*uniques +
    64*repeatCount;
    cerr << "Compression ratio: " << ((64*count)/((double)sizeInBits)) <<
    endl;
    map<uint64_t, int> runningLastIndexMap;
    uint64_t newSizeInBits = 0;
    for(int i = 0; i < count; i++){
        if(runningLastIndexMap.count(rawInput[i])){
            int tmp = i - runningLastIndexMap[rawInput[i]];
            int indexBits = 0;
            while(tmp){
                tmp >>= 2;
                indexBits += 2;
            }
            newSizeInBits += 1 + indexBits;
        }else{
            newSizeInBits +=65;
        }
        runningLastIndexMap[rawInput[i]] = i;

    }
    cerr << "New compression ratio: " << ((64*count)/((double)newSizeInBits))
    << endl;
    uint64_t repeatIndexBits = 0;
    for(int i = 0; i < count; i++){
        if(commonsReverseMap.count(rawInput[i])){
            int tmp = commonsReverseMap[rawInput[i]];
            int indexBits = 0;
            while(tmp){
                tmp >>=1;
                indexBits++;
            }
            repeatIndexBits += indexBits;
        }else{
        }
    }
    sizeInBits = count + repeatIndexBits + 64*uniques + 64*repeatCount;
    cerr << "another compression ratio: " << ((64*count)/((double)sizeInBits))
    << endl;
    return 0;
}
*/
bool compareBwt(int index1, int index2, vector<uint64_t> &data);
//quicksort fails
struct SortingNode{
    int prevIndex;
    int depth;
    map<uint64_t, SortingNode> children;
    bool isEof;
};
int reachedEnd = 0;
void cleanSortingNode(SortingNode &node, int depth){
    vector<SortingNode*> stack;
    stack.push_back(&node);
    while(stack.size()){
        //cerr << "INFO: " << "stack: " << stack.size() << " children: " << stack.back()->children.size() << endl;
        if(stack[stack.size()-1]->children.size())
            stack.push_back(&(stack[stack.size()-1]->children.begin()->second));
        else{
            //cerr << "reached end: " << reachedEnd++ << endl;
            stack.pop_back();
            if(stack.size())
                stack.back()->children.erase(stack.back()->children.begin());
        }
    }
}

void sortBwt(vector<int> &indexes, int first, int last, vector<uint64_t> &data){
    if(last == 1)
        return;
    SortingNode root;
    root.depth=0;
    root.prevIndex=indexes[first];
    if(data.size()-1 == root.prevIndex)
        root.isEof = true;
    else
        root.isEof = false;
    for(int i = first+1; i < last; i++){
        int index = indexes[i];
        uint64_t value = data[index];

        SortingNode* currentPtr = &root;
        while(true){
            if(currentPtr->children.size() == 0 && currentPtr->isEof == false){
                SortingNode newNode;
                newNode.depth = currentPtr->depth + 1;
                newNode.prevIndex = currentPtr->prevIndex + 1;
                if(newNode.prevIndex == data.size()-1)
                    newNode.isEof = true;
                else
                    newNode.isEof = false;
                currentPtr->children[data[currentPtr->prevIndex]] = newNode;
            }
            if(index==data.size()-1){
                currentPtr->prevIndex = data.size()-1;
                currentPtr->isEof=true;
                break;
            }
            if(currentPtr->children.count(value)){
                currentPtr = &(currentPtr->children[value]);
                index++;
                value = data[index];
            }else{
                SortingNode newNode;
                newNode.depth = currentPtr->depth+1;
                newNode.prevIndex = index+1;
                if(newNode.prevIndex == data.size()-1)
                    newNode.isEof = true;
                else
                    newNode.isEof = false;
                currentPtr->children[value] = newNode;
                break;
            }

        }
    }
    vector<int> newIndexes;
    vector<SortingNode*> stack;
    stack.push_back(&root);
    vector<map<uint64_t,SortingNode>::iterator> itStack;
    itStack.push_back(stack.back()->children.begin());
    while(stack.size()){
        if(stack.back()->children.end() != itStack.back()){
            stack.push_back(&(itStack.back()->second));
            itStack[itStack.size()-1]++;
            itStack.push_back(stack.back()->children.begin());
            
        }else{
            if(stack.back()->children.size() == 0){
                newIndexes.push_back(stack.back()->prevIndex - stack.back()->depth);
            }else if(stack.back()->isEof){
                newIndexes.push_back(stack.back()->prevIndex - stack.back()->depth);
            }
            stack.pop_back();
            itStack.pop_back();
        }
    }
    if(newIndexes.size() == indexes.size()){
    }else{
        cerr << "ERROR: check failed" << endl;
    }
    indexes = newIndexes;
    cleanSortingNode(root, 0);
}
/*
    cerr << "sortBwt called: " << first << ", " << last << endl;
    if(last - first < 2)
        return;
    int middleValue = indexes[first];
    bool topFree = false;
    int bottomPtr = first;
    int topPtr = last - 1;
    while(bottomPtr != topPtr){
        if(topFree){
            if(compareBwt(indexes[bottomPtr], middleValue, data)){
                bottomPtr++;
            }else{
                topFree = false;
                indexes[topPtr] = indexes[bottomPtr];
                topPtr--;
            }
        }else{
            if(compareBwt(middleValue, indexes[topPtr], data)){
                topPtr--;
            }else{
                topFree = true;
                indexes[bottomPtr] = indexes[topPtr];
                bottomPtr++;
            }
        }
    }
    indexes[bottomPtr] = middleValue;
    sortBwt(indexes, first, bottomPtr, data);
    sortBwt(indexes, bottomPtr+1, last, data);

}
*/
bool compareBwt(int index1, int index2, vector<uint64_t> &data){
    if(index1 == index2)
        cerr << "ERROR" << endl;
    if(index1 == data.size())
        return false;
    if(index2 == data.size())
        return true;
    if(data[index1] > data[index2])
        return false;
    if(data[index2] > data[index1])
        return true;
    if(data[index1] == data[index2]){
        //cerr << "excessive compares" << endl;
        return compareBwt(index1+1, index2+1, data);
    }
}
struct bwtSuffixNode{
    int depth;
    uint32_t treeLoc;
    int index;
    uint64_t value;
    int left;
    int right;
};
int bwt(vector<uint64_t> &data){
    vector<bwtSuffixNode> tree;
    tree.resize(data.size()+1);
    int root = data.size();
    tree[root].depth = 0;
    tree[root].treeLoc = 0;
    tree[root].index = root;
    tree[root].value = (uint64_t)(-1);
    tree[root].left = -1;
    tree[root].right = -1;
    bool isLeft;
    for(int currNode = root-1; currNode >= 0; currNode--){
        int nodePtr = root;
        uint64_t value = data[currNode];
        int nextNodePtr = root;
        while(nextNodePtr != -1){
            nodePtr = nextNodePtr;
            uint64_t cmpValue = tree[nodePtr].value;
            if(value > cmpValue){
                nextNodePtr = tree[nodePtr].right;
                isLeft = false;
            }else if(value < cmpValue){
                nextNodePtr = tree[nodePtr].left;
                isLeft = true;
            }else if(value == cmpValue){
                //TODO: compare currNode+1 with NodePtr+1
                int currPlusOne = currNode+1;
                int nodePtrPlusOne = nodePtr+1;
                if(tree[currPlusOne].depth < tree[nodePtrPlusOne].depth){
                    int minDepth = tree[currPlusOne].depth;
                    uint32_t m1 = ((uint32_t)-1)<<(32-minDepth);
                    uint32_t m2 = ((uint32_t)-1)<<(32-minDepth-1);
                    if((tree[currPlusOne].treeLoc & m1) < (tree[nodePtrPlusOne].treeLoc & m2)){
                        nextNodePtr = tree[nodePtr].left;
                        isLeft = true;
                    }else{
                        nextNodePtr = tree[nodePtr].right;
                        isLeft = false;
                    }
                }else{
                    int minDepth = tree[nodePtrPlusOne].depth;
                    uint32_t m1 = ((uint32_t)-1)<<(32-minDepth-1);
                    uint32_t m2 = ((uint32_t)-1)<<(32-minDepth);
                    if((tree[currPlusOne].treeLoc & m1) <= (tree[nodePtrPlusOne].treeLoc & m2)){
                        nextNodePtr = tree[nodePtr].left;
                        isLeft = true;
                    }else{
                        nextNodePtr = tree[nodePtr].right;
                        isLeft = false;
                    }
                }
            }
        }
        if(isLeft){
            tree[nodePtr].left = currNode;
            tree[currNode].treeLoc = tree[nodePtr].treeLoc;
        }else{
            tree[nodePtr].right = currNode;
            tree[currNode].treeLoc = tree[nodePtr].treeLoc | (1<<(31-tree[nodePtr].depth));
        }
        tree[currNode].depth = tree[nodePtr].depth+1;
        tree[currNode].index = currNode;
        tree[currNode].value = value;
        tree[currNode].left = -1;
        tree[currNode].right = -1;
    }
    //tree finished
    vector<uint64_t> newData;
    vector<int> stack;
    vector<int> isLeftRightDoneStack;
    stack.push_back(root);
    isLeftRightDoneStack.push_back(0);
    int key;
    while(stack.size()){
        if(stack.back() == -1){
            stack.pop_back();
            isLeftRightDoneStack.pop_back();
        }else if(isLeftRightDoneStack.back() == 0){
            isLeftRightDoneStack[isLeftRightDoneStack.size()-1]++;
            isLeftRightDoneStack.push_back(0);
            stack.push_back(tree[stack.back()].left);
        }else if(isLeftRightDoneStack.back() == 1){
            if(stack.back() == 0){
                key = newData.size();
                newData.push_back((uint64_t)-1);
                cerr << "Key found" << endl;
            }else
                newData.push_back(data[stack.back()-1]);
            isLeftRightDoneStack[isLeftRightDoneStack.size()-1]++;
            isLeftRightDoneStack.push_back(0);
            stack.push_back(tree[stack.back()].right);
        }else{
            stack.pop_back();
            isLeftRightDoneStack.pop_back();
        }
    }
    cerr << "newData size: " << newData.size() << endl;
    cerr << "data size: " << data.size() << endl;
    data = newData;
    return key;
}
    /*
    map<uint64_t, vector<int> > firstColumn;
    for(int i = 0; i < data.size(); i++){
        firstColumn[data[i]].push_back(i);
    }
    cerr << "first column size: " << firstColumn.size() << endl;
    int counter = 0;
    for(map<uint64_t, vector<int> >::iterator i = firstColumn.begin(); i != firstColumn.end(); i++){
        if(counter%1000 == 999)
            cerr << "iteration: " << counter << endl;
        cerr << "before odd" << endl;
        sortBwt(i->second, 0, i->second.size(), data);
        k
        cerr << "odd" << endl;
        counter++;
    }
    //TODO: add everything to map
    //TODO: sorting
}*/
