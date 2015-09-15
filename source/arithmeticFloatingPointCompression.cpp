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
#include "arithmeticFloatingPointCompression.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstring>
#include <tardis.h>
#include "common.h"
#include "sortByRepeats.h"
#include <exception>

using namespace std;

void arithmeticCompression(double* values, uint64_t count, uint64_t arithmeticCodeSize, uint64_t* arithmeticCodes, uint8_t* arithmeticCodeLengths);
void arithmeticCompressionPart2(double* values, uint64_t count, uint64_t arithmeticCodeSize, uint64_t* arithmeticCodes, uint8_t* arithmeticCodeLengths, uint16_t* arithmeticCodesStream, uint8_t* &leftOverStream, uint64_t &leftOverSize);

compressedArithmeticFloatingPointPackage_t* compressArithmeticFloatingPoint(double* data, uint64_t count){
    compressedArithmeticFloatingPoint_t cfp;
    return compressArithmeticFloatingPoint(data, count, &cfp);}

void decompressArithmeticFloatingPoint(double* dest, compressedArithmeticFloatingPointPackage_t* source){
    compressedArithmeticFloatingPoint_t cfp;
    cfp.count = source->count;
    cfp.codes = (uint64_t*)((uint64_t)source + source->codes);
    cfp.codeLengths = (uint8_t*)((uint64_t)source + source->codeLengths);
    cfp.codeStream = (uint16_t*)((uint64_t)source + source->codeStream);
    cfp.leftOverStream = (void*)((uint64_t)source + source->leftOverStream);
    cfp.leftOverSize = source->size - source->leftOverStream;
    decompressArithmeticFloatingPoint(dest, &cfp);
}

struct Node{
    int64_t left;
    int64_t right;
    uint64_t count;
    uint64_t value;
    uint8_t depth;
    bool included;
};

void clearNode(Node &a){a.left = -1; a.right = -1; a.count = 0; a.included = 0;}
void printNode(Node a){
    cout << "Node: " << a.count << " " << a.left << " " << a.right << " " << a.included << endl;
}
void printTree(vector<Node> a, uint64_t count){
    uint64_t includedCount = 0;
    for(uint64_t i = 0; i < count; i++){
        if(a[i].included)
            includedCount++;
        printNode(a[i]);
    }
    cout << "Included: " << includedCount << endl;
}

struct HeapSortNode{
    int64_t left;
    int64_t right;
    uint64_t count;
    uint64_t original;
};

void copyHeapSortTree(HeapSortNode* heapSortTree, vector<Node> tree, uint64_t count){
    for(uint64_t i = 0; i < count; i++){
        heapSortTree[i].left = tree[i].left;
        heapSortTree[i].right = tree[i].right;
        heapSortTree[i].count = tree[i].count;
        heapSortTree[i].original = i;
    }
}

uint64_t removeLargestNode(HeapSortNode* heapSortTree, uint64_t root){
    int64_t tmp = -1;
    if(heapSortTree[root].left == -1 && heapSortTree[root].right == -1){
        return -1;
    }else if(heapSortTree[root].left == -1){
        return heapSortTree[root].right;
    }else if(heapSortTree[root].right == -1){
        return heapSortTree[root].left;
    }else if(heapSortTree[heapSortTree[root].left].count > heapSortTree[heapSortTree[root].right].count){
        tmp = heapSortTree[root].right;
        root = heapSortTree[root].left;
        uint64_t newLeftNode = removeLargestNode(heapSortTree, root);
        heapSortTree[root].right = tmp;
        heapSortTree[root].left = newLeftNode;
        return root;
    }else{
        tmp = heapSortTree[root].left;
        root = heapSortTree[root].right;
        uint64_t newRightNode = removeLargestNode(heapSortTree, root);
        heapSortTree[root].left = tmp;
        heapSortTree[root].right = newRightNode;
        return root;
    }
}
void findLargestCount(vector<Node> tree, uint64_t currentNode, int64_t &largestCount, int64_t &largestCountNode);
void assignCodes(const vector<Node> &tree, uint64_t currentNode, uint64_t depth, uint64_t decode, uint64_t* arithmeticCodes, uint8_t* arithmeticCodeLengths, uint64_t &codeCount, uint64_t &maxDepth);
struct pointerNode{
    uint64_t pointer;
    int64_t left;
    int64_t right;
    int64_t parent;
    uint64_t treeDepth;
};
int treeSize(const vector<pointerNode> tree, int node){
    if(node == -1)
        return 0;
    int ret = 1;

    if(tree[node].left != -1)
        ret += treeSize(tree, tree[node].left);
    if(tree[node].right != -1)
        ret += treeSize(tree, tree[node].right);
    return ret;
}
//TODO: check parent Connections
int checkParenthood(const vector<pointerNode> tree, int node){
    if(node == -1)
        return 0;
    int parent = tree[node].parent;
    int left = tree[node].left;
    int right = tree[node].right;
    if(left == -1){
    }else if(tree[left].parent == node){
        checkParenthood(tree,left);
    }else{
        cerr << "ERROR: This child has forgetten the face of his father" << endl;
        return 1;
    }
    if(right == -1){
    }else if(tree[right].parent == node){
        checkParenthood(tree,right);
    }else{
        cerr << "ERROR: Gettel" << endl;
        return 1;
    }
    return 0;
}
void printPointerNode(pointerNode a){
    cerr << "pointer: " << a.pointer << endl;
    cerr << "left: " << a.left << endl;
    cerr << "right: " << a.right << endl;
    cerr << "parent: " << a.parent << endl;
    cerr << "treeDepth: " << a.treeDepth << endl;
}
void printPointerTree(const vector<pointerNode> tree){
    for(int i=0; i<tree.size(); i++){
        printPointerNode(tree[i]);
    }
}
void updateNodeDepth(vector<pointerNode> &utilityTree, int ret){
    if(ret == -1){
    }else if(utilityTree[ret].left == -1 && utilityTree[ret].right == -1){
        utilityTree[ret].treeDepth = 1;
    }else if(utilityTree[ret].left == -1){
        utilityTree[ret].treeDepth = utilityTree[utilityTree[ret].right].treeDepth + 1;
    }else if(utilityTree[ret].right == -1){
        utilityTree[ret].treeDepth = utilityTree[utilityTree[ret].left].treeDepth + 1;
    }else if(utilityTree[utilityTree[ret].right].treeDepth > utilityTree[utilityTree[ret].left].treeDepth){
        utilityTree[ret].treeDepth = utilityTree[utilityTree[ret].right].treeDepth + 1;
    }else{
        utilityTree[ret].treeDepth = utilityTree[utilityTree[ret].left].treeDepth + 1;
    }
}
int insertIntoUtility(vector<pointerNode> &utilityTree, int utilityRoot, int &freePtr, const vector<Node> &tree, int treeIndex){
    if(freePtr < 0 || freePtr >= utilityTree.size())
        cerr << "ERROR: free pointer" << endl;
    if(utilityRoot == -1){
        int newNode = freePtr;
        freePtr = utilityTree[freePtr].right;
        //TODO: finish
        utilityTree[newNode].left = -1;
        utilityTree[newNode].right = -1;
        utilityTree[newNode].parent = -1;
        utilityTree[newNode].pointer = treeIndex;
        utilityTree[newNode].treeDepth = 1;
        return newNode;
    }if(tree[utilityTree[utilityRoot].pointer].count <= tree[treeIndex].count){
        int newNode = freePtr;
        freePtr = utilityTree[freePtr].right;
        utilityTree[newNode].pointer = treeIndex;
        int parent = utilityTree[utilityRoot].parent;
        utilityTree[newNode].parent = parent;
        if(parent == -1){
        }else if(utilityTree[parent].left == utilityRoot){
            utilityTree[parent].left = newNode;
        }else if(utilityTree[parent].right == utilityRoot){
            utilityTree[parent].right = newNode;
        }else{
            cerr << "ERROR: Lost Parent" << endl;
            cerr << "utilityRoot: " << utilityRoot << endl;
            cerr << "parent: " << parent << endl;
            cerr << "newNode: " << newNode << endl;
        }
        utilityTree[newNode].left = -1;
        utilityTree[newNode].right = utilityRoot;
        utilityTree[utilityRoot].parent = newNode;
        updateNodeDepth(utilityTree, newNode);
        //TODO: update tree depth
        return newNode;
    }else if(utilityTree[utilityRoot].left == -1){
        int newNode = freePtr;
        freePtr = utilityTree[freePtr].right;
        utilityTree[newNode].pointer = treeIndex;
        utilityTree[newNode].left = -1;
        utilityTree[newNode].right = -1;
        utilityTree[newNode].parent = utilityRoot;
        utilityTree[utilityRoot].left = newNode;
        utilityTree[newNode].treeDepth = 1;
        //TODO: update tree depth
        return utilityRoot;
    }else if(utilityTree[utilityRoot].right == -1){
        int newNode = freePtr;
        freePtr = utilityTree[freePtr].right;
        utilityTree[newNode].pointer = treeIndex;
        utilityTree[newNode].left = -1;
        utilityTree[newNode].right = -1;
        utilityTree[newNode].parent = utilityRoot;
        utilityTree[utilityRoot].right = newNode;
        utilityTree[newNode].treeDepth = 1;
        //TODO: update tree depth
        return utilityRoot;
    }else if(utilityTree[utilityTree[utilityRoot].left].treeDepth < utilityTree[utilityTree[utilityRoot].right].treeDepth){
        utilityTree[utilityRoot].left = insertIntoUtility(utilityTree, utilityTree[utilityRoot].left, freePtr, tree, treeIndex);
        updateNodeDepth(utilityTree,utilityRoot);
    }else{
        utilityTree[utilityRoot].right = insertIntoUtility(utilityTree, utilityTree[utilityRoot].right, freePtr, tree, treeIndex);
        updateNodeDepth(utilityTree,utilityRoot);
    }
    return utilityRoot;
}
int removeFromUtility(vector<pointerNode> &utilityTree, int nodeToRemove, const vector<Node> &tree){
    if(nodeToRemove < 0 || nodeToRemove >= utilityTree.size())
        cerr << "ERROR: " << endl;
    int parent = utilityTree[nodeToRemove].parent;
    int left = utilityTree[nodeToRemove].left;
    int right = utilityTree[nodeToRemove].right;
    int ret = -1;
    if(utilityTree[nodeToRemove].left == -1){
        if(utilityTree[nodeToRemove].right == -1){
            ret = -1;
        }else{
            utilityTree[utilityTree[nodeToRemove].right].parent = -1;
            ret = utilityTree[nodeToRemove].right;
        }
    }else if(utilityTree[nodeToRemove].right == -1){
        ret = utilityTree[nodeToRemove].left;
        utilityTree[utilityTree[nodeToRemove].left].parent = -1;
    }else if(tree[utilityTree[utilityTree[nodeToRemove].left].pointer].count > tree[utilityTree[utilityTree[nodeToRemove].right].pointer].count){
        ret = left;
        utilityTree[left].left = removeFromUtility(utilityTree, left, tree);
        if(utilityTree[left].left != -1)
            utilityTree[utilityTree[left].left].parent = left;
        utilityTree[left].right = right;
        if(utilityTree[left].right != -1)
            utilityTree[utilityTree[left].right].parent = left;
        utilityTree[left].parent = -1;
    }else{
        ret = right;
        utilityTree[right].right = removeFromUtility(utilityTree, right, tree);
        if(utilityTree[right].right != -1)
            utilityTree[utilityTree[right].right].parent = right;
        utilityTree[right].left = left;
        if(utilityTree[right].left != -1)
            utilityTree[utilityTree[right].left].parent = right;
        utilityTree[right].parent = -1;
    }
    if(ret == -1){
    }else if(utilityTree[ret].left == -1 && utilityTree[ret].right == -1){
        utilityTree[ret].treeDepth = 1;
    }else if(utilityTree[ret].left == -1){
        utilityTree[ret].treeDepth = utilityTree[utilityTree[ret].right].treeDepth + 1;
    }else if(utilityTree[ret].right == -1){
        utilityTree[ret].treeDepth = utilityTree[utilityTree[ret].left].treeDepth + 1;
    }else if(utilityTree[utilityTree[ret].right].treeDepth > utilityTree[utilityTree[ret].left].treeDepth){
        utilityTree[ret].treeDepth = utilityTree[utilityTree[ret].right].treeDepth + 1;
    }else{
        utilityTree[ret].treeDepth = utilityTree[utilityTree[ret].left].treeDepth + 1;
    }
    return ret;
}
void assignValueAndDepth(vector<Node> &tree, int root, int depth, uint64_t value){
    if(root == -1)
        return;
    assignValueAndDepth(tree, tree[root].left, depth+1, value);
    assignValueAndDepth(tree, tree[root].right, depth+1, value | ((uint64_t)1<<63>>(depth)));
    tree[root].value = value;
    tree[root].depth = depth;
}
void setIncluded(vector<pointerNode> &maxUtilityTree, int node, vector<Node> &tree){
    if(node == -1)
        return;
    tree[maxUtilityTree[node].pointer].included = true;
    setIncluded(maxUtilityTree, maxUtilityTree[node].left, tree);
    setIncluded(maxUtilityTree, maxUtilityTree[node].right, tree);
}
uint64_t reverse(uint64_t value){
    uint64_t retValue = 0;
    uint64_t mask = 1;
    for(int i = 63; i >0; i -= 2){
        retValue |= value>>i & mask;
        mask <<= 1;
    }
    for(int i = 1; i < 64; i += 2){
        retValue |= value<<i & mask;
        mask <<= 1;
    }
    //TODO: left right right left
    return retValue;
}
uint32_t reverse(uint32_t value){
    uint32_t retValue = 0;
    uint32_t mask = 1;
    for(int i = 31; i >0; i -= 2){
        retValue |= value>>i & mask;
        mask <<= 1;
    }
    for(int i = 1; i < 32; i += 2){
        retValue |= value<<i & mask;
        mask <<= 1;
    }
    //TODO: left right right left
    return retValue;
}
//TODO: complete
JilesPackage_t* jilesCoding(double* values, uint64_t count, uint64_t maxCodeLength){
    maxCodeLength = 1;
    while(1<<maxCodeLength < count/100)
        maxCodeLength++;
    cerr << "Max code length: " << maxCodeLength << endl;
    cerr << "staring jiles coding" << endl;
    JilesPackage_t* Bogy = (JilesPackage_t*)malloc(sizeof(JilesPackage_t));
    Bogy->size = sizeof(JilesPackage_t);
    stealTardis();
    uint64_t* lValues = (uint64_t*)values;
    vector<Node> tree;
    tree.resize(2*count);
    cerr << "time .5" << endl;
    markTime();
    clearNode(tree[0]);
    uint64_t freePtr = 1;
    uint64_t currentPtr;
    for(uint64_t i = 0; i < count; i++){
        currentPtr = 0;
        uint64_t currentValue = lValues[i];
        for(int j = 0; j < 64; j++){
            tree[currentPtr].count++;
            bool leftRight = (currentValue >> (63-j)) & 1;
            if(leftRight){
                if(tree[currentPtr].right == -1){
                    tree[currentPtr].right = freePtr;
                    clearNode(tree[tree[currentPtr].right]);
                    freePtr++;
                    if(freePtr >= tree.size()){
                        tree.resize(tree.size()*2);
                    }
                }
                currentPtr= tree[currentPtr].right;
            }else{
                if(tree[currentPtr].left == -1){
                    tree[currentPtr].left= freePtr;
                    clearNode(tree[tree[currentPtr].left]);
                    freePtr++;
                    if(freePtr >= tree.size()){
                        tree.resize(tree.size()*2);
                    }
                }
                currentPtr = tree[currentPtr].left;
            }
        }
    }
    if(freePtr >= tree.size())
        cerr << "ERROR: go to hell" << endl;
    markTime();
    assignValueAndDepth(tree, 0, 0, 0);
    uint64_t* repeats = (uint64_t*)malloc(count*sizeof(uint64_t));
    cerr << "before sort" << endl;
    uint64_t repeatCount = sortByRepeats((uint64_t*)values, count, repeats, count, 3);
    cerr << "after sort" << endl;
    if(false)
        repeatCount = 0;
    repeats = (uint64_t*)realloc(repeats, repeatCount*sizeof(uint64_t));
    map <uint64_t, bool> mapOfIncludedRepeats;
    for(int i = 0; i < repeatCount; i++){
        mapOfIncludedRepeats[repeats[i]] = false;
    }
    int unincludedRepeats = repeatCount;
    cerr << "repeats: " << unincludedRepeats << endl;

    uint64_t codesUsedAtCurrentLevel = 0;
    uint64_t currentLevel = 0;
    vector<pointerNode> maxUtilityTree;
    maxUtilityTree.resize(1<<(maxCodeLength));
    for(int i = 0; i < 1<<(maxCodeLength);i++){
        maxUtilityTree[i].right = i+1;
    }
    int utilityFreePtr = 1;
    maxUtilityTree[0].pointer = 0;
    maxUtilityTree[0].left = -1;
    maxUtilityTree[0].right = -1;
    maxUtilityTree[0].parent = -1;
    maxUtilityTree[0].treeDepth = 1;
    codesUsedAtCurrentLevel = 1;
    currentLevel = 0;
    int utilityRoot = 0;
    int iteration = 0;
    int utilitySize = 1;
    int maxedOutCreditCards = 0;
    vector<Investor> suits;
    int eatenCake = 0;
    uint64_t suitCount = 0;
    int allTheWayRepeats = 0;
    bool lastIteration = false;
    if(utilitySize == treeSize(maxUtilityTree, utilityRoot)){
        //cerr << "Check passed" << endl;
    }else
        cerr << "ERROR: size mismatch" << endl;
    while(true){
        if(utilityRoot == -1){
            cerr << "straw" << endl;
            break;
        }
        if((utilitySize+eatenCake) == 1<<currentLevel){
            cerr << "iteration " << iteration << endl;
            if(currentLevel == maxCodeLength-1){
                cerr << "last it" << endl;
                //TODO: add together
                lastIteration = true;
                int tmp = utilitySize+(2*eatenCake)+unincludedRepeats;
                if(tmp >= 1<<(maxCodeLength)){
                    currentLevel++;
                    eatenCake *= 2;
                    cerr << "Early Break" << endl;
                    allTheWayRepeats = (1<<maxCodeLength) - utilitySize - eatenCake;
                    break;
                }
            }else if(currentLevel == maxCodeLength){
                cerr << "natural break" << endl;
                break;
            }
            currentLevel+=1;
            eatenCake *=2;
        }
        if(lastIteration){
            int tmp = utilitySize+eatenCake+unincludedRepeats;
            if(tmp >= 1<<(maxCodeLength)){
                cerr << "Last iteration break" << endl;
                allTheWayRepeats = (1<<maxCodeLength) - utilitySize - eatenCake;
                break;
            }
        }
        if((utilitySize+eatenCake) >= 1<<currentLevel){
            cerr << "ERROR: How does this happen?" << endl;
        }
        int left = tree[maxUtilityTree[utilityRoot].pointer].left;
        int right = tree[maxUtilityTree[utilityRoot].pointer].right;
        if(left == -1 && right == -1){
            eatenCake++;
            Investor pieceOfCake;
            pieceOfCake.value = tree[maxUtilityTree[utilityRoot].pointer].value; //TODO: fix
            pieceOfCake.depth = tree[maxUtilityTree[utilityRoot].pointer].depth; //TODO: fix
            pieceOfCake.level = currentLevel;
            pieceOfCake.suitNumber = suitCount++;
            suits.push_back(pieceOfCake);

            if(mapOfIncludedRepeats.count(pieceOfCake.value)){
                mapOfIncludedRepeats[pieceOfCake.value] = true;
                unincludedRepeats--;
            }
        }
        int oldRoot = utilityRoot;
        utilityRoot = removeFromUtility(maxUtilityTree, utilityRoot, tree);
        utilitySize--;
        maxUtilityTree[oldRoot].right = utilityFreePtr;
        utilityFreePtr = oldRoot;
        if(right != -1){
            utilityRoot = insertIntoUtility(maxUtilityTree, utilityRoot, utilityFreePtr, tree, right);
            utilitySize++;
        }
        if(left != -1){
            utilityRoot = insertIntoUtility(maxUtilityTree, utilityRoot, utilityFreePtr, tree, left);
            utilitySize++;
        }
        iteration++;
    }
    markTime();
    if(currentLevel != maxCodeLength)
        cerr << "ERROR: level doesn't match length" << endl;
    setIncluded(maxUtilityTree, utilityRoot, tree);
    cerr << "Last iteration: " << iteration << endl;
    cerr << "eaten cake: " << eatenCake << endl;
    cerr << "utility size: " << utilitySize << endl;
    cerr << "size: " << treeSize(maxUtilityTree, utilityRoot) << endl;
    for(int i = 0; i < utilitySize; i++){
        Investor pieceOfCake;
        pieceOfCake.value = tree[maxUtilityTree[utilityRoot].pointer].value; //TODO: fix
        pieceOfCake.depth = tree[maxUtilityTree[utilityRoot].pointer].depth; //TODO: fix
        pieceOfCake.level = currentLevel;
        pieceOfCake.suitNumber = suitCount++;
        suits.push_back(pieceOfCake);
        utilityRoot = removeFromUtility(maxUtilityTree, utilityRoot, tree);
    }
    //all the way repeats
    while(unincludedRepeats > allTheWayRepeats){
        allTheWayRepeats *= 2;
        currentLevel++;
    }
    int j = 0;
    int beforeRepeats = suits.size();
    for(int i = 0; i < unincludedRepeats; i++){
        while(mapOfIncludedRepeats[repeats[j]]){
            j++;
            if(!mapOfIncludedRepeats.count(repeats[j]))
                cerr << "ERROR: repeat value not mapped " << endl;
        }
        Investor pieceOfCake;
        pieceOfCake.value = repeats[j];
        pieceOfCake.depth = 64;
        pieceOfCake.level = currentLevel;
        pieceOfCake.suitNumber = suitCount++;
        suits.push_back(pieceOfCake);
        j++;
    }
    cerr << "unincluded repeats: " << unincludedRepeats << endl;
    cerr << "before repeats: " << beforeRepeats << endl;
    cerr << "suit count: " << suits.size() << endl;
    map<uint64_t, int> codeMap;
    map<uint64_t, int> codeMap64;
    for(int i = 0; i < suits.size(); i++){
        //cerr << hex << suits[i].value << " " << dec << suits[i].depth << " " << suits[i].level << endl;
        //cerr << *(double*)&(suits[i].value) << endl;
        if(suits[i].depth == 64){
            codeMap64[suits[i].value] = i;
        }else{
            codeMap[suits[i].value] = i;
        }
        //cerr << "suit value: " << hex << suits[i].value << endl;
    }
    cerr << "code map size: " << (codeMap.size() + codeMap64.size()) << endl;

    uint64_t* codes = (uint64_t*)malloc(128);
    int codeSize = 128;
    uint64_t* leftOver = (uint64_t*)malloc(128);
    int leftOverSize = 128;
    for(int i = 0; i < 128/8; i++)
        codes[i] = 0;
    for(int i = 0; i < leftOverSize/8; i++)
        leftOver[i] = 0;
    int currentCodeBit = 0;
    int currentLeftOverBit = 0;

    //TODO: reverse
    suitCount = 0;
    //cerr << "recallabrating suits" << endl;
    for(int i = 0; i < suits.size(); i++){
        suits[i].suitNumber =
        reverse(suitCount);
        suitCount += 1L<<(64-suits[i].level);
        //cerr << "suit: " << suits[i].suitNumber << endl;
        //cerr << "level: " << suits[i].level << endl;
    }
    if(unincludedRepeats != 0){
        cerr << "first repeats:" << endl;
        cerr << "value: " << hex << suits[beforeRepeats+1].value << endl;
        cerr << "depth: " << dec << suits[beforeRepeats+1].depth << endl;
        cerr << "level: " << dec << suits[beforeRepeats+1].level << endl;
        cerr << "suit number: " << hex <<
        reverse(suits[beforeRepeats+1].suitNumber) <<
        endl;
    }
    if(unincludedRepeats != 0){
        cerr << "first repeats:" << endl;
        cerr << "value: " << hex << suits[beforeRepeats].value << endl;
        cerr << "depth: " << dec << suits[beforeRepeats].depth << endl;
        cerr << "level: " << dec << suits[beforeRepeats].level << endl;
        cerr << "suit number: " << hex <<
        reverse(suits[beforeRepeats].suitNumber) <<
        endl;
    }
    cerr << "time 4" << endl;
    markTime();
    for(int i = 0; i < count; i++){
        //find depth
        int depth = 0;
        int currentNode = 0;
        uint64_t value = ((uint64_t*)values)[i];
        if(codeMap64.count(value)){
            depth = 64;
        }else{
            while(true){
                if(currentNode == -1){
                    cerr << "FUCK" << endl;
                    break;
                }
                if(depth == 64)
                    break;
                if(tree[currentNode].included){
                    break;
                }
                bool isRight = (value>>(63-depth)) & 1L;
                if(isRight){
                    currentNode = tree[currentNode].right;
                }else{
                    currentNode = tree[currentNode].left;
                }
                depth++;    
            }
        }     
        uint64_t lookup;
        if(depth == 64)
            lookup = value;
        else
            lookup = value & ((uint64_t)-1 << (64-depth));
        //find code
        //cerr << "lookup: " << lookup << endl;
        //cerr << "depth: " << depth << endl;
        map<uint64_t, int>::iterator it;
        if(depth == 64){
            it = codeMap64.find(lookup);
        }else{
            it = codeMap.find(lookup);
        }
        Investor currentInvestor = suits[(*it).second];
        uint64_t code = suits[(*it).second].suitNumber;
        //cerr << "suitNumber: " << code << endl;
        uint64_t codeWidth = currentInvestor.level;
        //cerr << "codeWidth: " << codeWidth << endl;
        uint64_t leftOverWidth = 64 - currentInvestor.depth;
        uint64_t leftOverPart = 0;
        //cerr << "leftOverWidth: " << dec << leftOverWidth << endl;
        if(leftOverWidth != 64)
            leftOverPart = value & (((uint64_t)-1)>>currentInvestor.depth);
        //cerr << "leftOverPart: " << hex << leftOverPart << endl;
        if(leftOverWidth == 0){
        }else if(currentLeftOverBit/64 == (currentLeftOverBit+leftOverWidth-1)/64){
            leftOver[currentLeftOverBit/64] |=
            leftOverPart<<(currentLeftOverBit%64);
        }else{
            leftOver[currentLeftOverBit/64] |=
            leftOverPart<<(currentLeftOverBit%64);
            leftOver[currentLeftOverBit/64+1] |=
            leftOverPart>>(64-currentLeftOverBit%64);
        }
        currentLeftOverBit += leftOverWidth;
        if(currentLeftOverBit +64 >  leftOverSize*8){
            leftOver = (uint64_t*)reallocAndZero(leftOver, leftOverSize, leftOverSize*2);
            leftOverSize *=2;
        }
        if(currentCodeBit/64 == (currentCodeBit+codeWidth-1)/64){
            //cerr << "code: " << code << endl;
            codes[currentCodeBit/64] |= code << (currentCodeBit%64);
        }else{
            codes[currentCodeBit/64] |= code << (currentCodeBit%64);
            codes[currentCodeBit/64+1] |= code >> (64 - currentCodeBit%64);
        }
        currentCodeBit += codeWidth;
        if(currentCodeBit + 64 > codeSize*8){
            codes = (uint64_t*)reallocAndZero(codes, codeSize, codeSize*2);
            codeSize *= 2;
        }
    }
    cerr << "time 5" << endl;
    markTime();
    //TODO: adjust suits
    leftOverSize = currentLeftOverBit/64*8+8;
    leftOver = (uint64_t*)realloc(leftOver, currentLeftOverBit/64*8+8);
    if(leftOver == NULL){
        cerr << "ERROR: cant allocate" << endl;
    }
    codeSize = currentCodeBit/64*8+8;
    codes = (uint64_t*)realloc(codes, currentCodeBit/64*8+8);
    if(codes == NULL){
        cerr << "ERROR: cant allocate" << endl;
    }
    //TODO: add place for repeats
    JilesPackage_t* retPackage =
    (JilesPackage_t*)malloc(sizeof(JilesPackage_t) +
    beforeRepeats*sizeof(Investor) + unincludedRepeats*sizeof(uint64_t) 
    + leftOverSize + codeSize);
    int packageSize =sizeof(JilesPackage_t) +
    beforeRepeats*sizeof(Investor) + unincludedRepeats*sizeof(uint64_t) 
    + leftOverSize + codeSize;
    retPackage->count = count;
    int currentLocation = sizeof(JilesPackage_t);
    retPackage->suits = currentLocation;
    currentLocation += beforeRepeats*sizeof(Investor);
    retPackage->finalRepeats = currentLocation;
    currentLocation += unincludedRepeats*sizeof(uint64_t);
    retPackage->codeStream = currentLocation;
    currentLocation += codeSize;
    retPackage->leftOverStream = currentLocation;
    currentLocation += leftOverSize;
    retPackage->size = currentLocation;
    Investor* packagedSuits = (Investor*)((uint8_t*)retPackage +
    retPackage->suits);
    for(int i = 0; i < beforeRepeats; i++){
        *packagedSuits = suits[i];
        packagedSuits++;
    }
    cerr << "HERE" << endl;
    for(int i = 0; i < unincludedRepeats; i++){
        ((uint64_t*)(((uint8_t*)retPackage) + retPackage->finalRepeats))[i] = suits[beforeRepeats + i].value;
    }
    cerr << "HERE 2" << endl;
    retPackage->lastLevel = currentLevel;

    memcpy((uint8_t*)retPackage + retPackage->codeStream, codes, codeSize);
    memcpy((uint8_t*)retPackage + retPackage->leftOverStream, leftOver,
    leftOverSize);
    //TODO: move codes
    //TODO: move left overs
    //TODO: return package
    //cerr << "returning" << endl;
    returnTardis();
    cerr << "fin" << endl;
    return retPackage;
}
void arithmeticCompression(double* values, uint64_t count, uint64_t arithmeticCodeSize, uint64_t* arithmeticCodes, uint8_t* arithmeticCodeLengths){
    /*uint64_t arithmeticCodeSize = 5;
    uint64_t *arithmeticCodes = new uint64_t[1<<arithmeticCodeSize];
    uint8_t* arithmeticCodeLengths = new uint8_t[1<<arithmeticCodeSize];*/
    uint64_t codeCount = 1<<arithmeticCodeSize;
    uint64_t* lValues = (uint64_t*)values;
    vector<Node> tree;
    uint64_t treeSize = 128;
    tree.resize(128);
    clearNode(tree[0]); 
    uint64_t freePtr = 1;
    uint64_t currentPtr;
    markTime();
    for(uint64_t i = 0; i < count; i++){
        currentPtr = 0;
        uint64_t currentValue = lValues[i];
        for(int j = 0; j < 64; j++){
            tree[currentPtr].count++;
            bool leftRight = (currentValue >> (63-j)) & 1;
            if(leftRight){
                if(tree[currentPtr].right == -1){
                    tree[currentPtr].right = freePtr;
                    clearNode(tree[tree[currentPtr].right]);
                    freePtr++;
                    if(freePtr >= treeSize){
                        treeSize = treeSize*2;
                        tree.resize(treeSize);
                    }
                }
                currentPtr= tree[currentPtr].right;
            }else{
                if(tree[currentPtr].left == -1){
                    tree[currentPtr].left= freePtr;
                    clearNode(tree[tree[currentPtr].left]);
                    freePtr++;
                    if(freePtr >= treeSize){
                        treeSize = treeSize*2;
                        tree.resize(treeSize);
                    }
                }
                currentPtr = tree[currentPtr].left;
            }
        }
    }
    markTime();
    //TODO: heap sort
    HeapSortNode* heapSortTree = new HeapSortNode[freePtr];
    copyHeapSortTree(heapSortTree, tree, freePtr);
    uint64_t root = 0;
    for(uint64_t i = 0; i < (1<<arithmeticCodeSize)-1;i++){
        tree[heapSortTree[root].original].included = true;
        uint64_t original = heapSortTree[root].original;
        if(tree[original].left == -1 && tree[original].right == -1){
        }else if(tree[original].left == -1 || tree[original].right == -1){
            i--;
        }
        root = removeLargestNode(heapSortTree, root);
        if(root == -1)
            break;
    }
    if(false)
    for(uint64_t i = 0; i < (1<<arithmeticCodeSize)-1; i++){
        int64_t largestCount = -1;
        int64_t largestCountNode = -1;
        uint64_t currentNode = 0;
        findLargestCount(tree, currentNode, largestCount, largestCountNode);
        if(largestCountNode != -1)
            tree[largestCountNode].included = true;
    }
    markTime();
    uint64_t tmp=0;
    uint64_t maxDepth =0;
    assignCodes(tree, 0, 0, 0, arithmeticCodes, arithmeticCodeLengths, tmp, maxDepth);
    cerr << "Max Depth: " << dec << maxDepth << endl;
    if(tmp != codeCount){
        while(tmp < codeCount){
            arithmeticCodes[tmp] = (uint64_t)-1;
            arithmeticCodeLengths[tmp] = 0;
            tmp++;
        }
    }
    markTime();
    return;
}

void findLargestCount(vector<Node> tree, uint64_t currentNode, int64_t &largestCount, int64_t &largestCountNode){
    if(tree[currentNode].included){
        if(tree[currentNode].left != -1)
            findLargestCount(tree, tree[currentNode].left, largestCount, largestCountNode);
        if(tree[currentNode].right != -1)
            findLargestCount(tree, tree[currentNode].right, largestCount, largestCountNode);
    }else if(int64_t(tree[currentNode].count) > largestCount){
        largestCount = tree[currentNode].count;
        largestCountNode = currentNode;
    }
}

void assignCodes(const vector<Node> &tree, uint64_t currentNode, uint64_t depth, uint64_t decode, uint64_t* arithmeticCodes, uint8_t* arithmeticCodeLengths, uint64_t &codeCount, uint64_t &maxDepth){
    if(maxDepth < depth)
        maxDepth = depth;
    //cerr << "assignCodes" << endl;
    if(tree[currentNode].included && depth != 64){
        if(tree[currentNode].left == -1){
            //arithmeticCodes[codeCount] = decode;
            //arithmeticCodeLengths[codeCount] = depth+1;
            //codeCount += 1;
            //cerr << "left empty" << endl;
            //cerr << hex << decode << endl;
            //cerr << dec << depth << endl;
        }else{
            assignCodes(tree, tree[currentNode].left, depth+1, decode, arithmeticCodes, arithmeticCodeLengths, codeCount, maxDepth);
        }
        if(tree[currentNode].right == -1){
            //arithmeticCodes[codeCount] = decode | ((uint64_t)1<<63>>(depth));
            //arithmeticCodeLengths[codeCount] = depth+1;
            //codeCount += 1;
            //cerr << "right empty" << endl;
            //cerr << hex << decode << endl;
            //cerr << dec << depth << endl;
        }else{
            assignCodes(tree, tree[currentNode].right, depth+1, decode| ((uint64_t)1<<63>>(depth)), arithmeticCodes, arithmeticCodeLengths, codeCount, maxDepth);
        }
    }else{
        arithmeticCodes[codeCount] = decode;
        arithmeticCodeLengths[codeCount] = depth;
        codeCount++;
    }
}
void arithmeticCompressionPart2(double* values, uint64_t count, uint64_t arithmeticCodeSize, uint64_t* arithmeticCodes, uint8_t* arithmeticCodeLengths, uint16_t* arithmeticCodesStream, uint8_t* &leftOverStream, uint64_t &leftOverSize){
    uint64_t currentLeftOverBit = 0;
    uint64_t* lValues = (uint64_t*)values;
    uint64_t* leftOverStreamLong = (uint64_t*)leftOverStream;
    for(uint64_t i = 0; i < count; i++){
        int code = 1<<arithmeticCodeSize;
        for(uint64_t j = 0; j < 1<<arithmeticCodeSize; j++){
            if(lValues[i] < arithmeticCodes[j]){
                code = j;
                break;
            }
        }
        code--;
        arithmeticCodesStream[i] = code;
        uint64_t currentValue = lValues[i];
        uint64_t leftOverLength = 64 - arithmeticCodeLengths[code];
        if(leftOverLength != 0)
            currentValue &= ((uint64_t)-1)>>(64-leftOverLength);
        if(leftOverLength==0){
        }else if((currentLeftOverBit / 64) == ((currentLeftOverBit+leftOverLength-1)/64)){
            leftOverStreamLong[currentLeftOverBit/64] |= (lValues[i]<<(currentLeftOverBit%64)) & (uint64_t(-1)>>(64-leftOverLength-(currentLeftOverBit%64)));
        }else{
            leftOverStreamLong[currentLeftOverBit/64] |= (lValues[i]<<(currentLeftOverBit%64));
            leftOverStreamLong[currentLeftOverBit/64+1] |= (currentValue>>(64-currentLeftOverBit%64));
        }
        currentLeftOverBit += leftOverLength;
        if(currentLeftOverBit/8+9 >= leftOverSize){
            leftOverStreamLong = (uint64_t*)reallocAndZero(leftOverStreamLong, leftOverSize, 2*leftOverSize);
            leftOverSize *=2;
        }
    }
    leftOverStreamLong = (uint64_t*)realloc(leftOverStreamLong, currentLeftOverBit/64*8 + 8);
    leftOverSize = currentLeftOverBit/64*8+8;
    leftOverStream = (uint8_t*)leftOverStreamLong;
}
compressedArithmeticFloatingPointPackage_t* compressArithmeticFloatingPoint(double* data, uint64_t count, compressedArithmeticFloatingPoint_t* cfp){
    uint64_t arithmeticCodeSize = 12;
    uint64_t codeCount = 1<<arithmeticCodeSize;
    uint64_t *arithmeticCodes = new uint64_t[1<<arithmeticCodeSize];
    uint8_t* arithmeticCodeLengths = new uint8_t[1<<arithmeticCodeSize];
    stealTardis();
    arithmeticCompression(data, count, arithmeticCodeSize, arithmeticCodes, arithmeticCodeLengths);
    markTime();
    //TODO: Jiles compression
    /*cerr << "Codes: " << endl;
    for(uint64_t i=0; i < codeCount; i++){
        cerr << dec << (int)arithmeticCodeLengths[i] << " " << hex << arithmeticCodes[i] << endl;
    }*/
    uint8_t* leftOverStream = (uint8_t*)malloc(128);
    for(int i=0; i< 128;i++)
        leftOverStream[i] = 0;
    uint16_t* arithmeticCodesStream = new uint16_t[count];
    uint64_t tmp = 128;
    arithmeticCompressionPart2(data, count, arithmeticCodeSize, arithmeticCodes, arithmeticCodeLengths, arithmeticCodesStream, leftOverStream, tmp);
    double firstValue;
    *(uint64_t*)&firstValue = arithmeticCodes[arithmeticCodesStream[0]] | (((uint64_t*)leftOverStream)[0] & ((uint64_t)-1>>(arithmeticCodeLengths[0])));
    double one = 1;
    //TODO: reduce leftOvers
    uint64_t size = sizeof(compressedArithmeticFloatingPointPackage_t) + (1<<arithmeticCodeSize)*sizeof(uint64_t) + (1<<arithmeticCodeSize)*sizeof(uint8_t) + tmp + count*sizeof(uint16_t);
    compressedArithmeticFloatingPointPackage_t* package = (compressedArithmeticFloatingPointPackage_t*)malloc(size);
    package->count = count;
    package->size = size;
    uint64_t pointer = sizeof(compressedArithmeticFloatingPointPackage_t);
    package->codes = pointer;
    pointer += codeCount * sizeof(uint64_t);
    package->codeLengths = pointer;
    pointer += codeCount * sizeof(uint8_t);
    package->codeStream = pointer;
    pointer += count * sizeof(uint16_t);
    package->leftOverStream = pointer;
    pointer += tmp;
    cfp->count = count;
    cfp->codes = (uint64_t*)((uint8_t*)package + package->codes);
    cfp->codeLengths = (uint8_t*)package + package->codeLengths;
    cfp->codeStream = (uint16_t*)((uint8_t*)package + package->codeStream);
    cfp->leftOverStream = (void*)((uint8_t*)package + package->leftOverStream);
    memcpy(cfp->codes, arithmeticCodes, codeCount*sizeof(uint64_t));
    memcpy(cfp->codeLengths, arithmeticCodeLengths, codeCount*sizeof(uint8_t));
    memcpy(cfp->codeStream, arithmeticCodesStream, count*sizeof(uint16_t));
    memcpy(cfp->leftOverStream, leftOverStream, tmp);
    cfp->leftOverSize = tmp;
    returnTardis();
    return package;
}
void decompressArithmeticFloatingPoint(double* dest, compressedArithmeticFloatingPoint_t* source){
    uint64_t count = source->count;
    uint64_t* codes = source->codes;
    uint8_t* codeLengths = source->codeLengths;
    uint16_t* codeStream = source->codeStream;
    uint64_t* leftOverStream = (uint64_t*)source->leftOverStream;
    uint64_t leftOverSize = source->leftOverSize;

    uint64_t leftOverP = 0;
    uint64_t leftOverPart = 0;
    for(uint64_t i = 0; i < count; i++){
        //cerr << "code number: " << codeStream[i] << endl;
        uint8_t length = codeLengths[codeStream[i]];
        uint8_t leftOverLength = 64-length;
        uint64_t code = codes[codeStream[i]];
        if(leftOverP/64 == (leftOverP+leftOverLength-1)/64){
            leftOverPart = (leftOverStream[leftOverP/64]>>(leftOverP%64) );// & (((uint64_t)-1)>>(leftOverLength));
        }else{
            leftOverPart = (leftOverStream[leftOverP/64]>>(leftOverP%64) );
            leftOverPart |= (leftOverStream[leftOverP/64+1]<<(64-leftOverP%64));// & /
        }
        if(leftOverLength == 0)
            leftOverPart = 0;
        else
            leftOverPart &=(((uint64_t)-1)>>(64-leftOverLength));
        leftOverPart |= code;
        ((uint64_t*)dest)[i] = leftOverPart;
        leftOverP += leftOverLength;
    }
}
JilesPackage_t* compressJilesFloatingPoint(double* data, uint64_t count){
    uint64_t maxEncodingLength = 16;
    JilesPackage_t* retPackage;
    retPackage = jilesCoding(data, count, maxEncodingLength);
    return retPackage;
}

double*
decompressJilesFloatingPoint(JilesPackage_t*
package){
    cerr << "HERE: a" << endl;
    //TODO: print RAW
    double* retValues = (double*)malloc(package->count * sizeof(double));
    //cerr << setw(2) << setfill('0');
    uint8_t* packageBytePointer = (uint8_t*)package;
    for(int i = package->size - 1; i >= 0; i--){
        //cerr << (int)packageBytePointer[i];
    }
    //cerr << endl;
    Investor* suitsPtr = (Investor*)(packageBytePointer + package->suits);
    int suitsSize = (package->finalRepeats - package->suits)/sizeof(Investor);
    vector<Investor> suits;
    for(int i = 0; i < suitsSize; i++){
        suits.push_back(suitsPtr[i]);
    }
    cerr << "HERE: b" << endl;
    int repeatsSize = (package->codeStream - package->finalRepeats)/sizeof(uint64_t);
    uint64_t* repeatPtr = (uint64_t*)(packageBytePointer + package->finalRepeats);
    int lastSuit = suits.size()-1;
    uint32_t suitCount = reverse(suits[suits.size()-1].suitNumber) + (1<<(32-suits[lastSuit].level));
    cerr << "HERE: c" << endl;
    for(int i = 0; i < repeatsSize; i++){
        Investor repeatInvestor;
        repeatInvestor.value = repeatPtr[i];
        repeatInvestor.depth = 64;
        repeatInvestor.level = package->lastLevel;
        repeatInvestor.suitNumber = reverse(suitCount);
        suitCount += 1<<(32-package->lastLevel);
        suits.push_back(repeatInvestor);
    }
    cerr << "first repeated suit: " << endl;
    if(lastSuit+1 != suits.size()){
        cerr << "value: " << hex << suits[lastSuit+1].value << endl;
        cerr << "depth: " << dec << suits[lastSuit+1].depth << endl;
        cerr << "level: " << dec << suits[lastSuit+1].level << endl;
        cerr << "suitnumber: " << hex <<
        reverse(suits[lastSuit+1].suitNumber) << endl;
    }

    cerr << "HERE: after FUCK" << endl;
    uint32_t currentCode;
    map<uint32_t, int> codeMap;
    for(int i = 0; i < suits.size(); i++){
        codeMap[reverse(suits[i].suitNumber)] = i;
    }
    for(map<uint32_t, int>::iterator it = codeMap.begin(); it != codeMap.end(); it++){
        //cerr << "codeMap: " << it->first << ", " << it->second << endl;
    }
    uint64_t currentCodeBit = 0;
    uint64_t currentLeftOverBit = 0;
    uint32_t* codes = (uint32_t*)(packageBytePointer + package->codeStream);
    uint64_t* leftOvers = (uint64_t*)(packageBytePointer + package->leftOverStream);
    //cerr << "leftOvers: " << hex << leftOvers[0] << endl;
    cerr << "HERE: decoding: " << endl;
    cerr << "count : " << package->count << endl;
    for(int i = 0; i < package->count; i++){
        if(i == 0)
            cerr << "one step a" << endl;
        if(currentCodeBit%32 == 0)
            currentCode = codes[currentCodeBit/32];
        else
            currentCode = codes[currentCodeBit/32]>>(currentCodeBit%32) | codes[currentCodeBit/32 + 1]<<(32-currentCodeBit%32);
        if(i == 0)
            cerr << "one step b" << endl;
        //cerr << "current code: " << currentCode << endl;
        currentCode = reverse(currentCode);
        //cerr << "current code: " << currentCode << endl;
        map<uint32_t, int>::iterator investorIndex = codeMap.upper_bound(currentCode);
        investorIndex--; //WOOT STL hacking
        //cerr << "Investor index: " << investorIndex->first << ", " << investorIndex->second << endl;
        Investor currentInvestor = suits[investorIndex->second];
        //cerr << "decoded value: " << hex << currentInvestor.value << endl;
        //cerr << "double: " << *(double*)&(currentInvestor.value) << endl;
        currentCodeBit += currentInvestor.level;
        if(currentCodeBit/8 > (package->leftOverStream - package->codeStream)){
            cerr << "ERROR: code going out of bounds at " << i << endl;
            break;
        }
        if(i == 0)
            cerr << "one step c" << endl;
        uint64_t value = currentInvestor.value;
        int depth = currentInvestor.depth;
        uint64_t leftOverValue = 0;
        if(i == 0)
            cerr << "one step c1" << endl;
        if(currentLeftOverBit%64 == 0){
            if(i == 0)
                cerr << "option a" << endl;
            leftOverValue = leftOvers[currentLeftOverBit/64];
        }else{
            if(i == 0)
                cerr << "option b" << endl;
            leftOverValue =
            leftOvers[currentLeftOverBit/64]>>(currentLeftOverBit%64) |
            leftOvers[currentLeftOverBit/64 + 1]<<(64-currentLeftOverBit%64);
        }
        if(i == 0)
            cerr << "one step d" << endl;
        currentLeftOverBit += 64-depth;
        if(currentLeftOverBit/8 > (package->size - package->leftOverStream)){
            cerr << "error: left overs going out of bounds at " << i << endl;
            break;
        }
        if(depth != 64)
            leftOverValue &= ((uint64_t)-1)>>(depth);
        else
            leftOverValue = 0;
        if(i == 0){
            cerr << "ERROR: mismatches: " << endl;
            cerr << hex << "code: " << currentCode << endl;
            cerr << dec << "depth: " << depth << endl;
            cerr << hex << "decoded: " << value << endl;
            cerr << "leftOver: " << leftOverValue << endl;
            cerr << dec << "code width: " << currentInvestor.level << endl;
            cerr << hex << "reversed suit num: " << reverse(currentInvestor.suitNumber) << endl;

        }
        value |= leftOverValue;
        if(i == 0)
            cerr << "one step e" << endl;
        ((uint64_t*)retValues)[i] = value;
    }
    return retValues;
}
