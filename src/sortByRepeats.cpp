#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <stdint.h>
#include <map>

using namespace std;

typedef struct{
    double value;
    int left;
    int right;
    int parent;
    int count;
    int smaller;
    int larger;
    int treeSize;
}valueNode;

struct frequency{
    double value;
    uint64_t count;
};

void clearValueNode(valueNode* node){
    node->value = 0;
    node->left = -1;
    node->right = -1;
    node->parent = -1;
    node->count = 0;
    node->smaller = -1;
    node->larger = -1;
    node->treeSize = 0;
    return;
}

bool countCompare(valueNode i, valueNode j){return (i.count > j.count);}

void setLeftChild(int parent, int child, valueNode* tree){
    tree[parent].left = child;
    if(child != -1)
        tree[child].parent = parent;
}
void setRightChild(int parent, int child, valueNode* tree){
    tree[parent].right = child;
    if(child != -1)
        tree[child].parent = parent;
}

bool cmpSubTree(int a, int b, valueNode* tree){
    if(b == -1)
        return false;
    if(a == -1)
        return true;
    if(tree[a].treeSize < tree[b].treeSize)
        return true;
}
/*
int cmpSubTree(int a, int b, valueNode* tree){
    if(b == -1)
        if(a== -1)
            return 0;
        else if(tree[a].treeCount > 2)
            return 2;
        else
            return 1;
    if(a == -1)
        if(tree[b].treeCount > 2)
            return -2;
        else
            return -1;
    if(tree[a].treeSize == tree[b].treeSize)
        return 0;
    else{}
}*/

void dfsMaxCount(valueNode* tree, uint64_t currNode, uint64_t* counts, double* values, uint64_t commonSize);

bool elementCountCmp(pair<uint64_t, int> a, pair<uint64_t, int> b){return a.second > b.second;}
uint64_t sortByRepeats(uint64_t* values, uint64_t nnz, uint64_t* common, uint64_t commonSize, int minRepeat = 2){
    vector<uint64_t> elements(values, values + nnz);
    sort(elements.begin(), elements.end());
    vector<pair<uint64_t,int> > elementCounts;
    elementCounts.push_back(make_pair(elements[0],1));
    int j = 0;
    for(int i = 1; i < elements.size(); i++){
        if(elements[i-1] == elements[i]){
            elementCounts.back().second++;
        }else{
            j++;
            elementCounts.push_back(make_pair(elements[i],1));
        }
    }
    sort(elementCounts.begin(), elementCounts.end(), elementCountCmp);
    if(commonSize > elementCounts.size())
        commonSize = elementCounts.size();
    while(commonSize > 0){
        if(elementCounts[commonSize-1].second > 1)
            break;
        else
            commonSize--;
    }
    for(int i = 0; i < commonSize; i++){
        common[i] = elementCounts[i].first;
    }
    return commonSize;
    //TODO: use stl
}
    /*
    valueNode* tree = (valueNode*)malloc(nnz * sizeof(valueNode));
    clearValueNode(tree);
    int root = 0;
    int currNode = 0;
    int nodeCount = 1;
    tree[0].value = values[0];
    tree[0].count = 1;
    tree[0].treeSize = 1;
    for(int i = 1; i < nnz; i++){
        currNode = root;
        //TODO: balance
        if(tree[currNode].left == -1){
            if(tree[currNode].right == -1){
            }else if(tree[tree[currNode].right].treeSize > 3){
                root = tree[currNode].right;
                if(cmpSubTree(tree[tree[currNode].right].left, tree[tree[currNode].right].right, tree)){
                    setRightChild(currNode, tree[tree[currNode].right].left, tree);
                    setLeftChild(root, currNode, tree);
                }else{
                    setRightChild(currNode, tree[tree[currNode].right].right, tree);
                    setRightChild(root, currNode, tree);
                }
                currNode = root;

            }
        }else if(tree[currNode].right == -1){
            if(tree[tree[currNode].left].treeSize > 3){
                root = tree[currNode].left;
                if(cmpSubTree(tree[tree[currNode].left].left, tree[tree[currNode].left].right, tree)){
                    setLeftChild(currNode, tree[tree[currNode].left].left, tree);
                    setLeftChild(root, currNode, tree);
                }else{
                    setLeftChild(currNode, tree[tree[currNode].left].right, tree);
                    setRightChild(root, currNode, tree);
                }
                currNode = root;
            }
        }else if(tree[tree[currNode].right].treeSize > 3 * tree[tree[currNode].left].treeSize){
            root = tree[currNode].right;
            if(cmpSubTree(tree[tree[currNode].right].left, tree[tree[currNode].right].right, tree)){
                setRightChild(currNode, tree[tree[currNode].right].left, tree);
                setLeftChild(root, currNode, tree);
            }else{
                setRightChild(currNode, tree[tree[currNode].right].right, tree);
                setRightChild(root, currNode, tree);
            }
            currNode = root;
        }else if(tree[tree[currNode].left].treeSize > 3 * tree[tree[currNode].right].treeSize){
            root = tree[currNode].left;
            if(cmpSubTree(tree[tree[currNode].left].left, tree[tree[currNode].left].right, tree)){
                setLeftChild(currNode, tree[tree[currNode].left].left, tree);
                setLeftChild(root, currNode, tree);
            }else{
                setLeftChild(currNode, tree[tree[currNode].left].right, tree);
                setRightChild(root, currNode, tree);
            }
            currNode = root;
        }
        while(true){
            //TODO: balance


            tree[currNode].treeSize++;
            if(tree[currNode].value == values[i]){
                tree[currNode].count += 1;
                break;
            }else if(values[i] < tree[currNode].value){
                if(tree[currNode].left == -1){
                    clearValueNode(&tree[nodeCount]);
                    tree[nodeCount].value = values[i];
                    tree[nodeCount].parent = currNode;
                    tree[nodeCount].count = 1;
                    tree[currNode].left = nodeCount;
                    tree[nodeCount].treeSize = 1;
                    nodeCount++;
                    break;
                }else{
                    currNode = tree[currNode].left;
                }
            }else if(values[i] > tree[currNode].value){
                if(tree[currNode].right == -1){
                    clearValueNode(&tree[nodeCount]);
                    tree[nodeCount].value = values[i];
                    tree[nodeCount].parent = currNode;
                    tree[nodeCount].count = 1;
                    tree[nodeCount].treeSize = 1;
                    tree[currNode].right = nodeCount;
                    nodeCount++;
                    break;
                }else{
                    currNode = tree[currNode].right;
                }
            }
            if(cmpSubTree(tree[currNode].left, tree[currNode].right, tree)){
            }
        }
    }
    cerr << "end counting" << endl;
    uint64_t *counts = new uint64_t[nodeCount];
    map<double, uint64_t> frequencyTree;
    cerr << "begin sorting" << endl;
    sort(tree,tree+nodeCount, countCompare);
    cerr << "debug: high: " <<  tree[0].count << " low: " << tree[nodeCount-1].count << endl;

    //dfsMaxCount(tree, root, counts, common, commonSize);
    cerr << "end sorting" << endl;
    cerr << "nodeCount: " << nodeCount << endl;
    cerr << "commonSize: " << commonSize << endl;
    for(uint64_t i=0; i < nodeCount;i++){
        if(i == commonSize)
            break;
        common[i] = tree[i].value;
        counts[i] = tree[i].count;
    }
    if(nodeCount < commonSize)
        commonSize = nodeCount;
    free(tree);

    while(counts[commonSize-1] < minRepeat){
        commonSize--;
        if(commonSize == 0)
            break;
    }
    delete[] counts;
    return commonSize;
    //TODO:
        uint64_t* depricatedFormat = (uint64_t*)malloc(2 * nnz * sizeof(uint64_t));
        for(uint64_t i = 0; i < nnz; i++){
            depricatedFormat[2 * i + 1] = values[i];
        }
        cerr << "entering find256mcvs" << endl;
        int ret = find256Mcv(mcvs, depricatedFormat, nnz);
        cerr << "exiting find256mcvs" << endl;
        free(depricatedFormat);
        *commonValueIndices = (int*)malloc(nnz * sizeof(int));
        for(uint64_t i = 0; i < nnz; i++){
            (*commonValueIndices)[i] = -1;
            for(int j = 0; j < 256; j++){
                if(values[i] == (*mcvs)[j]){
                    (*commonValueIndices)[i] = j;
                    break;
                }
            }
        }
    }
}
*/

void dfsMaxCount(valueNode* tree, uint64_t currNode, uint64_t* counts, double* values, uint64_t commonSize){
    if(tree[currNode].left != -1){
        dfsMaxCount(tree, tree[currNode].left, counts, values, commonSize);
    }
    if(tree[currNode].right != -1){
        dfsMaxCount(tree, tree[currNode].right, counts, values, commonSize);
    }
    int tmpCount = tree[currNode].count;
    double tmpVal = tree[currNode].value;
    int swapCount;
    double swapVal;
    for(int i = 0; i < commonSize; i++){
        if(counts[i] < tmpCount){
            swapCount = counts[i];
            counts[i] = tmpCount;
            tmpCount = swapCount;
            swapVal = values[i];
            values[i] = tmpVal;
            tmpVal = swapVal;
        }
    }
}

