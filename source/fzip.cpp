#include "fzip.hpp"

int main(int argc, char* argv[]){
    FzipOptions options(argc,argv);
    if(options.compress)
        fzipCompress(options);
    else
        fzipDecompress(options);
}
