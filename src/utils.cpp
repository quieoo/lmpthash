#include <vector>
#include <cstdint>
#include <string>
#include <cstdio>
#include <fstream>
#include <unordered_set>
#include <string.h>

const int page_size = 4096;

void parse_MSR_Cambridge(std::vector<uint64_t>&  uniq_lpn, std::vector<uint64_t>& lpns, std::string filename){
    // hash table for unique lpn
    std::unordered_set<uint64_t> ht;

    // open file with name "filename", and read by lines
    std::ifstream file(filename);
    std::string line;
    uint64_t timestamp, offset, size, t0;
    char trace_name[100];
    char op[100];
    int trace_id;
    uint64_t lpn;
    while (std::getline(file, line)) {
        sscanf(line.c_str(), "%lu,%100[^,],%d,%100[^,],%lu,%lu,%lu\n", &timestamp,trace_name,&trace_id,op,&offset,&size,&t0);
        for(int i=0;i<size/page_size;i++){
            lpn=offset/page_size+i;
            lpns.push_back(lpn);
            ht.insert(lpn);
        }
    }

    // get unique lpn
    for (auto it = ht.begin(); it != ht.end(); it++) {
        uniq_lpn.push_back(*it);
    }
    file.close();
    return;   
}

void output_access(std::vector<uint64_t>& access){
    // write access to file, each line is a access
    std::ofstream outfile;
    outfile.open("access.txt");
    for(auto it=access.begin();it!=access.end();it++){
        outfile<<*it<<std::endl;
    }
    outfile.close();
    return;
}

void locality(std::vector<uint64_t>& access, int interval, int window){
    // 对于每个Access，与LPN值相距在"interval"之内的上一个Access之间相隔的访问次数
    std::vector<int> precedor;
    precedor.resize(window);
    for(int i=0;i<access.size();i++){
        int in_window=0;
        for(int j=i-1;j>=0 && i-j<window;j--){
            if(access[j]<=access[i] && access[i]-access[j]<=interval){
                precedor[i-j]++;
                in_window=1;
                break;
            }
        }
        if(!in_window) precedor[0]++;
    }

    // output
    size_t sum=0;
    for(int i=1;i<window;i++){
        sum+=precedor[i];
        printf("%f\n", i, (double)sum/access.size());
    }
    return;
}

int main(int argc, char** argv){
    printf("Utils Usage\n");
    printf("----parse trace files----\n");
    printf("    <operation>: parse_csv, parse_output_csv\n");
    printf("    <filename>\n");
    printf("-------------------------\n");

    if(strcmp(argv[1], "parse_csv")==0){
        std::vector<uint64_t> uniq_lpn;
        std::vector<uint64_t> lpns;
        parse_MSR_Cambridge(uniq_lpn, lpns, std::string(argv[2]));

        printf("parse: %s, %lu lpns, %lu uniq lpns\n", argv[2], lpns.size(), uniq_lpn.size());
        // int window=1024*1024/8;
        int window=1000;
        int interval=10;
        locality(lpns, interval, window);
    }else if(strcmp(argv[1], "parse_output_csv")==0){
        std::vector<uint64_t> access;
        std::vector<uint64_t> lpns;
        parse_MSR_Cambridge(lpns, access, std::string(argv[2]));
        output_access(access);
    }else{
        printf("unknown operation\n");
    }

} 
