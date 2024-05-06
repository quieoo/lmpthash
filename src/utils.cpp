#include <vector>
#include <cstdint>
#include <string>
#include <cstdio>
#include <fstream>
#include <unordered_set>
#include <string.h>
#include <iostream>
#include <algorithm>
#include "lmpthash.hpp"

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

void num_segments(std::vector<uint64_t>& lpns){
    //sort lpns
    std::sort(lpns.begin(), lpns.end());
    // get the number of segments: each lpn is 1 bigger than the previous one
    int num_segments = 1;
    for(int i=1;i<lpns.size();i++){
        if(lpns[i]!=lpns[i-1]+1){
            num_segments++;
        }
    }
    printf("num lpns: %d\n", lpns.size());
    printf("num segments: %d\n", num_segments);
}


double mergeProfit(const std::vector<std::pair<int, int>>& segments, int m, int i) {
    // Calculate the profit from merging segments from m to i
    int total_length = 0;
    double profit = 0.0;
    for (int index = m; index <= i; ++index) {
        total_length += segments[index].second;
    }

    int start = segments[m].first;
    int end = start + total_length - 1;
    if (m == 0) {
        profit = total_length;  // Full length as profit if starts from first segment
    } else {
        int prev_end = segments[m-1].first + segments[m-1].second - 1;
        int gap = segments[m].first - prev_end - 1;
        if (gap > 0) {
            profit = static_cast<double>(total_length) / (gap + 1);
        } else {
            profit = total_length;  // Continuous case
        }
    }
    return profit;
}

double dpMaximizeProfit(int N, const std::vector<std::pair<int, int>>& segments, int P) {
    std::vector<double> dp(P + 1, 0.0);

    // Iterate over each segment
    for (int i = 1; i <= N; ++i) {
        // Update dp array from back to front
        for (int j = P; j > 0; --j) {
            double profit = 0.0;
            for (int m = i; m > 0; --m) {
                profit = mergeProfit(segments, m - 1, i - 1);
                // Update dp[j] if merging from m to i under j segments
                if (j > 1) {  // Ensure there are enough previous segments to form j-1 segments
                    dp[j] = std::max(dp[j], dp[j-1] + profit);
                }
            }
        }
    }
    return dp[P];
}

void test_dp() {
    int N = 5;  // Number of segments
    std::vector<std::pair<int, int>> segments = {{1, 3}, {5, 2}, {8, 3}, {12, 1}, {14, 2}};  // (F, A) pairs
    int P = 3;  // Desired number of segments
    std::cout << "Maximum profit: " << dpMaximizeProfit(N, segments, P) << std::endl;
}

void test_greedy() {
    MonoSegmentMerger<uint64_t> merger(0.5, 0.5, 0.5, 2);
    merger.mockkeys();
    merger.PrintSegs();
    merger.GreedyMerge();
    merger.PrintSegs();
}




int main(int argc, char** argv){
    printf("Utils Usage\n");
    printf("----parse trace files----\n");
    printf("    <operation>: parse_csv, parse_output_csv\n");
    printf("    <filename>\n");
    printf("----build segments----\n");
    printf("    <operation>: build_segs\n");
    printf("    <filename>\n");
    printf("-------------------------\n");

    if(strcmp(argv[1], "parse_csv")==0){
        std::vector<uint64_t> uniq_lpn;
        std::vector<uint64_t> lpns;
        parse_MSR_Cambridge(uniq_lpn, lpns, std::string(argv[2]));

        printf("parse: %s, %lu lpns, %lu uniq lpns\n", argv[2], lpns.size(), uniq_lpn.size());
        // int window=1024*1024/8;
        // int window=1000;
        // int interval=10;
        // locality(lpns, interval, window);
        // num_segments(uniq_lpn);

        MonoSegmentMerger<uint64_t> merger(0.5, 0.5, 0.5, 65536*9/10);
        // MonoSegmentMerger<uint64_t> merger(0.5, 0.5, 0.5, 3000);
        merger.LoadKeys(uniq_lpn);
        merger.GreedyMerge();
        merger.ScoreSegs();
    }else if(strcmp(argv[1], "parse_output_csv")==0){
        std::vector<uint64_t> access;
        std::vector<uint64_t> lpns;
        parse_MSR_Cambridge(lpns, access, std::string(argv[2]));
        output_access(access);
    }else if(strcmp(argv[1], "build_segs")==0){
        std::vector<uint64_t> uniq_lpn;
        std::vector<uint64_t> lpns;
        parse_MSR_Cambridge(uniq_lpn, lpns, std::string(argv[2]));
        printf("parse: %s, %lu lpns, %lu uniq lpns\n", argv[2], lpns.size(), uniq_lpn.size());
        MonoSegmentMerger<uint64_t> merger(0.5, 0.5, 0.5, 65536*9/10);
        // MonoSegmentMerger<uint64_t> merger(0.5, 0.5, 0.5, 3000);
        merger.LoadKeys(uniq_lpn);
        merger.GreedyMerge();
        merger.ScoreSegs();
    }
    else{
        printf("unknown operation\n");
    }

} 
