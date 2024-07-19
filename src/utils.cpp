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
#include <iomanip>





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
        printf("%d %f\n", i, (double)sum/access.size());
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


struct PhyAddr{
    uint8_t data[20];
    PhyAddr(uint64_t LogicAddr){
        for(int i=0;i<8;i++){
            data[7-i]=LogicAddr>>i*8;
        }
    }
    PhyAddr() = default;

    bool operator!=(const PhyAddr& rhs) const {
        return memcmp(data, rhs.data, 20)!=0;
    }

    void output(){
        for(int i=0;i<20;i++){
            printf("%02x ", data[i]);
        }
    }
};

int main(int argc, char** argv){
    printf("Utils Usage\n");
    printf("----parse trace files----\n");
    printf("    <operation>: parse_csv, parse_output_csv\n");
    printf("    <filename>\n");
    printf("----build segments----\n");
    printf("    <operation>: build_segs\n");
    printf("    <config_file_path>\n");
    printf("-------------------------\n");
    printf("----Parse Trace Query----\n");
    printf("    <operation>: parse_query\n");
    printf("    <config_file_path>\n");
    printf("-------------------------\n");
    printf("----parse femu trace----\n");
    printf("    <operation>: parse_femu\n");
    printf("    <config_file_path>\n");
    printf("-------------------------\n");
    printf("----parse_femu_last_n----\n");
    printf("    <operation>: parse_femu_last_n\n");
    printf("    <n>\n");
    printf("    <config_file_path>\n");
    

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
    }else if(strcmp(argv[1], "parse_femu")==0){
        std::vector<uint64_t> uniq_lpn;
        std::vector<uint64_t> lpns;
        parse_femu(uniq_lpn, lpns, std::string(argv[2]));
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
    }else if(strcmp(argv[1], "build_segs")==0){
        lmpthash_config cfg;
        cfg.load_config(std::string(argv[2]));

        std::vector<uint64_t> uniq_lpn;
        std::vector<uint64_t> lpns;
        parse_MSR_Cambridge(uniq_lpn, lpns, cfg.trace_path);
        printf("parse: %s, %lu lpns, %lu uniq lpns\n", cfg.trace_path.c_str(), lpns.size(), uniq_lpn.size());
        std::vector<PhyAddr> ppns;
        for(uint64_t i=0;i<uniq_lpn.size();i++){
            ppns.push_back(PhyAddr(uniq_lpn[i]));
        }

        LMPTHashBuilder<uint64_t, PhyAddr> builder(cfg);
        builder.Segmenting(uniq_lpn);
        builder.Learning();
        // builder.Bucketing();
        builder.Multi_Bucketing();

        builder.Tabling(uniq_lpn, ppns);

        builder.Verifing(uniq_lpn, ppns);
        builder.Cleaning();
    }else if(strcmp(argv[1], "parse_query")==0){
        lmpthash_config cfg;
        cfg.load_config(std::string(argv[2]));
        std::vector<uint64_t> uniq_lpn;
        std::vector<uint64_t> lpns;
        parse_MSR_Cambridge(uniq_lpn, lpns, cfg.trace_path);
        printf("parse: %s, %lu lpns, %lu uniq lpns\n", cfg.trace_path.c_str(), lpns.size(), uniq_lpn.size());
        output_query_to_file(lpns, 1000*10000);
    }else if (strcmp(argv[1], "model_lindex")==0){
        lmpthash_config cfg;
        cfg.load_config(std::string(argv[2]));

        std::vector<uint64_t> uniq_lpn;
        std::vector<uint64_t> lpns;
        parse_MSR_Cambridge(uniq_lpn, lpns, cfg.trace_path);
        printf("parse: %s, %lu lpns, %lu uniq lpns\n", cfg.trace_path.c_str(), lpns.size(), uniq_lpn.size());

        for(int ep=cfg.left_epsilon;ep<=cfg.right_epsilon;ep++){
            pgm::PGMIndex<uint64_t, 64,4,uint32_t> pgm(uniq_lpn, ep,ep);
            printf("%d \n", pgm.height());
        }
    }
    else if(strcmp(argv[1], "parse_femu_last_n")==0){
        int n=atoi(argv[2]);

        std::vector<uint64_t> uniq_lpn;
        std::vector<uint64_t> lpns;
        parse_femu(uniq_lpn, lpns, std::string(argv[3]), n);
        printf("parse: %s, %lu lpns, %lu uniq lpns\n", argv[2], lpns.size(), uniq_lpn.size());
        // int window=1024*1024/8;
        int window=1000;
        int interval=10;
        locality(lpns, interval, window);

    }
    else{
        printf("unknown operation\n");
    }
}
