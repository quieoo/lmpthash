#include <stdio.h>
#include <stdint.h>

typedef struct {
    uint64_t part1;
    uint64_t part2;
} A;

typedef uint64_t B;

B find_max_b_element(A* a, size_t a_len, size_t o, size_t x, B Y) {
    // 确保起始索引是64字节对齐的
    size_t start_index = (o * sizeof(A)) / 64 * 64 / sizeof(A);
    size_t skip_b_elements = x / 2;
    size_t aligned_start_index = start_index + (skip_b_elements / 4) * 4;

    B max_b = 0;
    int found = 0;

    for (size_t i = aligned_start_index; i < a_len; i += 4) {
        for (size_t j = 0; j < 4 && (i + j) < a_len; ++j) {
            B b1 = a[i + j].part1;
            B b2 = a[i + j].part2;

            if (b1 <= Y && (found == 0 || b1 > max_b)) {
                max_b = b1;
                found = 1;
            }
            if (b2 <= Y && (found == 0 || b2 > max_b)) {
                max_b = b2;
                found = 1;
            }
        }
    }

    return max_b;
}

int main() {
    // 示例数组
    A a[] = {
        {1, 2}, {3, 4}, {5, 6}, {7, 8},
        {9, 10}, {11, 12}, {13, 14}, {15, 16}
    };
    size_t a_len = sizeof(a) / sizeof(a[0]);
    size_t o = 1;  // 偏移量
    size_t x = 4;  // 需要跳过的B元素数量
    B Y = 10;      // 目标值

    B result = find_max_b_element(a, a_len, o, x, Y);
    printf("小于等于 %lu 的最大 B 元素是: %lu\n", Y, result);

    return 0;
}
