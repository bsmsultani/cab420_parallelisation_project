#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <chrono>
#include <iostream>
#include <omp.h>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <atomic>

int number_bacteria;
char** bacteria_name;
long M, M1, M2;
short code[27] = { 0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3};
#define encode(ch)		((ch >= 'A' && ch <= 'Z') ? code[ch-'A'] : -1)
#define LEN				6
#define AA_NUMBER		20
#define	EPSILON			1e-010


int NUM_CORES = 8; // Number of cores to use (user can change this value)


void Init()
{
    M2 = 1;
    for (int i=0; i<LEN-2; i++)
        M2 *= AA_NUMBER;
    M1 = M2 * AA_NUMBER;
    M  = M1 * AA_NUMBER;
}

class Bacteria
{
private:
    long* vector;
    long* second;
    long one_l[AA_NUMBER];
    long indexs;
    long total;
    long total_l;
    long complement;

    void InitVectors()
    {
        vector = new long [M];
        second = new long [M1];
        memset(vector, 0, M * sizeof(long));
        memset(second, 0, M1 * sizeof(long));
        memset(one_l, 0, AA_NUMBER * sizeof(long));
        total = 0;
        total_l = 0;
        complement = 0;
    }

    void init_buffer(char* buffer)
    {
        complement++;
        indexs = 0;
        for (int i=0; i<LEN-1; i++)
        {
            short enc = encode(buffer[i]);
            if (enc != -1) {
                one_l[enc]++;
                total_l++;
                indexs = indexs * AA_NUMBER + enc;
            }
        }
        second[indexs]++;
    }

    void cont_buffer(char ch)
    {
        short enc = encode(ch);
        if (enc != -1) {
            one_l[enc]++;
            total_l++;
            long index = indexs * AA_NUMBER + enc;
            if (index < M) {
                vector[index]++;
                total++;
                indexs = (indexs % M2) * AA_NUMBER + enc;
                if (indexs < M1) {
                    second[indexs]++;
                }
            }
        }
    }

public:
    long count;
    double* tv;
    long *ti;

    Bacteria(const char* filename)
    {
        FILE* bacteria_file = fopen(filename, "rb");
        if (bacteria_file == NULL)
        {
            fprintf(stderr, "Error: failed to open file %s\n", filename);
            exit(1);
        }
        InitVectors();
        fseek(bacteria_file, 0, SEEK_END);
        long file_size = ftell(bacteria_file);
        fseek(bacteria_file, 0, SEEK_SET);

        char* buffer = new char[file_size + 1];

        // Read file into buffer
        size_t read_size = fread(buffer, 1, file_size, bacteria_file);
        buffer[read_size] = '\0';

        fclose(bacteria_file);

        char* p = buffer;
        char* buffer_end = buffer + read_size;

        while (p < buffer_end)
        {
            char ch = *p++;
            if (ch == '>')
            {
                while (p < buffer_end && *p != '\n') p++;
                if (p < buffer_end) p++; // Skip '\n'
                // Read next LEN-1 characters into seq_buffer, skipping '\n'
                char seq_buffer[LEN-1];
                int i = 0;
                while (p < buffer_end && i < LEN-1)
                {
                    if (*p != '\n')
                    {
                        seq_buffer[i++] = *p;
                    }
                    p++;
                }
                if (i == LEN-1)
                {
                    init_buffer(seq_buffer);
                }
            }
            else if (ch != '\n')
            {
                cont_buffer(ch);
            }
        }

        delete[] buffer;

        long total_plus_complement = total + complement;
        double total_div_2 = total * 0.5;

        double one_l_div_total[AA_NUMBER];
        for (int i=0; i<AA_NUMBER; i++)
            one_l_div_total[i] = (total_l > 0) ? (double)one_l[i] / total_l : 0;

        double* second_div_total = new double[M1];
        for (int i=0; i<M1; i++)
            second_div_total[i] = (total_plus_complement > 0) ? (double)second[i] / total_plus_complement : 0;

        count = 0;
        double* t = new double[M];

        omp_set_num_threads(NUM_CORES);
        #pragma omp parallel for reduction(+:count)
        for (long i = 0; i < M; i++) {
            int i_mod_aa_number = i % AA_NUMBER;
            int i_div_aa_number = (i / AA_NUMBER) % M1;
            int i_mod_M1 = i % M1;
            int i_div_M1 = i / M1;

            double p1 = second_div_total[i_div_aa_number];
            double p2 = one_l_div_total[i_mod_aa_number];
            double p3 = second_div_total[i_mod_M1];
            double p4 = one_l_div_total[i_div_M1];

            double stochastic = (p1 * p2 + p3 * p4) * total_div_2;

            if (stochastic > EPSILON) {
                t[i] = (vector[i] - stochastic) / stochastic;
                count++;
            } else {
                t[i] = 0;
            }
        }

        delete[] second_div_total;
        delete[] vector;
        delete[] second;

        tv = new double[count];
        ti = new long[count];

        int pos = 0;
        for (long i=0; i<M; i++)
        {
            if (t[i] != 0)
            {
                tv[pos] = t[i];
                ti[pos] = i;
                pos++;
            }
        }
        delete[] t;
    }

    ~Bacteria()
    {
        delete[] tv;
        delete[] ti;
    }
};

void ReadInputFile(const char* input_name)
{
    FILE* input_file = fopen(input_name, "r");

    if (input_file == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s (Hint: check your working directory)\n", input_name);
        exit(1);
    }

    if (fscanf(input_file, "%d", &number_bacteria) != 1) {
        fprintf(stderr, "Error: failed to read number of bacteria from input file\n");
        fclose(input_file);
        exit(1);
    }

    bacteria_name = new char*[number_bacteria];
    for(long i=0; i<number_bacteria; i++)
    {
        char name[10];
        if (fscanf(input_file, "%9s", name) != 1) {
            fprintf(stderr, "Error: failed to read bacteria name from input file\n");
            fclose(input_file);
            exit(1);
        }
        bacteria_name[i] = new char[20];
        snprintf(bacteria_name[i], 20, "data/%s.faa", name);
    }
    fclose(input_file);
}

double CompareBacteria(Bacteria* b1, Bacteria* b2)
{
    double correlation = 0;
    double vector_len1=0;
    double vector_len2=0;
    long p1 = 0;
    long p2 = 0;

    while (p1 < b1->count && p2 < b2->count)
    {
        long n1 = b1->ti[p1];
        long n2 = b2->ti[p2];
        if (n1 < n2)
        {
            double t1 = b1->tv[p1];
            vector_len1 += (t1 * t1);
            p1++;
        }
        else if (n2 < n1)
        {
            double t2 = b2->tv[p2];
            vector_len2 += (t2 * t2);
            p2++;
        }
        else
        {
            double t1 = b1->tv[p1++];
            double t2 = b2->tv[p2++];
            vector_len1 += (t1 * t1);
            vector_len2 += (t2 * t2);
            correlation += t1 * t2;
        }
    }
    while (p1 < b1->count)
    {
        double t1 = b1->tv[p1++];
        vector_len1 += (t1 * t1);
    }
    while (p2 < b2->count)
    {
        double t2 = b2->tv[p2++];
        vector_len2 += (t2 * t2);
    }

    return (vector_len1 > 0 && vector_len2 > 0) ? correlation / (sqrt(vector_len1) * sqrt(vector_len2)) : 0;
}

class TaskQueue {
private:
    std::queue<std::pair<int, int>> tasks;
    std::mutex mtx;
    std::condition_variable cv;
    bool done = false;

public:
    void push_task(const std::pair<int, int>& task) {
        {
            std::lock_guard<std::mutex> lock(mtx);
            tasks.push(task);
        }
        cv.notify_one();
    }

    bool pop_task(std::pair<int, int>& task) {
        std::unique_lock<std::mutex> lock(mtx);
        while (tasks.empty() && !done) {
            cv.wait(lock);
        }
        if (!tasks.empty()) {
            task = tasks.front();
            tasks.pop();
            return true;
        } else {
            return false;
        }
    }

    void set_done() {
        {
            std::lock_guard<std::mutex> lock(mtx);
            done = true;
        }
        cv.notify_all();
    }
};

void CompareAllBacteria() {
    // Shared resources
    std::vector<Bacteria*> b;
    std::vector<int> b_indices; // Map from indices in b to original idx
    std::mutex b_mutex;
    TaskQueue task_queue;
    std::mutex output_mutex;

    // Number of threads
    int num_producers = NUM_CORES;
    int num_consumers = NUM_CORES;

    // Start consumer threads
    std::vector<std::thread> consumers;
    for (int i = 0; i < num_consumers; i++) {
        consumers.emplace_back([&task_queue, &b, &b_indices, &output_mutex, &b_mutex]() {
            std::pair<int, int> task;
            while (task_queue.pop_task(task)) {
                int i = task.first;
                int j = task.second;

                Bacteria* bi;
                Bacteria* bj;
                int idx_i, idx_j;

                // Lock b_mutex to safely access b and b_indices
                {
                    std::lock_guard<std::mutex> lock(b_mutex);
                    bi = b[i];
                    bj = b[j];
                    idx_i = b_indices[i];
                    idx_j = b_indices[j];
                }

                double correlation = CompareBacteria(bi, bj);

                // Store or print the result
                {
                    std::lock_guard<std::mutex> lock(output_mutex);
                    std::cout << std::setw(2) << idx_i << " " << std::setw(2) << idx_j << " -> "
                              << std::fixed << std::setprecision(20) << correlation << std::endl;
                }
            }
        });
    }

    // Start producer threads
    std::vector<std::thread> producers;
    std::atomic<int> next_bacteria_index(0);

    for (int i = 0; i < num_producers; i++) {
        producers.emplace_back([&]() {
            int idx;
            while ((idx = next_bacteria_index.fetch_add(1)) < number_bacteria) {
                Bacteria* bi = new Bacteria(bacteria_name[idx]);

                std::vector<std::pair<int, int>> new_tasks;

                int b_index;
                int existing_size;

                // Lock b to add bi and generate tasks
                {
                    std::lock_guard<std::mutex> lock(b_mutex);
                    b_index = b.size();
                    b.push_back(bi);
                    b_indices.push_back(idx); // Map b_index to idx
                    existing_size = b.size();
                }

                // Generate comparison tasks with all existing bacteria except itself
                for (int j = 0; j < existing_size; j++) {
                    if (j != b_index) {
                        // Push the task immediately
                        task_queue.push_task(std::make_pair(b_index, j));
                    }
                }
            }
        });
    }

    // Wait for producers to finish
    for (auto& t : producers) {
        t.join();
    }

    // Signal that no more tasks will be added
    task_queue.set_done();

    // Wait for consumers to finish
    for (auto& t : consumers) {
        t.join();
    }

    // Cleanup
    for (auto& bi : b) {
        delete bi;
    }
}

int main(int argc, char* argv[])
{
    time_t t1 = time(NULL);

    Init();
    ReadInputFile("list.txt");
    CompareAllBacteria();

    time_t t2 = time(NULL);
    printf("Time elapsed: %ld seconds\n", (long)(t2 - t1));

    for (int i = 0; i < number_bacteria; i++)
    {
        delete[] bacteria_name[i];
    }
    delete[] bacteria_name;

    return 0;
}
