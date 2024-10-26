#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <vector>

int number_bacteria;
char** bacteria_name;
long M, M1, M2;
short code[27] = { 0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3};
#define encode(ch)		((ch >= 'A' && ch <= 'Z') ? code[ch-'A'] : -1)
#define LEN				6
#define AA_NUMBER		20
#define	EPSILON			1e-010

void Init()
{
    M2 = 1;
    for (int i=0; i<LEN-2; i++)	// M2 = AA_NUMBER ^ (LEN-2);
        M2 *= AA_NUMBER; 
    M1 = M2 * AA_NUMBER;		// M1 = AA_NUMBER ^ (LEN-1);
    M  = M1 * AA_NUMBER;		// M  = AA_NUMBER ^ (LEN);
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
        FILE* bacteria_file = fopen(filename, "r");

        if (bacteria_file == NULL)
        {
            fprintf(stderr, "Error: failed to open file %s\n", filename);
            exit(1);
        }

        InitVectors();

        char ch;
        while ((ch = fgetc(bacteria_file)) != EOF)
        {
            if (ch == '>')
            {
                while (fgetc(bacteria_file) != '\n' && !feof(bacteria_file)); // skip rest of line

                char buffer[LEN-1];
                if (fread(buffer, sizeof(char), LEN-1, bacteria_file) == LEN-1) {
                    init_buffer(buffer);
                }
            }
            else if (ch != '\n')
                cont_buffer(ch);
        }

        long total_plus_complement = total + complement;
        double total_div_2 = total * 0.5;
        int i_mod_aa_number = 0;
        int i_div_aa_number = 0;
        long i_mod_M1 = 0;
        long i_div_M1 = 0;

        double one_l_div_total[AA_NUMBER];
        for (int i=0; i<AA_NUMBER; i++)
            one_l_div_total[i] = (total_l > 0) ? (double)one_l[i] / total_l : 0;
        
        double* second_div_total = new double[M1];
        for (int i=0; i<M1; i++)
            second_div_total[i] = (total_plus_complement > 0) ? (double)second[i] / total_plus_complement : 0;

        count = 0;
        double* t = new double[M];

        #pragma omp parallel for reduction(+:count)
        for (long i = 0; i < M; i++) {
            // Calculate indices directly from 'i'
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

        fclose(bacteria_file);
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
            p2++;
            vector_len2 += (t2 * t2);
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
        long n1 = b1->ti[p1];
        double t1 = b1->tv[p1++];
        vector_len1 += (t1 * t1);
    }
    while (p2 < b2->count)
    {
        long n2 = b2->ti[p2];
        double t2 = b2->tv[p2++];
        vector_len2 += (t2 * t2);
    }

    return (vector_len1 > 0 && vector_len2 > 0) ? correlation / (sqrt(vector_len1) * sqrt(vector_len2)) : 0;
}

// Shared data structures
std::vector<Bacteria*> bacteria_list;
std::mutex bacteria_list_mutex;

std::queue<std::pair<int,int>> task_queue;
std::mutex task_queue_mutex;
std::condition_variable task_queue_cv;

bool loading_done = false;

void producer(int start, int end)
{
    for(int i = start; i < end; i++)
    {
        printf("Loading bacteria %d of %d\n", i+1, number_bacteria);
        Bacteria* b = new Bacteria(bacteria_name[i]);

        // Lock the bacteria list mutex to update the list
        {
            std::lock_guard<std::mutex> lock(bacteria_list_mutex);
            bacteria_list[i] = b;

            // Generate comparison tasks with previously loaded bacteria
            for(int j = 0; j < i; j++)
            {
                if (bacteria_list[j] != nullptr)
                {
                    // Lock task queue mutex to add the task
                    {
                        std::lock_guard<std::mutex> task_lock(task_queue_mutex);
                        task_queue.push(std::make_pair(j, i));
                        task_queue_cv.notify_one();
                    }
                }
            }
        }
    }
}

void consumer()
{
    while(true)
    {
        std::pair<int,int> task;
        {
            std::unique_lock<std::mutex> lock(task_queue_mutex);
            task_queue_cv.wait(lock, []{ return !task_queue.empty() || loading_done; });
            if(task_queue.empty())
            {
                if(loading_done)
                    break;
                else
                    continue;
            }
            task = task_queue.front();
            task_queue.pop();
        }

        int i = task.first;
        int j = task.second;

        Bacteria* bi = nullptr;
        Bacteria* bj = nullptr;

        {
            // No need to lock bacteria_list_mutex here since the list is only appended to
            bi = bacteria_list[i];
            bj = bacteria_list[j];
        }

        if (bi != nullptr && bj != nullptr)
        {
            printf("Comparing bacteria %d and %d\n", i, j);
            double correlation = CompareBacteria(bi, bj);
            printf("Result between %d and %d: %.20lf\n", i, j, correlation);
        }
        else
        {
            fprintf(stderr, "Error: Bacteria %d or %d not loaded yet\n", i, j);
        }
    }
}

void CompareAllBacteria()
{
    int num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 2; // Default to 2 if unable to detect

    // Initialize bacteria_list with size number_bacteria
    bacteria_list.resize(number_bacteria, nullptr);

    // Start producer threads
    std::vector<std::thread> producer_threads;

    int bacteria_per_thread = (number_bacteria + num_threads - 1) / num_threads;

    for (int t = 0; t < num_threads; t++)
    {
        int start = t * bacteria_per_thread;
        int end = std::min(start + bacteria_per_thread, number_bacteria);
        producer_threads.emplace_back(producer, start, end);
    }

    // Start consumer threads
    std::vector<std::thread> consumer_threads;
    for (int t = 0; t < num_threads; t++)
    {
        consumer_threads.emplace_back(consumer);
    }

    // Wait for producers to finish
    for (auto& t : producer_threads)
        t.join();

    // Set loading_done to true and notify consumers
    {
        std::lock_guard<std::mutex> lock(task_queue_mutex);
        loading_done = true;
    }
    task_queue_cv.notify_all();

    // Wait for consumers to finish
    for (auto& t : consumer_threads)
        t.join();

    // Clean up
    for (int i = 0; i < number_bacteria; i++)
    {
        delete bacteria_list[i];
    }
}

int main(int argc,char * argv[])
{
    time_t t1 = time(NULL);

    Init();
    ReadInputFile("list.txt");
    CompareAllBacteria();

    time_t t2 = time(NULL);
    printf("time elapsed: %ld seconds\n", (long)(t2 - t1)); 

    for(int i=0; i<number_bacteria; i++)
    {
        delete[] bacteria_name[i];
    }
    delete[] bacteria_name;

    return 0;
}
