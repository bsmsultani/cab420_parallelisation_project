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

void CompareAllBacteria()
{
    Bacteria** b = new Bacteria*[number_bacteria];
    #pragma omp parallel for
    for (int i = 0; i < number_bacteria; i++)
    {
        printf("Loading bacteria %d of %d\n", i + 1, number_bacteria);
        b[i] = new Bacteria(bacteria_name[i]);
    }

    std::vector<std::string> results;

    #pragma omp parallel for
    for (int i = 0; i < number_bacteria - 1; i++)
    {
        #pragma omp parallel for
        for (int j = i + 1; j < number_bacteria; j++)
        {
            double correlation = CompareBacteria(b[i], b[j]);

            std::ostringstream result;
            result << std::setw(2) << i << " " << std::setw(2) << j << " -> ";
            result << std::fixed << std::setprecision(20) << correlation;

            #pragma omp critical
            results.push_back(result.str());
        }
    }

    for (int i = 0; i < number_bacteria; i++)
    {
        delete b[i];
    }
    delete[] b;

    for (const auto& result : results)
    {
        std::cout << result << std::endl;
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
