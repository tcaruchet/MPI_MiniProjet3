#define N_DIM 4

#include <mpi.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h> //memcpy

#pragma region "MatrixCourseGiven"
////////////////////////////////
// Matrix Course Given by Françoise Baude
////////////////////////////////

/**
 * @brief global variable to know the :
 * rank of the process among all the processors
 * size of of the processors
 */
int rank;
int size;

/**
 * @brief N : size of the vector and the matrix
 *        vector : the pointer that will hold the vector
 *        blocks_flatten_matrix : pointer holding the portion of the matrix
 *        blocks_result : pointer to the result of the vector * blocks_flatten_matrix
 *        swap : swap pointer
 */
int N, L, p;
int *vector = NULL;
int *blocks_flatten_matrix = NULL;
int *blocks_result = NULL;
int *swap = NULL;

struct matrix2D
{
        int **tab;
        int rows;
        int columns;
};

struct matrix2D *m = NULL;
int get(struct matrix2D *m, int row, int column)
{
        return m->tab[row][column];
}

void set(struct matrix2D *m, int row, int column, int val)
{
        m->tab[row][column] = val;
}

struct matrix2D *allocateMatrix(int rows, int columns)
{
        struct matrix2D *m = (struct matrix2D *)malloc(sizeof(struct matrix2D));
        m->rows = rows;
        m->columns = columns;
        m->tab = (int **)malloc(sizeof(int *) * rows);
        for (int i = 0; i < rows; i++)
        {
                m->tab[i] = (int *)malloc(sizeof(int) * columns);
        }
        return m;
}

void printMatrix(struct matrix2D m)
{
        printf("Matrix %d x %d\n", m.rows, m.columns);
        int v;
        for (int i = 0; i < m.rows; i++)
        {
                for (int j = 0; j < m.columns; j++)
                {
                        v = get(&m, i, j);
                        if (v == INT_MAX)
                        {
                                printf("i ");
                        }
                        else
                        {
                                printf("%d ", v);
                        }
                }
                printf("\n");
        }
        printf("\n");
}

int getSize(const char *filename)
{
        FILE *file;
        file = fopen(filename, "r");
        int c;
        int counter = 0;
        // we first read one single line to find the size of the matrix
        while ((c = fgetc(file)) != EOF)
        {
                // printf("%c ",c);
                if (c == ' ')
                {
                        counter++;
                }
                if (c == '\n')
                {
                        break;
                }
        }
        return counter;
}

struct matrix2D *readMatrix(const char *filename, const char TAG)
{
        FILE *file;
        file = fopen(filename, "r");
        struct matrix2D *m;
        fseek(file, 0, SEEK_SET);
        int val = 0;
        // Now read the file to build the matrix
        if (TAG == 'M')
        {
                m = allocateMatrix(N, N);
                for (int i = 0; i < m->rows; i++)
                {
                        for (int j = 0; j < m->columns; j++)
                        {
                                fscanf(file, "%d", &val);
                                set(m, i, j, val);
                        }
                }
        }
        // Now Read Vector which is a matrix M x 1 and as we know that N = M so we have already the vector size
        // Now read the file to build the matrix
        if (TAG == 'V')
        {
                m = allocateMatrix(N, 1);
                for (int i = 0; i < m->rows; i++)
                {
                        fscanf(file, "%d", &val);
                        set(m, i, 0, val);
                }
        }
        return m;
}

struct matrix2D *multseq(struct matrix2D *M, struct matrix2D *V)
{
        struct matrix2D *m = allocateMatrix(M->rows, V->columns);
        // multiply A by a vector V, result is in W, initialized to 0
        for (int i = 0; i < M->rows; i++)
        {
                set(m, i, 0, 0);
                for (int j = 0; j < V->rows; j++)
                {
                        set(m, i, 0, get(m, i, 0) + get(V, j, 0) * get(M, i, j));
                }
        }
        return m;
}

int *matrix_to_vector(struct matrix2D *m)
{
        for (int i = 0; i < m->rows; i++)
        {
                vector[i] = get(m, i, 0);
        }
        return vector;
}

void writefile()
{
        FILE *fptr;
        char name[50];
        sprintf(name, "./result_N=%d_X_P=%d.txt", N, size);
        fptr = fopen(name, "w");
        for (int i = 0; i < N; i++)
        {
                fprintf(fptr, "%d ", blocks_result[i]);
                if (i != N - 1)
                {
                        fprintf(fptr, "\n");
                }
        }
        fclose(fptr);
}

void print_array(const int *ARRAY, int size)
{
        putchar('[');
        for (int i = 0; i < size; i++)
        {
                if (i != size - 1)
                {
                        printf("%d,", ARRAY[i]);
                }
                else
                {
                        printf("%d", ARRAY[i]);
                }
        }
        putchar(']');
        putchar('\n');
}

int *flatten_matrix(struct matrix2D *m, const int SIZE)
{
        int *tmp = malloc(sizeof(int) * SIZE);
        for (int i = 0; i < m->rows; i++)
        {
                for (int j = 0; j < m->columns; j++)
                {
                        tmp[i * (N) + j] = get(m, i, j);
                }
        }
        return tmp;
}
void save(int begin)
{
        int *blocks_tmp = malloc(sizeof(int) * L);
        for (int i = 0; i < L; i++)
        {
                blocks_tmp[i] = blocks_flatten_matrix[begin++];
        }
        free(blocks_flatten_matrix);
        blocks_flatten_matrix = blocks_tmp;
}

/**
 * @brief printing the heap of each process
 *
 */
void print_heap()
{
        printf("Heap (rank : %d) {\n", rank);
        if (m != NULL)
        {
                printMatrix(*m);
        }
        if (vector != NULL)
        {
                printf("Vector >> ");
                print_array(vector, N);
        }
        if (blocks_flatten_matrix != NULL)
        {
                printf("Blocks >> ");
                print_array(blocks_flatten_matrix, L);
        }
        printf("}\n");
}

void free_mem()
{
        MPI_Barrier(MPI_COMM_WORLD);
        free(vector);
        vector = NULL;
        free(blocks_flatten_matrix);
        blocks_flatten_matrix = NULL;
        if (m != NULL)
        {
                for (int i = 0; i < N; i++)
                {
                        free(m->tab[i]);
                }
                m = NULL;
        }
        free(blocks_result);
        blocks_result = NULL;
}

void translate()
{
        // translate the first elements to the last
        for (int i = 0; i < N / size; i++)
        {
                blocks_result[i + (N - (N / size))] = blocks_result[i];
        }
}

#pragma endregion

/**
 * @brief Multiplication linéaire : vector * blocks_flatten_matrix
 */
void inlineMultplication()
{
        int taille = N / size;
        if (rank == 0)
        {
                taille = N;
        }
        blocks_result = malloc(sizeof(int) * (taille));
        for (int i = 0; i < N / size; i++)
        {
                blocks_result[i] = 0;
                for (int j = 0; j < N; j++)
                        blocks_result[i] = blocks_result[i] + vector[j] * blocks_flatten_matrix[(N * i) + j];
        }
}

/**
 * @brief Broadcast the result of the multiplication of the vector * blocks_flatten_matrix in a ring (naive way)
 * @details This application dispatches these values to all the processes in the same
 * communicator. Other processes just receive the dispatched value meant for them.
 **/
void broadcasting(const char *vectorName)
{
        int tag = 0;
        // Il faut prêter attention à la valeur du tag, car sinon, les processus bloqueront (recv).
        // Il faut donc savoir quand on est dans la phase de broadcast.
        MPI_Status status;
        vector = malloc(sizeof(int) * N);
        if (rank == 0)
        {                                                         // Emmeteur
                struct matrix2D *v = readMatrix(vectorName, 'V'); // Lecture du vecteur qui est une matrice N x 1 / V : vector
                vector = matrix_to_vector(v);
                MPI_Send(vector, N, MPI_INT, (rank + 1) % size, tag, MPI_COMM_WORLD);
        }
        else if (rank == size - 1)
        {
                MPI_Probe(((rank - 1) + size) % size, tag, MPI_COMM_WORLD, &status);
                int number;                               // Nombre d'entiers à circuler
                MPI_Get_count(&status, MPI_INT, &number); // dynamiquement découverts
                MPI_Recv(vector, number, MPI_INT, ((rank - 1) + size) % size, tag, MPI_COMM_WORLD, &status);
        }
        else
        { // Transmission simple du précedent reçu vers le suivant
                MPI_Probe(((rank - 1) + size) % size, tag, MPI_COMM_WORLD, &status);
                int number;                               // Nombre d'entiers à circuler
                MPI_Get_count(&status, MPI_INT, &number); // dynamiquement découverts
                MPI_Recv(vector, number, MPI_INT, ((rank - 1) + size) % size, tag, MPI_COMM_WORLD, &status);
                MPI_Send(vector, number, MPI_INT, ((rank + 1) + size) % size, tag, MPI_COMM_WORLD);
        }
        printf("-------------PROCESS %d >> @Broadcast done------------\n", rank);
        printf("(p%d) vector : ", rank);
        print_array(vector, N);
}

void broadcastSize(const char *arg)
{
        int tag = 0;
        // Il faut prêter attention à la valeur du tag, car sinon, les processus bloqueront (recv).
        // Il faut donc savoir quand on est dans la phase de broadcast.
        MPI_Status status;
        if (rank == 0)
        { // Emmeteur
                N = getSize(arg);
                MPI_Send(&N, 1, MPI_INT, (rank + 1) % size, tag, MPI_COMM_WORLD);
        }
        else if (rank == size - 1)
        {
                MPI_Probe(((rank - 1) + size) % size, tag, MPI_COMM_WORLD, &status);
                int number;                               // Nombre d'entiers à faire circuler
                MPI_Get_count(&status, MPI_INT, &number); // dynamiquement découverts
                MPI_Recv(&N, number, MPI_INT, ((rank - 1) + size) % size, tag, MPI_COMM_WORLD, &status);
        }
        else
        { // Transmission simple du précedent reçu vers le suivant
                MPI_Probe(((rank - 1) + size) % size, tag, MPI_COMM_WORLD, &status);
                int number;                               // Nombre d'entiers à faire circuler
                MPI_Get_count(&status, MPI_INT, &number); // et qu'on découvre dynamiquement
                MPI_Recv(&N, number, MPI_INT, ((rank - 1) + size) % size, tag, MPI_COMM_WORLD, &status);
                MPI_Send(&N, number, MPI_INT, ((rank + 1) + size) % size, tag, MPI_COMM_WORLD);
        }
        L = N * N / size; // block size 8 = 4*4/2
        printf("+++++++++++++++ (p%d) N=%d , L=%d +++++++++++++++++\n", rank, N, L);
}

void scatter_mpi(const char *matrix)
{
        int tag = 0;
        MPI_Status status;

        if (rank == 0)
        {                                                         // Initialisation de la matrice
                m = readMatrix(matrix, 'M');                      // Le fichier est une matrice N x N
                blocks_flatten_matrix = flatten_matrix(m, N * N); // Applatissement de la matrice
                // Formattage des données à envoyer
                for (int i = 0; i < size - 1; i++)
                        MPI_Send(&blocks_flatten_matrix[L * i], L, MPI_INT, 1, tag, MPI_COMM_WORLD);
                save(L * (size - 1));
        }
        else if (rank == (size - 1))
        { // Finalisation : Fin du Scatter
                MPI_Probe(((rank - 1) + size) % size, tag, MPI_COMM_WORLD, &status);
                int number;                               // Nombre d'entiers à faire circuler
                MPI_Get_count(&status, MPI_INT, &number); // et qu'on découvre dynamiquement
                blocks_flatten_matrix = malloc(sizeof(int) * number);
                MPI_Recv(blocks_flatten_matrix, number, MPI_INT, ((rank - 1) + size) % size, tag, MPI_COMM_WORLD, &status);
        }
        else
        {       // Phase de transfert
                // Transmission simple du précedent reçu vers le suivant
                MPI_Probe(((rank - 1) + size) % size, tag, MPI_COMM_WORLD, &status);
                int number;                               // Nombre d'entiers à faire circuler
                MPI_Get_count(&status, MPI_INT, &number); // dynamiquement découverts
                blocks_flatten_matrix = malloc(sizeof(int) * number);
                for (int i = 0; i < (size - 1) - rank; i++)
                {
                        MPI_Recv(blocks_flatten_matrix, number, MPI_INT, ((rank - 1) + size) % size, tag, MPI_COMM_WORLD, &status);
                        MPI_Send(blocks_flatten_matrix, number, MPI_INT, ((rank + 1) + size) % size, tag, MPI_COMM_WORLD);
                }
                MPI_Recv(blocks_flatten_matrix, number, MPI_INT, ((rank - 1) + size) % size, tag, MPI_COMM_WORLD, &status);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        printf("(p%d) SAVED : ", rank);
        print_array(blocks_flatten_matrix, L);
        inlineMultplication();
        printf("(p%d) result : ", rank);
        print_array(blocks_result, N / size);
}

void gather_mpi()
{
        int tag = 0;
        MPI_Status status;
        int block = N / size;

        if (rank == 0)
        {
                int counter = 0;
                translate();
                for (int i = 0; i < size - 1; i++)
                {
                        swap = malloc(sizeof(int) * block);
                        MPI_Recv(swap, block, MPI_INT, ((rank - 1) + size) % size, tag, MPI_COMM_WORLD, &status);
                        printf("match nadi : ");
                        print_array(swap, block);
                        for (int i = 0; i < block; i++)
                        {
                                blocks_result[counter++] = swap[i];
                        }
                        free(swap);
                }
                // DOnnée finale
                printf("result = ");
                print_array(blocks_result, N);
                writefile();
        }
        else
        {
                // envoyer les données à l'emetteur
                MPI_Send(blocks_result, block, MPI_INT, ((rank + 1) + size) % size, tag, MPI_COMM_WORLD);
                // phase de transfert
                swap = malloc(sizeof(int) * block);
                for (int i = 0; i < rank - 1; i++)
                {
                        MPI_Recv(swap, block, MPI_INT, ((rank - 1) + size) % size, tag, MPI_COMM_WORLD, &status);
                        MPI_Send(swap, block, MPI_INT, ((rank + 1) + size) % size, tag, MPI_COMM_WORLD);
                }
                free(swap);
        }
}

void read_mat_from_file(char *filename, int n_rows, int n_cols, double *matrix_data)
{
        FILE *fp;
        fp = fopen(filename, "r");
        if (fp == NULL)
        {
                printf("Error opening file!\n");
                exit(1);
        }
        int i, j;
        for (i = 0; i < n_rows; i++)
        {
                for (j = 0; j < n_cols; j++)
                {
                        fscanf(fp, "%lf", &matrix_data[i * n_cols + j]);
                }
        }
        fclose(fp);
}

int main(int argc, char *argv[])
{
        if (argc != 3)
        {
                printf("Wrong number of arguments.\n");
                return 1;
        }

#pragma region "Inutile"

        // variables string
        char *file_matrix_name = argv[5];
        char *file_vector_name = argv[6];

        int i, j, k, l;
        int n = N_DIM;
        int p = size;
        int n_lines = n / p;

        double matrix_data[N_DIM][N_DIM]; // matrix
        double vector_data[N_DIM];        // vector
        double result[N_DIM] = {0.0};     // final result holder

        if (rank == 0)
        {
                read_mat_from_file(file_matrix_name, N_DIM, N_DIM, (double *)matrix_data); // Populating the Matrix
                read_mat_from_file(file_vector_name, N_DIM, 1, vector_data);               // Populating the Vector
        }
#pragma endregion

        // Initialise MPI and check its completion
        // return MPI_SUCESS
        MPI_Init(&argc, &argv);
        // get the number of processes
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        // get the rank of the process
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        // we know that the matrix is square so we only read the size once (M * N and that M = N)
        // and then we can use it to allocate the matrix
        broadcastSize(argv[1]);
        MPI_Barrier(MPI_COMM_WORLD); // check every process has the right size N and blockSize L
        broadcasting(argv[2]);       // broadcasting vector to all the processes
        MPI_Barrier(MPI_COMM_WORLD); // check if every process has the right vector v
        scatter_mpi(argv[1]);        // scattering the matrix in the processes who receive the blocks
        gather_mpi();                // gathering the result in the process 0
        free_mem();                  // desallocation of the matrix and the vector
        MPI_Finalize();
        return 0;
}