#include <stdio.h>
#include <time.h>
#include <stdlib.h>


void writeEyeMatrix(const int N){
      FILE *fptr;
      fptr = fopen("./matrixEye.txt","w");
      for(int i=0;i<N;i++){
          for(int j = 0; j< N;j++){
            if(i==j){
                fprintf(fptr,"%d ",1);
            }else{
                fprintf(fptr,"%d ",0);
            }
            if(i< N-1 &&j==N-1){fprintf(fptr,"\n");}
          }
      }
      fclose(fptr);
}

void writeRandMatrix(const int N){
      srand(time(NULL));
      FILE *fptr;
      fptr = fopen("./RandMatrix.txt","w");
      for(int i=0;i<N;i++){
          for(int j = 0; j< N;j++){
            fprintf(fptr,"%d ",rand()%10);
            if(i< N-1 &&j==N-1){fprintf(fptr,"\n");}
          }
      }
      fclose(fptr);
}

void writeVector(const int N){
    srand(time(NULL));
      FILE *fptr;
      fptr = fopen("./vector.txt","w");
      for(int i=0;i<N;i++){
            fprintf(fptr,"%d ",rand()%2);
            if(i!=N-1){fprintf(fptr,"\n");}
      }
      fclose(fptr);
}

int main(int argc,char** argv){
    if(argc!=2){
        printf("wring arguments\n");
        return 0;
    }
    writeEyeMatrix(atoi(argv[1]));
    writeRandMatrix(atoi(argv[1]));
    writeVector(atoi(argv[1]));
}