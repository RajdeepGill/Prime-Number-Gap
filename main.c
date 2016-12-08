#include <stdio.h>
#include "mpi.h"
#include <gmp.h>

main(int argc, char** argv){
    int my_rank;   /* My process rank           */
    int p;         /* The number of processes   */
    unsigned long int MAX_VALUE = 1000000000;
    int source;    /* Process sending integral  */
    int dest = 0;  /* All messages go to 0      */
    int tag = 0;
    unsigned long int range;
    int N = 5; // size of MPI_LONG_INT array to send to proc 0
    unsigned long int data[N];
    MPI_Status  status;
    double elapsed_time;
    
    // these are the range
    mpz_t startval; // first value to look for prime
    mpz_t endval; // last value to look for prime
    
    // these are for communication
    mpz_t startprime; // store first prime found in range
    mpz_t endprime; // store last prime found in range
    
    // these are for final results
    mpz_t prime1; // first prime realizing start of gap
    mpz_t prime2; // second prime realizing end of gap
    mpz_t gap; // difference between local_prime1 and 2
    
    // these are temp vars
    mpz_t curprime; // current prime found
    mpz_t prevprime; // previously found prime
    
    //START PARELLEL CODE
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    /* Find out how many processes are being used */
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    range = MAX_VALUE/p;
    // initialize all mpz_t variables as longs
    mpz_init_set_ui(startval, my_rank*range);
    mpz_init_set_ui(endval, (my_rank+1)*range);
    mpz_init2(prevprime, 32);
    mpz_init2(curprime, 32);
    mpz_init2(prime1, 32);
    mpz_init2(prime2, 32);
    mpz_init2(startprime, 32);
    mpz_init2(endprime, 32);
    mpz_init_set(count, startval); // initialize count to first number in interval
    mpz_init_set_ui(gap, 0);
    
    // set startprime to first prime found in range
    mpz_nextprime(startprime, startval);
    mpz_set(curprime, startprime);
    mpz_set(prevprime, curprime);
    // find another prime for curprime
    mpz_nextprime(curprime, prevprime);
    
    // tempgap holds gap found right now. will be compared to stored gap
    mpz_t tempgap;
    mpz_init2(tempgap, 32);
    // calculate each prime in range and set gap accordingly
    
    //Set timer
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();
    
    while (mpz_cmp(curprime, endval) < 0){
        // calculate gap between cur and prev primes
        mpz_sub(tempgap, curprime, prevprime);
        // check if temp gap is greater than gap
        if ( mpz_cmp(tempgap, gap) > 0) {
            mpz_set(gap, tempgap);
            mpz_set(prime1, prevprime);
            mpz_set(prime2, curprime);
            //gmp_printf("found new gap: %Zd\n", gap);
            //gmp_printf("prime1: %Zd, prime2: %Zd\n", prime1, prime2);
        }
        // calc next prime
        mpz_set(prevprime, curprime);
        mpz_nextprime(curprime, prevprime);
    }
    mpz_set(endprime, curprime);
    
    // communicate startprime and endprime to process 0
    if (my_rank == 0){
        // variables used to calculate the largest gap
        unsigned long int finalgap = mpz_get_ui(gap);
        unsigned long int finalprime1 = mpz_get_ui(prime1);
        unsigned long int finalprime2 = mpz_get_ui(prime2);
        unsigned long int procdata[p][2]; // slot for each process, holding 2 pieces of data each
        procdata[0][0] = mpz_get_ui(startprime); // store startprime
        procdata[0][1] = mpz_get_ui(endprime); // store endprime
        
        for (source = 1; source < p; source++){
            MPI_Recv(&data, 5, MPI_UNSIGNED_LONG, source, tag, MPI_COMM_WORLD, &status);
            
            //printf("Recieved from %d: Prime 1 = %ld, Prime 2 = %ld, Gap = %ld, Startprime = %ld, Endprime = %ld\n", source, data[0], data[1], data[2], data[3], data[4]);
            // check if the received process gap is larger than the current largest gap
            if (data[2] > finalgap){
                finalprime1 = data[0];
                finalprime2 = data[1];
                finalgap = data[2];
            }
            procdata[source][0] = data[3]; // store startprime
            procdata[source][1] = data[4]; // store endprime
        }
        
        // check if interprocess gaps are larger than current gap
        int i;
        for (i = 0; i < p - 1; i++){
            if( (procdata[i+1][0] - procdata[i][1]) > finalgap){
                finalprime1 = procdata[i][1];
                finalprime2 = procdata[i+1][0];
                finalgap = procdata[i+1][0] - procdata[i][1];
                //printf("new gap found from interprocess: %ld\n", finalgap);
            }
        }
        
        //Update timer
        elapsed_time += MPI_Wtime();
        printf("Prime 1: %ld\n",finalprime1);
        printf("Prime 2: %ld\n",finalprime2);
        printf("Gap: %ld\n",finalgap);
        printf("Elasped Time: %f\n",elapsed_time);
    }else{
        // fill data array before sending
        data[0] = mpz_get_ui(prime1);
        data[1] = mpz_get_ui(prime2);
        data[2] = mpz_get_ui(gap);
        data[3] = mpz_get_ui(startprime);
        data[4] = mpz_get_ui(endprime);
        //printf("%d sending: Prime 1 = %ld, Prime 2 = %ld, Gap = %ld, Startprime = %ld, Endprime = %ld\n", my_rank ,data[0], data[1], data[2], data[3], data[4]);
        MPI_Send(&data, 5, MPI_UNSIGNED_LONG, dest, tag, MPI_COMM_WORLD);
    }
    MPI_Finalize();
}