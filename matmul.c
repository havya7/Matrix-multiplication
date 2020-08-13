/*
 * Program: Speedup calculation of matrix multiplication with
 *          multi-processing and multi-threading.
 * Author:  K.Havya
 * Roll# :  CS18BTECH11022
 */

#include <stdlib.h> /* for exit, atoi */
#include <stdio.h>  /* for fprintf */
#include <errno.h>  /* for error code eg. E2BIG */
#include <getopt.h> /* for getopt */
#include <assert.h> /* for assert */
#include <time.h>
#include <pthread.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

/*
 * Forward declarations
 */
# define z 8
void usage(int argc, char *argv[]);
void input_matrix(int *mat, int nrows, int ncols);
void output_matrix(int *mat, int nrows, int ncols);
void init_matrix(int *mat, int nrows, int ncols);

int *A, *B, *C;
int crows, ccols;
int arows, acols, brows, bcols;
char interactive = 0;

struct node{         // structure for sending arguments to multi-thread function
	int a;
	int b;
};

unsigned long long int start,end;

unsigned long long 
properclock()    // function for time
{
    struct timespec spec;
    clock_gettime(CLOCK_REALTIME, &spec);

    return (spec.tv_sec*1000000)+(spec.tv_nsec/1000);
}

unsigned long long 
single_thread_mm()
{
	int i,j,k;
	unsigned long long ans;

	A = (int *)malloc(arows*acols*sizeof(int));    /*memory allocation from heap*/
	B = (int *)malloc(brows*bcols*sizeof(int));
	C = (int *)malloc(arows*bcols*sizeof(int));

	if(interactive == 1)
	{
		printf("Enter A:\n");
		fflush(stdout);
		input_matrix(A,arows,acols);
		printf("Enter B:\n");
		fflush(stdout);
		input_matrix(B,brows,bcols);
	}

	else
	{
		init_matrix(A,arows,acols);
		init_matrix(B,brows,bcols);
	}

	start = properclock();         /*starts measuring time*/

	for(i=0;i<arows;i++)
	{
		for(j=0;j<bcols;j++)
		{
			C[i*bcols+j] = 0;

			for(k=0;k<acols;k++)
			{	
				C[i*bcols+j] += (A[acols*i+k])*(B[k*bcols+j]);
			}

		}

	}

	end = properclock();     /*starts measuring time*/

	if(interactive == 1)
	{
		printf("Result:\n");
		fflush(stdout);
		output_matrix(C,arows,bcols);
	}

	ans = (end - start);
	return ans;
}

unsigned long long 
multi_process_mm()
{
	int shmid;
	int *ptr;
	unsigned long long ans;

	shmid = shmget(IPC_PRIVATE,(arows*acols+brows*bcols+arows*bcols)*sizeof(int),IPC_CREAT | 0666);

	ptr = (int *) shmat(shmid, NULL, 0);

    	if(interactive == 1)
   	{
 		printf("Enter A:\n");
 		fflush(stdout);
 		input_matrix(ptr,arows,acols);
 		printf("Enter B:\n");
 		fflush(stdout);
 		input_matrix(ptr+(arows*acols),brows,bcols);
 	}

 	else
 	{
 		init_matrix(ptr,arows,acols);
		init_matrix(ptr+(arows*acols),brows,bcols);
 	}

 	start = properclock();

	
	pid_t cpid = fork();   //creates a child process

	if(cpid == 0)   //child process
	{
		for(int i=0;i<arows/2;i++)  // calculates the first half rows of the resultant matrix 
		{
			for(int j=0;j<bcols;j++)
			{
				*(ptr+(arows*acols)+(brows*bcols)+(i*bcols+j)) = 0;

	 			for(int k=0;k<brows;k++)
	 			{
	 				*(ptr+(arows*acols)+(brows*bcols)+(i*bcols+j)) += 
						(*(ptr+i*acols+k))*(*(ptr+(arows*acols)+(k*bcols+j)));
				}
			}	
		}
		_exit(1);   // for killing child process
	}

	else   // parent process
	{

		for(int i=arows/2;i<arows;i++) // calculates the other half of the resultant matrix
		{
			for(int j=0;j<bcols;j++)
			{
				*(ptr+(arows*acols)+(brows*bcols)+(i*bcols+j)) = 0;

	 			for(int k=0;k<brows;k++)
				{
	 				*(ptr+(arows*acols)+(brows*bcols)+(i*bcols+j)) += 
						(*(ptr+i*acols+k))*(*(ptr+(arows*acols)+(k*bcols+j)));
				}
			}
	
		}
		wait(NULL);  //waits for the child process to die
	}
	

	end = properclock();

	if(interactive == 1)
	{
		printf("Result:\n");
		fflush(stdout);
		output_matrix(ptr+(arows*acols)+(brows*bcols),arows,bcols);
	}

	shmdt((void *) ptr);
	shmctl(shmid, IPC_RMID, NULL);

	ans = (end - start);

	return ans;
}

void *mulfunc(void *arg)
{
	struct node *ptr1 = (struct node *)arg;
	int r1 = ptr1->a;  // gives the no. of rows of the resultant matrix which should be calculated
	int r2 = ptr1->b;  // gives the starting row of the resultant matrix from which it should be calculated
	int j,k,i,r;

	for(i=0;i<r1;i++)
	{
		r = r2;
		for(j=0;j<bcols;j++)
		{
			C[r*bcols+j] = 0;

	 		for(k=0;k<brows;k++)
	 		{
	 			C[r*bcols+j] += A[r*acols+k]*B[k*bcols+j];
	 		}
		}
		r++;
	}

}

unsigned long long 
multi_thread_mm()
{
	pthread_t thr[z];
	int i;
	struct node  s[z+1];
	unsigned long long ans;

	A = (int *)malloc(arows*acols*sizeof(int));
	B = (int *)malloc(brows*bcols*sizeof(int));
	C = (int *)malloc(arows*bcols*sizeof(int));

	if(interactive == 1)
	{
		printf("Enter A:\n");
		fflush(stdout);
		input_matrix(A,arows,acols);
		printf("Enter B:\n");
		fflush(stdout);
		input_matrix(B,brows,bcols);
	}

	else
	{
		init_matrix(A,arows,acols);
		init_matrix(B,brows,bcols);
	}

	int j = arows/z;
	int k = arows%z;

	for(i=0;i<z;i++)
	{
		
		if(k>0)
		{
			s[i].a = j+1;
			k--;
			if(i == 0)
			{
				s[i].b = 0;
			}
			else
			{
				s[i].b = s[i-1].b + s[i-1].a;
			}
		}
		else
		{	
			s[i].a = j;
			
			if(i == 0)
			{
				s[i].b = 0;
			}
			else
			{
				s[i].b = s[i-1].b + s[i-1].a;
			}
		}
		
	}

	start = properclock();

	for(i=0;i<z;i++)
	{
		pthread_create(&thr[i],NULL,mulfunc,(void *)&s[i]);  // creates 8 threads
	}

	for(i=0;i<z;i++)
	{
		pthread_join(thr[i],NULL);
	}

	end = properclock();

	if(interactive == 1)
	{
		printf("Result:\n");
		fflush(stdout);
		output_matrix(C,arows,bcols);
	}

	ans = (end - start);

	return ans;
}

int 
main(int argc, char *argv[])
{
	int c;

	/* Loop through each option (and its's arguments) and populate variables */
	while (1) 
	{
		int this_option_optind = optind ? optind : 1;
		int option_index = 0;
		static struct option long_options[] = {
			{"help",	no_argument,		0, 'h'},
			{"ar",		required_argument,	0, '1'},
			{"ac",		required_argument,	0, '2'},
			{"br",		required_argument,	0, '3'},
			{"bc",		required_argument,	0, '4'},
			{"interactive",	no_argument, 		0, '5'},
			{0,		0,			0,  0 }
		};

		c = getopt_long(argc, argv, "h1:2:3:4:5", long_options, &option_index);

		if (c == -1)
			break;

		switch (c) 
		{
		case 0:
			fprintf(stdout, "option %s", long_options[option_index].name);
			if (optarg)
				fprintf(stdout, " with arg %s", optarg);
				fprintf(stdout, "\n");
			break;

		case '1':
			arows = atoi(optarg);
			break;

		case '2':
			acols = atoi(optarg);
			break;

		case '3':
			brows = atoi(optarg);
			break;

		case '4':
			bcols = atoi(optarg);
			break;

		case '5':
			interactive = 1;
			break;

		case 'h':
		case '?':
			usage(argc, argv);

		default:
			fprintf(stdout, "?? getopt returned character code 0%o ??\n", c);
			usage(argc, argv);
		}
	}

	if (optind != argc) 
	{
		fprintf(stderr, "Unexpected arguments\n");
		usage(argc, argv);
	}

	if(acols == brows)
	{
		unsigned long long time_single, time_multi_process, time_multi_thread;

		/* Add your code here */
		/* TODO */
		time_single = single_thread_mm();
		time_multi_process = multi_process_mm();
		time_multi_thread = multi_thread_mm();

		fprintf(stdout, "Time taken for single threaded: %llu us\n",
			time_single);
		fflush(stdout);

		fprintf(stdout, "Time taken for multi process: %llu us\n",
				time_multi_process);
		fflush(stdout);

		fprintf(stdout, "Time taken for multi threaded: %llu us\n",
				time_multi_thread);
		flush(stdout);

		fprintf(stdout, "Speedup for multi process : %4.2f x\n",
				(double)time_single/time_multi_process);
		fflush(stdout);

		fprintf(stdout, "Speedup for multi threaded : %4.2f x\n",
				(double)time_single/time_multi_thread); 
		fflush(stdout);

		exit(EXIT_SUCCESS);
	}

	else
	{
		exit(EXIT_FAILURE);
	}
}

/*
 * Show usage of the program
 */
void 
usage(int argc, char *argv[])
{
	fprintf(stderr, "Usage:\n");
	fprintf(stderr, "%s --ar <rows_in_A>  --ac <cols_in_A>"
			" --br <rows_in_B>  --bc <cols_in_B>"
			" [--interactive]",
			argv[0]);
	exit(EXIT_FAILURE);
}

/*
 * Input a given 2D matrix
 */
void 
input_matrix(int *mat, int rows, int cols)
{
	for (int i=0; i<rows; i++) 
	{
		for (int j=0; j<cols; j++) 
		{
			fscanf(stdin, "%d", mat+(i*cols+j));
		}
	}
}

/*
 * Output a given 2D matrix
 */
void 
output_matrix(int *mat, int rows, int cols)
{
	for (int i=0; i<rows; i++) 
	{
		for (int j=0; j<cols; j++) 
		{
			fprintf(stdout, "%d ", *(mat+(i*cols+j)));
			fflush(stdout);
		}

		fprintf(stdout, "\n");
		fflush(stdout);
	}
}

void 
init_matrix(int *mat,int rows,int cols)
{
	for (int i=0; i<rows; i++) 
	{
		for (int j=0; j<cols; j++) 
		{
			*(mat+(i*cols+j)) = rand() ;
		}
	}
}
