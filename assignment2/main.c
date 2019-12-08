#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "bitmap.h"
#include <math.h>
#include <string.h>

#define XSIZE 2560 // Size of before image
#define YSIZE 2048

int main(int argc, char **argv) {
		int ierr,num_procs, my_id;
    ierr = MPI_Init(&argc, &argv);
    ierr=MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    ierr=MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
		uchar *r;
		int extra=0;
		int coloumns=0;
		uchar *temp;
		//if the number of processes are 1 then it will only run on that one core
		if(num_procs==1)
		{
			uchar *image = calloc(XSIZE * YSIZE * 3, 1); // Three uchars per pixel (RGB)
			uchar *newimage = calloc(4*XSIZE * YSIZE * 3, 1);
			readbmp("before.bmp", image);
    	alter(image,newimage,XSIZE * YSIZE * 3);
			savebmp("after.bmp", newimage,2560*2,2048*2);
			free(image);
			free(newimage);
		}
		//otherwise when there are several cores the rank 0 process will first read the bmp and then send pieces to the other processes
		//the other processes will work on their part and return a modified array which the rank 0 process is waiting for.
		//by combining the different pieces together and saving them the image is scaled and altered
		else
		{
			if(YSIZE%(num_procs-1)==0)
			{
				coloumns=YSIZE/(num_procs-1);
			}
			else
			{
				coloumns=(int)floor(YSIZE/(num_procs-1));
				extra=YSIZE%(num_procs-1)+coloumns;
			}

	    if(my_id==0)
	    {
					uchar *image = calloc(XSIZE * YSIZE * 3, 1);
	        readbmp("before.bmp", image);
	        for(int i=1;i<num_procs;i++)
	        {
	            int size=0;
							if(extra!=0&&i==num_procs-1)
							{
								size=extra*XSIZE*3;
							}
							else
							{
								size=XSIZE*3*coloumns;
							}
	            temp = calloc(size, 1);

						  for(int t=0;t<size;t++)
							{
								temp[t]=image[t+XSIZE*3*coloumns*(i-1)];
							}
							MPI_Send(temp, size, MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD);
							free(temp);
	        }
					free(image);
					image=calloc(XSIZE*YSIZE*3*4,1);
					for(int i=1;i<num_procs;i++)
					{
							int rlength=0;
							if(i==num_procs-1&&extra!=0)
							{
								rlength=extra*XSIZE*3*4;
							}
							else
							{
								rlength=XSIZE*3*coloumns*4;
							}
							r=calloc(rlength,1);
							MPI_Recv(r,rlength, MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							for(int t=0;t<rlength;t++)
							{

									image[t+(i-1)*XSIZE*3*coloumns*4]=r[t];
							}

							free(r);
					}

	        savebmp("after.bmp", image,2560*2,2048*2);
					free(image);

	    }
	    else{
					int n_p=0;
					if(extra!=0&&my_id==num_procs-1)
					{
							n_p = extra*XSIZE*3;
					}
					else
					{
						n_p = XSIZE*3*coloumns;
					}
					temp = calloc(n_p,1);
	        MPI_Recv(temp, n_p, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				  r=calloc(n_p*4,1);
				  alter(temp,r,n_p);
				  MPI_Send(r, 4*n_p, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
					free(temp);
					free(r);
	    }
		}

	ierr = MPI_Finalize();
	return 0;
}
