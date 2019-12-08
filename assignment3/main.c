#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <stdlib.h>
#include "libs/bitmap.h"
#include <mpi.h>
#include <time.h>
#include <math.h>

// Convolutional Kernel Examples, each with dimension 3,
// gaussian kernel with dimension 5
// If you apply another kernel, remember not only to exchange
// the kernel but also the kernelFactor and the correct dimension.

int const sobelYKernel[] = {-1, -2, -1,
                             0,  0,  0,
                             1,  2,  1};
float const sobelYKernelFactor = (float) 1.0;

int const sobelXKernel[] = {-1, -0, -1,
                            -2,  0, -2,
                            -1,  0, -1 , 0};
float const sobelXKernelFactor = (float) 1.0;


int const laplacian1Kernel[] = {  -1,  -4,  -1,
                                 -4,  20,  -4,
                                 -1,  -4,  -1};

float const laplacian1KernelFactor = (float) 1.0;

int const laplacian2Kernel[] = { 0,  1,  0,
                                 1, -4,  1,
                                 0,  1,  0};
float const laplacian2KernelFactor = (float) 1.0;

int const laplacian3Kernel[] = { -1,  -1,  -1,
                                  -1,   8,  -1,
                                  -1,  -1,  -1};
float const laplacian3KernelFactor = (float) 1.0;


//Bonus Kernel:

int const gaussianKernel[] = { 1,  4,  6,  4, 1,
                               4, 16, 24, 16, 4,
                               6, 24, 36, 24, 6,
                               4, 16, 24, 16, 4,
                               1,  4,  6,  4, 1 };

float const gaussianKernelFactor = (float) 1.0 / 256.0;


// Helper function to swap bmpImageChannel pointers

void swapImageChannel(bmpImageChannel **one, bmpImageChannel **two) {
  bmpImageChannel *helper = *two;
  *two = *one;
  *one = helper;
}
//This method goes trough the 2d unsigned array and converts it into a 1d array
void loadChannelToBuffer(unsigned char **in,unsigned int width, unsigned int height, unsigned char *out)
{
  for (unsigned int y = 0; y < height; y++) {
    for (unsigned int x = 0; x < width; x++) {

        out[y*width+x]=in[y][x];
    }
  }
}
//This makes a 1d array into a 2d dimensional array
void loadBufferToChannel(unsigned char **in,unsigned int width, unsigned int height, unsigned char *out)
{
  for (unsigned int y = 0; y < height; y++) {
    for (unsigned int x = 0; x < width; x++) {
        in[y][x]=out[y*width+x];
    }
  }
}
// Apply convolutional kernel on image data
//The applyKernel method needs to see the ghos cells that are included as ghostTop and ghostBot
void applyKernel(unsigned char **out, unsigned char **in, unsigned int width, unsigned int height, int *kernel, unsigned int kernelDim, float kernelFactor, unsigned char *ghostTop, unsigned char *ghostBot, int my_id, int num_procs) {
  unsigned int const kernelCenter = (kernelDim / 2);
  for (unsigned int y = 0; y < height; y++) {
      for (unsigned int x = 0; x < width; x++) {
        int aggregate = 0;
        for (unsigned int ky = 0; ky < kernelDim; ky++) {
          int nky = kernelDim - 1 - ky;
          for (unsigned int kx = 0; kx < kernelDim; kx++) {
            int nkx = kernelDim - 1 - kx;

            int yy = y + (ky - kernelCenter);
            int xx = x + (kx - kernelCenter);
            if (xx >= 0 && xx < (int) width && yy >=0 && yy < (int) height)
            {
              aggregate += in[yy][xx] * kernel[nky * kernelDim + nkx];
            }
            //When an area outside of the regular image is checked, it uses the ghostBot value instead if its not process 0 and yy is less than 0
            else if(yy<0&&my_id!=0)
            {
                aggregate += ghostBot[xx] * kernel[nky * kernelDim + nkx];
            }
            //If yy is larger than the height and its not the last process then the ghostTop array is checked
            else if(yy==(int) height&&my_id!=num_procs-1)
            {
                aggregate += ghostTop[xx] * kernel[nky * kernelDim + nkx];
            }
          }
        }
        aggregate *= kernelFactor;

          if (aggregate > 0) {
            out[y][x] = (aggregate > 255) ? 255 : aggregate;
          } else {
            out[y][x] = 0;
          }

      }
    }
}


void help(char const *exec, char const opt, char const *optarg) {
    FILE *out = stdout;
    if (opt != 0) {
        out = stderr;
        if (optarg) {
            fprintf(out, "Invalid parameter - %c %s\n", opt, optarg);
        } else {
            fprintf(out, "Invalid parameter - %c\n", opt);
        }
    }
    fprintf(out, "%s [options] <input-bmp> <output-bmp>\n", exec);
    fprintf(out, "\n");
    fprintf(out, "Options:\n");
    fprintf(out, "  -i, --iterations <iterations>    number of iterations (1)\n");

    fprintf(out, "\n");
    fprintf(out, "Example: %s in.bmp out.bmp -i 10000\n", exec);
}

int main(int argc, char **argv) {
  /*
    Parameter parsing, don't change this!
   */
   /*
    Starts MPI and decleares some variables that are going to be used,
    remaining is the number of heights not divisiable by the number of processes remaining=height%num_procs
    int sum is used to give the scatterv sendcounts correct values,
    then some common empty pointer are also decleared these are only actually given values in the process 0
   */
  int remaining=0;
  int sum = 0;
  int sendheight=0;
  int ierr,num_procs, my_id;
  ierr = MPI_Init(&argc, &argv);
  ierr=MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  ierr=MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  bmpImageChannel *imageChannel;
  unsigned char *imageChannelBuffer;
  bmpImage *image;
  unsigned int iterations = 1;
  char *output = NULL;
  char *input = NULL;
  int ret = 0;
  int height, width;
  if(my_id==0)
  {




  static struct option const long_options[] =  {
      {"help",       no_argument,       0, 'h'},
      {"iterations", required_argument, 0, 'i'},
      {0, 0, 0, 0}
  };

  static char const * short_options = "hi:";
  {
    char *endptr;
    int c;
    int option_index = 0;
    while ((c = getopt_long(argc, argv, short_options, long_options, &option_index)) != -1) {
      switch (c) {
      case 'h':
        help(argv[0],0, NULL);
        goto graceful_exit;
      case 'i':
        iterations = strtol(optarg, &endptr, 10);
        if (endptr == optarg) {
          help(argv[0], c, optarg);
          goto error_exit;
        }
        break;
      default:
        abort();
      }
    }
  }

  if (argc <= (optind+1)) {
    help(argv[0],' ',"Not enough arugments");
    goto error_exit;
  }
  input = calloc(strlen(argv[optind]) + 1, sizeof(char));
  strncpy(input, argv[optind], strlen(argv[optind]));
  optind++;

  output = calloc(strlen(argv[optind]) + 1, sizeof(char));
  strncpy(output, argv[optind], strlen(argv[optind]));
  optind++;

  /*
    End of Parameter parsing!
   */

  /*
    Create the BMP image and load it from disk.
   */
   image = newBmpImage(0,0);
    if (image == NULL) {
      fprintf(stderr, "Could not allocate new image!\n");
    }

    if (loadBmpImage(image, input) != 0) {
      fprintf(stderr, "Could not load bmp image '%s'!\n", input);
      freeBmpImage(image);
      goto error_exit;
    }

    // Create a single color channel image. It is easier to work just with one color
    imageChannel = newBmpImageChannel(image->width, image->height);
    if (imageChannel == NULL) {
      fprintf(stderr, "Could not allocate new image channel!\n");
      freeBmpImage(image);
      goto error_exit;
    }

    // Extract from the loaded image an average over all colors - nothing else than
    // a black and white representation
    // extractImageChannel and mapImageChannel need the images to be in the exact
    // same dimensions!
    // Other prepared extraction functions are extractRed, extractGreen, extractBlue
    if(extractImageChannel(imageChannel, image, extractAverage) != 0) {
      fprintf(stderr, "Could not extract image channel!\n");
      freeBmpImage(image);
      freeBmpImageChannel(imageChannel);
      goto error_exit;
    }
    height=imageChannel->height;
    width=imageChannel->width;
    //The imageChannelBuffer is the array being scattered between the processes, by calling loadChannelToBuffer with the imageChannel->data, width, height and the imageChannelBuffer
    //The content in the image is transformed into a sendable 1d array
    imageChannelBuffer = calloc(sizeof(unsigned char)*width*height,1);
    loadChannelToBuffer(imageChannel->data,width,height,imageChannelBuffer);
}
  //This shares the widht, height and iterations information from process 0 to the other processes
  int *sendbuffer=calloc(sizeof(int)*3,1);
  if(my_id==0)
  {
    sendbuffer[0]=width;
    sendbuffer[1]=height;
    sendbuffer[2]=iterations;
    for(int i=0; i<num_procs;i++)
    {
      MPI_Send(sendbuffer,3,MPI_INT,i,0, MPI_COMM_WORLD);
    }
  }
  else
  {
    MPI_Recv(sendbuffer,3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    width=sendbuffer[0];
    height=sendbuffer[1];
    iterations=sendbuffer[2];
  }
  //sendbuffer is freed
  free(sendbuffer);



  //To use scatterv sendcounts and displs are needed, sendcounts are the number of data per process and displs are the displacement of the data
  int *sendcounts=malloc(sizeof(int)*num_procs);
  int *displs=malloc(sizeof(int)*num_procs);

  //the value sendheight is calculated this is the height of each imagepart that each process have
  if(height%num_procs==0)
  {
    sendheight=height/num_procs;
  }
  else
  {
    remaining=height%num_procs;
    sendheight=(height-remaining)/num_procs;
  }

  //Here sum is used to make sure that the displacement displs are set correctly and sendcounts are given the value of sendheight*width
  for (int i = 0; i < num_procs; i++) {
      sendcounts[i] = sendheight*width;
      displs[i] = sum;
      sum += sendcounts[i];
  }

  //The last process has a slightly larger image part and this is done by adding the remaining*width to the current sendcount value
  sendcounts[num_procs-1]+=remaining*width;

  //This variable is the number of elements per image array
  int workingsize=0;

  if(my_id!=num_procs-1)
  {
      workingsize=sendheight*width;
  }
  else
  {
    workingsize=(sendheight+remaining)*width;
  }
  //The workingheight is the height for each imagepart, the last process might have a higher workingheight than the others depending on the vale of remaining
  int workingheight=workingsize/width;

  //The reciverbuffer is used to recive the 1d array information from the process 0, and later used to send back the information with a MPI_Gatherv
  //As with the original code there are two imagechannels that are being used in the loop, the difference is that workingImageChannel will contain the image data
  //But is still empty
  unsigned char *reciverbuffer = calloc(sizeof(unsigned char)*workingsize,1);
  bmpImageChannel *workingImageChannel = newBmpImageChannel(width, workingheight);
  bmpImageChannel *processImageChannel = newBmpImageChannel(width, workingheight);
  //Mpi scatter is used to send the imageChannelBuffer data to all the processes with the values from sendcounts displs, reciverbuffer and workingsize
  MPI_Scatterv(imageChannelBuffer, sendcounts, displs, MPI_UNSIGNED_CHAR, reciverbuffer, workingsize, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
  //Then the workingImageChannel->data is loaded with the information from reciverbuffer that turns the 1d array into a 2d array
  loadBufferToChannel(workingImageChannel->data,width,workingheight,reciverbuffer);

  //Here we do the actual computation!
  // imageChannel->data is a 2-dimensional array of unsigned char which is accessed row first ([y][x])

  //These four buffers are used for sending a process outer imagedata and reciving the outer imagedata of other processes
  //They are all 1d array with the size of the image width
  unsigned char *ghostSendTop = calloc(width,1);
  unsigned char *ghostSendBot = calloc(width,1);
  unsigned char *ghostBufferTop= calloc(width,1);
  unsigned char *ghostBufferBottom=calloc(width,1);

  for (unsigned int i = 0; i < iterations; i ++) {
    //This fills in the information from the image into the send buffers
    for(int x=0;x<width;x++)
    {
      ghostSendTop[x]=workingImageChannel->data[workingheight-1][x];
      ghostSendBot[x]=workingImageChannel->data[0][x];
    }
    //ready to send
    /*
      When communicating the buffers the processes are split into odd and even processes
      The even numbered processes are sending first and then reciving, and the odd are reciving and sending

      The sending and reciving is locked in on direction at a time. First sending and reciving is done upwards,
      in the perspective of the even numbered processes (process 0 sending and reciving with process 1)

      After the upwards communication the downwards communication is done, (process 2 sending and reciving with process 1)
      The last process cant send or recive upwards
      And process 0 cant send or recive downwards
    */
    if(my_id%2==0)
    {
      if(my_id!=num_procs-1)
      {
        MPI_Send(ghostSendTop,width,MPI_UNSIGNED_CHAR,my_id+1,0, MPI_COMM_WORLD);
        MPI_Recv(ghostBufferTop,width, MPI_UNSIGNED_CHAR, my_id+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      if(my_id!=0)
      {
        MPI_Send(ghostSendBot,width,MPI_UNSIGNED_CHAR,my_id-1,0, MPI_COMM_WORLD);
        MPI_Recv(ghostBufferBottom,width, MPI_UNSIGNED_CHAR, my_id-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }
    else
    {
      MPI_Recv(ghostBufferBottom,width, MPI_UNSIGNED_CHAR, my_id-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(ghostSendBot,width,MPI_UNSIGNED_CHAR,my_id-1,0, MPI_COMM_WORLD);
      //After the top border has been exchanged next is bottom
      if(my_id!=num_procs-1)
      {
        MPI_Recv(ghostBufferTop,width, MPI_UNSIGNED_CHAR, my_id+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(ghostSendTop,width,MPI_UNSIGNED_CHAR,my_id+1,0, MPI_COMM_WORLD);
      }
    }
    applyKernel(processImageChannel->data,
                workingImageChannel->data,
                width,
                workingheight,
                (int *)laplacian1Kernel, 3, laplacian1KernelFactor,
                ghostBufferTop,ghostBufferBottom,my_id,num_procs
              );

    swapImageChannel(&processImageChannel,&workingImageChannel);



  }
  //The 4 buffers used to send information between processes are freed
  free(ghostSendBot);
  free(ghostSendTop);
  free(ghostBufferTop);
  free(ghostBufferBottom);
  //The 0 process allocates space to store the finished image
  //The values from workingImageChannel are loaded into the 1d array reciverbuffer
  loadChannelToBuffer(workingImageChannel->data,width,workingheight,reciverbuffer);
  //The bmpImageChannels are freed after use
  freeBmpImageChannel(processImageChannel);
  freeBmpImageChannel(workingImageChannel);
  //MPI gather is called gathering all the processes reciverbuffer information and putting it on process 0 imageChannelBuffer
  MPI_Gatherv(reciverbuffer, workingsize, MPI_UNSIGNED_CHAR, imageChannelBuffer, sendcounts, displs, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
  //reciverbuffer is freed
  free(reciverbuffer);
  //Mpi finalize is called by all process to ensure safe exit as only process 0 is actually exiting the program
  MPI_Finalize();
  if(my_id==0)
  {
    //the information from imageChannelBuffer is stored in the 2d array imageChannel->data
    loadBufferToChannel(imageChannel->data,width,height,imageChannelBuffer);
    //imageChannelBuffer is freed
    free(imageChannelBuffer);

  // Map our single color image back to a normal BMP image with 3 color channels
  // mapEqual puts the color value on all three channels the same way
  // other mapping functions are mapRed, mapGreen, mapBlue
  if (mapImageChannel(image, imageChannel, mapEqual) != 0) {
    fprintf(stderr, "Could not map image channel!\n");
    freeBmpImage(image);
    freeBmpImageChannel(imageChannel);
    goto error_exit;
  }
  freeBmpImageChannel(imageChannel);

  //Write the image back to disk
  if (saveBmpImage(image, output) != 0) {
    fprintf(stderr, "Could not save output to '%s'!\n", output);
    freeBmpImage(image);
    goto error_exit;
  };

graceful_exit:
  ret = 0;
error_exit:
  if (input)
    free(input);
  if (output)
    free(output);
  return ret;
}
};
