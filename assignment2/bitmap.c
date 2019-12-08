#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "bitmap.h"

// save 24-bits bmp file, buffer must be in bmp format: upside-down
void savebmp(char *name,uchar *buffer,int x,int y) {
	FILE *f=fopen(name,"wb");
	if(!f) {
		printf("Error writing image to disk.\n");
		return;
	}
	unsigned int size=x*y*3+54;
	uchar header[54]={'B','M',size&255,(size>>8)&255,(size>>16)&255,size>>24,0,
                    0,0,0,54,0,0,0,40,0,0,0,x&255,x>>8,0,0,y&255,y>>8,0,0,1,0,24,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	fwrite(header,1,54,f);
	fwrite(buffer,1,x*y*3,f);
	fclose(f);
}

// read bmp file and store image in contiguous array
void readbmp(char* filename, uchar* array) {
	FILE* img = fopen(filename, "rb");   //read the file
	uchar header[54];
	fread(header, sizeof(uchar), 54, img); // read the 54-byte header

  // extract image height and width from header
	int width = *(int*)&header[18];
	int height = *(int*)&header[22];
	int padding=0;
	while ((width*3+padding) % 4!=0) padding++;

	int widthnew=width*3+padding;
	uchar* data = calloc(widthnew, sizeof(uchar));

	for (int i=0; i<height; i++ ) {
		fread( data, sizeof(uchar), widthnew, img);
		for (int j=0; j<width*3; j+=3) {
			array[3 * i * width + j + 0] = data[j+0];
			array[3 * i * width + j + 1] = data[j+1];
			array[3 * i * width + j + 2] = data[j+2];
		}
	}
	fclose(img); //close the file
}
void alter(uchar* array,uchar* image,int size)
{
    /*
    Starts by declaring a two times as large empty uchar array that can be used to fill inn the information from the image
    the variables n and t are used to write double horizontal lines
    */
		int width=2560;
    int n=0;
    int t=0;

    /*
    A for loop looping over the array the array is constructed, the array is a quarter the size of our new image
    for each iteration three values are read from the array and adjusted the values are then placed twelve different
    places, each "rgb block" is placed in the relative positions (0,0) (0,1) (1,0) and (1,1) from the array
    ie if the array was in 2d it would look like this array[x,y]=image[x,y] image[x+1,y] image[x,y+1] image[x+1,y+1]
    there are some adjustments to do 2d operations on the 1d array
    */
    for(int i=0;i<size;i+=3)
    {

            image[2*i+2*width*3*n]=256-array[i];
            image[2*i+1+2*width*3*n]=256-array[i+1];
            image[2*i+2+2*width*3*n]=256-array[i+2];

            image[2*i+3+2*width*3*n]=256-array[i];
            image[2*i+4+2*width*3*n]=256-array[i+1];
            image[2*i+5+2*width*3*n]=256-array[i+2];

            image[2*i+2*width*3*(n+1)]=256-array[i];
            image[2*i+1+2*width*3*(n+1)]=256-array[i+1];
            image[2*i+2+2*width*3*(n+1)]=256-array[i+2];

            image[2*i+3+2*width*3*(n+1)]=256-array[i];
            image[2*i+4+2*width*3*(n+1)]=256-array[i+1];
            image[2*i+5+2*width*3*(n+1)]=256-array[i+2];

            t++;
            if(t==2560)
            {
                n+=1;
                t=0;
            }
    }
}
