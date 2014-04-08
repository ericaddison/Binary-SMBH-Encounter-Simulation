/*
	Eric Addison
	Utah State University
	Summer 2009

	File and string routines header
*/


#include <stdio.h>
#include <stdlib.h>		//is where atof() is
#include <string.h>
#include <math.h>

void get_basename(char * filename, char * base);
void filename_index(char * base, char * label, int j);
FILE* open_file(char * filename, char * ch);
FILE* open_file_no_ask(char * filename, char * ch);
