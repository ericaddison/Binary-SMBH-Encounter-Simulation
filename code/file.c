/*
	Eric Addison
	Utah State University
	Summer 2009

	File and string routines
*/

#include "file.h"
#include "kb.h"

//this routine will take in a file name and output the basename, i.e. strip off the extension

void get_basename(char * filename, char * base)
{

	int j=0;

	while( filename[j] != '.' )
	{
		base[j] = filename[j];
		j++;

		if(filename[j] == '\0')
		{
			printf("\nInvalid Filename\n");
			printf("\nPress any key to continue...\n");
			getch();
			return;
		}
	}

	return;
}



/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

//takes in the destination string base, the label of the series label, and the
//current index j
//writes "label00j" to base
//i.e. if label = "star" and j = 6, then written to base is: base = "star006"
void filename_index(char * base, char * label, int j)
{

		int size = sizeof(label);
		char num[3];


		strncpy(base,label,size);

		if(j < 10)
		{
			base[size] = '0';
			base[size+1] = '0';
			sprintf(num,"%d",j);
			base[size+2] = num[0];
		}
		else if(j < 100)
		{
			base[size] = '0';
			sprintf(num,"%d",j);
			base[size+1] = num[0];
			base[size+2] = num[1];
		}
		else if(j < 1000)
		{
			sprintf(num,"%d",j);
			base[size] = num[0];
			base[size+1] = num[1];
			base[size+2] = num[2];
		}
		else
			printf("\nIndex out of bounds\n");

		base[size+3] = '\0';


	return;
}


//attempts to open file
//ch holds r,w,or whatever fopen paramter

FILE* open_file(char * filename, char * ch)
{
    FILE *fp;
	printf("\nEnter file name: ");
	scanf("%s",filename);

	//check for file
	if( (fp = fopen(filename, ch)) == NULL)
	{
		printf("\nCannot open file: %s\n",filename);
		anykey();
		return NULL;
	}

	return fp;

}

FILE* open_file_no_ask(char * filename, char * op)
{
    FILE *fp;

	//check for file
	if( (fp = fopen(filename, op)) == NULL)
	{
		printf("\nCannot open file: %s\n",filename);
		anykey();
		return NULL;
	}

	return fp;

}



