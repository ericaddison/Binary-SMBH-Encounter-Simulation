#include <stdio.h>
#include "BEMRI.h"

int getParamsFromFile(char * fname, iParams * iPars);

int args(int argc,char * argv[],iParams *iPars, flags * f, int * N, double * frac, long * seed, double * PbN, char * fname, double * ebh, int * ipsval, double * Hrat)
{

	f->cflag=0; f->sflag=0; f->tflag=0; f->ndflag=0; f->Nflag=0; f->soflag=0; f->PbNflag=0; f->th0flag=0; f->bs2flag = 0; f->fnameflag = 0; f->betaflag = 0; f->ipsflag = 0; f->ehflag = 0; 
	f->Htflag = 0; f->Hrat = 0; f->Xflag = 0; f->geoflag = 0; f->angflag = 0; f->gamflag = 0;
	f->fParflag = 0;

if(argc > 1)
{
	printf("\nFlags entered: \n");

	for(int ii = 1; ii < argc; ii+=1)
	{
		if(argv[ii][0] == '-' && isalpha(argv[ii][1]) != 0) 	//look for flags
		{
			if(strcmp(argv[ii],"-s")==0)		//if flag is -s for seed
			{
				errno = 0;
				*seed = strtol(argv[ii+1],NULL,0);
				if(errno != 0)
				{
					fprintf(stderr,"\nError in -s seed entry: %s\n",strerror(errno));
					return errno;
				}
				f->sflag = 1;
				printf("-s: Manual seed = %ld\n",*seed);
			}
			else if( strcmp(argv[ii],"-c") == 0 )	//if flag is -c for completion
			{
				if(f->tflag == 1)
				{
					fprintf(stderr,"\nError: Cannot use both -c and -t\n");
					return 1;
				}
				if(f->ndflag == 1)
				{
					fprintf(stderr,"\nError: Cannot use both -c and -nd\n");
					return 1;
				}
				f->cflag = 1;
				*frac = 1e300;
				printf("-c: Run simulations to completion\n");
			}
			else if( strcmp(argv[ii],"-t") == 0 )	//if flag is -t for tmax fraction
			{
				if(f->cflag == 1)
				{
					fprintf(stderr,"\nError in -t tmax entry: Cannot use both -c and -t\n");
					return 1;
				}
				if(argc == ii+1)
				{
					fprintf(stderr,"\nError in -t tmax entry: Invalid Entry\n");
					return 1;
				}
				*frac = strtod(argv[ii+1],NULL);
				if(*frac <= 0)
				{
					fprintf(stderr,"\nError in -t tmax entry: Invalid Entry\n");
					return 1;
				}
				f->tflag = 1;
				printf("-t: tmax = %.2f*Ph\n",*frac);
			}
			else if( strcmp(argv[ii],"-nd") == 0 )	//if flag is -nd for do not stop when disrupted
			{
				if(f->cflag == 1)
				{
					fprintf(stderr,"\nError: Cannot use both -c and -nd\n");
					return 1;
				}
				f->ndflag = 1;
				printf("-nd: Do not stop simulation if BERMI disrupts\n");
			}
			else if( strcmp(argv[ii],"-N") == 0 )	//if flag is -N for specifying number of runs
			{
				errno = 0;
				*N = strtol(argv[ii+1],NULL,0);
				if(errno != 0)
				{
					fprintf(stderr,"\nError in -N entry: %s\n",strerror(errno));
					return errno;
				}
				else if( f->fParflag)
				{
					fprintf(stderr,"\nError in -N entry: cannot use both -N and -fpars\n");
					return -1;
				}
				printf("-N: Run %d simulations -- no position or quadrupole output\n",*N);
				f->Nflag = 1;
			}
			else if( strcmp(argv[ii],"-so") == 0 )	//if flag is -so
			{
				printf("-so: Suppress all on-screen output\n");
				f->soflag = 1;
			}
			else if( strcmp(argv[ii],"-PbN") == 0 )	//if flag is -PbN for forcing Pb = fraction of tnode
			{
				errno = 0;
				*PbN = strtod(argv[ii+1],NULL);
				if(errno != 0)
				{
					fprintf(stderr,"\nError in -PbN entry: %s\n",strerror(errno));
					return errno;
				}
				printf("-Pnode: Pb = tnode * %.2f",*PbN);
				f->PbNflag = 1;
			}
			else if( strcmp(argv[ii],"-th0") == 0 )	//if flag is -th0 for setting theta0 = [value]
			{
				errno = 0;
				iPars->th0 = strtod(argv[ii+1],NULL);
				if(errno != 0)
				{
					fprintf(stderr,"\nError in -th0 entry: %s\n",strerror(errno));
					return errno;
				}
				printf("-th0: theta0 = %.2f",iPars->th0);
				f->th0flag = 1;
			}
			else if( strcmp(argv[ii],"-RK4") == 0 )	//if flag is -bs2 for using the Burlisch-Stoer integrator
			{
				printf("-RK4: use 4th order Runge-Kutta integrator\n");
				f->bs2flag = 1;
			}
			else if( strcmp(argv[ii],"-fname") == 0 ) //if -fname flag is used for specifying output filename
			{
				strcpy(fname,argv[ii+1]);
				printf("-fname: use output filename %s instead of BEMRI.dat\n",fname);
				f->fnameflag = 1;

			}
			else if( strcmp(argv[ii],"-beta") == 0 ) //if -beta flag is used to specify penetration ratio r_t/r_p
			{
				errno = 0;
				iPars->beta = strtod(argv[ii+1],NULL);
				iPars->gamma = 1.0/iPars->beta;
				if(errno != 0)
				{
					fprintf(stderr,"\nError in -beta entry: %s\n",strerror(errno));
					return errno;
				}
				printf("-beta: beta = %.2f\n",iPars->beta);
				f->betaflag = 1;

			}
			else if( strcmp(argv[ii],"-gam") == 0 ) //if -gam flag is used to specify penetration ratio r_p/r_t
			{
				errno = 0;
				iPars->gamma = 1.0/strtod(argv[ii+1],NULL);
				iPars->beta = 1.0/iPars->gamma;
				if(errno != 0)
				{
					fprintf(stderr,"\nError in -gam entry: %s\n",strerror(errno));
					return errno;
				}
				printf("-gam: gamma = %.2f\n",iPars->gamma);
				f->gamflag = 1;
				
			}
			else if( strcmp(argv[ii],"-ips") == 0 ) //if -ips flag is used to ask for initial phase sampling
			{
				errno = 0;
				*ipsval = (int)strtod(argv[ii+1],NULL);
				if(errno != 0)
				{
					fprintf(stderr,"\nError in -ips entry: %s\n",strerror(errno));
					return errno;
				}

				printf("-ips: sample evenly over %d initial phases\n",*ipsval);
				f->ipsflag = 1;

			}
			else if( strcmp(argv[ii],"-eh") == 0 ) //if -eh flag is used set BH orbit eccentricity
			{
				errno = 0;
				*ebh = strtod(argv[ii+1],NULL);
				if(errno != 0)
				{
					fprintf(stderr,"\nError in -eh entry: %s\n",strerror(errno));
					return errno;
				}
				printf("-eh: eh = %.4f\n",*ebh);
				f->ehflag = 1;

			}
			else if( strcmp(argv[ii],"-Htest") == 0 ) //if -Htest flag is entered to do Heggie test
			{
				printf("-Htest: Test Heggie values\n");
				f->Htflag = 1;
				
			}
			else if( strcmp(argv[ii],"-Hrat") == 0 ) //if -Hrat flag is entered to specify rp/a ratio in Heggie test
			{				
				errno = 0;
				*Hrat = strtod(argv[ii+1],NULL);
				if(errno != 0)
				{
					fprintf(stderr,"\nError in -Hrat entry: %s\n",strerror(errno));
					return errno;
				}
				printf("-Hrat: Use value of %.3g for rp/a ratio in Heggie test\n",*Hrat);
				f->Hrat = 1;
				
			}
			else if( strcmp(argv[ii],"-Xgrid") == 0 ) //if -X flag is used to specify running sims on Xgrid
			{
				printf("-Xgrid: Running on Xgrid\n");
				f->Xflag = 1;
			}		
			else if( strcmp(argv[ii],"-geo") == 0 ) //if -geo flag is use geometricized units
			{
				printf("-geo: Use geometricized units -- G = c = 1\n");
				f->geoflag = 1;
			}
			else if( strcmp(argv[ii],"-ang") == 0 ) //if -ang flag means randomize inclination and lan
			{
				printf("-ang: Randomized binary orientation angles\n");
				f->angflag = 1;
			}
			else if( strcmp(argv[ii],"-fpars") == 0 ) // if -fpars flag - use parameters from file in string
			{
				printf("-fpars: use parameters from filename %s\n",argv[ii+1]);
				printf("\toverriding any other parameter specifying flags\n");
				printf("\tsetting number of runs N=1\n");
				*N = 1;		// set for a single run
				if( getParamsFromFile(argv[ii+1],iPars))
					return 2;
				f->fParflag = 1;
			}
		}
	}
}
else
	printf("\nNo flags entered\n");

	printf("\n");
	return 0;
}

// get input parameters from file filename. Any unspecified params are set to default values.
int getParamsFromFile(char * fname, iParams * iPars)
{

	char check[80], str[80], ch;
	char ar[][10] = {"gamma","beta","inc","lan","th0"};
	FILE *fp;

	//check for file
	if( (fp = fopen(fname, "r")) == NULL)
	{
		fprintf(stderr,"\nError in getParamsFromFile: Cannot open file: %s\n",fname);
		return 1;
	}


	//check for proper header
	fscanf(fp,"%s",check);

	if(strcmp(check,"BEMRI_params") != 0)
	{
		fprintf(stderr,"Error in getParamsFromFile: Invalid BEMRI paramter file: %s\n",fname);
		fclose(fp);
		return 2;
	}

//file is OK, set default values
iPars->gamma = 1.0;
iPars->beta = 1.0;
iPars->incA = 0.0;
iPars->lanA = 0.0;
iPars->th0 = 0.0;


// start loading data
	while(fscanf(fp,"%s",str)==1)
	{
		//check for comment by //
		if(str[0] == '/' && str[1] == '/')
		{
			//printf("\ncomment found");
			//if comment is found, stall until finding a new line
			while( (ch = getc(fp)) != '\n')
			{}
		}

	//check for parameter identifier
		for(int j=0;j<5;j++)
		{
			if(strcmp(str,ar[j]) == 0)			//if parameter identifier found
			{
			    fscanf(fp,"%s",str);            //parameter was found, so scan the value into str

			    switch(j)
                {
                    case 0:                         //gamma
                        iPars->gamma = atof(str);
			iPars->beta = 1.0/iPars->gamma;
			printf("\tinput gamma = %.3f\n",iPars->gamma);
                    break;

                    case 1:                         //beta
                        iPars->beta= atof(str);
			iPars->gamma = 1.0/iPars->beta;
			printf("\tinput beta = %.3f\n",iPars->beta);
                    break;

                    case 2:                         //inc
			iPars->incA = atof(str);
			printf("\tinput inc = %.3f\n",iPars->incA);
                    break;

                    case 3:                         //lan
			iPars->lanA = atof(str);
			printf("\tinput lan = %.3f\n",iPars->lanA);
                    break;

                    case 4:                         //th0
			iPars->th0 = atof(str);
			printf("\tinput theta0 = %.3f\n",iPars->th0);
                    break;

                }
		}
		}
	}


	printf("\nData successfully loaded from file: %s\n\n",fname);

	fclose(fp);
	return 0;

}
