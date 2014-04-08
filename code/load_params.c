/********************************************************************
*********************************************************************

    BEMRI Tidal Disruption Simulation

    Author: Eric Addison
    Initial Date: 8 Feb 2010
    Affiliation: USU Astrophysics

    File: load_params.c
    Purpose: Routine to load parameters from file.

*********************************************************************
*********************************************************************/

#include "BEMRI.h"


void load_params(struct BEMRI * b)
{

	char filename[80], check[80], str[80], ch;
	char ar[][12] = {"m1","m2","m3","e_b","e_h","P_b","P_h","ia","aop","lan","a_b","a_h"};
	int used_a = 0, used_P = 0;
	FILE *fp;

	//printf("\nEnter file name: ");
	//scanf("%s",filename);
    strcpy(filename,"params");

	//check for file
	if( (fp = fopen(filename, "r")) == NULL)
	{
		printf("\nCannot open file: %s\n",filename);
		anykey();
		return;
	}


	//check for proper header
	fscanf(fp,"%s",check);

	if(strcmp(check,"BEMRI_params") != 0)
	{
		printf("Invalid BEMRI paramter file: %s\n",filename);
		anykey();
		fclose(fp);
		return;
	}

//store basename of the file
	//get_basename(filename,base);

//file is OK, start loading data

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
		for(int j=0;j<12;j++)
		{
			if(strcmp(str,ar[j]) == 0)			//if parameter identifier found
			{
			    fscanf(fp,"%s",str);            //parameter was found, so scan the value into str

			    switch(j)
                {
                    case 0:                         //m1
                        m1 = atof(str)*MSUN;
                    break;

                    case 1:                         //m2
                        m2 = atof(str)*MSUN;
                    break;

                    case 2:                         //m3
                        m3 = atof(str)*MSUN;
                    break;

                    case 3:                         //e_b
                        eb = atof(str);
                    break;

                    case 4:                         //e_h
                        eh = atof(str);
                    break;

                    case 5:                          //P_b
                        Pb = atof(str);
                        used_P = 1;
                    break;

                    case 6:                         //P_h
                        Ph = atof(str);
                    break;

                    case 7:                         //ia
                        ia = atof(str)*PI/180;
                    break;

                    case 8:                         //aop
                        aop = atof(str)*PI/180;
                    break;

                    case 9:                         //lan
                        lan = atof(str)*PI/180;
                    break;

                    case 10:                         //a_b
                        ab = atof(str);
                        used_a = 1;
                    break;
                    case 11:                         //a_h
                        ah = atof(str);
                    break;
                }
			}
		}
	}


    if(used_P)
    {
    //calculate semi-major axes from Kepler III
        ab = pow((G*(m1+m2)*Pb*Pb/(4*PI*PI)),0.3333333333);

        //for now, change to P_tid, tidal disrution period
        //double rp = ab*(1+eb);
        //Ph = 2*PI*pow( (2*rp*rp*rp/ ( (1-eh)*(1-eh)*(1-eh) * G * m1)),0.5)*2 ;
        //printf("\nClose to disruption radius, Ph = %.3e\n",Ph);

        ah = pow((G*(m1+m2+m3)*Ph*Ph/(4*PI*PI)),0.3333333333);
    }
    else if(used_a)
    {
    //calculate Orbital periods from Kepler III
        Pb = 2*PI*sqrt(ab*ab*ab/(G*(m1+m2)));
        Ph = 2*PI*sqrt(ah*ah*ah/(G*(m1+m2+m3)));
    }


    b->binary_b.M = m1+m2;
    b->binary_b.mu = m1*m2/(m1+m2);
    b->binary_h.mass2 = m1+m2;
    b->binary_h.M = m1+m2+m3;
	b->binary_h.mu = m3*(m1+m2)/(m1+m2+m3);


printf("\nData successfully loaded from file: %s\n\n",filename);
anykey();

fclose(fp);
return;

}
