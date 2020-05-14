#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif


/* This array represents the amino acids as they are defined in the genetic_codes array. The Amino Acids are numbered from 0 to 21 in the following order; */
/* Stop = X, Phe = F, Trp = W, Tyr = Y, His = H, Met = M, Leu = L, Ile = I, Val = V, Pro = P, Cys = C, Ala = A, Gly = G, Thr = T, Ser = S, Gln = Q, Asn = N, Lys = K, Arg = R, Glu = E, Asp = D, Gap = -.										   */								
char amino_acids[22] = {
  'X', 'F', 'W', 'Y', 'H', 'M', 'L', 'I', 'V', 'P', 'C', 'A', 'G', 'T', 'S', 'Q', 'N', 'K', 'R', 'E', 'D', '-'};

/* this array contains the amino acid categories as defined in "Prediction of protein secondary structure and active sites using the alignment of homologous sequences",
    M.J> Zvelebil, G.J. Barton, W.R. Taylor, and Sternberg, (J. Mol. Biol., 195, p957 - 1987)
    Wherever a particular amino acid belongs to a category, there is a 1, otherwise a 0 (or TRUE and FALSE), The order of the amino acids is the same as the array
    amino_acids and so stop and gaps are included even though they don't belong to any category, this is to ensure compatibility
    There are nine categories listed and are in the following order:
    hydrophobic, positive, negative, polar, charged, small, tiny, aromatic, aliphatic */
int AA_categories[22][9] = {
	0,0,0,0,0,0,0,0,0,  /*X */
	1,0,0,0,0,0,0,1,0,	/*F */
	1,0,0,1,0,0,0,1,0,	/*W */
	1,0,0,1,0,0,0,1,0,	/*Y */
	1,1,0,1,1,0,0,1,0,	/*H */
	1,0,0,0,0,0,0,0,0,	/*M */
	1,0,0,0,0,0,0,0,1,	/*L */
	1,0,0,0,0,0,0,0,1,	/*I */
	1,0,0,0,0,1,0,0,1,	/*V */
	0,0,0,0,0,1,0,0,0,	/*P */
	1,0,0,0,0,1,0,0,0,	/*C */
	1,0,0,0,0,1,1,0,0,	/*A */
	1,0,0,0,0,1,1,0,0,	/*G */
	1,0,0,1,0,1,0,0,0,	/*T */
	0,0,0,1,0,1,1,0,0,	/*S */
	0,0,0,1,0,0,0,0,0,	/*Q */
	0,0,0,1,0,1,0,0,0,	/*N */
	1,1,0,1,1,0,0,0,0,	/*K */
	0,1,0,1,1,0,0,0,0,	/*R */ 
	0,0,1,1,1,0,0,0,0,	/*E */
	0,0,1,1,1,1,0,0,0,	/*D */
	0,0,0,0,0,0,0,0,0  }; /*- */
    



int main(int argc, char *argv[])
	{
	FILE *file = '\0', *outfile = '\0', *infile2 = '\0', *nstatesfile = '\0', *meanidfile = '\0', *tmpout = '\0';
	char command[1000], c = '\0', codon[3], string[1000], **alignment ='\0', **names = '\0', *chars = '\0';
	int found=FALSE, first=0, second=0, con_method = -1, all=TRUE, bysite=FALSE, bytaxa=FALSE, error=FALSE, nchars =0, nseqs=0, maxnamelen=0, i=0, j=0, x=0, y=0, z=0, nummissing=0, num_undercutoff=0,  *include='\0', *states='\0', excluded_sites=0, *num_times_masked = '\0', numparsuninf = 0;
	float *nstates = '\0', *meanid = '\0', con_cutoff = 0, euc_sum=0, euc_count=0, tmpcount=0, value=0, overallmeanid=0;

	string[0] = '\0';
	if(argc < 3)
		{
		printf(" the usage of this program is:\n\n\tstrip_parsinform sequence.file <bysite|bytaxa|all> <gt|lt> <cutoff> <nstates|aacats|tuples|seqdiff> <seqdiff file name>\n\n\t'sequence.file' the alignment that is to be stripped in fasta format\n\t'Options:\n\t\t\"bysite\": Strip parsimony uninformative sites from the alignment\n\t\t\"bytaxa\": Mask parsimony uninformative taxa within sites\n\t\t\"all\": Strip parsimony uninformative sites and mask uninformative taxa in remaining sites (default)\n\n\t<gt|lt>: specifies whether the sites that are less then (lt) or greater then (gt) the conservation cutoff should be stripped\n\tConservation cutoff: An option which removes sites according to the relative conservation.\n\n\t<nstates|aacats> refer to the method to determine the conservation at each site.\n\t\tnstates is a value representing the conservation based on the number of states; \n\t\taacats uses the AA categories to determine the average difference in biochemical proerty represented by the different states at this site.\n\n\ttuples identifies all sites which have instances of characters that appear only twice (ie tuples)\n\n\tmeanid causes the program to use the meanID calculated externally to be used to rank sites by the hetergeneity\n\n\t<mean ID file> Is a file containing the names of each taxon in the sequence file followed by the mean percentage Identity of the taxa to all others\n\n"); 
		exit(1);
		}
	
	if(argc > 3)
		{
		con_cutoff = atof(argv[4]);
		/*if( con_cutoff < 0 || con_cutoff > 1)
			{
			printf("ERROR: conservation cutoff should be between 0 and 1\n");
			error = TRUE;
			}
	*/	}

	con_method = 0; /* by default uses the number of states "nstates" to calculate the conservation at each site*/
	if(argc > 5)
		{
		if(strcmp(argv[5], "nstates") == 0)
			{
			con_method=0;
			}
		else
			{
			if(strcmp(argv[5], "aacats") == 0)
				{
				con_method=1;
				}
			else
				{
				if(strcmp(argv[5], "tuples") == 0)
					{
					con_method=2;
					}
				else
					{
					if(strcmp(argv[5], "seqdiff") == 0)
						{
						con_method=3;
						if(argc == 7)
							{
							if((meanidfile = fopen(argv[6], "r")) == '\0')		/* check to see if the file is there */
								{								/* Open the file */
								printf("Error: Cannot open file %s\n", argv[1]);
								error=TRUE;
								}
							}
						else
							{
							printf("ERROR: you must provide a file containing the mean percentage Identities of each taxon if using the meanid method\n");
							error=TRUE;
							}
						}
					else
						{
						printf("ERROR: option %s does not exist, please use \"nstates\" or \"aacats\"\n");
						error=TRUE;
						}
					}
				}
			}
		}
		
	
	
	
	if(strcmp(argv[2], "all") == 0)
		{
		all = TRUE;
		}
	else
		{
		if(strcmp(argv[2], "bysite") == 0)
			{
			bysite = TRUE;
			all=FALSE;
			}
		else
			{
			if(strcmp(argv[2], "bytaxa") == 0)
				{
				bytaxa=TRUE;
				all=FALSE;
				}
			else
				{
				printf("ERROR: option '%s' does not exist, please use either 'bytaxa', 'bysite' or 'all'\n", argv[2]);
				exit(1);
				}
			}
		}

		

	if(!error)
		{
		 /*Open the sequence file  */
		 if((file = fopen(argv[1], "r")) == '\0')		/* check to see if the file is there */
				{								/* Open the file */
				printf("Error: Cannot open file %s\n", argv[1]);
				exit(1);
				}

		command[0] = '\0';
		strcpy(command, argv[1]);
		strcat(command, ".strip");
		outfile = fopen(command, "w");
		i=0; j=0;

		

		/* scan through the alignment checking that they are all the same length and counting the size of the necessary matrices */
		while(!feof(file))
			{
			c=getc(file);
			if(c == '>')
				{
				nseqs++;
				if( nchars == 0 ) nchars = i;		
				if(nchars != i)
					{
					printf("Error: The sequences are not all the same length! %d != %d @ seqnum: %d\n", i, nchars, nseqs );
					fclose(file);
					exit(1);
					}
				i=0; j=0;
				while(!feof(file) && (c=getc(file)) != '\n' && c != '\r')j++;
				if(j>maxnamelen) maxnamelen=j;
				}
			if(!feof(file) && c != '\n' && c != '\r' && c != '\t' && c != ' ' && c !='>') i++;
			}
		if(nchars != i && nchars != 0)
			{ 
			printf("Error: The sequences are not all the same length! %d != %d @ seqnum: %d\n", i, nchars, nseqs );
			fclose(file);
			exit(1);
			}
		rewind(file);

		/* create the array to hold the alignment and names*/
		printf("number of sequences %d\nOriginal alignment length = %d\nMax name length = %d\n", nseqs, nchars, maxnamelen);
		names=malloc(nseqs*sizeof(char *));
		alignment=malloc(nseqs*sizeof(char*));
		nstates=malloc(nchars*sizeof(float));
		include=malloc(nchars*sizeof(int));
		states=malloc(nseqs*sizeof(int));
		chars=malloc(nseqs*sizeof(char));
		num_times_masked=malloc(nseqs*sizeof(int));
		meanid=malloc(nseqs*sizeof(float));
		for(i=0; i<nseqs; i++)
			{
			meanid[i] = 0;
			states[i] = 0;
			num_times_masked[i] = 0;
			chars[i] = '\0';
			names[i]=malloc((maxnamelen+1)*sizeof(char));
			names[i][0] = '\0';
			alignment[i] = malloc((nchars+1)*sizeof(char));
			alignment[i][0] = '\0';
			}
		for(i=0; i<nchars; i++) 
			{
			include[i] = TRUE;
			nstates[i] = 0;
			}
		
		/*read throuugh the alignment file again, reading in everything */
		i=-1; j=0;
		while(!feof(file))
			{
			c=getc(file);
			if(c == '>')
				{
				if(i!=-1) alignment[i][j] = '\0';
				i++;
				j=0;
				while(!feof(file) && (c=getc(file)) != '\n' && c != '\r')
					{
					names[i][j] = c;
					j++;
					}
				names[i][j] = '\0';
				j=0;
				}
			if(c != '\n' && c != '\r' && c != '\t' && c != ' ')
				{
				alignment[i][j] = c;
				j++;
				}
			}
		
		/* if using the meanid method, read in the mean ids calculated externally */
		if(con_method==3)
			{
			j=0;
			while(!feof(meanidfile))
				{
				string[0] = '\0';
				fscanf(meanidfile, "%s\t%f\n", string, &value);
				/*printf("%s\t%f\n", string, value); */
				found=FALSE;
				for(i=0; i<nseqs; i++)
					{
					if(strcmp(string, names[i]) == 0)
						{
						meanid[i] = value;
						found=TRUE;
						i=nseqs;
						}
					}
				if(found==FALSE)
					{
					printf("ERROR: Taxa %s in meanID file cannot be found in the alignment... exiting\n", string);
					error=TRUE;
					}
				j++;
				}
			if(j<nseqs)
				{
				printf("ERROR: there were fewer taxa listed in the meanid file than in the alignment file....exiting\n");
				error=TRUE;
				}
			}
						
				
		if(con_method==3)
			{
			/* calculate the overall average mean % ID for all taxa */
			overallmeanid=0;
			for(j=0; j<nseqs; j++)
				overallmeanid+=meanid[j];

			overallmeanid=overallmeanid/nseqs;
			printf("overallmeanid = %f\n", overallmeanid); 
			}
			
		
		if(error == FALSE)
			{
			/* Now scan throuigh the alignment excluding sites as specified by the user */
			
			for(i=0; i<nchars; i++)
				{
				for(j=0;j<nseqs; j++)
					{
					states[j]=0;
					chars[j] = '\0';
					}
				/* count the number of each character in this position of the alignment */
				nummissing=0;
				for(j=0; j<nseqs; j++)
					{
					x=0;
					while(x < nseqs && chars[x] != '\0' && chars[x] != alignment[j][i]) x++;
					if(chars[x] == '\0')
						{
						chars[x] = alignment[j][i];
						states[x]++;
						}
					else
						states[x]++;
						
					/* count the number of sequences with missing/indel data at this site */
					if(alignment[j][i] == '-' || alignment[j][i] == '?' )
						{
						nummissing++;
						}
					
					}
				
				/* Depending on the option chosen by the user, identify the appropritate sites etc */
				if(bytaxa || all ) /* mask sites with a '?' which are singletons at this position in the alignment */
					{
					for(j=0; j<nseqs; j++)
						{
						if(states[j] == 1 && chars[j] != '-' && chars[j] != '?') /* this is a singleton, so mask it and change the details in the states array */
							{
							/* search for this character in this opsition in the alinment and mask it*/
							for(x=0; x<nseqs; x++)
								{
								if(alignment[x][i] == chars[j])
									{
									alignment[x][i] = '?';
									num_times_masked[x]++;
									x = nseqs;
									states[j] = 0;
									}
								}
							}
						}
					}
				
				
				/* calculate the conservation sttistic for this position*/
				if(con_method == 0)
					{
					/* count the numbers of different characters at this position (not counting gaps or missing data) */
					x=0;
					while(x< nseqs && chars[x] != '\0')
						{
						if(chars[x] != '-' && chars[x] != '?') nstates[i]++;
						x++;
						}
					if((nseqs-nummissing) > 0)
						{
						nstates[i] = 1-(nstates[i]/(nseqs-nummissing));	
						}
					else
						nstates[i] =0;
					}
				if(con_method == 1)
					{
					/* use the difference in the biochemical proerties of the AA at this site to calculate the conservation */
					/* we will calculatwe the euclidean distance between all possibile non-gap states to calculate the average distance for this site*/
					x=0;
					euc_sum=0;
					euc_count=0;
					while(x< nseqs && chars[x] != '\0')
						{
						/*identify this AA character */
						for(first=0; first<22; first++)
							{	
							if(amino_acids[first] == chars[x]) break;
							}
						if(chars[x] != '-' && chars[x] != '?' && first != 22)
							{
							y=x+1;
							while(y< nseqs && chars[y] != '\0')
								{
								/*identify this AA character */
								for(second=0; second<22; second++)
									{	
									if(amino_acids[second] == chars[y]) break;
									}
								if(chars[y] != '-' && chars[y] != '?' && second != 22)
									{
									tmpcount=0;
									for(z=0; z<9; z++) /* gop through all 9 AA categories */
										{
										tmpcount+=pow((AA_categories[first][z] - AA_categories[second][z]),2);
										}
									euc_sum+=sqrt(tmpcount);
									euc_count++;
									}
								y++;
								}
							}
						x++;
						}
					if(euc_sum > 0)
						{
						nstates[i] = euc_sum/euc_count;
						}
					else
						nstates[i] = 0;
					}
					
				if(con_method == 2) /* identify all sites with tuples */
					{
					found=FALSE;
					for(x=0; x<nseqs; x++)
						{
						if(states[x] != '?' && states[x] != '-' && states[x] <= 1 && states[x] > 0) found=TRUE;
						}
						
					if(found)
						nstates[i] = 1;
					else
						nstates[i] = 0;
					
					
					
					}
				
				if(con_method == 3) /* use species differences as part of calculation of site heterogeneity */
					{

					/* for each character calculate the average differences of all the taxa with that character and sum for all characters*/
					
					nstates[i] =0;
					y=0;
					while(y< nseqs && chars[y] != '\0')/* scan through the saved characters */
						{
						z=0;
						tmpcount=0;
						if(chars[y] != '-' && chars[y] != '?')
							{
							for(j=0; j<nseqs; j++)
								{
								if(alignment[j][i] == chars[y])
									{
									z++;
									tmpcount+=meanid[j];
									}
								}
							if(z>0)
								{
								nstates[i]=nstates[i]+(tmpcount/z); /* calculate the average % ID for all taxa containing this character and add to the total so far*/
								}
							}
						y++;
						}
					}
					
					
				

				if(bysite || all) /* identify all sites that are parsimony uninformative */
					{
					x=0;
					for(j=0; j<nseqs; j++)
						{
						if(chars[j] != '-' && chars[j] != '?' && states[j] > 1 ) x++;
						}
					if(x<2) /* if there are fewer than 2 characters with a least two instances, this site is not parsimony informative */
						{
						include[i] = FALSE;
						numparsuninf++;
						}
					}
				}
			
			if(argc>3) /* include sites based on parsimony informativeness and also based on rate calculation of some kind */
				{
				/* Now print out the parsimony infomrative alignment */
				if(strcmp(argv[3], "lt") == 0 ) /* we are to strip the sites that are less than the specified cutoff */
					{
					for(j=0; j<nseqs; j++)
						{
						fprintf(outfile, ">%s\n", names[j]);
						for(i=0; i<nchars; i++)
							{
							if( include[i] && nstates[i] >= con_cutoff)
								{
								fprintf(outfile, "%c", alignment[j][i]);
								}
							else
								{
								excluded_sites++;
								if(nstates[i] < con_cutoff) num_undercutoff++;
								}
							}
						fprintf(outfile, "\n");
						}
					}
				else /* we are to strip the sites that are greate rthan the specified cutoff */
					{
					for(j=0; j<nseqs; j++)
						{
						fprintf(outfile, ">%s\n", names[j]);
						for(i=0; i<nchars; i++)
							{
							if( include[i] && nstates[i] <= con_cutoff)
								{
								fprintf(outfile, "%c", alignment[j][i]);
								}
							else
								{
								excluded_sites++;
								if(nstates[i] > con_cutoff) num_undercutoff++;
								}
							}
						fprintf(outfile, "\n");
						}
					}
				}
			else  /* Just include sites based on whether they are parsimony informative  */
				{
				for(j=0; j<nseqs; j++)
					{
					fprintf(outfile, ">%s\n", names[j]);
					for(i=0; i<nchars; i++)
						{
						if( include[i] )
							{
							fprintf(outfile, "%c", alignment[j][i]);
							}
						else
							{
							excluded_sites++;
							}
						}
					fprintf(outfile, "\n");
					}	
				}
			
			/* now print out the information of number of states at each site */
			strcpy(command, argv[1]);
			strcat(command, ".nstates");
			nstatesfile = fopen(command, "w");
			for(i=0; i<nchars; i++)
				{
				fprintf(nstatesfile, "%f", nstates[i]);
				if(!include[i]) fprintf(nstatesfile, "\tUninformative");
				fprintf(nstatesfile, "\n");
				
				}
			fclose(nstatesfile);
			nstatesfile='\0';

			
			
			if(bysite || all)printf("Number of parsimony uniformative sites excluded = %d\n", numparsuninf);
			if(argc>3)printf("Number of sites %s cutoff of %f = %d\n", argv[3], con_cutoff, num_undercutoff/nseqs);
			printf("New Alignment length: %d\n", nchars-(excluded_sites/nseqs));
			printf("Site Masking report:\n\tName\tNum masked Sites\n");
			for(i=0; i<nseqs; i++)
				printf("\t%s\t%d\n", names[i], num_times_masked[i]);
			
			}


			
		/* clean up memory */
		
		for(i=0; i<nseqs; i++)
			{
			free(names[i]);
			names[i] = '\0';
			free(alignment[i]);
			alignment[i] = '\0';
			}
		free(names);
		free(alignment);
		free(include);
		free(nstates);
		free(states);
		free(chars);
		free(num_times_masked);
		free(meanid);


		fclose(file);
		fclose(outfile);
		if(con_method==3)fclose(meanidfile);
		}

	}

