#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif


int main(int argc, char *argv[])
	{
	FILE *file = '\0', *outfile = '\0', *infile2 = '\0', *nstatesfile = '\0';
	char command[1000], c = '\0', codon[3], string[1000], **alignment ='\0', **names = '\0', *chars = '\0';
	int chosentaxa = -1, all=TRUE, bysite=FALSE, bytaxa=FALSE, error=FALSE, nchars =0, nseqs=0, maxnamelen=0, i=0, j=0, x=0, y=0, nummissing=0, num_undercutoff=0,  *include='\0', *states='\0', excluded_sites=0, *num_times_masked = '\0';
	float *nstates = '\0', con_cutoff = 0;

	string[0] = '\0';
	if(argc < 3)
		{
		printf(" the usage of this program is:\n\n\tstrip_sitesunknown sequence.file taxa\n\n\t'sequence.file' the alignment that is to be stripped in fasta format\ntaxa, is the taxa to be used to identify the sites to be stripped by (if it has a '-' or '?'\n\n"); 
		exit(1);
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
		strcat(command, ".StripBy");
		strcat(command, argv[2]);
		
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
		for(i=0; i<nseqs; i++)
			{
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
				if(strcmp(argv[2], names[i]) == 0) chosentaxa=i;
				j=0;
				}
			if(c != '\n' && c != '\r' && c != '\t' && c != ' ')
				{
				alignment[i][j] = c;
				j++;
				}
			}
		
		
		if(chosentaxa == -1) 
			{
			printf("ERROR: taxa %s not found in alignment\n", argv[2]);
			exit(1);
			}
		outfile = fopen(command, "w");
		/* Now scan throuigh the alignment excluding sites as specified by the user */
		
		for(i=0; i<nchars; i++)
			{
			if(alignment[chosentaxa][i] == '-' || alignment[chosentaxa][i] == '?' ) include[i] = FALSE;
			}
		
		
		/* Now print out the parsimony infomrative alignment */
		for(j=0; j<nseqs; j++)
			{
			fprintf(outfile, ">%s\n", names[j]);
			for(i=0; i<nchars; i++)
				{
				if(include[i] )
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
		
		printf("New Alignment length: %d\n", nchars-(excluded_sites/nseqs));

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


		fclose(file);
		fclose(outfile);
		}

	}

