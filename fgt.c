#include <stdlib.h>
#include <stdio.h>

typedef struct {
	char 			filename[128];
	unsigned long 	allele;
} Sample;

int uniquePairs(unsigned long p1, unsigned long p2)
{
	unsigned char uniques[4] = { 0, 0, 0, 0 };

	register unsigned long bitMask = 1UL;
	for (unsigned short i = 0; i < 24; i++)
	{
		register unsigned long p1_bit = p1 & bitMask;
		register unsigned long p2_bit = p2 & bitMask;

		if (p1_bit && p2_bit) /*if both are 1*/
			uniques[0] = 1;
		else if (p1_bit)	/*if p1 is one and p2 is 0*/
			uniques[1] = 1;
		else if (p2_bit)	/*if p2 is 1 and p1 is 0*/
			uniques[2] = 1;
		else				/*both are 0 (all the others are false)*/
			uniques[3] = 1;

		bitMask <<= 1;
	}

	return uniques[0] + uniques[1] + uniques[2] + uniques[3];
}

int main(int argc, char* argv[])
{
    char const* const fileName = argv[1]; /* should check that argc > 1 */
	FILE* file = fopen(fileName, "r"); /* should check the result */
	char line[256];
	char filename[128];
	char numbers[24];

	Sample *loci;

	loci = calloc(200000, sizeof(Sample));


	unsigned long ns = 0;
	while (fgets(line, sizeof(line), file)) 
	{
		/* note that fgets don't strip the terminating \n, checking its
		   presence would allow to handle lines longer that sizeof(line) */

		sscanf(line, 
			"%s " 
			"%c %c %c %c "
			"%c %c %c %c "
			"%c %c %c %c "
			"%c %c %c %c "
			"%c %c %c %c "
			"%c %c %c %c", 
			loci[ns].filename, 
			&numbers[0],
			&numbers[1],
			&numbers[2],
			&numbers[3],
			&numbers[4],
			&numbers[5],
			&numbers[6],
			&numbers[7],
			&numbers[8],
			&numbers[9],
			&numbers[10],
			&numbers[11],
			&numbers[12],
			&numbers[13],
			&numbers[14],
			&numbers[15],
			&numbers[16],
			&numbers[17],
			&numbers[18],
			&numbers[19],
			&numbers[20],
			&numbers[21],
			&numbers[22],
			&numbers[23]);

		/* convert numbers into an integer */
		unsigned long n = 0;

		loci[ns].allele = 0;

		for (short i = 0; i < 24; i++)
		{
			loci[ns].allele <<= 1;
			loci[ns].allele |= numbers[i] == '0' ? 0UL : 1UL;
		}

		/*
		printf("[%s] - %s - 0x%x\n", line, 
				loci[ns].filename, loci[ns].allele); 
		*/
		ns++;
	}

	/* may check feof here to make a difference between 
	   eof and io failure -- network
	   timeout for instance */

	fclose(file);

	/* loci are in loci[] ... run a double loop to compare */
	/* ns holds the number of loci read */
	
	unsigned long uniques = 0;
	fprintf(stderr, "Processing %ld loci\n", ns);
	for (unsigned long int row = 0; row < ns - 1; row++)
	{
		if (row % 100 == 0)
			fprintf(stderr, ".");

		for (unsigned long int cmp_row = row + 1; cmp_row < ns; cmp_row++)
		{
			if (uniquePairs(loci[row].allele, loci[cmp_row].allele) == 4)
			{
				uniques++;
				printf("%s %s\n", 
						loci[row].filename, loci[cmp_row].filename);
			}
		}
	}

	fprintf(stderr, "\nFound %ld allele matches\n", uniques);

    return 0;

}
