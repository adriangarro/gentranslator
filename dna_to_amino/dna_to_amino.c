/* 
	--------------------------------------------------------------
	dna_to_amino.c
	
	DNA sequence to Amino acid sequences.
	
	Notes: 1. The file's name that contains the DNA sequence 
	          must be not longer of 79 chars.
	       2. The DNA sequence must contain a maximum 
	          of 4095 chars.
	
	(c) E. Adrian Garro S, Giaccomo Ubaldo P. 
		Costa Rica Institute of Technology. 
		February 2017.
	--------------------------------------------------------------
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#define NELEMS(arr) (sizeof(arr) / sizeof((arr)[0]))

static bool sorted = false;

/* Genetic code */
struct code {
	const char *codon;
	const char aminoacid;
};

static struct code codes[] = {
	{"UUU", 'F'},	/* F - Phenylalanine */
	{"UUC", 'F'},
	{"UUA", 'L'},	/* L - Leucine */ 
	{"UUG", 'L'}, 
	{"CUU", 'L'}, 
	{"CUC", 'L'}, 
	{"CUA", 'L'}, 
	{"CUG", 'L'},
	{"AUU", 'I'},	/* I - Isoleucine */ 
	{"AUC", 'I'}, 
	{"AUA", 'I'},
	{"AUG", 'M'},	/* M - Methionine */
	{"GUU", 'V'},	/* V - Valine */ 
	{"GUC", 'V'},
	{"GUA", 'V'}, 
	{"GUG", 'V'},
	{"UCU", 'S'},	/* S - Serina */
	{"UCC", 'S'}, 
	{"UCA", 'S'}, 
	{"UCG", 'S'}, 
	{"AGU", 'S'}, 
	{"AGC", 'S'},
	{"CCU", 'P'},	/* P - Proline */
	{"CCC", 'P'},
	{"CCA", 'P'}, 
	{"CCG", 'P'},
	{"ACU", 'T'},	/* T - Threonine */ 
	{"ACC", 'T'},
	{"ACA", 'T'}, 
	{"ACG", 'T'},
	{"GCU", 'A'},	/* A - Alanina */ 
	{"GCC", 'A'}, 
	{"GCA", 'A'}, 
	{"GCG", 'A'},
	{"UAU", 'Y'},	/* Y - Tirosina */ 
	{"UAC", 'Y'},
	{"UAA", '.'},	/* . - Ocre Stop */
	{"UAG", '.'},	/* . - Amber Stop */
	{"CAU", 'H'},	/* H - Histidine */ 
	{"CAC", 'H'},
	{"CAA", 'Q'},	/* Q - Glutamine */
	{"CAG", 'Q'},
	{"AAU", 'N'},	/* N - Asparagine */ 
	{"AAC", 'N'},
	{"AAA", 'K'},	/* K- Lysine */ 
	{"AAG", 'K'},
	{"GAU", 'D'},	/* D - Aspartic acid */ 
	{"GAC", 'D'},
	{"GAA", 'E'},	/* E - Glutamic acid */ 
	{"GAG", 'E'},
	{"UGU", 'C'},	/* C - Cysteine */ 
	{"UGC", 'C'},
	{"UGA", '.'},	/* . - Opal Stop */
	{"UGG", 'W'},	/* W - Tryptophan */
	{"CGU", 'R'},	/* R - Arginine */ 
	{"CGC", 'R'},
	{"CGA", 'R'}, 
	{"CGG", 'R'}, 
	{"AGA", 'R'}, 
	{"AGG", 'R'},
	{"GGU", 'G'},	/* G - Glycine */ 
	{"GGC", 'G'},
	{"GGA", 'G'},
	{"GGG", 'G'}
};

static int compare(const void *p1, const void *p2)
{
	return strcmp(*((const char **)p1), *((const char **)p2));
}

char get_aminoacid(const char *codon)
{
	if(!sorted) {
		qsort(codes, NELEMS(codes), sizeof(*codes), compare);
		sorted = true;
	}
	struct code *code = bsearch(&codon, codes, NELEMS(codes), 
		sizeof(*codes), compare);
	return code ? code->aminoacid : '\0';
}

bool is_DNA(char *sequence)
{
	int seq_len = strlen(sequence);
	for (int i = 0; i < seq_len; ++i) {
		char c = sequence[i];
		if (c !='A' && c != 'T' && c != 'G' && c !='C'){
			return false;
		}
	}
	return true;
}

void to_RNA(char *sequence)
{
	int seq_len = strlen(sequence);
	for (int i = 0; i < seq_len; ++i) {
		if (sequence[i] == 'T') {
			sequence[i] = 'U';
		}
	}
}

void RNA_to_aminoacids(char *sequence)
{
	int seq_len = strlen(sequence);
	/* make sure to read only codon triplets */
	while (seq_len % 3 != 0) {
		seq_len -= 1;
	}
	seq_len -= 3;
	int j = 0;
	for (int i = 0; i <= seq_len; i += 3) {
		char codon[4];
		strncpy(codon, sequence + i, 3);
		codon[3] = '\0';
		sequence[j] = get_aminoacid(codon);
		++j;
	}
	sequence[j] = '\0';
}

void translate_file_with_DNA() 
{
	puts("**DNA-to-Amino-Acids**");
	/* scans file with DNA */
	printf("Please enter the path of file to scan: ");
	char source_file_name[80];
	scanf("%79s", source_file_name);
	/* scans DNA sequence */
	FILE *source_file_ptr;
	if (!(source_file_ptr = fopen(source_file_name, "r"))) {
		fputs("translate_file_with_DNA() error: file open failed.", stderr);
		fclose(source_file_ptr);	
		return;        
	}
	/* scans text until newline */ 
	char translation_buffer[4096];
	fscanf(source_file_ptr,"%4095[^\n]", translation_buffer);
	fclose(source_file_ptr);
	printf("Data from the file: %s\n", translation_buffer);
	if (!is_DNA(translation_buffer)) {
		fputs("translate_file_with_DNA() error: this is not DNA.", stderr);
		return;
	}
	to_RNA(translation_buffer);
	printf("RNA: %s\n", translation_buffer);
	RNA_to_aminoacids(translation_buffer);
	printf("Aminoacids: %s\n", translation_buffer); 
}

int main(int argc, char *argv[]) 
{	
	translate_file_with_DNA(); 
}
