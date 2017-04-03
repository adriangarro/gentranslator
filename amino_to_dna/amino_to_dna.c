/* 
    --------------------------------------------------------------
    amino_to_dna.c
    
    Amino acid sequences to DNA sequence.
    
    Notes: 1. The file's name that contains the aminoacids sequence 
              must be not longer of 79 chars.
              
           2. The aminoacids sequence must contain a maximum 
              of 1364 chars.
    
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
    if (!sorted) {
        qsort(codes, NELEMS(codes), sizeof(*codes), compare);
        sorted = true;
    }
    struct code *code = bsearch(&codon, codes, NELEMS(codes), 
        sizeof(*codes), compare);
    return code ? code->aminoacid : '\0';
}

bool are_aminoacids(char *sequence)
{
    int seq_len = strlen(sequence);
    for (int i = 0; i < seq_len; ++i) {
        char c = sequence[i];
        if (c !='F' && c != 'L' && c != 'I' && c !='M' &&
            c != 'V' && c != 'S' && c != 'P' && c != 'T' && 
            c != 'A' && c != 'Y' && c != '.' && c != 'H' &&
            c != 'Q' && c != 'N' && c != 'K' && c != 'D' &&
            c != 'E' && c != 'C' && c != 'W' && c != 'R' && c != 'G') {
            return false;
        }
    }
    return true;
}

void RNA_to_aminoacids(char *RNA_sequence, char *aminoacids_sequence)
{
    int seq_len = strlen(RNA_sequence);
    /* make sure to read only codon triplets */
    while (seq_len % 3 != 0) {
        seq_len -= 1;
    }
    seq_len -= 3;
    int j = 0;
    for (int i = 0; i <= seq_len; i += 3) {
        char codon[4];
        strncpy(codon, RNA_sequence + i, 3);
        codon[3] = '\0';
        aminoacids_sequence[j] = get_aminoacid(codon);
        ++j;
    }
    aminoacids_sequence[j] = '\0';
}

void permutations_aux(FILE *target_file_ptr, char *set, 
    char *prefix, int set_len, int k) 
{
    /* Base case: if k is 0 then print prefix */
    if (k == 0) {
        fprintf(target_file_ptr, "%s\n", prefix);
        return;
    }
    /* One by one add all characters from set and recursively */ 
    /* call for k equals to k-1 */
    for (int i = 0; i < set_len; ++i) {
        /* Next character of input added */
        /* Only permutations of 4095 lenght are allowed */
        char new_prefix[4096]; 
        strcpy(new_prefix, prefix); 
        int len = strlen(new_prefix);
        new_prefix[len] = set[i];;
        new_prefix[len+1] = '\0';
        /* k is decreased because a new character is added */
        permutations_aux(target_file_ptr, set, new_prefix, set_len, k-1); 
    }
}

void permutations(char set[], char *target_file_name, int k)
{
    int set_len = strlen(set);
    FILE *file_ptr = fopen(target_file_name, "w");
    if (file_ptr == NULL) {
        fputs("permutations error: file open failed.", stderr);
        return;
    }
    permutations_aux(file_ptr, set, "", set_len, k);
    fclose(file_ptr);
}

void generate_RNA(char *target_file_name, int k) 
{
    char nitrogenous_bases[] = "ACUG";
    permutations(nitrogenous_bases, target_file_name, k); 
}

void get_valid_RNA(char *source_file_name, char *target_file_name, char *aminoacids)
{
    FILE* source_file_ptr = fopen(source_file_name, "r");
    FILE* target_file_ptr = fopen(target_file_name, "w");
    char current_RNA[4096];
    char current_aminoacids[1365];
    while (fgets(current_RNA, sizeof(current_RNA), source_file_ptr)) {
        RNA_to_aminoacids(current_RNA, current_aminoacids);
        if (strcmp(aminoacids, current_aminoacids) == 0) {
            fputs(current_RNA, target_file_ptr);
        }
    }
    fclose(source_file_ptr);
    fclose(target_file_ptr);
}

void RNA_to_DNA(char *source_file_name, char *target_file_name) 
{
    char sequence[4096];
    FILE* source_file_ptr = fopen(source_file_name, "r");
    if (source_file_ptr == NULL) {
        fputs("RNA_to_DNA error: file open failed.", stderr);
        return;
    }
    FILE* target_file_ptr = fopen(target_file_name, "w");
    if (target_file_ptr == NULL) {
        fputs("RNA_to_DNA error: file open failed.", stderr);
        return;
    }
    puts("Possible DNA sequences:");
    while (fgets(sequence, sizeof(sequence), source_file_ptr)) {
        int seq_len = strlen(sequence);
        for (int i = 0; i < seq_len; ++i) {
            if (sequence[i] == 'U') {
                sequence[i] = 'T';
            }
        }
        printf("%s", sequence);
        fputs(sequence, target_file_ptr);
    }
    fclose(source_file_ptr);
    fclose(target_file_ptr);
} 

void translate_file_with_aminoacids() 
{
    puts("**Amino-Acids-to-DNA**");
    /* scans file with aminoacids */
    printf("Please enter the path of file to scan: ");
    char source_file_name[80];
    scanf("%79s", source_file_name);
    /* scans aminoacids sequence */
    FILE *source_file_ptr;
    if (!(source_file_ptr = fopen(source_file_name, "r"))) {
        fputs("translate_file_with_aminoacids() error: file open failed.", stderr);
        fclose(source_file_ptr);	
        return;        
    }
    /* scans text until newline */
    char aminoacids_buffer[1365]; 
    fscanf(source_file_ptr,"%1364[^\n]", aminoacids_buffer);
    fclose(source_file_ptr);
    printf("Data from the file: %s\n", aminoacids_buffer);
    if (!are_aminoacids(aminoacids_buffer)) {
        fputs("translate_file_with_aminoacids() error: this is not aminoacids.", stderr);
        return;
    }
    generate_RNA("output_RNA.txt", strlen(aminoacids_buffer) * 3);
    get_valid_RNA("output_RNA.txt", "output_valid_RNA.txt", aminoacids_buffer);
    RNA_to_DNA("output_valid_RNA.txt", "output_DNA.txt");
}

int main()
{
    translate_file_with_aminoacids();
}
