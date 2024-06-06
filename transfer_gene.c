/**
 * @file    script.c
 * @brief   Identify horizontal gene transfer (HGT) genes
 *
 * @author  hanfc
 * @date    2024/06/05
 * @license MIT License
 *
 * @version 1.1
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE_LENGTH 1024
#define MAX_GENE_NAME_LENGTH 100

// Define macros for min and max
#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

// Define a struct to hold gene information
typedef struct {
    char name[MAX_GENE_NAME_LENGTH];
    int start;
    int end;
    int length;
    int strand;
} Gene;

// Define a struct to hold BLASTN alignment information
typedef struct {
    char query[MAX_GENE_NAME_LENGTH];
    char subject[MAX_GENE_NAME_LENGTH];
    float identity;
    int alignment_length;
    int q_start;
    int q_end;
    int s_start;
    int s_end;
} Blastn;

// Function prototypes
void print_usage(const char *program_name);
void parse_arguments(int argc, char *argv[], char **transfer_file, char **location_file, char **output_file);
void read_genes(const char *filename, Gene **genes, int *gene_count);
void read_blastn(const char *filename, Blastn **alignments, int *alignment_count);
void find_transfer_genes(Gene *genes, int gene_count, Blastn *alignments, int alignment_count, const char *output_file);
void free_memory(Gene *genes, Blastn *alignments);

void print_usage(const char *program_name) {
    fprintf(stderr, "Usage: %s -t <transfer_file> -l <location_file> -o <output_file>\n", program_name);
    fprintf(stderr, "Required options:\n");
    fprintf(stderr, "   -t, --transfer  Gene transfer file\n");
    fprintf(stderr, "   -l, --location  Gene location file\n");
    fprintf(stderr, "   -o, --output    Output file\n");
    fprintf(stderr, "Optional options:\n");
    fprintf(stderr, "   -h, --help      Display this help message\n");
}

void parse_arguments(int argc, char *argv[], char **transfer_file, char **location_file, char **output_file) {
    int i;
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--transfer") == 0) {
            if (i + 1 < argc) {
                *transfer_file = argv[++i];
            } else {
                fprintf(stderr, "Error: Missing transfer file argument\n");
                print_usage(argv[0]);
                exit(EXIT_FAILURE);
            }
        } else if (strcmp(argv[i], "-l") == 0 || strcmp(argv[i], "--location") == 0) {
            if (i + 1 < argc) {
                *location_file = argv[++i];
            } else {
                fprintf(stderr, "Error: Missing location file argument\n");
                print_usage(argv[0]);
                exit(EXIT_FAILURE);
            }
        } else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
            if (i + 1 < argc) {
                *output_file = argv[++i];
            } else {
                fprintf(stderr, "Error: Missing output file argument\n");
                print_usage(argv[0]);
                exit(EXIT_FAILURE);
            }
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            exit(EXIT_SUCCESS);
        } else {
            fprintf(stderr, "Error: Invalid option '%s'\n", argv[i]);
            print_usage(argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    if (*transfer_file == NULL || *location_file == NULL || *output_file == NULL) {
        fprintf(stderr, "Error: Please provide all required arguments\n");
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }
}

void read_genes(const char *filename, Gene **genes, int *gene_count) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error opening gene file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    char line[MAX_LINE_LENGTH];
    *gene_count = 0;

    while (fgets(line, sizeof(line), file)) {
        if (line[0] != '#' && strncmp(line, "Gene", 4) && strlen(line) > 1) {
            (*gene_count)++;
        }
    }

    rewind(file);

    *genes = (Gene *)malloc((*gene_count) * sizeof(Gene));
    if (*genes == NULL) {
        fprintf(stderr, "Error allocating memory for genes\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    int i = 0;
    while (fgets(line, sizeof(line), file)) {
        if (line[0] != '#' && strncmp(line, "Gene", 4) && strlen(line) > 1) {
            sscanf(
                line, "%s %d %d %d %d", 
                (*genes)[i].name, 
                &(*genes)[i].start, 
                &(*genes)[i].end, 
                &(*genes)[i].length, 
                &(*genes)[i].strand
            );
            i++;
        }
    }

    fclose(file);
}

void read_blastn(const char *filename, Blastn **alignments, int *alignment_count) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error opening BLASTN file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    char line[MAX_LINE_LENGTH];
    *alignment_count = 0;

    while (fgets(line, sizeof(line), file)) {
        if (line[0] != '#' && strlen(line) > 1) {
            (*alignment_count)++;
        }
    }

    rewind(file);

    *alignments = (Blastn *)malloc((*alignment_count) * sizeof(Blastn));
    if (*alignments == NULL) {
        fprintf(stderr, "Error allocating memory for alignments\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    int i = 0;
    while (fgets(line, sizeof(line), file)) {
        if (line[0] != '#' && strlen(line) > 1) {
            sscanf(
                line, "%s %s %f %d %*d %*d %d %d %d %d", 
                (*alignments)[i].query, 
                (*alignments)[i].subject, 
                &(*alignments)[i].identity, 
                &(*alignments)[i].alignment_length, 
                &(*alignments)[i].q_start, 
                &(*alignments)[i].q_end, 
                &(*alignments)[i].s_start, 
                &(*alignments)[i].s_end
            );
            i++;
        }
    }

    fclose(file);
}

void find_transfer_genes(Gene *genes, int gene_count, Blastn *alignments, int alignment_count, const char *output_file) {
    FILE *file = fopen(output_file, "w");
    if (file == NULL) {
        fprintf(stderr, "Error opening output file: %s\n", output_file);
        exit(EXIT_FAILURE);
    }
    fprintf(file, "No\tCp\tMt\tIdentity\tlength\tq.start\tq.end\ts.start\ts.end\tHGT gene\n");
    for (int j = 0; j < alignment_count; j++) {
        char hgt_gene[MAX_GENE_NAME_LENGTH * 100] = ""; // Increase the buffer size as needed
        char  incomplete_gene[MAX_GENE_NAME_LENGTH * 100] = ""; // Increase the buffer size as needed
        for (int i = 0; i < gene_count; i++) {
            if (min(genes[i].start, genes[i].end) >= min(alignments[j].q_start, alignments[j].q_end) &&
                max(genes[i].start, genes[i].end) <= max(alignments[j].q_start, alignments[j].q_end)) {
                strcat(hgt_gene, genes[i].name);
                strcat(hgt_gene, "* ");
            }else if (genes[i].start > min(alignments[j].q_start, alignments[j].q_end) && genes[i].start < max(alignments[j].q_start, alignments[j].q_end) ||
                genes[i].end > min(alignments[j].q_start, alignments[j].q_end) && genes[i].end < max(alignments[j].q_start, alignments[j].q_end)) {
                    printf("ceshi\n");
                    printf("%s %d %d %d %d\n",genes[i].name, genes[i].start, genes[i].end, alignments[j].q_start, alignments[j].q_end);
                    strcat(incomplete_gene, genes[i].name);
                    strcat(incomplete_gene, " ");
            }
            
            
        }
        fprintf(file, "%d\t%s\t%s\t%.2f\t%d\t%d\t%d\t%d\t%d\t%s%s\n",
                j + 1,
                "Cp",
                "Mt",
                alignments[j].identity,
                alignments[j].alignment_length,
                alignments[j].q_start,
                alignments[j].q_end,
                alignments[j].s_start,
                alignments[j].s_end,
                hgt_gene,
                incomplete_gene);
    }

    fclose(file);
}

void free_memory(Gene *genes, Blastn *alignments) {
   free(genes);
   free(alignments);
}

int main(int argc, char *argv[]) {
   char *transfer_file = NULL;
   char *location_file = NULL;
   char *output_file = NULL;

   parse_arguments(argc, argv, &transfer_file, &location_file, &output_file);

   Gene *genes = NULL;
   int gene_count = 0;
   Blastn *alignments = NULL;
   int alignment_count = 0;

   read_genes(location_file, &genes, &gene_count);
   read_blastn(transfer_file, &alignments, &alignment_count);

   find_transfer_genes(genes, gene_count, alignments, alignment_count, output_file);

   free_memory(genes, alignments);

   return 0;
}