/**
 * @file    /get_seq.c
 * @brief   script description
 * 
 * @author  hanfc
 * @date    2024/06/06
 * @license MIT License
 * 
 * @version 1.0
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>


#define MAX_LINE_LEN 1024
#define MAX_GENE_LEN 100
#define MAX_LOCATION_LEN 100
// #define MAX_SEQUENCE_LEN 1000000
#define MIN_SEQUENCE_LEN 10000


// Function prototypes
void print_usage(const char *prog_name) {
    fprintf(stdout, "Usage:%s -g <genbank_file> -a\n", prog_name);
    fprintf(stdout, "Required options:\n");
    fprintf(stdout, "   -g, --genbank  Intput genbank file\n");

    fprintf(stdout, "Optional options:\n");
    fprintf(stdout, "   -pre, --prefix  Prefix of the output\n");
    fprintf(stdout, "   -a, --all    Flag to output all annotations\n");
    fprintf(stdout, "   -f, --faa    Flag to output fasta\n");
    fprintf(stdout, "   -p, --pep    Flag to output pep\n");
    fprintf(stdout, "   -c, --cds    Flag to output cds\n");
    fprintf(stdout, "   -t, --trn    Flag to output trn\n");
    fprintf(stdout, "   -r, --rrn    Flag to output rrn\n");
    fprintf(stdout, "   -o, --output The output path\n");
    fprintf(stdout, "   -h, --help      Display this help message\n");
}


void parse_arguments(int argc, char *argv[], char *genbank_file, char *prefix, int *all_flag, int *faa_flag, int *pep_flag, int *cds_flag, int *trn_flag, int *rrn_flag, char *output_file) {
    int i;
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--genbank") == 0) {
            if (i + 1 < argc) {
                // genbank_file = argv[++i];
                char *genbank_path = argv[++i];
                strcpy(genbank_file, genbank_path);
            } else {
                fprintf(stderr, "Error: Missing genbank file argument\n");
                print_usage(argv[0]);
                exit(EXIT_FAILURE);
            }
        } else if (strcmp(argv[i], "-pre") == 0 || strcmp(argv[i], "--prefix") == 0) {
            if (i + 1 < argc) {
                // *prefix = argv[++i];
                char *prefix_path = argv[++i];
                strcpy(prefix, prefix_path);
            }
        } else if (strcmp(argv[i], "-a") == 0 || strcmp(argv[i], "--all") == 0) {
                *all_flag = 1, *faa_flag = 1, *pep_flag = 1, *cds_flag = 1, *trn_flag = 1, *rrn_flag = 1;
        } else if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--faa") == 0) {
                *faa_flag = 1;
        } else if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--pep") == 0) {
                *pep_flag = 1;
        } else if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--cds") == 0) {
                *cds_flag = 1;
        } else if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--trn") == 0) {
                *trn_flag = 1;
        } else if (strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "--rrn") == 0) {
                *rrn_flag = 1;
        } else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
            if (i + 1 < argc) {
                char *output_path = argv[++i];
                strcpy(output_file, output_path);
                // output_file = argv[++i];
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

    if (strlen(genbank_file) == 0) {
        fprintf(stderr, "Error: Please provide all required arguments\n");
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    } else {
        char *ext = strrchr(genbank_file, '.'); // get the extension of the file name
        if (strcmp(ext, ".gb") != 0) { // check if the extension is ".gb"
            fprintf(stderr, "Error: Genbank file must have a extension (.gb)\n");
            print_usage(argv[0]);
            exit(EXIT_FAILURE);
        } else {
            size_t bnmlen = ext - genbank_file; // get the length of the base name
            char *bnm = (char *)malloc(bnmlen +  1);

            if (bnm == NULL) 
            {
                fprintf(stderr, "Error: Failed to allocate memory for base name\n");
                exit(EXIT_FAILURE);
            }
            strncpy(bnm, genbank_file, bnmlen);
            bnm[bnmlen] = '\0';

            if (strlen(prefix) == 0) {
                // prefix = strrchr(bnm, '/') + 1;
                strcpy(prefix, strrchr(bnm, '/') + 1);
            }

            if (strlen(output_file) == 0)
            {
                char *last_prefix = strrchr(bnm, '/') + 1;
                size_t outlen = last_prefix - bnm;
                strncpy(output_file, bnm, outlen);
                output_file[outlen] = '\0';
                printf("The output file is: %s", output_file);
            } 
            free(bnm);
        }
    }

    if (*all_flag == 0 && *faa_flag == 0 && *pep_flag == 0 && *cds_flag == 0 && *trn_flag == 0 && *rrn_flag == 0) {
        *all_flag = 1;
    }
    
}


// define a struct to store the annotation information
typedef struct {
    char *gene;
    char *location;
    char *sequence;
} Cds;

// define a struct to store all the annotations
typedef struct {
    char *gene;
    char *sequence;
} Faa;

// define a struct to store all the annotations
typedef struct {
    char *gene;
    char *sequence;
} Pep;

//
typedef struct {
    char *gene;
    char *location;
    char *sequence;
} Rrn;

//
typedef struct {
    char *gene;
    char *location;
    char *sequence;
} Trn;


// Function to get the reverse complement of a sequence
char* reverse_complement(char *seq) {
    char *rc_seq = (char*)malloc(strlen(seq) + 1);
    int len = strlen(seq);
    for (int i = len - 1, m = 0; i >= 0; --i, ++m) {
        switch (seq[i]) {
            case 'A': case 'a': rc_seq[m] = 'T'; break;
            case 'T': case 't': rc_seq[m] = 'A'; break;
            case 'C': case 'c': rc_seq[m] = 'G'; break;
            case 'G': case 'g': rc_seq[m] = 'C'; break;
        }
    }
    rc_seq[len] = '\0';
    return rc_seq;
}

// Function to extract a substring
char* subseq(const char* str, int start, int end) {
    start--;
    end--;
    int len = end - start + 1;
    char *sub = (char*)malloc((len + 1) * sizeof(char));
    if (end >= strlen(str)) {
        fprintf(stderr, "Error: Invalid location '%d..%d'\n", start, end);
        exit(EXIT_FAILURE);
    }
    strncpy(sub, str + start, len);
    sub[len] = '\0';
    return sub;
}


char* extract_sequence(char *seq, char *location) {
    int loclen = strlen(location);
    int join_flag = 0;
    int com_flag = 0;
    int cj_flag = 0;


    char *tk_loc = malloc(loclen + 1);
    memset(tk_loc, 0, loclen + 1);

    for (int i = 0, j = 0; i < loclen; i++) {
        if (location[i] == '>' || location[i] == '<') {
            continue;
        } else if (com_flag == 0 && join_flag == 0 && cj_flag == 0 && location[i] != 'j' && location[i] != 'c') {
            if (location[i] == ',') {
                tk_loc[j] = ' ';
                j++;
            } else if (location[i] == ')') {
                continue;
            } else {
                tk_loc[j] = location[i];
                j++;
            }
        } else if (location[i] == 'c' && location[i+11] == 'j') {
            cj_flag = 1;
            tk_loc[j] = 'c';
            i += 15;
            j++;
        } else if (cj_flag == 1) {
            if (location[i] == ')') {
                if (i < loclen -2) {
                    tk_loc[j] = ' ';
                    cj_flag = 0;
                    j++;
                } else {
                    break;
                }
            } else {
                tk_loc[j] = location[i];
                j++;
            }
        } else if (location[i] == 'j') {
            i += 4;
        } else if (location[i] == 'c' && location[i+11] != 'j') {
            com_flag = 1;
            tk_loc[j] = '-';
            i += 10;
            j++;
        } else if (com_flag == 1) {
            if (location[i] == ')') {
                if (i < loclen -2) {
                    com_flag = 0;
                    tk_loc[j] = ' ';
                    j++;
                } else {
                    break;
                }
            } else {
                tk_loc[j] = location[i];
                j++;
            }
        } 
        tk_loc[j] = '\0';
    }
    char *tk_subseq = malloc(MIN_SEQUENCE_LEN);
    memset(tk_subseq, 0, MIN_SEQUENCE_LEN);


    if (strchr(tk_loc, ' ') == NULL) {
        if (strchr(tk_loc, ',') == NULL) {
            if (strstr(tk_loc, "-")) {
                char cp_loc[strlen(tk_loc)];
                strcpy(cp_loc, tk_loc + 1);
                int l_loc = 0, r_loc = 0;
                sscanf(cp_loc, "%d..%d", &l_loc, &r_loc);
                char *rc_subseq = subseq(seq, l_loc, r_loc);
                // char *tk_subseq = malloc(strlen(rc_subseq) + 1);
                tk_subseq = reverse_complement(rc_subseq);
                // free(rc_subseq);
            } else {
                int l_loc = 0, r_loc = 0;
                sscanf(tk_loc, "%d..%d", &l_loc, &r_loc);
                tk_subseq = subseq(seq, l_loc, r_loc);
            }
        } else {
            char cm_subseq[MIN_SEQUENCE_LEN];
            memset(cm_subseq, 0, MIN_SEQUENCE_LEN);
            char cp_loc[strlen(tk_loc)];
            strcpy(cp_loc, tk_loc + 1);
            int l_loc = 0, r_loc = 0;
            char *tokendot = strtok(cp_loc, ",");
            while (tokendot != NULL) {
                sscanf(tokendot, "%d..%d", &l_loc, &r_loc);
                char *temp_subseq = subseq(seq, l_loc, r_loc);
                strcat(cm_subseq, temp_subseq);
                // free(temp_subseq);
                tokendot = strtok(NULL, ",");
            }
            tk_subseq =  reverse_complement(cm_subseq);
        }
    } else {
        char *tokenspace = strtok(tk_loc, " ");
        while (tokenspace != NULL) {
            int l_loc = 0, r_loc = 0;
            // char *temp_subseq = NULL;
            if (strstr(tokenspace, "-")) {
                char cp_space[strlen(tokenspace)];
                strcpy(cp_space, tokenspace + 1);
                sscanf(cp_space, "%d..%d", &l_loc, &r_loc);
                char *temp_subseq = subseq(seq, l_loc, r_loc);
                char *rc_temp_subseq = reverse_complement(temp_subseq);
                strcat(tk_subseq, rc_temp_subseq);
                // free(temp_subseq);
                // free(rc_temp_subseq);
            } else {
                sscanf(tokenspace, "%d..%d", &l_loc, &r_loc);
                char *temp_subseq = subseq(seq, l_loc, r_loc);
                strcat(tk_subseq, temp_subseq);
                // free(temp_subseq);
            }
            tokenspace = strtok(NULL, " ");
        }
    }
    return tk_subseq;
}


void extract_annotation(int *cds_count, int *rna_count, int *trn_count, Cds **cds_list, Faa **faa, Pep **pep_list, Rrn **rrna_list, Trn **trna_list, char *genbank_file) {

    FILE *gbk = fopen(genbank_file, "r");
    if (!gbk) {
        fprintf(stderr, "Error: Failed to open genbank file '%s'\n", genbank_file);
        exit(EXIT_FAILURE);
    }


    char line[MAX_LINE_LEN];

    *pep_list = malloc(sizeof(Pep) * 100);
    *cds_list = malloc(sizeof(Cds) * 100);
    *rrna_list = malloc(sizeof(Rrn) * 100);
    *trna_list = malloc(sizeof(Trn) * 100);
    if (*pep_list == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for pep_list\n");
        fclose(gbk);
        exit(EXIT_FAILURE);
    }
    *faa = malloc(sizeof(Faa));

    int cds_flag = 0;
    int pep_flag = 0;
    int rrn_flag = 0;
    int trn_flag = 0;
    int faa_flag = 0;
    int cds_loc_flag = 0;
    int trn_loc_flag = 0;
    int rrn_loc_flag = 0;
    int cds_seq = 0;
    int pep_seq = 0;
    int rrn_seq = 0;
    int trn_seq = 0;
    int cds_len = 0;
    int seq_flag = 0;
    char *organ = NULL;
    char *acces = NULL;

    char temp_seq[MIN_SEQUENCE_LEN] = "";
    char temp_loc[MAX_LOCATION_LEN] = "";
    char gene_id[MAX_GENE_LEN] = "";


    while (fgets(line, MAX_LINE_LEN, gbk)) {
        if (strstr(line, "ORGANISM")) {
            organ = (char *)malloc(strlen(line + 12) + 1); // allocate memory for the organism name
            if (organ == NULL)
            {
                fprintf(stderr, "Error: Failed to allocate memory for organism name\n");
                fclose(gbk);
                exit(EXIT_FAILURE);
            }
            strcpy(organ, line + 12);
            printf("The organism is: %s", organ);
        }
        if(strstr(line, "ACCESSION")) {
            acces = (char *)malloc(strlen(line + 12) + 1); // allocate memory for the accession
            if (acces == NULL)
            {
                fprintf(stderr, "Error: Failed to allocate memory for accession\n");
                fclose(gbk);
                exit(EXIT_FAILURE);
            }
            strcpy(acces, line + 12);
            printf("The accession is: %s", acces);
        }

        char seq1[60] = "";
        char seq2[60] = "";
        char seq3[60] = "";
        char seq4[60] = "";
        char seq5[60] = "";
        char seq6[60] = "";
        char temp_faa[60 * 6] = "";
        if (strstr(line, "ORIGIN")) {
            (*faa)->sequence = malloc(1);
            faa_flag = 1;
        } else if (faa_flag == 1) {
            sscanf(line, "        %*d %s %s %s %s %s %s", seq1, seq2, seq3, seq4, seq5, seq6);
            sprintf(temp_faa, "%s%s%s%s%s%s", seq1, seq2, seq3, seq4, seq5, seq6);
            (*faa)->sequence = realloc((*faa)->sequence, strlen((*faa)->sequence) + strlen(temp_faa) + 1);
            if ((*faa)->sequence == NULL)
            {
                fprintf(stderr, "Error: gb file format error or incomplete sequence.\n");
                fclose(gbk);
                exit(EXIT_FAILURE);
            }

            int i = 0;
            while (i < strlen(temp_faa))
            {
                temp_faa[i] = toupper(temp_faa[i]);
                i++;
            }
            
            strcat((*faa)->sequence, temp_faa);
        }
    }

    if ((*faa)->sequence == NULL) {
        fprintf(stderr, "Error: gb file format error or incomplete sequence.\n");
        fclose(gbk);
        exit(EXIT_FAILURE);
    }

    if (organ != NULL) {
        (*faa)->gene = malloc(strlen(organ) + 2);
    } else {
        organ = malloc(6);
        strcpy(organ, "Chr1\n");
        (*faa)->gene = malloc(strlen(organ) + 2);
    }
    sprintf((*faa)->gene, "%s", organ);


    rewind(gbk); // rewind the file to the beginning

    while (fgets(line, MAX_LINE_LEN, gbk)) {

        if (strstr(line, "     CDS             ")) {
            cds_flag = 1;
            sscanf(line, "     CDS             %s", temp_loc);
            if (line[strlen(line) - 2] == ',') {
                cds_loc_flag = 1;
                (*cds_list)[*cds_count].location = malloc(strlen(temp_loc) + 1);
                if ((*cds_list)[*cds_count].location == NULL)
                {
                    fprintf(stderr, "Error: Failed to allocate memory for cds location\n");
                    fclose(gbk);
                    exit(EXIT_FAILURE);
                }                
                strcpy((*cds_list)[*cds_count].location, temp_loc);
            } else {
                cds_loc_flag = 0;
                (*cds_list)[*cds_count].location = malloc(strlen(temp_loc) + 1);

                if ((*cds_list)[*cds_count].location == NULL)
                {
                    fprintf(stderr, "Error: Failed to allocate memory for cds location\n");
                    fclose(gbk);
                    exit(EXIT_FAILURE);
                }                
                strcpy((*cds_list)[*cds_count].location, temp_loc);
                char *tk_seq = extract_sequence((*faa)->sequence, (*cds_list)[*cds_count].location);
                (*cds_list)[*cds_count].sequence = malloc(strlen(tk_seq) + 1);
                strcpy((*cds_list)[*cds_count].sequence, tk_seq);

            }
        } else if (cds_loc_flag == 1 && cds_flag == 1) {
            sscanf(line, "                     %s", temp_loc);
            if (line[strlen(line) - 2] == ',') {
                (*cds_list)[*cds_count].location = realloc((*cds_list)[*cds_count].location, strlen((*cds_list)[*cds_count].location) + strlen(temp_loc) + 1);
                if ((*cds_list)[*cds_count].location == NULL)
                {
                    fprintf(stderr, "Error: Failed to reallocate memory for cds location\n");
                    fclose(gbk);
                    exit(EXIT_FAILURE);
                }
                strcat((*cds_list)[*cds_count].location, temp_loc);
            } else {
                cds_loc_flag = 0;
                (*cds_list)[*cds_count].location = realloc((*cds_list)[*cds_count].location, strlen((*cds_list)[*cds_count].location) + strlen(temp_loc) + 1);
                if ((*cds_list)[*cds_count].location == NULL)
                {
                    fprintf(stderr, "Error: Failed to reallocate memory for cds location\n");
                    fclose(gbk);
                    exit(EXIT_FAILURE);
                }
                strcat((*cds_list)[*cds_count].location, temp_loc);
                char *tk_seq = extract_sequence((*faa)->sequence, (*cds_list)[*cds_count].location);
                (*cds_list)[*cds_count].sequence = malloc(strlen(tk_seq) + 1);
                strcpy((*cds_list)[*cds_count].sequence, tk_seq);
            }
        } else if (cds_flag == 1 && strstr(line, "/gene="))
        {
            // gene_id = (char *)malloc(MAX_GENE_LEN);
            sscanf(line, "                     /gene=\"%[^\"]", gene_id);
            (*pep_list)[*cds_count].gene = malloc(strlen(gene_id) + 1);
            (*cds_list)[*cds_count].gene = malloc(strlen(gene_id) + 1);
            
            // *pep_list = realloc(*pep_list, sizeof(Pep) * ((*cds_count) + 1)); 
            if ((*pep_list)[*cds_count].gene == NULL) 
            {
                fprintf(stderr, "Error: Failed to reallocate memory for pep_list gene\n");
                fclose(gbk);
                exit(EXIT_FAILURE);
            }
            strcpy((*pep_list)[*cds_count].gene, gene_id);
            strcpy((*cds_list)[*cds_count].gene, gene_id);

        } else if (cds_flag == 1 && strstr(line, "/translation"))
        {
            if (line[strlen(line) - 2] == '\"')
            {
                seq_flag = 0;
                cds_flag = 0;
                sscanf(line, "                     /translation=\"%[^\"]\"", temp_seq);
                if (strlen(temp_seq) >= MIN_SEQUENCE_LEN) 
                {
                    fprintf(stderr, "Error: Sequence length exceeds limit\n");
                    break;
                }
                (*pep_list)[*cds_count].sequence = malloc(strlen(temp_seq) + 1);

                if ((*pep_list)[*cds_count].sequence == NULL)
                {
                    fprintf(stderr, "Error: Failed to reallocate memory for pep_list sequence\n");
                    fclose(gbk);
                    exit(EXIT_FAILURE);
                }

                strcpy((*pep_list)[*cds_count].sequence, temp_seq);
                (*cds_count)++;

            } else {
                seq_flag = 1;
                cds_flag = 1;
                sscanf(line, "                     /translation=\"%[^\n]", temp_seq);
                if (strlen(temp_seq) >= MIN_SEQUENCE_LEN) 
                {
                    fprintf(stderr, "Error: Sequence length exceeds limit\n");
                    break;
                }

                (*pep_list)[*cds_count].sequence = malloc(strlen(temp_seq) + 1);

                if ((*pep_list)[*cds_count].sequence == NULL)
                {
                    fprintf(stderr, "Error: Failed to reallocate memory for pep_list sequence\n");
                    fclose(gbk);
                    exit(EXIT_FAILURE);
                }                
                strcpy((*pep_list)[*cds_count].sequence, temp_seq);
            }
        } else if (cds_flag == 1 && seq_flag == 1) {
            if (line[strlen(line) - 2] == '\"')
            {
                seq_flag = 0;
                cds_flag = 0;
                sscanf(line, "                     %[^\"]", temp_seq);
                if ((strlen((*pep_list)[*cds_count].sequence) + strlen(temp_seq)) >= MIN_SEQUENCE_LEN) 
                {
                    fprintf(stderr, "Error: Sequence length exceeds limit\n");
                    break;
                }

                (*pep_list)[*cds_count].sequence = realloc((*pep_list)[*cds_count].sequence, strlen((*pep_list)[*cds_count].sequence) + strlen(temp_seq) + 1);
                if ((*pep_list)[*cds_count].sequence == NULL)
                {
                    fprintf(stderr, "Error: Failed to reallocate memory for pep_list sequence\n");
                    fclose(gbk);
                    exit(EXIT_FAILURE);
                }
                strcat((*pep_list)[*cds_count].sequence, temp_seq);
                (*cds_count)++;
            } else {
                seq_flag = 1;
                cds_flag = 1;
                sscanf(line, "                     %[^\n]", temp_seq);
                if ((strlen((*pep_list)[*cds_count].sequence) + strlen(temp_seq)) >= MIN_SEQUENCE_LEN) 
                {
                    fprintf(stderr, "Error: Sequence length exceeds limit\n");
                    break;
                }
                (*pep_list)[*cds_count].sequence = realloc((*pep_list)[*cds_count].sequence, strlen((*pep_list)[*cds_count].sequence) + strlen(temp_seq) + 1);
                if ((*pep_list)[*cds_count].sequence == NULL)
                {
                    fprintf(stderr, "Error: Failed to reallocate memory for pep_list sequence\n");
                    fclose(gbk);
                    exit(EXIT_FAILURE);
                }
                strcat((*pep_list)[*cds_count].sequence, temp_seq);
            }
        } else if (strstr(line, "     rRNA            ")) {
            rrn_flag = 1;
            sscanf(line, "     rRNA             %s", temp_loc);
            if (line[strlen(line) - 2] == ',') {
                rrn_loc_flag = 1;
                (*rrna_list)[*rna_count].location = malloc(strlen(temp_loc) + 1);
                if ((*rrna_list)[*rna_count].location == NULL)
                {
                    fprintf(stderr, "Error: Failed to allocate memory for rRNA location\n");
                    fclose(gbk);
                    exit(EXIT_FAILURE);
                }                
                strcpy((*rrna_list)[*rna_count].location, temp_loc);
            } else {
                rrn_loc_flag = 0;
                (*rrna_list)[*rna_count].location = malloc(strlen(temp_loc) + 1);

                if ((*rrna_list)[*rna_count].location == NULL)
                {
                    fprintf(stderr, "Error: Failed to allocate memory for rRNA location\n");
                    fclose(gbk);
                    exit(EXIT_FAILURE);
                }                
                strcpy((*rrna_list)[*rna_count].location, temp_loc);
                char *tk_seq = extract_sequence((*faa)->sequence, (*rrna_list)[*rna_count].location);
                (*rrna_list)[*rna_count].sequence = malloc(strlen(tk_seq) + 1);
                strcpy((*rrna_list)[*rna_count].sequence, tk_seq);

            }
        } else if (rrn_loc_flag == 1 && rrn_flag == 1) {
            sscanf(line, "                     %s", temp_loc);
            if (line[strlen(line) - 2] == ',') {
                (*rrna_list)[*rna_count].location = realloc((*rrna_list)[*rna_count].location, strlen((*rrna_list)[*rna_count].location) + strlen(temp_loc) + 1);
                if ((*rrna_list)[*rna_count].location == NULL)
                {
                    fprintf(stderr, "Error: Failed to reallocate memory for rRNA location\n");
                    fclose(gbk);
                    exit(EXIT_FAILURE);
                }
                strcat((*rrna_list)[*rna_count].location, temp_loc);
            } else {
                rrn_loc_flag = 0;
                (*rrna_list)[*rna_count].location = realloc((*rrna_list)[*rna_count].location, strlen((*rrna_list)[*rna_count].location) + strlen(temp_loc) + 1);
                if ((*rrna_list)[*rna_count].location == NULL)
                {
                    fprintf(stderr, "Error: Failed to reallocate memory for rRNA location\n");
                    fclose(gbk);
                    exit(EXIT_FAILURE);
                }
                strcat((*rrna_list)[*rna_count].location, temp_loc);
                char *tk_seq = extract_sequence((*faa)->sequence, (*rrna_list)[*rna_count].location);
                (*rrna_list)[*rna_count].sequence = malloc(strlen(tk_seq) + 1);
                strcpy((*rrna_list)[*rna_count].sequence, tk_seq);
            }
        } else if (rrn_flag == 1 && strstr(line, "/gene=")) {
            // gene_id = (char *)malloc(MAX_GENE_LEN);
            sscanf(line, "                     /gene=\"%[^\"]", gene_id);

            (*rrna_list)[*rna_count].gene = malloc(strlen(gene_id) + 1);
            
            if ((*rrna_list)[*rna_count].gene == NULL) 
            {
                fprintf(stderr, "Error: Failed to reallocate memory for pep_list gene\n");
                fclose(gbk);
                exit(EXIT_FAILURE);
            }
            strcpy((*rrna_list)[*rna_count].gene, gene_id);

            (*rna_count)++;
            rrn_flag = 0;
        } else if (strstr(line, "     tRNA            ")) {
            trn_flag = 1;
            sscanf(line, "     tRNA             %s", temp_loc);
            if (line[strlen(line) - 2] == ',') {
                trn_loc_flag = 1;
                (*trna_list)[*trn_count].location = malloc(strlen(temp_loc) + 1);
                if ((*trna_list)[*trn_count].location == NULL)
                {
                    fprintf(stderr, "Error: Failed to allocate memory for tRNA location\n");
                    fclose(gbk);
                    exit(EXIT_FAILURE);
                }                
                strcpy((*trna_list)[*trn_count].location, temp_loc);
            } else {
                trn_loc_flag = 0;
                (*trna_list)[*trn_count].location = malloc(strlen(temp_loc) + 1);

                if ((*trna_list)[*trn_count].location == NULL)
                {
                    fprintf(stderr, "Error: Failed to allocate memory for tRNA location\n");
                    fclose(gbk);
                    exit(EXIT_FAILURE);
                }                
                strcpy((*trna_list)[*trn_count].location, temp_loc);
                char *tk_seq = extract_sequence((*faa)->sequence, (*trna_list)[*trn_count].location);
                (*trna_list)[*trn_count].sequence = malloc(strlen(tk_seq) + 1);
                strcpy((*trna_list)[*trn_count].sequence, tk_seq);

            }
        } else if (trn_loc_flag == 1 && trn_flag == 1) {
            sscanf(line, "                     %s", temp_loc);
            if (line[strlen(line) - 2] == ',') {
                (*trna_list)[*trn_count].location = realloc((*trna_list)[*trn_count].location, strlen((*trna_list)[*trn_count].location) + strlen(temp_loc) + 1);
                if ((*trna_list)[*trn_count].location == NULL)
                {
                    fprintf(stderr, "Error: Failed to reallocate memory for tRNA location\n");
                    fclose(gbk);
                    exit(EXIT_FAILURE);
                }
                strcat((*trna_list)[*trn_count].location, temp_loc);
            } else {
                trn_loc_flag = 0;
                (*trna_list)[*trn_count].location = realloc((*trna_list)[*trn_count].location, strlen((*trna_list)[*trn_count].location) + strlen(temp_loc) + 1);
                if ((*trna_list)[*trn_count].location == NULL)
                {
                    fprintf(stderr, "Error: Failed to reallocate memory for tRNA location\n");
                    fclose(gbk);
                    exit(EXIT_FAILURE);
                }
                strcat((*trna_list)[*trn_count].location, temp_loc);
                char *tk_seq = extract_sequence((*faa)->sequence, (*trna_list)[*trn_count].location);
                (*trna_list)[*trn_count].sequence = malloc(strlen(tk_seq) + 1);
                strcpy((*trna_list)[*trn_count].sequence, tk_seq);
            }
        } else if (trn_flag == 1 && strstr(line, "/gene=")) {
            // gene_id = (char *)malloc(MAX_GENE_LEN);
            sscanf(line, "                     /gene=\"%[^\"]", gene_id);
            (*trna_list)[*trn_count].gene = malloc(strlen(gene_id) + 1);
            
            if ((*trna_list)[*trn_count].gene == NULL) 
            {
                fprintf(stderr, "Error: Failed to reallocate memory for pep_list gene\n");
                fclose(gbk);
                exit(EXIT_FAILURE);
            }
            strcpy((*trna_list)[*trn_count].gene, gene_id);
            (*trn_count)++;
            trn_flag = 0;
        }
    }
    fclose(gbk);
    free(organ);
    free(acces);
}



int main(int argc, char *argv[]) {

    char *genbank_file = malloc(1024);
    char *prefix = malloc(1024);
    int all_flag = 0;
    int faa_flag = 0;
    int pep_flag = 0;
    int cds_flag = 0;
    int trn_flag = 0;
    int rrn_flag = 0;
    char *output_file = malloc(1024);
    
    parse_arguments(argc, argv, genbank_file, prefix, &all_flag, &faa_flag, &pep_flag, &cds_flag, &trn_flag, &rrn_flag, output_file);
    printf("[Info] The genbank file: %s\n", genbank_file);
    printf("[Info] The prefix: %s\n", prefix);
    printf("[Info] The output path: %s\n", output_file);

    if (access(output_file, F_OK) == -1) {
        printf("[error] Output Path does not exist.\n");
        exit(1);
    }
    

    // Define arrays to store annotations
    Cds *cds_list = NULL;
    Faa *faa = NULL;
    Pep *pep_list = NULL;
    Rrn *rrn_list = NULL;
    Trn *trn_list = NULL;
    int cds_count = 0;
    int rrn_count = 0;
    int trn_count = 0;

    // if (!cds_list || !faa || !pep_list || !rrn_list || !trn_list) {
    //     fprintf(stderr, "Error: Failed to allocate memory\n");
    //     return 1;
    // }


    extract_annotation(&cds_count, &rrn_count, &trn_count, &cds_list, &faa, &pep_list, &rrn_list, &trn_list, genbank_file);

    // Print the annotations
    if (cds_flag == 1) {
        char *output_cds = malloc(1024);
        strcpy(output_cds, output_file);
        strcat(output_cds, "/");
        strcat(output_cds, prefix);
        strcat(output_cds, ".cds");
        FILE *fcds = fopen(output_cds, "w");
        for (int i = 0; i < cds_count; i++)
        {
            fprintf(fcds, ">%s\n", cds_list[i].gene);
            fprintf(fcds, "%s\n", cds_list[i].sequence);
        }
        fclose(fcds);
    }

    if (rrn_flag == 1) {
        char *output_rrn = malloc(1024);
        strcpy(output_rrn, output_file);
        strcat(output_rrn, "/");
        strcat(output_rrn, prefix);
        strcat(output_rrn, ".rrn");
        FILE *frrn = fopen(output_rrn, "w");
        for (int i = 0; i < rrn_count; i++)
        {
            fprintf(frrn, ">%s\n", rrn_list[i].gene);
            fprintf(frrn, "%s\n", rrn_list[i].sequence);
        }
        fclose(frrn);
    }


    if (trn_flag == 1) {
        char *output_trn = malloc(1024);
        strcpy(output_trn, output_file);
        strcat(output_trn, "/");
        strcat(output_trn, prefix);
        strcat(output_trn, ".trn");
        FILE *ftrn = fopen(output_trn, "w");
        for (int i = 0; i < trn_count; i++)
        {
            fprintf(ftrn, ">%s\n", trn_list[i].gene);
            fprintf(ftrn, "%s\n", trn_list[i].sequence);
        }
        fclose(ftrn);
    }

    
    if (pep_flag == 1) {
        char *output_pep = malloc(1024);
        strcpy(output_pep, output_file);
        strcat(output_pep, "/");
        strcat(output_pep, prefix);
        strcat(output_pep, ".pep");
        FILE *fpep = fopen(output_pep, "w");
        for (int i = 0; i < cds_count; i++)
        {
            fprintf(fpep, ">%s\n", pep_list[i].gene);
            fprintf(fpep, "%s\n", pep_list[i].sequence);
        }
        fclose(fpep);
    }

    if (faa_flag == 1) {
        char *output_faa = malloc(1024);
        strcpy(output_faa, output_file);
        strcat(output_faa, "/");
        strcat(output_faa, prefix);
        strcat(output_faa, ".faa");
        FILE *ffaa = fopen(output_faa, "w");
        fprintf(ffaa, ">%s", faa->gene);
        fprintf(ffaa, "%s\n", faa->sequence);
        fclose(ffaa);
    }

    // Clean up memory
    if (cds_list) free(cds_list);
    if (faa) free(faa);
    if (pep_list) free(pep_list);
    if (rrn_list) free(rrn_list);
    if (trn_list) free(trn_list);
    if (genbank_file) free(genbank_file);
    if (prefix) free(prefix);

    // Clean up dynamically allocated strings
    if (output_file) free(output_file);
    return 0;

}
