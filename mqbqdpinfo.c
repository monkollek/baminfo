/* This program demonstrates how to generate pileup from multiple BAMs
 * simutaneously, to achieve random access and to use the BED interface.
 * To compile this program separately, you may:
 *
 *   gcc -o mqbqdpinfo mqbqdpinfo.c -I../htslib ../htslib/libhts.a -lz -lpthread ../samtools/bedidx.o
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <sys/stat.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "samtools.h"

typedef struct {     // auxiliary data structure
	samFile *fp;     // the file handle
	bam_hdr_t *hdr;  // the file header
	hts_itr_t *iter; // NULL if a region not specified
	int min_mapQ, min_len; // mapQ filter; length filter
} aux_t;

void *bed_read(const char *fn); // read a BED or position list file
void bed_destroy(void *_h);     // destroy the BED data structure
int bed_overlap(const void *_h, const char *chr, int beg, int end); // test if chr:beg-end overlaps


// This function reads a BAM alignment from one BAM file.
static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
	aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
    int ret;
    while (1)
    {
        ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
        if ( ret<0 ) break;
        if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
        if ( (int)b->core.qual < aux->min_mapQ ) continue;
        if ( aux->min_len && bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)) < aux->min_len ) continue;
        break;
    }
	return ret;
}

//int read_file_list(const char *file_list,int *n,char **argv[]);

#define MAX_PATH_LEN 1024
int read_file_list(const char *file_list,int *n,char **argv[])
{
    char buf[MAX_PATH_LEN];
    int len, nfiles = 0;
    char **files = NULL;
    struct stat sb;

    *n = 0;
    *argv = NULL;

    FILE *fh = fopen(file_list,"r");
    if ( !fh )
    {
        fprintf(stderr,"%s: %s\n", file_list,strerror(errno));
        return 1;
    }

    files = calloc(nfiles,sizeof(char*));
    nfiles = 0;
    while ( fgets(buf,MAX_PATH_LEN,fh) )
    {
        // allow empty lines and trailing spaces
        len = strlen(buf);
        while ( len>0 && isspace(buf[len-1]) ) len--;
        if ( !len ) continue;

        // check sanity of the file list
        buf[len] = 0;
        if (stat(buf, &sb) != 0)
        {
            // no such file, check if it is safe to print its name
            int i, safe_to_print = 1;
            for (i=0; i<len; i++)
                if (!isprint(buf[i])) { safe_to_print = 0; break; }
            if ( safe_to_print )
                fprintf(stderr,"The file list \"%s\" appears broken, could not locate: %s\n", file_list,buf);
            else
                fprintf(stderr,"Does the file \"%s\" really contain a list of files and do all exist?\n", file_list);
            return 1;
        }

        nfiles++;
        files = realloc(files,nfiles*sizeof(char*));
        files[nfiles-1] = strdup(buf);
    }
    fclose(fh);
    if ( !nfiles )
    {
        fprintf(stderr,"No files read from %s\n", file_list);
        return 1;
    }
    *argv = files;
    *n    = nfiles;
    return 0;
}
#undef MAX_PATH_LEN

char translate_base(int enc_base){
	if(enc_base == 1){
		return 'A';
	}
	else if(enc_base == 2){
		return 'C';
	}
	else if(enc_base == 4){
		return 'G';
	}
	else if(enc_base == 8){
		return 'T';
	}
	else{
		return '?';
	}


}


int main(int argc, char *argv[])
{
	int i, n, tid, beg, end, pos, *n_plp, baseQ = 0, mapQ = 0, min_len = 0, status = EXIT_SUCCESS, nfiles;
	const bam_pileup1_t **plp;
	char *reg = 0; // specified region
	void *bed = 0; // BED data structure
    char *file_list = NULL, **fn = NULL;
	bam_hdr_t *h = NULL; // BAM header of the 1st input
	aux_t **data;
	bam_mplp_t mplp;
	faidx_t *fai;
	char *ref = 0;
	int len;


	// parse the command line
	while ((n = getopt(argc, argv, "r:b:q:Q:l:f:")) >= 0) {
		switch (n) {
			case 'l': min_len = atoi(optarg); break; // minimum query length
			case 'r': reg = strdup(optarg); break;   // parsing a region requires a BAM header
			case 'b':
				bed = bed_read(optarg); // BED or position list file can be parsed now
				if (!bed) { //print_error_errno("Could not read file \"%s\"", optarg); 
							return 1; }
				break;
			case 'q': baseQ = atoi(optarg); break;   // base quality threshold
			case 'Q': mapQ = atoi(optarg); break;    // mapping quality threshold
			case 'f': file_list = optarg; break;
		}
	}
	if (optind == argc && !file_list) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: samtools depth [options] in1.bam [in2.bam [...]]\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "   -b <bed>            list of positions or regions\n");
        fprintf(stderr, "   -f <list>           list of input BAM filenames, one per line [null]\n");
        fprintf(stderr, "   -l <int>            minQLen\n");
        fprintf(stderr, "   -q <int>            base quality threshold\n");
        fprintf(stderr, "   -Q <int>            mapping quality threshold\n");
        fprintf(stderr, "   -r <chr:from-to>    region\n");
        fprintf(stderr, "\n");
		return 1;
	}

	fai = fai_load("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta");
	ref = fai_fetch(fai,reg, &len);

	// initialize the auxiliary data structures
    if (file_list) 
    {
        if ( read_file_list(file_list,&nfiles,&fn) ) return 1;
        n = nfiles;
        argv = fn;
        optind = 0;
    }
    else
        n = argc - optind; // the number of BAMs on the command line
	data = calloc(n, sizeof(aux_t*)); // data[i] for the i-th input
	beg = 0; end = 1<<30;  // set the default region
	for (i = 0; i < n; ++i) {
		data[i] = calloc(1, sizeof(aux_t));
		data[i]->fp = sam_open(argv[optind+i], "r"); // open BAM
		if (data[i]->fp == NULL) {
			//print_error_errno("Could not open \"%s\"", argv[optind+i]);
			status = EXIT_FAILURE;
			goto depth_end;
		}
		data[i]->min_mapQ = mapQ;                    // set the mapQ filter
		data[i]->min_len  = min_len;                 // set the qlen filter
		data[i]->hdr = sam_hdr_read(data[i]->fp);    // read the BAM header
		if (reg) { // if a region is specified
			hts_idx_t *idx = sam_index_load(data[i]->fp, argv[optind+i]);  // load the index
			if (idx == NULL) {
				//print_error("can't load index for \"%s\"", argv[optind+i]);
				status = EXIT_FAILURE;
				goto depth_end;
			}
			data[i]->iter = sam_itr_querys(idx, data[i]->hdr, reg); // set the iterator
			hts_idx_destroy(idx); // the index is not needed any more; free the memory
			if (data[i]->iter == NULL) {
				//print_error("can't parse region \"%s\"", reg);
				status = EXIT_FAILURE;
				goto depth_end;
			}
		}
	}

	h = data[0]->hdr; // easy access to the header of the 1st BAM
	if (reg) {
		beg = data[0]->iter->beg; // and to the parsed region coordinates
		end = data[0]->iter->end;
	}

	// the core multi-pileup loop
	mplp = bam_mplp_init(n, read_bam, (void**)data); // initialization
	n_plp = calloc(n, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
	plp = calloc(n, sizeof(bam_pileup1_t*)); // plp[i] points to the array of covering reads (internal in mplp)

	bam_mplp_set_maxcnt(mplp, 250); // cap to 300
	bam_mplp_init_overlaps(mplp); // don't double count bases when reads overlap

	int mq_total, bq_total, total_bases;
	while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) { // come to the next covered position
		if (pos < beg || pos >= end) continue; // out of range; skip
		if (bed && bed_overlap(bed, h->target_name[tid], pos, pos + 1) == 0) continue; // not in BED; skip
		total_bases = 0;
		//fputs(h->target_name[tid], stdout); printf("\t%d", pos+1); // a customized printf() would be faster
		
		// for each bam
		for (i = 0; i < n; ++i) { // base level filters have to go here
			int j, m = 0;
			mq_total = 0;
			bq_total = 0;
			
			// for each read that overlaps position
			for (j = 0; j < n_plp[i]; ++j) {
				const bam_pileup1_t *p = plp[i] + j; // DON'T modfity plp[][] unless you really know
				if (p->is_del || p->is_refskip){ 
					++m; // having dels or refskips at tid:pos
					continue;
				}
				//else if (bam_get_qual(p->b)[p->qpos] < baseQ) ++m; // low base quality

				bq_total += bam_get_qual(p->b)[p->qpos];
				mq_total += ((p->b)->core).qual;

			}

			total_bases += n_plp[i] - m;
			
			/*			
			if(n_plp[i] - m < 250){
				total_bases += n_plp[i] - m;
			}
			else{
				total_bases += 250;
			}
			*/
		}

		float mean_cov = (n == 0) ? 0 : (float) total_bases/n;
		float mean_bq = (total_bases == 0) ? 0 : (float) bq_total/total_bases;
		float mean_mq = (total_bases == 0) ? 0 : (float) mq_total/total_bases;


		printf("%s\t%d\t%.2f\t%.2f\t%.2f\n",h->target_name[tid],pos+1,mean_cov,mean_bq,mean_mq);


		//putchar('\n');
	}
	free(n_plp); free(plp);
	bam_mplp_destroy(mplp);

depth_end:
	for (i = 0; i < n && data[i]; ++i) {
		bam_hdr_destroy(data[i]->hdr);
		sam_close(data[i]->fp);
		if (data[i]->iter) hts_itr_destroy(data[i]->iter);
		free(data[i]);
	}
	free(data); free(reg);
	if (bed) bed_destroy(bed);
    if ( file_list )
    {
        for (i=0; i<n; i++) free(fn[i]);
        free(fn);
    }
	return status;
}

