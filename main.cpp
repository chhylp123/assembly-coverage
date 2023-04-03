#include <stdio.h>
#include <stdlib.h>
#include "cov.h"
#include "ketopt.h"


static ko_longopt_t long_options[] = {
	{ "cpafs",       ko_required_argument, 300 },
	{ "upafs",       ko_required_argument, 301 },
	{ "cgfas",     ko_required_argument, 302 },
	{ "ugfas",      ko_required_argument, 303 },
	{ "Flen",      ko_required_argument, 304 },
	{ "Fasm",      ko_required_argument, 305 },
	{ "Fcon",      ko_required_argument, 306 },
	{ "FmaskBed",      ko_required_argument, 307 },
	{ "FmaskPaf",      ko_required_argument, 308 },
	{ "FmaskCut",      ko_required_argument, 309 },
	{ "Fmaskasm",      ko_required_argument, 310 },
	{ "Fmerge",      ko_required_argument, 311 },
	{ 0, 0, 0 }
};

static void print_usage(FILE *fp, const cov_opt_t *opt)
{
	fprintf(fp, "Usage: cov_cal [options] <in.paf>\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  -c STR      chromosome id [NULL]\n");
	fprintf(fp, "  -s INT      interval beg pos [%lu]\n", opt->chrS);
	fprintf(fp, "  -e INT      interval end pos [%lu]\n", opt->chrE);
	fprintf(fp, "  -x INT      coverage threshold [%lu]\n", opt->covX);
	fprintf(fp, "  -m INT      alignment length threshold [%lu]\n", opt->minLen);
	fprintf(fp, "  -A INT      match score [%ld]\n", opt->match_score);
	fprintf(fp, "  -B INT      mismatch penalty [%ld]\n", opt->mis_penalty);
	fprintf(fp, "  -b INT      minimum size of miss-assembly size\n");
	fprintf(fp, "  -u          print uncovered regions [disable]\n");
	fprintf(fp, "  -a          print covered regions [disable]\n");
	fprintf(fp, "  -p          only consider primary alignment in paf file\n");
	fprintf(fp, "  -1 STR      paf file to hap1\n");
	fprintf(fp, "  -2 STR      paf file to hap2\n");
	fprintf(fp, "  -q STR      gfa file to query\n");
	fprintf(fp, "  -r STR      gfa file to trio reference\n");
	fprintf(fp, "  -E STR      yak hamming log\n");
	fprintf(fp, "  -S INT      window size for hamming evaluation\n");
	fprintf(fp, "  -H 	       regard one contig as a whole window; should be 'cov_cal paf1:yak_log1 paf2:yak_log2'\n");
	fprintf(fp, "  -C STR      centromere bed files\n");
	fprintf(fp, "  -P 	       directly output plotting file, instead of intermediate file\n");
	fprintf(fp, "Detect genes:\n");
	fprintf(fp, "  --cpafs/-5  [paf1,paf2,...]\n");
	fprintf(fp, "        	   paf files of gene to contig\n");
	fprintf(fp, "  --upafs/-6  [paf1,paf2,...]\n");
	fprintf(fp, "        	   paf files of gene to unitig\n");
	fprintf(fp, "  --cgfas/-7  [gfa1,gfa2,...]\n");
	fprintf(fp, "        	   gfa files of contigs\n");
	fprintf(fp, "  --ugfas/-7  [gfa1,gfa2,...]\n");
	fprintf(fp, "        	   gfa files of unitigs\n");
	fprintf(fp, "Auxiliary functions:\n");
	fprintf(fp, "  -N          randomly replace bases with Ns\n");
	fprintf(fp, "Merged SVs Functions:\n");
	fprintf(fp, "  -M          print core-SVs/pan-SVs plots\n");
	fprintf(fp, "T2T detection:\n");
	fprintf(fp, "  -T          print T2T contigs from minigraph paf\n");
	fprintf(fp, "  -t INT      print chains longer than INT [%ld]\n", opt->t2tlen);
	fprintf(fp, "Flagger process:\n");
	fprintf(fp, "  -F          STR\n");
	fprintf(fp, "        	   flagger results in bed format\n");
	fprintf(fp, "  --Flen 	   INT\n");
	fprintf(fp, "        	   only consider contigs/scaffolds longer than INT [%ld]\n", opt->flen);
	fprintf(fp, "  --Fasm 	   STR\n");
	fprintf(fp, "        	   assembly file to detect gap coordinates\n");
	fprintf(fp, "  --Fcon 	   STR\n");
	fprintf(fp, "        	   only consider within the bed file\n");
	fprintf(fp, "  --FmaskBed  STR\n");
	fprintf(fp, "        	   the bed file of another flagger result; mask with --FmaskPaf\n");
	fprintf(fp, "  --FmaskPaf  STR\n");
	fprintf(fp, "        	   the overlap file to another flagger result; work with --FmaskBed\n");
	fprintf(fp, "  --FmaskCut  FLOAT\n");
	fprintf(fp, "        	   ignore flagger errors if the rate comparing with masked regions less than FLOAT [%f]\n", opt->fcut);
	fprintf(fp, "  --Fmaskasm  STR\n");
	fprintf(fp, "        	   assembly file of another flagger result; work with --FmaskBed\n");
	// fprintf(fp, "  --Fmerge    FLOAT");
	// fprintf(fp, "        	   merge multiple flagger regions if FLOAT bases are errors [%f]\n", opt->t2tlen);
}

static void init_opt(cov_opt_t *opt)
{
	opt->fin_bed = NULL;
	opt->fin_asm = NULL;
	opt->fin_con = NULL;
	opt->fmask_bed = NULL;
	opt->fmask_paf = NULL;
	opt->fmask_asm = NULL;
	opt->flen = 0;
	opt->fcut = 1.1;

	opt->in = NULL;
	opt->t_in_1 = NULL;
	opt->t_in_2 = NULL;
	opt->yak_hh = NULL;
	opt->chr = NULL;
	opt->chrS = 0;
	opt->chrE = (uint64_t)-1;
	opt->covX = 1;
	opt->minLen = 100000;
	opt->is_uncov = 0;
	opt->is_cov = 0;
	opt->p_only = 0;
	opt->mis_assembLen = 1000000;
	opt->mis_penalty = 3;
	opt->match_score = 1;
	opt->q_gfa = NULL;
	opt->r_gfa = NULL;
	opt->sec_rate = 0.333333;
	opt->hhStep = 100000;
	opt->N_rate = -1;
	opt->is_m_sv = 0;
	opt->cpafs = NULL;
	opt->upafs = NULL;
	opt->cgfas = NULL;
	opt->ugfas = NULL;
	opt->cen = NULL;
	opt->plot = 0;
	opt->t2t = 0;
	opt->t2tlen = INT64_MAX;
}

void get_hic_enzymes(char *argv, enzyme** x, int check_name)
{
	// fprintf(stderr, "%s\n", argv);
    int i, k, pre_i, len = strlen(argv);
    (*x) = (enzyme*)calloc(1, sizeof(enzyme));
    if(len == 0)
    {
        (*x)->n = 0; (*x)->l = NULL; (*x)->a = NULL;
        return;
    }


    (*x)->n = 1;
    for (i = pre_i = 0; i < len; i++)
    {
        if(argv[i] == ',')
        {
            (*x)->n++;
            continue;
        } 
        
        if(check_name)
        {
            if(argv[i] != 'A' && argv[i] != 'C' && argv[i] != 'G' && argv[i] != 'T' &&
                argv[i] != 'a' && argv[i] != 'c' && argv[i] != 'g' && argv[i] != 't' && 
                argv[i] != 'N' && argv[i] != 'n')
            {
                (*x)->n = 0;
                (*x)->l = NULL;
                (*x)->a = NULL;
                return;
            }
        }
        
    }
    (*x)->l = (int*)calloc((*x)->n, sizeof(int));
    (*x)->a = (char**)calloc((*x)->n, sizeof(char*));

    for (i = pre_i = k = 0; i < len; i++)
    {
        if(argv[i] == ',')
        {
            (*x)->l[k] = i - pre_i;
            (*x)->a[k] = (char*)malloc(sizeof(char)*((*x)->l[k]+1));
            memcpy((*x)->a[k], argv + pre_i, (*x)->l[k]);
            (*x)->a[k][(*x)->l[k]] = '\0';
            pre_i = i + 1;
            k++;
        }
    }

    (*x)->l[k] = i - pre_i;
    (*x)->a[k] = (char*)malloc(sizeof(char)*((*x)->l[k]+1));
    memcpy((*x)->a[k], argv + pre_i, (*x)->l[k]);
    (*x)->a[k][(*x)->l[k]] = '\0';
}

int main(int argc, char *argv[])
{
	// int i, ret;
	int c;
	ketopt_t o = KETOPT_INIT;
	cov_opt_t opt;
	init_opt(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "c:s:e:x:m:ua1:2:pb:B:q:r:E:S:H5:6:7:8N:C:MPTt:F:", long_options)) >= 0) {
		if (c == 'c') opt.chr = o.arg;
		else if (c == 's') opt.chrS = atoll(o.arg);
		else if (c == 'e') opt.chrE = atoll(o.arg);
		else if (c == 'x') opt.covX = atoll(o.arg);
		else if (c == 'm') opt.minLen = atoll(o.arg);
		else if (c == 'u') opt.is_uncov = 1;
		else if (c == 'a') opt.is_cov = 1;
		else if (c == '1') opt.t_in_1 = o.arg;
		else if (c == '2') opt.t_in_2 = o.arg;
		else if (c == 'p') opt.p_only = 1;
		else if (c == 'b') opt.mis_assembLen = atoll(o.arg);
		else if (c == 'B') opt.mis_penalty = atoll(o.arg);
		else if (c == 'A') opt.match_score = atoll(o.arg);
		else if (c == 'q') opt.q_gfa = o.arg;
		else if (c == 'r') opt.r_gfa = o.arg;
		else if (c == 'E') opt.yak_hh = o.arg;
		else if (c == 'S') opt.hhStep = atoll(o.arg);
		else if (c == 'H') opt.hhStep = (uint64_t)-1;
		else if (c == 'N') opt.N_rate = atof(o.arg);
		else if (c == 'M') opt.is_m_sv = 1;
		else if (c == 'C') opt.cen = o.arg;
		else if (c == 'P') opt.plot = 1;
		else if (c == 'T') opt.t2t = 1;
		else if (c == 't') opt.t2tlen = atoll(o.arg);
		else if (c == '5' || c == 300) get_hic_enzymes(o.arg, &(opt.cpafs), 0);
        else if (c == '6' || c == 301) get_hic_enzymes(o.arg, &(opt.upafs), 0);
        else if (c == '7' || c == 302) get_hic_enzymes(o.arg, &(opt.cgfas), 0);
		else if (c == '8' || c == 303) get_hic_enzymes(o.arg, &(opt.ugfas), 0);
		else if (c == 'F') opt.fin_bed = o.arg;
		else if(c == 304) opt.flen = atoll(o.arg);
		else if(c == 305) opt.fin_asm = o.arg;
		else if(c == 306) opt.fin_con = o.arg;
		else if(c == 307) opt.fmask_bed = o.arg;
		else if(c == 308) opt.fmask_paf = o.arg;
		else if(c == 309) opt.fcut = atof(o.arg);
		else if(c == 310) opt.fmask_asm = o.arg;
	}
	if(opt.fin_bed) {
		flagger_filter_main(&opt, opt.fin_bed, opt.fin_asm, opt.fin_con, opt.fmask_bed, opt.fmask_asm, 
		opt.fmask_paf, opt.flen, opt.fcut);
		return 0;
	}
	if (o.ind == argc)
	{
		if(opt.cpafs && opt.upafs && opt.cgfas && opt.ugfas) {
			detect_gene_main(&opt);
			return 1;
		}
		print_usage(stderr, &opt);
		return 0;
	}
	opt.in = argv[o.ind];
	if(opt.t2t) t2t_main(&opt, argv[o.ind]);
	else if(opt.q_gfa && opt.r_gfa) trio_cov_main_gfa(&opt);
	else if(opt.t_in_1 && opt.t_in_2) trio_cov_main(&opt);
	else if(opt.hhStep == (uint64_t)-1) contig_bin_main(&opt, argv + o.ind, argc - o.ind);
	else if(opt.yak_hh) chr_bin_main(&opt);
	else if(opt.N_rate >= 0) replace_Ns(&opt, argv + o.ind, argc - o.ind);
	else if(opt.is_m_sv) plot_merge_Svs(&opt, argv[o.ind]);
	// else cov_bp_check(&opt);
	else cov_main(&opt);

    return 1;
}
