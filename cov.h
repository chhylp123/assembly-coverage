#ifndef __RCUT__
#define __RCUT__
#include <stdio.h>
#include <stdint.h>

#define CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))
#define MAX(x, y) (((x) >= (y))?(x):(y))
#define MIN(x, y) (((x) <= (y))?(x):(y))
#define OVL(s_0, e_0, s_1, e_1) ((MIN((e_0), (e_1)) > MAX((s_0), (s_1)))? MIN((e_0), (e_1)) - MAX((s_0), (s_1)):0)

typedef struct{
    int *l, n; 
    char **a;
}enzyme;

typedef struct {
	char *in, *t_in_1, *t_in_2;
    char *q_gfa, *r_gfa, *yak_hh;
    char *chr, *cen;
    enzyme *cpafs, *upafs, *cgfas, *ugfas;
    uint64_t chrS, chrE, covX, minLen, mis_assembLen, hhStep;
    int64_t mis_penalty, match_score, plot, t2t, t2tlen;
    uint8_t is_uncov, is_cov, p_only, is_m_sv;
    double sec_rate, N_rate;
    char *fin_bed, *fin_asm, *fin_con, *fmask_bed, *fmask_asm, *fmask_paf;
    int64_t flen;
    double fcut;
}cov_opt_t;

void cov_main(const cov_opt_t *opt);
void chr_bin_main(const cov_opt_t *opt);
void trio_cov_main(const cov_opt_t *opt);
void trio_cov_main_gfa(const cov_opt_t *opt);
void cov_bp_check(const cov_opt_t *opt);
void contig_bin_main(const cov_opt_t *opt, char **pafs, uint32_t pafn);
void detect_gene_main(const cov_opt_t *opt);
void replace_Ns(const cov_opt_t *opt, char **fs, uint32_t fn);
void plot_merge_Svs(const cov_opt_t *opt, char *fn);
void t2t_main(const cov_opt_t *opt, char *paf);
void flagger_filter_main(const cov_opt_t *opt, char *fin_bed, char *fin_assemb, char *fin_con, 
char *fmask_bed, char *fmask_assemb, char *fmask_paf, int64_t flen, double fcut);

#endif