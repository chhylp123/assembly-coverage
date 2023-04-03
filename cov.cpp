#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <stdlib.h>
#include <zlib.h>
#include <ctype.h>
#include <time.h>
#include "cov.h"
#include "ksort.h"
#include "kvec.h"
#include "khashl.h"
#include "kseq.h"
#include "assert.h"
KSEQ_INIT(gzFile, gzread)
KHASHL_MAP_INIT(, nidx_t, nidx, const char*, int32_t, kh_hash_str, kh_eq_str)
KSORT_INIT_GENERIC(uint64_t)

///Match
#define MF 0   
///Insertion
#define IF 1
///Deletion
#define DF 2
///Softclip
#define SF 3
///Hardclip
#define HF 4

typedef struct { size_t n, m; uint64_t *a; } asg64_v;
typedef struct { size_t n, m; uint32_t *a; } asg32_v;

typedef struct {
    size_t n, m;
    uint64_t *a;
}kv_uint64_t_warp;

typedef struct {
    size_t n, m;
    uint32_t *a;
}kv_uint32_t_warp;


typedef struct {
    size_t n, m;
    int8_t *a;
}kv_int8_t_warp;

typedef struct {
    size_t n, m;
    char *a;
    FILE *fp;
}rc_t;

typedef struct {
	char *name;
} nm_t;

typedef struct {
    kvec_t(nm_t) n;
	void *h_name;
    kvec_t(uint8_t) tf;
    kvec_t(uint64_t) l;
} q_idx_t;


typedef struct {
    uint64_t qn, qs, qe, ql;
    uint64_t tn, ts, te, tl;
    uint64_t mlen, mq, cidx, clen;
    uint8_t tf;
    uint8_t is_p, is_i; 
    uint8_t rev;
    int64_t err;
} interval_t;

typedef struct {
    size_t n, m;
    interval_t *a;
    kvec_t(uint32_t) cg;
}ist_t;


typedef struct {
    uint64_t qn, qs, qe;
    uint8_t type;
} bed_t;

typedef struct {
    size_t n, m;
    bed_t *a;
    kvec_t(uint64_t) idx;
}vbed_t;

typedef struct {
    uint64_t s, e, l, id, pl;
    uint64_t chr_id, chr_len;
    char *chr;
} ref_in;

typedef struct {
    size_t n, m;
    ref_in *a;
}ref_in_t;

typedef struct {
    uint64_t cov[2], cov_a;
    uint64_t chr_id, tlen, qid;
    char *chr;
} ref_cov;

typedef struct {
    size_t n, m;
    ref_cov *a;
}ref_cov_t;

#define bed_qn_key(x) ((x).qn)
KRADIX_SORT_INIT(bed_qn, bed_t, bed_qn_key, member_size(bed_t, qn))
#define bed_qs_key(x) ((x).qs)
KRADIX_SORT_INIT(bed_qs, bed_t, bed_qs_key, member_size(bed_t, qs))

#define ist_tn_key(x) ((x).tn)
KRADIX_SORT_INIT(ist_tn, interval_t, ist_tn_key, member_size(interval_t, tn))
#define ist_mq_key(x) ((x).mq)
KRADIX_SORT_INIT(ist_mq, interval_t, ist_mq_key, member_size(interval_t, mq))
#define ist_ts_key(x) ((x).ts)
KRADIX_SORT_INIT(ist_ts, interval_t, ist_ts_key, member_size(interval_t, ts))
#define ist_qn_key(x) ((x).qn)
KRADIX_SORT_INIT(ist_qn, interval_t, ist_qn_key, member_size(interval_t, qn))
#define ist_tqn_key(x) (((x).tn<<32)|(x).qn)
KRADIX_SORT_INIT(ist_tqn, interval_t, ist_tqn_key, 8)
#define ist_qs_key(x) ((x).qs)
KRADIX_SORT_INIT(ist_qs, interval_t, ist_qs_key, member_size(interval_t, qs))
#define ist_qe_key(x) ((x).qe)
KRADIX_SORT_INIT(ist_qe, interval_t, ist_qe_key, member_size(interval_t, qe))

#define ist_err_key(x) ((x).err)
KRADIX_SORT_INIT(ist_err, interval_t, ist_err_key, member_size(interval_t, err))
#define ist_mlen_key(x) ((x).mlen)
KRADIX_SORT_INIT(ist_mlen, interval_t, ist_mlen_key, member_size(interval_t, mlen))
#define ref_cov_chr_id_key(x) ((x).chr_id)
KRADIX_SORT_INIT(ref_cov_chr_id, ref_cov, ref_cov_chr_id_key, member_size(ref_cov, chr_id))
#define ref_cov_chr_cov_key(x) ((uint64_t)((uint64_t)(!!((x).cov[0]<=(x).cov[1]))<<63) | (uint64_t)(MAX((x).cov[0], (x).cov[1]) - MIN((x).cov[0], (x).cov[1])))
KRADIX_SORT_INIT(ref_cov_chr_cov, ref_cov, ref_cov_chr_cov_key, 8)

typedef struct {
	uint64_t s, e, dp;
} ma_sub_t;

#define ma_dp_key(x) ((x).dp)
KRADIX_SORT_INIT(ma, ma_sub_t, ma_dp_key, member_size(ma_sub_t, dp))

#define ma_s_key(x) ((x).s)
KRADIX_SORT_INIT(s, ma_sub_t, ma_s_key, member_size(ma_sub_t, s))

typedef struct {
    size_t n, m;
	ma_sub_t *a;
    kvec_t(uint64_t) idx;
} kv_warp_ma_sub_t;

typedef struct {
	int64_t hL[2], len, w, id;
} hapLen_t;

typedef struct {
	kvec_t(ma_sub_t) Ns;
    kvec_t(char) seq;
    uint64_t id;
} assemb_t;

typedef struct {
    size_t n, m;
	assemb_t *a;
} vec_assemb_t;

#define hapLen_key(x) ((x).w)
KRADIX_SORT_INIT(hapLen, hapLen_t, hapLen_key, member_size(hapLen_t, w))

typedef struct {
	uint32_t len:31, circ:1; // len: length of the unitig; circ: circular if non-zero
	uint32_t start, end; // start: starting vertex in the string graph; end: ending vertex
	uint32_t m, n; // number of reads
	uint64_t *a; // list of reads
	char *s; // unitig sequence is not null
} ma_utg_t;

typedef struct {
	uint64_t ul;
	uint32_t v;
	uint32_t ol:31, del:1;
	uint8_t strong;
	uint8_t el;
	uint8_t no_l_indel;
} asg_arc_t;

typedef struct {
	uint32_t len:31, del:1;
	uint8_t c;
} asg_seq_t;

typedef struct {
	uint32_t m_arc, n_arc:31, is_srt:1;
	asg_arc_t *arc;
	uint32_t m_seq, n_seq:31, is_symm:1;
	uint32_t r_seq;

	asg_seq_t *seq;
	uint64_t *idx;

	uint8_t* seq_vis;

	uint32_t n_F_seq;
	ma_utg_t* F_seq;
} asg_t;

typedef struct { size_t n, m; ma_utg_t *a;} ma_utg_v;

typedef struct {
	ma_utg_v u;
	asg_t *g;
	kvec_t(uint64_t) occ;
} ma_ug_t;

int rc_t_init(rc_t *kt, const char *fn)
{
    memset(kt, 0, sizeof(rc_t));
    kt->fp = fopen(fn, "r");
    if (!(kt->fp)) {
        return 0;
    }
    return 1;
}

void rc_t_des(rc_t *kt)
{
    kv_destroy(*kt);
    fclose(kt->fp);
}

int getLine(rc_t *kt)
{
    char t;
    kt->n = 0;
    while (1)
    {
        fread(&t, 1, 1, kt->fp);
        if(feof(kt->fp)) break;
        kv_push(char, *kt, t);
        if(t == '\n') break;
    }

    if(kt->n > 0)
    {
        t = '\0';
        kv_push(char, *kt, t);
        return 1;
    } 
    
    if(feof(kt->fp)) return 0;

    t = '\0';
    kv_push(char, *kt, t);
    return 1;
}

char *nidx_strdup(const char *src)
{
	int32_t len;
	char *dst;
	len = strlen(src);
	MALLOC(dst, len + 1);
	memcpy(dst, src, len + 1);
	return dst;
}

uint32_t nidx_append(q_idx_t *g, const char *name, uint8_t tf)
{
	nidx_t *h = (nidx_t*)g->h_name;
	int absent;
	khint_t k;
	k = nidx_put(h, name, &absent);
	if (absent) {
		nm_t *t = NULL;
        kv_pushp(nm_t, g->n, &t);
        kv_push(uint8_t, g->tf, tf);
        kv_push(uint64_t, g->l, (uint64_t)-1);
		kh_key(h, k) = t->name = nidx_strdup(name);
		kh_val(h, k) = g->n.n - 1;
	}
    g->tf.a[kh_val(h, k)] |= tf;
	return kh_val(h, k);
}

uint32_t nidx_query(q_idx_t *g, const char *name)
{
	nidx_t *h = (nidx_t*)g->h_name;
    khint_t n_buckets = (h->keys? 1U<<h->bits : 0U);
	khint_t k = nidx_get(h, name);
    if(k < n_buckets) return kh_val(h, k);
    return (uint32_t)-1;
}

q_idx_t *q_idx_t_init()
{
    q_idx_t *qI = NULL;
    CALLOC(qI, 1);
    qI->h_name = nidx_init();
    kv_init(qI->tf);
    kv_init(qI->l);
    return qI;
}

void q_idx_t_des(q_idx_t **qI)
{
    uint64_t i;
    nidx_destroy((nidx_t*)(*qI)->h_name);
    for (i = 0; i < (*qI)->n.n; ++i) free((*qI)->n.a[i].name);
    kv_destroy((*qI)->n); 
    kv_destroy((*qI)->tf); 
    kv_destroy((*qI)->l);
    free((*qI));
}

int cmp_ref_cov(const void * a, const void * b)
{
    return strcmp((*(ref_cov*)a).chr, (*(ref_cov*)b).chr);
}

void push_cigar_paf(char *in, asg32_v *res)
{
    uint64_t occ = strlen(in), k, l, cl, cn, clen; char cc;
    res->n = 0;
    if(!occ) return;
    for (l = 0, k = 1; (k < occ) && (in[k] != '\n'); k++) {
        if(/**k == occ ||**/ (!isdigit(in[k]))) {
            // if(!(k > l)) {
            //     fprintf(stderr, "%s\n", in);
            //     fprintf(stderr, "k::%lu, l::%lu, occ::%lu, in[k]::%c, in[k]::%u\n", k, l, occ, in[k], (uint32_t)in[k]);
            // }
            assert(k > l);
            cc = in[k]; in[k] = '\0';
            cl = atoll(in + l); in[k] = cc;
            if(in[k] == 'M') {
                cn = MF;
            } else if(in[k] == 'I') {
                cn = IF;
            } else if(in[k] == 'D') {
                cn = DF;
            } else if(in[k] == 'S') {
                cn = SF;
            } else if(in[k] == 'H') {
                cn = HF;
            } else {
                fprintf(stderr, "illegal cigar base::%c", in[k]);
                exit(1);
            }
            if(cn!=4) {///hard-clip does not affect coordinate
                while (cl) {
                    if(cl <= (0x3fffffff)) {
                        clen = cl; cl = 0;
                    } else {
                        clen = (0x3fffffff); cl -= (0x3fffffff);
                    }
                    clen <<= 2; clen |= cn;
                    kv_push(uint32_t, *res, clen);
                }
            }
            l = k+1;
        }
    }
}

uint64_t push_ist_t_cigar(ist_t *l, interval_t *t, asg32_v *ciagr)
{
    int64_t s = 0, e = ((int64_t)(ciagr->n))-1, scut, ecut, ns, ne;
    for (scut = 0; (s <= e) && ((ciagr->a[s]&3) == SF); s++) scut += (ciagr->a[s]>>2);
    for (ecut = 0; (e >= s) && ((ciagr->a[s]&3) == SF); e--) ecut += (ciagr->a[s]>>2);
    // assert(((int64_t)(t->qe - t->qs)) >= (e - s));
    ns = t->qs + scut; ne = ((int64_t)t->qe) - ecut;
    assert(ns <= (int64_t)t->ql && ne <= (int64_t)t->ql && ns >= 0 && ne >= 0);
    t->qs = ns; t->qe = ne; t->cidx = l->cg.n; t->clen = e + 1 - s;
    kv_resize(uint32_t, l->cg, l->cg.n + t->clen);
    // fprintf(stderr, "[M::%s::] t->cidx::%lu, t->clen::%lu, l->cg.n::%lu, l->cg.m::%lu, ciagr->n::%lu\n", __func__, 
    // t->cidx, t->clen, (uint64_t)l->cg.n, (uint64_t)l->cg.m, (uint64_t)ciagr->n);
    memcpy(l->cg.a + t->cidx, ciagr->a + s, t->clen*sizeof((*(l->cg.a))));
    l->cg.n += t->clen;
    // fprintf(stderr, "[M::%s::] t->cidx::%lu, t->clen::%lu, l->cg.n::%lu\n", __func__, t->cidx, t->clen, (uint64_t)l->cg.n);
    return 1;
}

uint32_t load_paf(const cov_opt_t *opt, const char* fn, q_idx_t *qI, q_idx_t *tI, ist_t *l, uint8_t tf, uint8_t filter_def, uint8_t append_hapid, uint8_t keep_cigar)
{
    rc_t kt; memset(&kt, 0, sizeof(kt));
    uint64_t i, ts, te, tl, ql, qs, qe, mq;
    int64_t err;
    uint8_t is_p, is_i, rev;
    char *pch = NULL;
    char *qn = NULL;
    char *tn = NULL;
    char *def = NULL;
    interval_t *t = NULL;
    asg32_v ciagr; kv_init(ciagr);
    
    kvec_t(char) qnA; kv_init(qnA);
    if(!rc_t_init(&kt, fn)) {
        rc_t_des(&kt);
        fprintf(stderr, "[M::%s::] ==> Cnnot open paf %s\n", __func__, fn);
        kv_destroy(qnA);
        return 0;
    }

    while (getLine(&kt)) {
        i = 0; is_p = 0; is_i = 0; err = 0; rev = 0; mq = 0; ciagr.n = 0;
        ts = te = tl = ql = qs = qe = (uint64_t)-1; def = NULL;
        pch = strtok(kt.a, "\t");
        while (pch != NULL)
        {
            if(i == 0) qn = pch;
            if(i == 1) ql = atoll(pch);
            if(i == 2) qs = atoll(pch);
            if(i == 3) qe = atoll(pch);
            if(i == 4 && pch[0] == '-') rev = 1;
            if(i == 5) tn = pch;
            if(i == 6) tl = atoll(pch);
            if(i == 7) ts = atoll(pch);
            if(i == 8) te = atoll(pch);
            if(i == 11) mq = atoll(pch);
            if(strlen(pch) >= 6) {
                if(!memcmp(pch, "tp:A:P", 6))
                {
                    is_p = 1;
                }
                else if(!memcmp(pch, "tp:A:I", 6))
                {
                    is_i = 1;
                }
                else if(!memcmp(pch, "NM:i:", 5))
                {
                    err = atoll(pch+5);
                }
                else if(filter_def && !memcmp(pch, "de:f:", 5))
                {
                    def = pch;
                }
            }
            if(keep_cigar && strlen(pch) > 5) {
                if(!memcmp(pch, "cg:Z:", 5)) push_cigar_paf(pch+5, &ciagr);
            }

            pch = strtok (NULL, "\t");
            i++;
        }
        if(opt->chr && strcmp(tn, opt->chr)) continue;
        if(OVL(ts, te, opt->chrS, opt->chrE) == 0) continue;
        if(te - ts < opt->minLen) continue;
        if(qe - qs < opt->minLen) continue;
        if(te <= ts || qe <= qs) continue;
        if(filter_def && def == NULL) continue;
        kv_pushp(interval_t, *l, &t);
        if(append_hapid) {
            kv_resize(char, qnA, strlen(qn)+128); sprintf(qnA.a, "%s_hap%u", qn, tf);
            t->qn = nidx_append(qI, qnA.a, tf);
        }
        else {
            t->qn = nidx_append(qI, qn, tf);
        }
        t->ql = ql; t->qs = qs; t->qe = qe;
        t->tn = nidx_append(tI, tn, tf);
        t->tl = tl; t->ts = ts; t->te = te;
        t->tf = tf; t->is_p = is_p; t->is_i = is_i;
        t->err = err; t->cidx = t->clen = (uint64_t)-1;
        t->rev = rev; t->mq = mq;
        if(ciagr.n) push_ist_t_cigar(l, t, &ciagr);
        if((t->qe <= t->qs) || (t->te <= t->ts)) l->n--;
    }

    rc_t_des(&kt); kv_destroy(qnA); kv_destroy(ciagr);
    fprintf(stderr, "[M::%s::] ==> Loaded %s\n", __func__, fn);
    return 1;
}

uint32_t load_mul_paf(const cov_opt_t *opt, const enzyme *fn, q_idx_t *qI, q_idx_t *tI, ist_t *l, uint8_t filter_def)
{
    uint64_t i;
    for (i = 0; i < (uint64_t)fn->n; i++) {
        if(!load_paf(opt, fn->a[i], qI, tI, l, i, filter_def, 0, 0)) {
            fprintf(stderr, "Wrong file name: %s\n", fn->a[i]);
            return 0;
        }
    }
    return 1;
}

uint64_t get_ol(char *cigar)
{
    int64_t ov = 0;
    char *p = cigar, *q = NULL, *r = NULL;
    ov = strtol(p, &r, 10);
    if (isupper(*r)) { // CIGAR
        ov = 0; q = p;
        do {
            long l;
            l = strtol(q, &q, 10);
            if (*q == 'M' || *q == 'D' || *q == 'N') ov += l;
            // if (*q == 'M' || *q == 'D' || *q == 'N') ov += l;
            // if (*q == 'M' || *q == 'I' || *q == 'S') ow += l;
            ++q;
        } while (isdigit(*q));
        return ov;
    }
    return (uint64_t)-1;
}

asg_t *asg_init(void)
{
	return (asg_t*)calloc(1, sizeof(asg_t));
}

void asg_destroy(asg_t *g)
{
	if (g == 0) return;
	free(g->seq); free(g->idx); free(g->arc); free(g->seq_vis);
     
    if(g->n_F_seq > 0 && g->F_seq)
    {
        uint32_t i = 0;
        for (i = 0; i < g->n_F_seq; i++)
        {
            if(g->F_seq[i].a) free(g->F_seq[i].a);
            if(g->F_seq[i].s) free(g->F_seq[i].s);
        }

        free(g->F_seq);
    }
    

    free(g);
}


void asg_seq_set(asg_t *g, int sid, int len, int del)
{
	///just malloc size
	if (sid >= (int)g->m_seq) {
		g->m_seq = sid + 1;
		kv_roundup32(g->m_seq);
		g->seq = (asg_seq_t*)realloc(g->seq, g->m_seq * sizeof(asg_seq_t));
	}

	if (sid >= g->n_seq) g->n_seq = sid + 1;
	
    g->seq[sid].del = !!del;
    g->seq[sid].len = len;
}

static inline asg_arc_t *asg_arc_pushp(asg_t *g)
{
	if (g->n_arc == g->m_arc) {
		g->m_arc = g->m_arc? g->m_arc<<1 : 16;
		g->arc = (asg_arc_t*)realloc(g->arc, g->m_arc * sizeof(asg_arc_t));
	}
	return &g->arc[g->n_arc++];
}

#define asg_arc_len(arc) ((uint32_t)(arc).ul)
#define asg_arc_n(g, v) ((uint32_t)(g)->idx[(v)])
#define asg_arc_a(g, v) (&(g)->arc[(g)->idx[(v)]>>32])

#define asg_arc_key(a) ((a).ul)
KRADIX_SORT_INIT(asg, asg_arc_t, asg_arc_key, 8)
#define ug_u_key(a) ((uint32_t)(a))
KRADIX_SORT_INIT(ug_u, uint64_t, ug_u_key, 8)

void asg_arc_sort(asg_t *g)
{
	radix_sort_asg(g->arc, g->arc + g->n_arc);
}

// hard remove arcs marked as "del"
void asg_arc_rm(asg_t *g)
{
	/**
	p->ul: |____________31__________|__________1___________|______________32_____________|
	                    qns            direction of overlap       length of this node (not overlap length)
	p->v : |___________31___________|__________1___________|
				        tns              relative strand between query and target
	p->ol: overlap length
	**/
	uint32_t e, n;
	///just clean arc requiring: 1. arc it self must be available 2. both the query and target are available
	for (e = n = 0; e < g->n_arc; ++e) {
		//u and v is the read id
		uint32_t u = g->arc[e].ul>>32, v = g->arc[e].v;
		if (!g->arc[e].del && !g->seq[u>>1].del && !g->seq[v>>1].del)
			g->arc[n++] = g->arc[e];
	}
	if (n < g->n_arc) { // arc index is out of sync
		if (g->idx) free(g->idx);
		g->idx = 0;
	}
	g->n_arc = n;
}

uint64_t *asg_arc_index_core(size_t max_seq, size_t n, const asg_arc_t *a)
{
	size_t i, last;
	uint64_t *idx;
	idx = (uint64_t*)calloc(max_seq * 2, 8);


    /**
		 * ul: |____________31__________|__________1___________|______________32_____________|
	                       qns            direction of overlap       length of this node (not overlap length)
	**/
    ///so if we use high 32-bit, we store the index of each qn with two direction
	for (i = 1, last = 0; i <= n; ++i)
		if (i == n || a[i-1].ul>>32 != a[i].ul>>32)
			idx[a[i-1].ul>>32] = (uint64_t)last<<32 | (i - last), last = i;

    
	return idx;
}

void asg_arc_index(asg_t *g)
{
	if (g->idx) free(g->idx);
	g->idx = asg_arc_index_core(g->n_seq, g->n_arc, g->arc);
}

void asg_cleanup(asg_t *g)
{
    ///remove overlaps, instead of reads
    ///remove edges with del, and free idx
	asg_arc_rm(g);
	if (!g->is_srt) {
		/**
		 * sort by ul, that is, sort by qns + direction
		 * ul: |____________31__________|__________1___________|______________32_____________|
	                       qns            direction of overlap       length of this node (not overlap length)
		**/
		asg_arc_sort(g);
		g->is_srt = 1;
	}
	///index the overlaps in graph with query id
	if (g->idx == 0) asg_arc_index(g);
}

void usg_seq_set(ma_ug_t *ug, int sid, uint64_t len, int del)
{
    uint64_t i, mm;
	///just malloc size
	if (sid >= (int)ug->g->m_seq) {
        mm = ug->g->m_seq;
		ug->g->m_seq = sid + 1;
		kv_roundup32(ug->g->m_seq);
		ug->g->seq = (asg_seq_t*)realloc(ug->g->seq, ug->g->m_seq * sizeof(asg_seq_t));
        for (i = mm; i < ug->g->m_seq; i++) {
            ug->g->seq[i].len = ((uint32_t)-1)>>1; ug->g->seq[i].del = 1; ug->g->seq[i].c = 0;
        }
	}
    if(sid >= (int)ug->u.m) {
        mm = ug->u.m;
        ug->u.m = sid + 1;
        kv_roundup32(ug->u.m);
        ug->u.a = (ma_utg_t*)realloc(ug->u.a, ug->u.m * sizeof(ma_utg_t));
        memset(ug->u.a + mm, 0, (ug->u.m - mm)*sizeof(ma_utg_t));
    }

	if(sid >= ug->g->n_seq) ug->g->n_seq = sid + 1;
    if(sid >= (int32_t)ug->u.n) ug->u.n = sid + 1;
	
    ug->g->seq[sid].del = !!del;
    if(len != (uint64_t)-1) ug->g->seq[sid].len = len, ug->u.a[sid].len = len;
}

// delete multi-arcs
/**
 * remove edges like:   v has two out-edges to w
**/
int asg_arc_del_multi(asg_t *g)
{
	//the number of nodes are number of read times 2
	uint32_t *cnt, n_vtx = g->n_seq * 2, n_multi = 0, v;
	cnt = (uint32_t*)calloc(n_vtx, 4);
	for (v = 0; v < n_vtx; ++v) {
		///out-nodes of v
		asg_arc_t *av = asg_arc_a(g, v);
		int32_t i, nv = asg_arc_n(g, v);
		///if v just have one out-node, there is no muti-edge
		if (nv < 2) continue;
		for (i = nv - 1; i >= 0; --i) ++cnt[av[i].v];
		for (i = nv - 1; i >= 0; --i)
			if (--cnt[av[i].v] != 0)
				av[i].del = 1, ++n_multi;
	}
	free(cnt);
	if (n_multi) asg_cleanup(g);
    
    fprintf(stderr, "[M::%s] removed %d multi-arcs\n", __func__, n_multi);
	return n_multi;
}

// remove asymmetric arcs: u->v is present, but v'->u' not
int asg_arc_del_asymm(asg_t *g)
{
	uint32_t e, n_asymm = 0;
	///g->n_arc is the number of overlaps 
	for (e = 0; e < g->n_arc; ++e) {
		uint32_t v = g->arc[e].v^1, u = g->arc[e].ul>>32^1;
		uint32_t i, nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		for (i = 0; i < nv; ++i)
			if (av[i].v == u) break;
		if (i == nv) g->arc[e].del = 1, ++n_asymm;
	}
	if (n_asymm) asg_cleanup(g);

	fprintf(stderr, "[M::%s] removed %d asymmetric arcs\n", __func__, n_asymm);
	return n_asymm;
}

void asg_symm(asg_t *g)
{
	asg_arc_del_multi(g);
	asg_arc_del_asymm(g);
	g->is_symm = 1;
}

void ma_ug_destroy(ma_ug_t *ug)
{
	uint32_t i;
	if (ug == 0) return;
	for (i = 0; i < ug->u.n; ++i) {
		free(ug->u.a[i].a);
		free(ug->u.a[i].s);
	}
	free(ug->u.a);
	asg_destroy(ug->g);
    kv_destroy(ug->occ);
	free(ug);
}

ma_ug_t *load_gfa(enzyme* fn, q_idx_t *uidx, q_idx_t *ridx, asg_t *rg)
{
    rc_t kt; memset(&kt, 0, sizeof(kt));
    uint64_t i, k;
    uint64_t rs, src_len, src_d, des_d, rd, r_len, rf, eol, src_id, des_id, r_id;
    char *pch = NULL;
    char *src_n = NULL;
    char *des_n = NULL;
    char *rn = NULL;
    char ff;
    ma_ug_t *ug; ma_utg_t *p; asg_arc_t *e; 
	ug = (ma_ug_t*)calloc(1, sizeof(ma_ug_t));
	ug->g = asg_init();

    for (k = 0; k < (uint64_t)fn->n; k++) {
        // fprintf(stderr, "+[M::%s::] gfa: %s, ug->u.n: %u\n", __func__, fn->a[k], (uint32_t)ug->u.n);
        if(!rc_t_init(&kt, fn->a[k]))
        {
            rc_t_des(&kt);
            fprintf(stderr, "[M::%s::] ==> Cnnot open paf %s\n", __func__, fn->a[k]);
            return 0;
        }
        while (getLine(&kt)) {
            ff = kt.a[0]; i = 0;
            if(ff != 'S' && ff != 'A' && ff != 'L') continue;
            pch = strtok(kt.a, "\t");
            src_n = des_n = rn = NULL; rs = src_d = des_d = rd = rf = (uint64_t)-1;
            src_len = r_len = eol = src_id = des_id = r_id = (uint64_t)-1;
            while (pch != NULL) {
                if(i == 1) {
                    src_n = pch;
                }
                else {
                    if(ff == 'S') {
                        if(strlen(pch) >= 6 && !memcmp(pch, "LN:i:", 5)) src_len = atoll(pch+5);
                    }
                    else if(ff == 'A') {
                        if(i == 2) rs = atoll(pch);
                        if(i == 3) rd = (pch[0]=='+'?0:(pch[0]=='-'?1:(uint64_t)-1));
                        if(i == 4) rn = pch;
                        if(i == 6) r_len = atoll(pch);
                        if(strlen(pch) >= 6 && !memcmp(pch, "HG:A:", 5)) {
                            rf = (pch[5]== 'p'?1:(pch[5]== 'm'?2:0));
                        }
                    }
                    else if(ff == 'L') {
                        if(i == 2) src_d = (pch[0]=='+'?0:(pch[0]=='-'?1:(uint64_t)-1));
                        if(i == 3) des_n = pch;
                        if(i == 4) des_d = (pch[0]=='+'?0:(pch[0]=='-'?1:(uint64_t)-1));
                        if(i == 5) eol = get_ol(pch);
                    }
                }
                pch = strtok (NULL, "\t");
                i++;
            }

            if(src_n) {
                src_id = nidx_append(uidx, src_n, k);
                usg_seq_set(ug, src_id, src_len, 0);
                // fprintf(stderr, "src_n: %s, src_id: %lu\n", src_n, src_id);
            }
            if(des_n) {
                des_id = nidx_append(uidx, des_n, k);
                usg_seq_set(ug, des_id, (uint64_t)-1, 0);
                // fprintf(stderr, "des_n: %s, des_id: %lu\n", des_n, des_id);
            }

            if(rn) {
                r_id = nidx_append(ridx, rn, rf);
                p = &ug->u.a[src_id];
                kv_push(uint64_t, *p, (r_id<<33)|((rd&1)<<32)|rs);
                asg_seq_set(rg, r_id, r_len, 0); rg->seq[r_id].c = 0;
            }

            if(ff == 'L') {
                e = asg_arc_pushp(ug->g);
                memset(e, 0, sizeof(*e));
                e->ol = eol; e->v = (des_id<<1)|(des_d&1);
                e->ul = (src_id<<33)|((src_d&1)<<32);
            }
        }
       rc_t_des(&kt);
       // fprintf(stderr, "-[M::%s::] gfa: %s, ug->u.n: %u\n", __func__, fn->a[k], (uint32_t)ug->u.n);
    }

    for (i = 0; i < ug->g->n_arc; i++) {
        e = &(ug->g->arc[i]);
        e->ul += (ug->g->seq[e->ul>>33].len - e->ol);
    }

    asg_cleanup(ug->g);
    asg_symm(ug->g);

    // fprintf(stderr, "ug->u.n: %u\n\n", (uint32_t)ug->u.n);
    
    for (i = 0; i < ug->u.n; i++) {
        p = &ug->u.a[i];
        radix_sort_ug_u(p->a, p->a + p->n);
        for (k = 0; k + 1 < p->n; k++) {
            p->a[k] = ((p->a[k]>>32)<<32) + ((uint32_t)p->a[k+1] - (uint32_t)p->a[k]);
        }
        p->a[k] = ((p->a[k]>>32)<<32) + (p->len - (uint32_t)p->a[k]);
        uint32_t v = (i<<1);
        asg_arc_t *av = asg_arc_a(ug->g, v);
        uint32_t nv = asg_arc_n(ug->g, v);
        for (k = 0, e = NULL; k < nv; k++) {
            if(av[k].del) continue;
            if(!e) {
                e = &(av[k]);
            }
            else {
                e = NULL;
                break;
            }            
        }

        if(e && e->v == v) {
            v = (i<<1);
            av = asg_arc_a(ug->g, v);
            nv = asg_arc_n(ug->g, v);
            for (k = 0, e = NULL; k < nv; k++) {
                if(av[k].del) continue;
                if(!e) {
                    e = &(av[k]);
                }
                else {
                    e = NULL;
                    break;
                }            
            }
        }

        if(e && e->v == v) p->circ = 1;
        if(!p->circ) {
            p->start = p->a[0]>>32;
            p->end = (p->a[p->n-1]>>32)^1;
        }
        else {
            p->start = p->end = UINT32_MAX;
        }
    }

    return ug;
}


uint32_t load_yak_hh_S(const cov_opt_t *opt, const char* fn, q_idx_t *tI, kv_uint64_t_warp *tI_hh, uint8_t tf, uint8_t append_hapid)
{
    rc_t kt; memset(&kt, 0, sizeof(kt));
    uint64_t i, pat, mat, m;
    char *pch = NULL;
    char *tn = NULL;
    kvec_t(char) tnA; kv_init(tnA);
    if(!rc_t_init(&kt, fn))
    {
        rc_t_des(&kt);
        fprintf(stderr, "[M::%s::] ==> Cannot open paf %s\n", __func__, fn);
        kv_destroy(tnA);
        return 0;
    }


    while (getLine(&kt))
    {
        if(kt.n < 2) continue;
        if(kt.a[0] != 'S' || kt.a[1] != '\t') continue;
        i = 0; pat = mat = (uint64_t)-1;
        pch = strtok(kt.a, "\t");
        while (pch != NULL)
        {
            if(i == 1) tn = pch;
            if(i == 2) pat = atoll(pch);
            if(i == 3) mat = atoll(pch);
            pch = strtok (NULL, "\t");
            i++;
        }
        if(pat == (uint64_t)-1 || mat == (uint64_t)-1) continue;

        if(append_hapid){
            kv_resize(char, tnA, strlen(tn)+128); sprintf(tnA.a, "%s_hap%u", tn, tf);
            i = nidx_append(tI, tnA.a, tf);
        } else {
            i = nidx_append(tI, tn, tf);
        }

        if(tI_hh->m <= i) {
            m = tI_hh->m;
            kv_resize(uint64_t, *tI_hh, i+1);
            memset(tI_hh->a+m, -1, sizeof(uint64_t)*(tI_hh->m-m));
        }
        if(tI_hh->n <= i) tI_hh->n = i+1;
        tI_hh->a[i] = (pat<<32)|mat;
    }

    rc_t_des(&kt); kv_destroy(tnA);
    return 1;
}


uint32_t load_yak_hh(const cov_opt_t *opt, const char* fn, q_idx_t *tI, ist_t *l, uint8_t tf)
{
    rc_t kt; memset(&kt, 0, sizeof(kt));
    uint64_t i, ts, te, tt, tc;
    char *pch = NULL;
    char *tn = NULL;
    interval_t *t = NULL;
    if(!rc_t_init(&kt, fn))
    {
        rc_t_des(&kt);
        fprintf(stderr, "[M::%s::] ==> Cannot open paf %s\n", __func__, fn);
        return 0;
    }

    while (getLine(&kt))
    {
        if(kt.a[0] != 'F') continue;
        i = 0; tt = 0; tc = 0;
        ts = te = (uint64_t)-1;
        pch = strtok(kt.a, "\t");
        while (pch != NULL)
        {
            if(i == 1) tn = pch;
            if(i == 2) tt = atoll(pch);
            if(i == 3) ts = atoll(pch);
            if(i == 4) te = atoll(pch);
            if(i == 5) tc = atoll(pch);
            pch = strtok (NULL, "\t");
            i++;
        }
        kv_pushp(interval_t, *l, &t);
        memset(t, 0, sizeof(*t));
        t->tn = nidx_append(tI, tn, tf);
        t->mq = tt; t->tl = tc; 
        t->ts = ts; t->te = te;
    }

    rc_t_des(&kt);
    return 1;
}

void cov_main(const cov_opt_t *opt)
{
    q_idx_t *qI = q_idx_t_init();
    q_idx_t *tI = q_idx_t_init();
    ist_t l; memset(&l, 0, sizeof(l));///all overlaps
    uint64_t i, c_tl, k, m, t_c_tl = 0, aligned_tl = 0;
    kvec_t(uint64_t) b; kv_init(b);
    kvec_t(ma_sub_t) res; kv_init(res);
    kvec_t(ma_sub_t) contigs; kv_init(contigs);
    ma_sub_t *p = NULL;
    int64_t dp, old_dp, start = 0, min_dp = opt->covX;
    if(!load_paf(opt, opt->in, qI, tI, &l, 1, 0, 0, 0)) goto des_cov_main;///load all overlaps
    radix_sort_ist_tn(l.a, l.a+l.n);
    for (k = 1, m = 0; k <= l.n; ++k)
    {
        if (k == l.n || l.a[k].tn != l.a[m].tn) 
        {
            
            b.n = 0; 
            for (i = m; i < k; i++)
            {
                kv_push(uint64_t, b, l.a[i].ts<<1);
                kv_push(uint64_t, b, l.a[i].te<<1|1);
            }
            ks_introsort_uint64_t(b.n, b.a);


            res.n = 0;
            for (i = 0, dp = 0, start = 0; i < b.n; ++i) 
            {
                old_dp = dp;
                ///if a[j] is qe
                if (b.a[i]&1) 
                {
                    --dp;
                }
                else
                {
                    ++dp;
                } 


                if (old_dp < min_dp && dp >= min_dp) ///old_dp < dp, b.a[j] is qs
                { 
                    ///case 2, a[j] is qs
                    start = b.a[i]>>1;
                } 
                else if (old_dp >= min_dp && dp < min_dp) ///old_dp > min_dp, b.a[j] is qe
                {
                    if(OVL((uint64_t)start, b.a[i]>>1, opt->chrS, opt->chrE) == 0) continue;

                    kv_pushp(ma_sub_t, res, &p);
                    p->s = MAX((uint64_t)start, opt->chrS); p->e = MIN(b.a[i]>>1, opt->chrE);
                    p->dp = old_dp;
                }
            }

            for (i = c_tl = 0; i < res.n; ++i)
            {
                c_tl += res.a[i].e - res.a[i].s;
                // fprintf(stderr, "i-%lu, [%lu, %lu), dp-%lu\n", i, res.a[i].s, res.a[i].e, res.a[i].dp);
            } 

            

            
            t_c_tl += c_tl;
            aligned_tl += l.a[m].tl;
            fprintf(stderr, "\n[%s] len: %lu, covered len: %lu\n", tI->n.a[l.a[m].tn].name, l.a[m].tl, c_tl);
            if(opt->is_cov)
            {
                for (i = 0; i < res.n; ++i)
                {
                    fprintf(stderr, "[#-covered] %s:%lu-%lu\n", tI->n.a[l.a[m].tn].name, res.a[i].s, res.a[i].e);
                }
            }

            if(opt->is_uncov)
            {
                uint64_t st, ed;
                for (i = 0, st = 0, ed = l.a[m].tl; i < res.n; ++i)
                {
                    ed = res.a[i].s;
                    if(ed - st > 0) fprintf(stderr, "[*-uncovered] %s:%lu-%lu\n", tI->n.a[l.a[m].tn].name, st, ed);
                    st = res.a[i].e;
                }

                ed = l.a[m].tl;
                if(ed - st > 0) fprintf(stderr, "[*-uncovered] %s:%lu-%lu\n", tI->n.a[l.a[m].tn].name, st, ed);
            }
            

            kv_pushp(ma_sub_t, contigs, &p);
            p->s = m; p->e = k; p->dp = l.a[m].tl - c_tl;
            m = k;
        }
    }
    fprintf(stderr, "*****3*****\n");

    radix_sort_ma(contigs.a, contigs.a + contigs.n);
    for (i = 0; i < contigs.n; i++)
    {
        p = &(contigs.a[i]); 
        fprintf(stderr, "[%s] len: %lu, covered len: %lu\n", 
                    tI->n.a[l.a[p->s].tn].name, l.a[p->s].tl, l.a[p->s].tl - p->dp);
    }
    fprintf(stderr, "*****4*****\n");
    
    fprintf(stderr, "\n[M::%s] t_c_tl: %lu, aligned_tl: %lu\n", __func__, t_c_tl, aligned_tl);


    

    des_cov_main:
    q_idx_t_des(&qI);
    q_idx_t_des(&tI);
    kv_destroy(l);
    kv_destroy(b);
    kv_destroy(res);
    kv_destroy(contigs);
}

int cmp_mq_t(const void * a, const void * b)
{
    if((*(ma_sub_t*)a).dp != (*(ma_sub_t*)b).dp){
        return (*(ma_sub_t*)a).dp > (*(ma_sub_t*)b).dp? -1:1;
    } else {
        if((*(ma_sub_t*)a).s == (*(ma_sub_t*)b).s) return 0;
        return (*(ma_sub_t*)a).s > (*(ma_sub_t*)b).s? -1 : 1;
    }
}

uint64_t get_eovlp_len(ma_sub_t *a, uint64_t a_n, uint64_t idx, ma_sub_t *vis) 
{
    uint64_t k, sidx = vis[idx].s, eidx = vis[idx].e;
    uint64_t sp = a[sidx].dp>>1, ep = a[eidx].dp>>1, ovlp = 0;
    int64_t dp, old_dp, start = 0;
    if(a[sidx].s != idx || a[eidx].s != idx) fprintf(stderr, "ERROR-debug1\n");
    if(vis[idx].dp) return ep + 1 - sp;
    start = a[sidx].dp>>1; dp = a[sidx].e;
    for (k = sidx; k <= eidx; k++) {
        if(vis[a[k].s].dp == 0) continue;
        old_dp = dp;
        if (a[k].dp&1) --dp;
        else ++dp;

        if (old_dp < 1 && dp >= 1) { ///ts 
            start = a[k].dp>>1;
        } else if (old_dp >= 1 && dp < 1) {///te
            ovlp += OVL((uint64_t)start, (a[k].dp>>1)+1, sp, ep+1);
        }
    }
    if(dp > 0) ovlp += ep + 1 - start;
    return ovlp;
}

void set_eovlp(ma_sub_t *a, uint64_t a_n, uint64_t idx, ma_sub_t *vis)
{
    if(vis[idx].dp) return;
    uint64_t k, sidx = vis[idx].s, eidx = vis[idx].e;
    vis[idx].dp = 1;
    for (k = sidx; k < eidx; k++) a[k].e++;
}

uint64_t extract_ist_t(ma_sub_t *a, uint64_t a_n, ma_sub_t *vis, interval_t *ref, ist_t *res)
{
    uint64_t k, mq, mmq, mmq_i, s, e, sn = res->n, oo, to = 0;
    int64_t dp, old_dp, sk, m, lk;
    interval_t *p = NULL;
    for (k = 0, dp = 0, sk = 0, lk = -1; k < a_n; k++) {
        if(vis[a[k].s].dp == 0) continue;
        old_dp = dp;
        if (a[k].dp&1) --dp;
        else ++dp;
        if(old_dp == 0 || lk < 0) {
            lk = k;
            continue;
        }
        s = (a[lk].dp>>1) + (a[lk].dp&1); 
        e = (a[k].dp>>1) + (a[k].dp&1); 
        if(s >= e) {
            lk = k;
            continue;
        }
        if(old_dp > 1) to += (e -s);
        sk = 0, mmq = 0, mmq_i = (uint64_t)-1;
        for (m = a[k].s; m < (int64_t)a_n; m++) {
            if(vis[m].dp == 0) continue;
            if(ref[m].ts >= e) break;
            oo = OVL(s, e, ref[m].ts, ref[m].te);
            if(!oo) continue;
            if(oo != (e-s)) fprintf(stderr, "ERROR-debug5\n");
            mq = (ref[a[m].s].te - ref[a[m].s].ts)*ref[a[m].s].mq;
            if(mmq_i == (uint64_t)-1 || mmq < mq) {
                mmq_i = m; mmq = mq;
            }
            sk++;
            if(sk >= old_dp) break;
        }
        if (sk < old_dp){
            for (m = (int64_t)(a[k].s)-1; m >= 0; m--) {
                if(vis[m].dp == 0) continue;
                oo = OVL(s, e, ref[m].ts, ref[m].te);
                if(!oo) continue;
                if(oo != (e-s)) fprintf(stderr, "ERROR-debug5\n");
                mq = (ref[a[m].s].te - ref[a[m].s].ts)*ref[a[m].s].mq;
                if(mmq_i == (uint64_t)-1 || mmq < mq) {
                    mmq_i = m; mmq = mq;
                }
                sk++;
                if(sk >= old_dp) break;
            }
        }
        
        if(sk != old_dp) {
            fprintf(stderr, "\nERROR-debug2, sk: %ld, old_dp: %ld, dp: %ld, s: %lu, e: %lu, k: %ld, lk: %ld\n", 
                            sk, old_dp, dp, s, e, k, lk);
            // fprintf(stderr, "ref-i: %lu\n", a[k].s);
            // sk = 0, mmq = 0, mmq_i = (uint64_t)-1;
            // for (m = a[k].s; m < (int64_t)a_n; m++) {
            //     if(vis[m].dp == 0) continue;
            //     fprintf(stderr, "+ts: %lu, +te: %lu\n", ref[m].ts, ref[m].te);
            //     if(ref[m].ts >= e) break;
            //     oo = OVL(s, e, ref[m].ts, ref[m].te);                
            //     sk++;
            //     if(sk >= old_dp) break;
            // }
            // if (sk < old_dp){
            //     for (m = (int64_t)(a[k].s)-1; m >= 0; m--) {
            //         if(vis[m].dp == 0) continue;
            //         fprintf(stderr, "-ts: %lu, -te: %lu\n", ref[m].ts, ref[m].te);
            //         oo = OVL(s, e, ref[m].ts, ref[m].te);
            //         sk++;
            //         if(sk >= old_dp) break;
            //     }
            // }
        }
        if(res->n > sn && res->a[res->n-1].tn == mmq_i && res->a[res->n-1].te == s) {
            res->a[res->n-1].te = e;
        }
        else {
            kv_pushp(interval_t, *res, &p);
            p->tn = mmq_i; p->ts = s; p->te = e;
        }   
        lk = k;
        // if (old_dp < dp && dp >= 1) { ///ts 
        //     start = a[k].dp>>1; sk = k;
        // } else if (old_dp > dp && old_dp >= 1) {///te
        //     for (m = sk, mmq = 0, mmq_i = (uint64_t)-1; m <= k; m++) {
        //         mq = (ref[a[m].s].te - ref[a[m].s].ts)*ref[a[m].s].mq;
        //         if(mmq_i == (uint64_t)-1 || mmq < mq) {
        //             mmq_i = m; mmq = mq;
        //         }
        //     }
        //     s = start; e = (a[k].dp>>1)+1;
        //     if(res->n > sn && res->a[res->n-1].tn == a[mmq_i].s && res->a[res->n-1].te == s) {
        //         res->a[res->n-1].te = e;
        //     }
        //     else {
        //         kv_pushp(interval_t, *res, &p);
        //         p->tn = a[mmq_i].s; p->ts = s; p->te = e;
        //     }   
        // }
    }

    for (k = sn; k < res->n; k++) {
        m = res->a[k].tn; s = res->a[k].ts; e = res->a[k].te;
        res->a[k] = ref[m]; 
        if(!ref[m].rev) {
            res->a[k].qs = (((double)(s - ref[m].ts))/((double)(ref[m].te - ref[m].ts)))*(ref[m].qe - ref[m].qs) + ref[m].qs;
            res->a[k].qe = (((double)(e - ref[m].ts))/((double)(ref[m].te - ref[m].ts)))*(ref[m].qe - ref[m].qs) + ref[m].qs;
        } else {
            res->a[k].qs = ref[m].qe - (((double)(e - ref[m].ts))/((double)(ref[m].te - ref[m].ts)))*(ref[m].qe - ref[m].qs);
            res->a[k].qe = ref[m].qe - (((double)(s - ref[m].ts))/((double)(ref[m].te - ref[m].ts)))*(ref[m].qe - ref[m].qs);
        }
        
        res->a[k].ts = s; res->a[k].te = e;
    }

    return to;
}

uint64_t get_subflag(uint64_t sp, uint64_t ep, uint64_t tn, int64_t idx, ist_t *qE, uint64_t *qidx, double *ff)
{
    ff[0] = ff[1] = ff[2] = 0;
    interval_t *a = qE->a + (qidx[tn]>>32);
    uint64_t a_n = (uint32_t)qidx[tn], oo;
    int64_t i;
    if(idx < 0 || idx >= (int64_t)a_n) idx = 0;
    for (i = 0; i < (int64_t)a_n; i++) {
        if(a[i].ts >= ep) break;
        oo = OVL(sp, ep, a[i].ts, a[i].te);
        ff[a[i].mq] += (oo/(a[i].te - a[i].ts))*a[i].tl;
    }
    for (i = idx-1; i >= 0; i--) {
        if(sp >= a[i].te) break;
        oo = OVL(sp, ep, a[i].ts, a[i].te);
        ff[a[i].mq] += (oo/(a[i].te - a[i].ts))*a[i].tl;
    }
    return i;
}

uint64_t get_type_occ(uint64_t sp, uint64_t ep, interval_t *a, uint64_t a_n, int64_t idx, ist_t *qE, uint64_t *qidx, double *rf)
{
    int64_t i, ii = 0;
    uint64_t ss, ee, qs, qe;
    double ff[3]; rf[0] = rf[1] = rf[2] = 0;
    if(idx < 0 || idx >= (int64_t)a_n) idx = 0;
    for (i = idx; i < (int64_t)a_n; i++) {
        if(a[i].ts >= ep) break;
        ss = MAX(sp, a[i].ts); ee = MIN(ep, a[i].te);
        if(ee <= ss) continue;
        if(!a[i].rev) {
            qs = (((double)(ss - a[i].ts))/((double)(a[i].te - a[i].ts)))*(a[i].qe - a[i].qs) + a[i].qs;
            qe = (((double)(ee - a[i].ts))/((double)(a[i].te - a[i].ts)))*(a[i].qe - a[i].qs) + a[i].qs;
        } else {
            qs = a[i].qe - (((double)(ee - a[i].ts))/((double)(a[i].te - a[i].ts)))*(a[i].qe - a[i].qs);
            qe = a[i].qe - (((double)(ss - a[i].ts))/((double)(a[i].te - a[i].ts)))*(a[i].qe - a[i].qs);
        }
        ii = get_subflag(qs, qe, a[i].qn, ii, qE, qidx, ff);
        rf[0] += ff[0]; rf[1] += ff[1]; rf[2] += ff[2];
    }
    for (i = idx-1; i >= 0; i--) {
        if(sp >= a[i].te) break;
        ss = MAX(sp, a[i].ts); ee = MIN(ep, a[i].te);
        if(ee <= ss) continue;
        if(!a[i].rev) {
            qs = (((double)(ss - a[i].ts))/((double)(a[i].te - a[i].ts)))*(a[i].qe - a[i].qs) + a[i].qs;
            qe = (((double)(ee - a[i].ts))/((double)(a[i].te - a[i].ts)))*(a[i].qe - a[i].qs) + a[i].qs;
        } else {
            qs = a[i].qe - (((double)(ee - a[i].ts))/((double)(a[i].te - a[i].ts)))*(a[i].qe - a[i].qs);
            qe = a[i].qe - (((double)(ss - a[i].ts))/((double)(a[i].te - a[i].ts)))*(a[i].qe - a[i].qs);
        }
        ii = get_subflag(qs, qe, a[i].qn, ii, qE, qidx, ff);
        rf[0] += ff[0]; rf[1] += ff[1]; rf[2] += ff[2];
    }

    return i;
}

void adjust_tid(q_idx_t *tI, ist_t *tbin)
{
    nidx_t *h = (nidx_t*)tI->h_name;
    khint_t k;
    const char *tname = NULL;
    uint32_t tid, tnl, i, kk;
    uint64_t *idx = NULL, t_occ = 0; MALLOC(idx, tI->n.n); memset(idx, -1, tI->n.n*8);
    for (k = 0; k < kh_end(h); ++k) {
        if (kh_exist(h, k)) {
            tname = kh_key(h, k);
            tid = kh_val(h, k);   
            tnl = strlen(tname);
            kk = (uint32_t)-1;
            if(tnl > 3 && tname[0] == 'c' && tname[1] == 'h' && tname[2] == 'r') {
                for (i = 3; i < tnl; i++) {
                    if(!isdigit(tname[i])) break;
                }
                if(i==tnl) kk = atoll(tname+3);
            }
            idx[t_occ] = kk; idx[t_occ] <<= 32; idx[t_occ] |= tid;
            t_occ++;
        }
    }
    ks_introsort_uint64_t(t_occ, idx);


    uint64_t *oidx = NULL, z, m; CALLOC(oidx, tI->n.n);
    for (z = 1, m = 0; z <= tbin->n; ++z) {
        if (z == tbin->n || tbin->a[z].tn != tbin->a[m].tn) {
            if(z < tbin->n && tbin->a[m].tn >= tbin->a[z].tn) break;
            oidx[tbin->a[m].tn] = (m<<32)|(z-m);
            m = z;
        }
    }
    if(z <= tbin->n) fprintf(stderr, "ERROR-sort-1, tbin->n: %u, z: %lu\n", (uint32_t)tbin->n, z);
    interval_t *a = NULL; MALLOC(a, tbin->n);
    for (i = m = 0; i < t_occ; i++) {
        memcpy(a + m, tbin->a + (oidx[(uint32_t)idx[i]]>>32), 
                                                        ((uint32_t)(oidx[(uint32_t)idx[i]]))*sizeof(interval_t));
        m += (uint32_t)(oidx[(uint32_t)idx[i]]);
    }
    if(m != tbin->n) fprintf(stderr, "ERROR-sort-2\n");

    free(tbin->a); tbin->m = tbin->n; tbin->a = a;
    free(idx); free(oidx);
}

void print_hamming(q_idx_t *qI, q_idx_t *tI, ist_t *tbin, const cov_opt_t *opt)
{
    // fprintf(stderr, "\n[M::%s]\n", __func__);
    fprintf(stderr, "chr,beg,end,mid,hap1,hap2,occ\n");
    uint64_t i, k, m, t, *qidx, s = opt->hhStep;
    double rf[3];
    ist_t qE; memset(&qE, 0, sizeof(qE));
    if(!load_yak_hh(opt, opt->yak_hh, qI, &qE, 1)) return;
    CALLOC(qidx, qI->n.n);
    radix_sort_ist_tn(qE.a, qE.a+qE.n);
    for (k = 1, m = 0; k <= qE.n; ++k) {
        if (k == qE.n || qE.a[k].tn != qE.a[m].tn) {
            radix_sort_ist_ts(qE.a+m, qE.a+k);
            qidx[qE.a[m].tn] = (m<<32)|(k-m);
            m = k;
        }
    }

    adjust_tid(tI, tbin);
    for (k = 1, m = 0; k <= tbin->n; ++k) {
        if (k == tbin->n || tbin->a[k].tn != tbin->a[m].tn) {
            ///one chr
            for (t = 1, i = 0; t < tbin->a[m].tl; t += s) {
                i = get_type_occ(t, MIN(t+s, tbin->a[m].tl), tbin->a+m, k-m, i, &qE, qidx, rf);
                rf[0] = rf[1] + rf[2];
                if(rf[0] == 0) continue;
                if(rf[0] > 0) rf[1] /= rf[0], rf[2] /= rf[0];
                rf[1] *= 100; rf[2] *=100;
                // fprintf(stderr, "[M::%s]\t%s\t%lu\t%lu\t%f\t%f\n", __func__, tI->n.a[tbin->a[m].tn].name, 
                // t, MIN(t+s, tbin->a[m].tl), rf[1], rf[2]);
                fprintf(stderr, "%s,%lu,%lu,%lu,%f,%f,%f\n", tI->n.a[tbin->a[m].tn].name, 
                t, MIN(t+s, tbin->a[m].tl), (t + MIN(t+s, tbin->a[m].tl))>>1, rf[1], rf[2], rf[0]);
            }
           /**
            for (t = 1, i = 0; t < tbin->a[m].tl; t += 500) {
                i = get_type_occ(t, MIN(t+s, tbin->a[m].tl), tbin->a+m, k-m, i, &qE, qidx, rf);
                rf[0] = rf[1] + rf[2];
                if(rf[0] == 0) continue;
                if(rf[0] > 0) rf[1] /= rf[0], rf[2] /= rf[0];
                rf[1] *= 100; rf[2] *=100;
                // fprintf(stderr, "[M::%s]\t%s\t%lu\t%lu\t%f\t%f\n", __func__, tI->n.a[tbin->a[m].tn].name, 
                // t, MIN(t+s, tbin->a[m].tl), rf[1], rf[2]);
                fprintf(stderr, "%s,%lu,%lu,%lu,%f,%f\n", tI->n.a[tbin->a[m].tn].name, 
                t, MIN(t+s, tbin->a[m].tl), (t + MIN(t+s, tbin->a[m].tl))>>1, rf[1], -rf[2]);
            }
            **/
            // for (i = m, no = 0; i < k; i++) {
            //     fprintf(stderr, "%s\t%lu\t%lu\t%lu\t%c\t%s\t%lu\t%lu\t%lu\n", 
            //     tI->n.a[tl.a[i].tn].name, tl.a[i].tl, tl.a[i].ts, tl.a[i].te, rr[tl.a[i].rev], 
            //     qI->n.a[tl.a[i].qn].name, tl.a[i].ql, tl.a[i].qs, tl.a[i].qe);
            //     no += tl.a[i].te - tl.a[i].ts;
            // }
            // i = tl.a[m].tn;
            // fprintf(stderr, "[M::%s]\t%s\t%lu(tl)\t%lu(ol-%f)\t%lu(nol-%f)\n", __func__, 
            // tI->n.a[i].name, lchr[i], ochr[i], (((double)ochr[i])/((double)lchr[i]))*100, 
            // lchr[i]-no, (((double)lchr[i]-no)/((double)lchr[i]))*100);
            m = k;
        }
    }

    kv_destroy(qE);
}

void chr_bin_main(const cov_opt_t *opt)
{
    q_idx_t *qI = q_idx_t_init();
    q_idx_t *tI = q_idx_t_init();
    ist_t l; memset(&l, 0, sizeof(l));///all overlaps
    ist_t tl; memset(&tl, 0, sizeof(tl));
    uint64_t i, k, m, oo, *ochr = NULL, *lchr = NULL/**, no**/;
    kvec_t(ma_sub_t) b_idx; kv_init(b_idx);
    kvec_t(ma_sub_t) b_srt; kv_init(b_srt);
    kvec_t(ma_sub_t) b_sc; kv_init(b_sc);
    ma_sub_t *p = NULL;
    // char rr[2] = {'+', '-'};
    if(!load_paf(opt, opt->in, qI, tI, &l, 1, 0, 0, 0)) goto des_cov_main;///load all overlaps
    CALLOC(ochr, tI->n.n); CALLOC(lchr, tI->n.n);
    if(opt->p_only){
        for (k = m = 0; k <= l.n; ++k){
            if(!l.a[k].is_p) continue; 
            if(m!=k) l.a[m] = l.a[k];
            m++;
        }
        l.n = m;
    }
    
    radix_sort_ist_tn(l.a, l.a+l.n);
    for (k = 1, m = 0; k <= l.n; ++k) {
        if (k == l.n || l.a[k].tn != l.a[m].tn) {
            radix_sort_ist_ts(l.a+m, l.a+k);
            kv_resize(ma_sub_t, b_idx, k-m); b_idx.n = 0;
            for (i = m; i < k; i++) {
                kv_pushp(ma_sub_t, b_idx, &p);
                p->dp = 0; p->s = p->e = (uint64_t)-1;
            }

            kv_resize(ma_sub_t, b_srt, (k-m)<<1); b_srt.n = 0;
            for (i = m; i < k; i++){
                kv_pushp(ma_sub_t, b_srt, &p);
                p->dp = l.a[i].ts<<1; p->s = i-m; p->e = 0;
                kv_pushp(ma_sub_t, b_srt, &p);
                p->dp = (l.a[i].te-1)<<1|1; p->s = i-m; p->e = 0;
            }
            radix_sort_ma(b_srt.a, b_srt.a + b_srt.n);
            for (i = 0; i < b_srt.n; i++) {
                if (b_srt.a[i].dp&1) b_idx.a[b_srt.a[i].s].e = i;
                else b_idx.a[b_srt.a[i].s].s = i;
            }
            ///debug
            for (i = 0; i < b_idx.n; i++) {
                if(b_idx.a[i].s == (uint64_t)-1 || b_idx.a[i].e == (uint64_t)-1) {
                    fprintf(stderr, "ERROR-debug0\n");
                }
            }

            kv_resize(ma_sub_t, b_sc, k-m); b_sc.n = 0;
            for (i = m; i < k; i++) {
                kv_pushp(ma_sub_t, b_sc, &p);
                p->dp = l.a[i].mq * (l.a[i].te - l.a[i].ts);
                p->s = l.a[i].te - l.a[i].ts; p->e = i-m; 
            }
            qsort(b_sc.a, b_sc.n, sizeof(ma_sub_t), cmp_mq_t);
            for (i = 0; i < b_sc.n; i++) {
                oo = get_eovlp_len(b_srt.a, b_srt.n, b_sc.a[i].e, b_idx.a);
                if(oo < 10 || oo <= b_sc.a[i].s*opt->sec_rate) {
                    set_eovlp(b_srt.a, b_srt.n, b_sc.a[i].e, b_idx.a);
                }
            }

            ochr[l.a[m].tn] = extract_ist_t(b_srt.a, b_srt.n, b_idx.a, l.a+m, &tl);
            lchr[l.a[m].tn] = l.a[m].tl;
            m = k;
        }
    }
    
    
    // for (k = 1, m = 0; k <= tl.n; ++k) {
    //     if (k == tl.n || tl.a[k].tn != tl.a[m].tn) {
    //         for (i = m, no = 0; i < k; i++) {
    //             fprintf(stderr, "%s\t%lu\t%lu\t%lu\t%c\t%s\t%lu\t%lu\t%lu\n", 
    //             tI->n.a[tl.a[i].tn].name, tl.a[i].tl, tl.a[i].ts, tl.a[i].te, rr[tl.a[i].rev], 
    //             qI->n.a[tl.a[i].qn].name, tl.a[i].ql, tl.a[i].qs, tl.a[i].qe);
    //             no += tl.a[i].te - tl.a[i].ts;
    //         }
    //         i = tl.a[m].tn;
    //         fprintf(stderr, "[M::%s]\t%s\t%lu(tl)\t%lu(ol-%f)\t%lu(nol-%f)\n", __func__, 
    //         tI->n.a[i].name, lchr[i], ochr[i], (((double)ochr[i])/((double)lchr[i]))*100, 
    //         lchr[i]-no, (((double)lchr[i]-no)/((double)lchr[i]))*100);
    //         m = k;
    //     }
    // }

    print_hamming(qI, tI, &tl, opt);
    
    // fprintf(stderr, "\n[M::%s] t_c_tl: %lu, aligned_tl: %lu\n", __func__, t_c_tl, aligned_tl);

    des_cov_main:
    q_idx_t_des(&qI);
    q_idx_t_des(&tI);
    kv_destroy(l);
    kv_destroy(tl);
    kv_destroy(b_idx);
    kv_destroy(b_srt);
    kv_destroy(b_sc);
}

void replace_Ns(const cov_opt_t *opt, char **fs, uint32_t fn)
{
    uint32_t i, k;
    int ret, r, tr = opt->N_rate*128;
    srand(time(NULL));
    for (i = 0; i < fn; i++){
        gzFile fp; kseq_t *ks = NULL;
        if ((fp = gzopen(fs[i], "r")) == 0) {
            fprintf(stderr, "[M::%s::] ==> Cannot open %s\n", __func__, fs[i]);
            return;
        }
        ks = kseq_init(fp);
        while ((ret = kseq_read(ks)) >= 0) {
            fprintf(stderr, ">%s\n", ks->name.s);
            if(tr > 0) {
                for (k = 0; k < ks->seq.l; k++) {
                    r = rand(); 
                    if((r&127) < tr) ks->seq.s[k] = 'N';
                }
            }
            fprintf(stderr, "%s\n", ks->seq.s);
        }
        kseq_destroy(ks);
        gzclose(fp);
    }

}

void plot_merge_Svs(const cov_opt_t *opt, char *fn)
{
    rc_t kt; memset(&kt, 0, sizeof(kt));
    uint64_t i, k, m, l, z, zl, vn, sn, *sn_idx = NULL, *ta, tn, occ;
    uint8_t *vn_vec = NULL;
    char *pch = NULL, *a = NULL;
    kvec_t(uint64_t) res;
    kvec_t(uint64_t) sn_srt;
    if(!rc_t_init(&kt, fn)) {
        rc_t_des(&kt);
        fprintf(stderr, "[M::%s::] ==> Cannot open paf %s\n", __func__, fn);
        return;
    }
    kv_init(res); kv_init(sn_srt); vn = sn = 0;
    while (getLine(&kt))
    {
        if(kt.a[0] == '#') continue;
        i = 0; 
        pch = strtok(kt.a, "\t");
        while (pch != NULL)
        {
            if(i == 7) {
                l = strlen(pch);
                // fprintf(stdout, "i=7->%s\n", pch);
                for (k = 0, m = 0; k <= l; ++k) {
                    if (k == l || pch[k] == ';') {
                        // fprintf(stdout, "m:%lu, k:%lu, %.*s\n", m, k, (int)(k + 1 - m), pch+m);
                        if(k - m > 9 && memcmp(pch+m, "SUPP_VEC=", 9) == 0) {
                            // fprintf(stdout, "%.*s\n", (int)(k + 1 - m), pch+m);
                            a = pch + m + 9;
                            // fprintf(stdout, "%.*s\n", (int)(k - m - 9), a);
                            for (z = 0, zl = k - m - 9; z < zl; z++) {
                                if(a[z] != '0' && a[z] != '1') break;
                                if(a[z] == '1') kv_push(uint64_t, res, (z<<32)|vn);
                                if(sn < z + 1) sn = z + 1;
                            }
                            break;
                        }
                        m = k+1;
                    }
                }
            }

            pch = strtok (NULL, "\t");
            i++;
        }
        vn++;
        if(vn > 0 && vn%100000 == 0) fprintf(stdout, "[M::%s::] ==> # SVs: %lu\n", __func__, vn);
    }

    fprintf(stdout, "[M::%s::] ==> # SVs: %lu; # samples: %lu\n", __func__, vn, sn);

    if(sn > 0) {
        CALLOC(sn_idx, sn); MALLOC(vn_vec, vn);
        ks_introsort_uint64_t(res.n, res.a);
        for (k = 1, m = 0; k <= res.n; ++k) {
            if (k == res.n || (res.a[k]>>32) != (res.a[m]>>32)) {
                sn_idx[res.a[m]>>32] = (m<<32) | (k-m);
                kv_push(uint64_t, sn_srt, (((uint64_t)((uint32_t)(-1) - ((uint32_t)(k-m))))<<32)|(res.a[m]>>32));
                m = k;
            }
        }
        ks_introsort_uint64_t(sn_srt.n, sn_srt.a);

        fprintf(stderr, "Sample,SVs,type\n");
        // fprintf(stderr, "[M::%s::] ==> pan-SVs\n", __func__);
        memset(vn_vec, 0, vn);
        for (i = occ = 0; i < sn_srt.n; i++) {
            if(i > 0 && ((uint32_t)(sn_idx[(uint32_t)sn_srt.a[i]])) > ((uint32_t)(sn_idx[(uint32_t)sn_srt.a[i-1]]))) {
                fprintf(stderr, "ERROR-sort\n");
            }
            ta = res.a + (sn_idx[(uint32_t)sn_srt.a[i]]>>32);
            tn = (uint32_t)(sn_idx[(uint32_t)sn_srt.a[i]]);
            for (k = 0; k < tn; k++) {
                if(vn_vec[(uint32_t)ta[k]] == 0) occ++;
                vn_vec[(uint32_t)ta[k]] = 1;
            }
            if(i > 0) fprintf(stderr, "%lu,%lu,Pan\n", i+1, occ);
        }

        // fprintf(stderr, "[M::%s::] ==> core-SVs\n", __func__);
        memset(vn_vec, 0, vn);
        for (i = 0; i < sn_srt.n; i++) {
            ta = res.a + (sn_idx[(uint32_t)sn_srt.a[i]]>>32);
            tn = (uint32_t)(sn_idx[(uint32_t)sn_srt.a[i]]);
            for (k = occ = 0; k < tn; k++) {
                vn_vec[(uint32_t)ta[k]]++;
                if(vn_vec[(uint32_t)ta[k]] == i+1) occ++;
            }
            if(i > 0) fprintf(stderr, "%lu,%lu,Core\n", i+1, occ);
        }
    }
    rc_t_des(&kt); kv_destroy(res); kv_destroy(sn_srt); free(sn_idx); free(vn_vec);
}

void label_rid_by_ovlp(ma_ug_t *ug, asg_t *rg, ist_t *paf, uint8_t *rf, uint8_t ff)
{
    uint64_t i, k, us, ue, rs, re;
    ma_utg_t *u; 
    interval_t *p = NULL;
    for (i = 0; i < paf->n; i++) {
        p = &(paf->a[i]);
        us = p->ts; ue = p->te; u = &(ug->u.a[p->tn]);
        for (k = rs = 0; k < u->n; k++) {
            re = rs + rg->seq[u->a[k]>>33].len;
            if(OVL(rs, re, us, ue)) rf[u->a[k]>>33] |= ff;
            if(rs >= ue) break;
            rs += (uint32_t)u->a[k];
        }
    }
}

uint64_t travese_nodes(asg_t *g, uint32_t src, uint64_t *feas, uint8_t *vis, kv_uint64_t_warp *res, 
kv_uint32_t_warp *vec)
{
    uint64_t k, v, o = 0, nv;
    asg_arc_t *av = NULL;
    if(res) res->n = 0; vec->n = 0;
    if(feas[src] == 0 || vis[src]) return 0;
    kv_push(uint32_t, *vec, src);

    while (vec->n) {
        v = vec->a[vec->n-1];
        vec->n--;
        if(feas[v] == 0 || vis[v]) continue;
        if(res) kv_push(uint64_t, *res, v);
        o += feas[v]; vis[v] = 1;

        av = asg_arc_a(g, v<<1);
        nv = asg_arc_n(g, v<<1);
        for (k = 0; k < nv; k++) {
            if(av[k].del) continue;
            if(feas[av[k].v>>1] == 0) continue;
            if(vis[av[k].v>>1]) continue;
            kv_push(uint32_t, *vec, av[k].v>>1);
        }

        av = asg_arc_a(g, (v<<1)+1);
        nv = asg_arc_n(g, (v<<1)+1);
        for (k = 0; k < nv; k++) {
            if(av[k].del) continue;
            if(feas[av[k].v>>1] == 0) continue;
            if(vis[av[k].v>>1]) continue;
            kv_push(uint32_t, *vec, av[k].v>>1);
        }
    }

    return o;
}

void print_label_paf(q_idx_t *uI, ma_ug_t *ug, asg_t *rg, ist_t *paf, uint8_t *rf, uint8_t ff, double rr)
{
    kvec_t(uint64_t) ol; kv_init(ol);
    kvec_t(uint64_t) out_ol; kv_init(out_ol);
    uint64_t i, k, m, h, us, ue, rs, re, o, a, num_f = 0, num_uf = 0;
    ma_utg_t *u; 
    interval_t *p = NULL;
    uint64_t *vis; CALLOC(vis, ug->g->n_seq);
    uint64_t *idx; CALLOC(idx, ug->g->n_seq);
    radix_sort_ist_tn(paf->a, paf->a+paf->n);

    for (h = 1, m = 0; h <= paf->n; ++h) {
        if (h == paf->n || paf->a[h].tn != paf->a[m].tn) {
            for (i = m; i < h; i++) {
                p = &(paf->a[i]);
                us = p->ts; ue = p->te; u = &(ug->u.a[p->tn]);
                for (k = rs = o = a = 0; k < u->n; k++) {
                    re = rs + rg->seq[u->a[k]>>33].len;
                    if(OVL(rs, re, us, ue)) {
                        a++;
                        if(rf[u->a[k]>>33]==ff) o++;
                    }
                    if(rs >= ue) break;
                    rs += (uint32_t)u->a[k];
                }

                // fprintf(stderr, "us: %lu, ue: %lu, a: %lu, o: %lu\n", us, ue, a, o);

                p->qs = a; p->qe = o;
                if(o <= a*rr) {
                    num_uf += ue -us;
                    // o = ue - us; o <<= 32; o += i;
                    // kv_push(uint64_t, ol, o);
                    vis[p->tn] += ue - us;
                    if(vis[p->tn] > ug->g->seq[p->tn].len) vis[p->tn] = ug->g->seq[p->tn].len;
                }
                else {
                    num_f += ue -us;
                }
            }
            if(vis[paf->a[m].tn] > 0) {
                o = vis[paf->a[m].tn]; o <<= 32; o += paf->a[m].tn;
                kv_push(uint64_t, ol, o);
            }
            idx[paf->a[m].tn] = (m<<32)|(h-m);
            m = h;
        }
    }

    ks_introsort_uint64_t(ol.n, ol.a);
    for (i = 0; i < (ol.n>>1); i++) {
        o = ol.a[i];
        ol.a[i] = ol.a[ol.n - i - 1];
        ol.a[ol.n - i - 1] = o;
    }
    
    uint8_t *vv; CALLOC(vv, ug->g->n_seq);
    kv_uint64_t_warp res; kv_init(res);
    kv_uint32_t_warp vec; kv_init(vec);
    char rrr[2] = {'+', '-'};
    for (i = 0; i < ol.n; i++) {
        if(vv[(uint32_t)(ol.a[i])] || vis[(uint32_t)(ol.a[i])] == 0) continue;
        o = travese_nodes(ug->g, (uint32_t)(ol.a[i]), vis, vv, NULL, &vec);
        if(o) {
            o <<= 32; o |= (uint32_t)(ol.a[i]);
            kv_push(uint64_t, out_ol, o);
        }
    }


    memset(vv, 0, ug->g->n_seq);
    ks_introsort_uint64_t(out_ol.n, out_ol.a);
    for (i = 0; i < (out_ol.n>>1); i++) {
        o = out_ol.a[i];
        out_ol.a[i] = out_ol.a[out_ol.n - i - 1];
        out_ol.a[out_ol.n - i - 1] = o;
    }

    for (i = 0; i < out_ol.n; i++) {
        o = travese_nodes(ug->g, (uint32_t)(out_ol.a[i]), vis, vv, &res, &vec);
        if(o) {
            fprintf(stderr, "\nS\t%s\to:%lu\n", uI->n.a[(uint32_t)(out_ol.a[i])].name, o);
            for (k = 0; k < res.n; k++) {
                for (m = 0; m < (uint32_t)(idx[res.a[k]]); m++) {
                    p = &(paf->a[(idx[res.a[k]]>>32) + m]);
                    fprintf(stderr, "A\t%s\t%lu\t%lu\t%lu\t%c\t%lu/%lu\n", 
                    uI->n.a[p->tn].name, p->tl, p->ts, p->te, rrr[p->rev], p->qe, p->qs);
                }
            }
        }
    }

    kv_destroy(ol); 
    free(vis); free(idx); free(vv); 
    kv_destroy(res); kv_destroy(vec); kv_destroy(out_ol);
}


void detect_gene_main(const cov_opt_t *opt)
{
    q_idx_t *qI = q_idx_t_init(); 
    q_idx_t *rI = q_idx_t_init(); 
    q_idx_t *uI = q_idx_t_init();
    q_idx_t *cI = q_idx_t_init();
    ist_t u_ov; memset(&u_ov, 0, sizeof(u_ov));
    ist_t c_ov; memset(&c_ov, 0, sizeof(c_ov)); 
    asg_t *rg = asg_init();
    ma_ug_t *u_ug = NULL, *c_ug = NULL;
    uint8_t *rf = NULL;

    if(!(u_ug = load_gfa(opt->ugfas, uI, rI, rg))) goto des_detect_gene_main;
    if(!(c_ug = load_gfa(opt->cgfas, cI, rI, rg))) goto des_detect_gene_main;
    asg_cleanup(rg); asg_symm(rg);

    fprintf(stderr, "[M::%s->u_ug] # nodes: %u, # edges: %u\n", __func__, u_ug->g->n_seq, u_ug->g->n_arc);
    fprintf(stderr, "[M::%s->c_ug] # nodes: %u, # edges: %u\n", __func__, c_ug->g->n_seq, c_ug->g->n_arc);
    fprintf(stderr, "[M::%s->rg] # nodes: %u, # edges: %u\n", __func__, rg->n_seq, rg->n_arc);
    

    if(!load_mul_paf(opt, opt->upafs, qI, uI, &u_ov, 0)) goto des_detect_gene_main;///load all overlaps
    if(!load_mul_paf(opt, opt->cpafs, qI, cI, &c_ov, 0)) goto des_detect_gene_main;///load all overlaps
    

    CALLOC(rf, rg->n_seq);
    label_rid_by_ovlp(c_ug, rg, &c_ov, rf, 1);
    label_rid_by_ovlp(u_ug, rg, &u_ov, rf, 2);
    print_label_paf(uI, u_ug, rg, &u_ov, rf, 1+2, 0.25);

    des_detect_gene_main:
    q_idx_t_des(&qI);
    q_idx_t_des(&rI);
    q_idx_t_des(&uI);
    q_idx_t_des(&cI);
    kv_destroy(u_ov);
    kv_destroy(c_ov);
    asg_destroy(rg);
    ma_ug_destroy(u_ug);
    ma_ug_destroy(c_ug);    
}

void print_sorted_ovlps(q_idx_t *qI, q_idx_t *tI, kv_uint64_t_warp *qI_hh, ist_t *l)
{
    uint64_t m;
    for (m = 0; m < l->n; m++) {
        if(qI_hh->a[l->a[m].qn] == (uint64_t)-1) fprintf(stderr, "ERROR+++%s\n", qI->n.a[l->a[m].qn].name);
        fprintf(stderr, "%s\t%lu\t%lu\t%u\t%lu\t%u\t%s\t%lu\t%lu\t%lu\n", qI->n.a[l->a[m].qn].name,  
        l->a[m].qs, l->a[m].qe, l->a[m].tf, qI_hh->a[l->a[m].qn]>>32, (uint32_t)qI_hh->a[l->a[m].qn], 
        tI->n.a[l->a[m].tn].name, l->a[m].ts, l->a[m].te, l->a[m].tl);        
    }
}

void plot_cen(char *fn, q_idx_t *tI, uint32_t *alph_idx, uint32_t n_chr, double shift_ratio)
{
    if(!fn) return;
    uint32_t i, cid, y0, sk; double y1, y2, x1, x2; char *pch = NULL;
    rc_t kt; memset(&kt, 0, sizeof(kt));
    if(!rc_t_init(&kt, fn)) {
        rc_t_des(&kt);
        fprintf(stderr, "[M::%s::] ==> Cannot open cen bed file %s\n", __func__, fn);
        return;
    }

    while (getLine(&kt)) {
        i = 0; y1 = y2 = x1 = x2 = sk = 0;
        pch = strtok(kt.a, "\t");
        while (pch != NULL) {
            if(i == 0) {
                cid = nidx_query(tI, pch);
                if(cid == ((uint32_t)-1)) {
                    sk = 1; break;
                }
                // if(!(cid != (uint32_t)-1)) {
                //     fprintf(stderr, "[M::%s::] pch::%s\n", __func__, pch);
                // }
                // assert(cid != (uint32_t)-1);
                cid = alph_idx[cid];
                assert(cid != (uint32_t)-1);
                y0 = n_chr - cid;
                y1 = (double)y0 - 0.4;
                y2 = (double)y0 + 0.4;
            }

            if(i == 1) {
                x1 = (double)(atoll(pch)) + shift_ratio;
            }

            if(i == 2) {
                x2 = (double)(atoll(pch)) + shift_ratio;
            }

            pch = strtok(NULL, "\t");
            i++;
        }
        if(sk) continue;
        fprintf(stderr, "set obj rect from %f, %f to %f, %f fc rgb \"#c0c0c0\" fs solid noborder\n", x1, y1, x2, y2);
    }

    rc_t_des(&kt);
}

void print_plot_dot(q_idx_t *qI, q_idx_t *tI, kv_uint64_t_warp *qI_hh, ist_t *l, char *cen, kv_warp_ma_sub_t *N_sites, double thres, double shift_ratio, uint32_t output_lengend, uint32_t set_gaps)
{
    uint64_t k, m, i, pi, c, z, z_s, z_e, min_e, max_s, q_pos, q_subl, t_subl; 
    uint32_t *alph_idx = NULL, n_chr/** = tI->n.n**/ = 0;
    double h, xran, yran, max_len = 0, shift;
    MALLOC(alph_idx, tI->n.n); memset(alph_idx, -1, sizeof(*(alph_idx))*tI->n.n);
    for (k = 1, m = 0; k <= l->n; ++k) {
        if (k == l->n || l->a[k].tn != l->a[m].tn) {
            max_len = MAX(max_len, (double)l->a[m].tl);
            alph_idx[l->a[m].tn] = n_chr;
            m = k; n_chr++;
        }
    }
    h = (((double)n_chr+1)*0.03)*2;
    yran = (double)n_chr + 1.2;
    xran = max_len*(1+shift_ratio);
    shift = max_len*shift_ratio;
    fprintf(stderr, "set t po eps co so \"Helvetica,16\"\n");
    fprintf(stderr, "set out \"test.eps\"\n");
    fprintf(stderr, "set size 2, %f\n", h);
    fprintf(stderr, "set multiplot\n\n");

    fprintf(stderr, "set origin 0, 0\n");
    fprintf(stderr, "set size 2, %f\n", h);
    fprintf(stderr, "set xran [0:%f]\n", xran);
    fprintf(stderr, "set yran [0:%f]\n", yran);
    fprintf(stderr, "set border 0\n");
    fprintf(stderr, "unset xtics\n");
    fprintf(stderr, "unset ytics\n");
    fprintf(stderr, "set bmargin 0\n");
    fprintf(stderr, "set tmargin 0\n");

    plot_cen(cen, tI, alph_idx, n_chr, shift);
    
    for (k = 1, m = 0, i = 0; k <= l->n; ++k) {
        if (k == l->n || l->a[k].tn != l->a[m].tn) {
            fprintf(stderr, "set label \"%s\" at 0.0, first %lu\n", tI->n.a[l->a[m].tn].name, (uint64_t)n_chr-i);
            m = k; i++;
        }
    }

    double x1, x2, y1, y2, r; uint32_t y0, nn, np, nm;
    const char *c0 = "\"#999999\"";
    // const char *cm = "\"#e41a1c\"";
    // const char *cp = "\"#377eb8\"";
    const char *pp;
    for (k = 1, m = 0, c = 0; k <= l->n; ++k) {
        if (k == l->n || l->a[k].tn != l->a[m].tn) {
            for (i = m; i < k; i++) {
                x1 = l->a[i].ts + shift;
	            x2 = l->a[i].te + shift;
                y0 = n_chr - c;
                y1 = l->a[i].tf == 1? (double)y0 + 0.00 : (double)y0 - 0.00;
                y2 = l->a[i].tf == 1? (double)y0 + 0.25 : (double)y0 - 0.25;
                if(qI_hh->a[l->a[i].qn] == (uint64_t)-1) {
                    fprintf(stderr, "ERROR\n");
                    exit(1);
                }
                pp = c0;
                np = qI_hh->a[l->a[i].qn]>>32; nm = (uint32_t)(qI_hh->a[l->a[i].qn]); 
                nn = np + nm;
                if(nn > 0) {
                    r = (double)np/(double)nn;
                    fprintf(stderr, "set obj rect from %f, %f to %f, %f fc rgb \"#%.2x%.2x%.2x\"\n", 
                    x1, y1, x2, y2, ((uint32_t)(255 * (1.0 - r))), 80, ((uint32_t)(255 * r)));
                } else {
                    fprintf(stderr, "set obj rect from %f, %f to %f, %f fc rgb %s\n", x1, y1, x2, y2, pp);
                }
                
            }
            m = k; c++;
        }
    }

    if(set_gaps){
        interval_t *x, *y;
        for (k = 1, m = 0, c = 0; k <= l->n; ++k) {
            if (k == l->n || l->a[k].tn != l->a[m].tn) {
                for (i = m, pi = (uint64_t)-1; i < k; i++) {
                    if(l->a[i].tf != 1) continue;
                    if(pi == (uint64_t)-1) {
                        pi = i;
                        continue;
                    }
                    x = &(l->a[i]); y = &(l->a[pi]); pi = i;
                    if(x->qn == y->qn) continue;
                    ///10000 is the offset
                    // if((y->qe + 10000 < y->ql) && (y->qe + y->ql*0.1 < y->ql)) continue;
                    // if((x->qs > 10000) && (x->qs > x->ql*0.1)) continue;
                    if(((y->qe + 10000 < y->ql) && (y->qe + y->ql*0.1 < y->ql)) && 
                                            ((x->qs > 10000) && (x->qs > x->ql*0.1))) continue;

                    if(y->te < x->ts) {
                        x1 = y->te + shift; x2 = x->ts + shift;
                    } else if(y->te > x->ts) {
                        x1 = x->ts + shift; x2 = y->te + shift;
                    } else {
                        x1 = y->te + shift; x2 = x1 + 32;
                    }

                    y0 = n_chr - c;

                    // y1 = l->a[i].tf == 1? (double)y0 + 0.25 : (double)y0 - 0.25;
                    // y2 = l->a[i].tf == 1? (double)y0 + 0.25 + 0.125 : (double)y0 - 0.25 - 0.125;
                    // fprintf(stderr, "set obj rect from %f, %f to %f, %f fc rgb \"#f09e2b\"\n", x1, y1, x2, y2);  
                    y1 = l->a[i].tf == 1? (double)y0 + 0.00 : (double)y0 - 0.00;      
                    y2 = l->a[i].tf == 1? (double)y0 + 0.25 : (double)y0 - 0.25;
                    fprintf(stderr, "set label \"\" at %f,%f point pointtype 3 pointsize 0.5\n", (x1+x2)/2, (y1+y2)/2);
                }


                for (i = m, pi = (uint64_t)-1; i < k; i++) {
                    if(l->a[i].tf == 1) continue;
                    if(pi == (uint64_t)-1) {
                        pi = i;
                        continue;
                    }
                    x = &(l->a[i]); y = &(l->a[pi]); pi = i;
                    if(x->qn == y->qn) continue;
                    ///10000 is the offset
                    // if((y->qe + 10000 < y->ql) && (y->qe + y->ql*0.1 < y->ql)) continue;
                    // if((x->qs > 10000) && (x->qs > x->ql*0.1)) continue;
                    if(((y->qe + 10000 < y->ql) && (y->qe + y->ql*0.1 < y->ql)) &&
                                            ((x->qs > 10000) && (x->qs > x->ql*0.1))) continue;

                    if(y->te < x->ts) {
                        x1 = y->te + shift; x2 = x->ts + shift;
                    } else if(y->te > x->ts) {
                        x1 = x->ts + shift; x2 = y->te + shift;
                    } else {
                        x1 = y->te + shift; x2 = x1 + 32;
                    }

                    y0 = n_chr - c;

                    // y1 = l->a[i].tf == 1? (double)y0 + 0.25 : (double)y0 - 0.25;
                    // y2 = l->a[i].tf == 1? (double)y0 + 0.25 + 0.15 : (double)y0 - 0.25 - 0.15;
                    // fprintf(stderr, "set obj rect from %f, %f to %f, %f fc rgb \"#f09e2b\"\n", x1, y1, x2, y2); 
                    y1 = l->a[i].tf == 1? (double)y0 + 0.00 : (double)y0 - 0.00;      
                    y2 = l->a[i].tf == 1? (double)y0 + 0.25 : (double)y0 - 0.25;
                    fprintf(stderr, "set label \"\" at %f,%f point pointtype 3 pointsize 0.5\n", (x1+x2)/2, (y1+y2)/2);       
                }
                
                m = k; c++;
            }
        }
    }

    if(N_sites->n > 0) {
        radix_sort_ma(N_sites->a, N_sites->a + N_sites->n);
        CALLOC(N_sites->idx.a, qI->n.n); N_sites->idx.n = qI->n.n;
        for (k = 1, m = 0; k <= N_sites->n; ++k) {
            if (k == N_sites->n || N_sites->a[k].dp != N_sites->a[m].dp) {
                if(k > m) {
                    N_sites->idx.a[N_sites->a[m].dp] = m; 
                    N_sites->idx.a[N_sites->a[m].dp] <<= 32;
                    N_sites->idx.a[N_sites->a[m].dp] += (k - m);
                    radix_sort_s(N_sites->a + m, N_sites->a + k);
                }                    
                m = k;
            }
        }


        for (k = 1, m = 0, c = 0; k <= l->n; ++k) {
            if (k == l->n || l->a[k].tn != l->a[m].tn) {
                for (i = m; i < k; i++) {
                    z_s = N_sites->idx.a[l->a[i].qn]>>32;
                    z_e = z_s + ((uint32_t)(N_sites->idx.a[l->a[i].qn]));
                    for (z = z_s; z < z_e; z++) {
                        assert(N_sites->a[z].dp == l->a[i].qn);
                        if(N_sites->a[z].s >= l->a[i].qe) break;
                        max_s = MAX(N_sites->a[z].s, l->a[i].qs);
                        min_e = MIN(N_sites->a[z].e, l->a[i].qe);
                        if(min_e <= max_s) continue;
                        q_pos = (max_s + min_e)/2; 
                        assert(q_pos >= l->a[i].qs && q_pos < l->a[i].qe);
                        q_pos -= l->a[i].qs;
                        q_subl = l->a[i].qe - l->a[i].qs; 
                        t_subl = l->a[i].te - l->a[i].ts;
                        if(l->a[i].rev) {
                            if(q_subl >= q_pos + 1) q_pos = q_subl - q_pos - 1;
                            else q_pos = 0;
                        }
                        x1 = l->a[i].ts + shift + ((((double)q_pos) / ((double)q_subl))*((double)t_subl));
                        if(x1 - shift < l->a[i].ts || x1 - shift >= l->a[i].te) {
                            fprintf(stdout, "sbsbsbsbsbsbsb\n");
                        }
                        y0 = n_chr - c;
                        y1 = l->a[i].tf == 1? (double)y0 + 0.33 : (double)y0 - 0.33;
                        fprintf(stderr, "set label \"\" at %f,%f point pointtype 12 pointsize 0.5\n", x1, y1);       
                    }
                }
                m = k; c++;
            }
        }
    }
    

    if(output_lengend) {
        fprintf(stderr, "set label \"maternal\" at screen 1.71,0.79\n");
        fprintf(stderr, "set label \"paternal\" at screen 1.71,0.52\n");
    }
    
    fprintf(stderr, "plot \"<echo '0 0'\" not w d\n");

    if(output_lengend) {
        fprintf(stderr, "set origin 1.3, 0.5\n");
        fprintf(stderr, "set size 0.5, 0.3\n");
        fprintf(stderr, "unset label\n");
        fprintf(stderr, "unset obj\n");
        fprintf(stderr, "set palette defined (0 \"#0040FF\", 1 \"#FF4000\");\n");
        fprintf(stderr, "unset key\n");
        fprintf(stderr, "unset xtics\n");
        fprintf(stderr, "unset ytics\n");
        fprintf(stderr, "unset ztics\n");
        fprintf(stderr, "set colorbox\n");
        fprintf(stderr, "unset cbtics\n");
        fprintf(stderr, "plot \"<echo 0 0 0\" w d lc palette z\n");
    }
}

void collect_N_sites(const char* fn, q_idx_t *tI, uint64_t append_hapid, uint8_t tf, kv_warp_ma_sub_t *N_sites)
{
    gzFile fp; kseq_t *ks = NULL; kvec_t(char) idA; kv_init(idA); 
    char *nid; uint32_t cid, s, e, k; int32_t ret; ma_sub_t *x;
    if ((fp = gzopen(fn, "r")) == 0) {
        fprintf(stderr, "[M::%s::] ==> Cannot open %s\n", __func__, fn);
        kv_destroy(idA);
        return;
    }
    ks = kseq_init(fp);
    while ((ret = kseq_read(ks)) >= 0) {
        nid = NULL;
        if(append_hapid) {
            idA.n = 0; kv_resize(char, idA, ks->name.l+128);
            sprintf(idA.a, "%s_hap%u", ks->name.s, tf);
            nid = idA.a;
        } else {
            nid = ks->name.s;
        }
        cid = nidx_query(tI, nid);
        if(cid == (uint32_t)-1) continue;
        for (k = 0, s = e = (uint32_t)-1; k < ks->seq.l; k++) {
            if(ks->seq.s[k] == 'N') {
                if(s == (uint32_t)-1) {
                    s = k; e = k + 1;
                } else {
                    e = k + 1;
                }                
            } else {
                if(s != (uint32_t)-1) {
                    kv_pushp(ma_sub_t, *N_sites, &x);
                    x->dp = cid; x->s = s; x->e = e;
                }
                s = e = (uint32_t)-1;
            }
        }

        if(s != (uint32_t)-1 && e > s) {
            kv_pushp(ma_sub_t, *N_sites, &x);
            x->dp = cid; x->s = s; x->e = e;
        }
    }
    kseq_destroy(ks);
    gzclose(fp);
    kv_destroy(idA);
}


void contig_bin_main(const cov_opt_t *opt, char **pafs, uint32_t pafn)
{
    q_idx_t *qI = q_idx_t_init();///contig idx
    q_idx_t *tI = q_idx_t_init();///chr idx
    kv_uint64_t_warp qI_hh; kv_init(qI_hh);
    ist_t l; memset(&l, 0, sizeof(l));///all overlaps
    uint64_t i, k, m, ss, tt; 
    kvec_t(uint64_t) n_occ; kv_init(n_occ);
    kv_warp_ma_sub_t N_sites; memset(&N_sites, 0, sizeof(N_sites));

    for (i = 0; i < pafn; i++) {
        m = strlen(pafs[i]); n_occ.n = 0;
        for (k = ss = 0; k < m; k++) {
            if(pafs[i][k] == ':') {
                if(k > ss) {
                    tt = ss; tt <<= 32; tt += k;
                    kv_push(uint64_t, n_occ, tt);
                }
                ss = k + 1;
            }
        }
        if(n_occ.n > 0 && k > ss) {
            tt = ss; tt <<= 32; tt += k;
            kv_push(uint64_t, n_occ, tt);
        }
        if(n_occ.n < 2) {
            fprintf(stderr, "Wrong file name: %s\n", pafs[i]);
            exit(1);
        }

        for (k = 0; k < n_occ.n; k++) {
            assert(pafs[i][(uint32_t)(n_occ.a[k])] == ':' || pafs[i][(uint32_t)(n_occ.a[k])] == '\0');
            pafs[i][(uint32_t)(n_occ.a[k])] = '\0';
        }
        
        // pafs[i][k] = '\0';
        if((!load_paf(opt, pafs[i] + (n_occ.a[0]>>32), qI, tI, &l, i+1, 1, 1, 0)) || (!load_yak_hh_S(opt, pafs[i]+ (n_occ.a[1]>>32), qI, &qI_hh, i+1, 1))) {
            fprintf(stderr, "Wrong file name: %s\n", pafs[i]);///load all overlaps
            exit(1);
        }

        if(n_occ.n > 2) collect_N_sites(pafs[i]+ (n_occ.a[2]>>32), qI, 1, i + 1, &N_sites);
        // pafs[i][k] = ':';
        for (k = 0; k < n_occ.n; k++) {
            assert(pafs[i][(uint32_t)(n_occ.a[k])] == '\0');
            if(((uint32_t)(n_occ.a[k])) != m) {
                pafs[i][(uint32_t)(n_occ.a[k])] = ':';
            }
        }
    }
    if(qI_hh.n != qI->n.n) fprintf(stderr, "ERROR-ss\n");///all contigs should have binning results

    for (k = m = 0; k <= l.n; ++k){
        if(l.a[k].is_p == 0 && l.a[k].is_i == 0) continue; 
        if(m!=k) l.a[m] = l.a[k];
        m++;
    }
    l.n = m;
    //sort overlaps by chr_id:chr_pos
    radix_sort_ist_tn(l.a, l.a+l.n);
    for (k = 1, m = 0; k <= l.n; ++k) {
        if (k == l.n || l.a[k].tn != l.a[m].tn) {
            radix_sort_ist_ts(l.a+m, l.a+k);
            m = k;
        }
    }
    adjust_tid(tI, &l);///sort alphabetically by chr_id
    if(!(opt->plot)) {
        print_sorted_ovlps(qI, tI, &qI_hh, &l);
    } else {
        print_plot_dot(qI, tI, &qI_hh, &l, opt->cen, &N_sites, 0.9, 0.05, 1, 1);
    }
    
    

    q_idx_t_des(&qI);
    q_idx_t_des(&tI);
    kv_destroy(qI_hh);
    kv_destroy(l);
    kv_destroy(n_occ);
    kv_destroy(N_sites);
}


int cmp_hapLen_t(const void * a, const void * b)
{
    if((*(hapLen_t*)a).w == (*(hapLen_t*)b).w) return 0;
    return (*(hapLen_t*)a).w > (*(hapLen_t*)b).w? 1:-1;
}

void cov_bp_check(const cov_opt_t *opt)
{
    q_idx_t *qI = q_idx_t_init();
    q_idx_t *tI = q_idx_t_init();
    ist_t l; memset(&l, 0, sizeof(l));
    uint64_t i, c_tl, k, m, tn, qn;
    kvec_t(uint64_t) b; kv_init(b);
    kvec_t(ma_sub_t) res; kv_init(res);
    kvec_t(ma_sub_t) contigs; kv_init(contigs);
    kvec_t(hapLen_t) contig_cov; kv_init(contig_cov);
    ma_sub_t *p = NULL;
    int64_t dp, old_dp, start = 0, min_dp = opt->covX;
    if(!load_paf(opt, opt->in, qI, tI, &l, 1, 0, 0, 0)) goto des_cov_main;
    contig_cov.n = contig_cov.m = tI->n.n;
    CALLOC(contig_cov.a, contig_cov.n);
    radix_sort_ist_tqn(l.a, l.a+l.n);
    for (k = 1, m = 0; k <= l.n; ++k)
    {
        if (k == l.n || l.a[k].tn != l.a[m].tn || l.a[k].qn != l.a[m].qn) 
        {
            
            tn = l.a[m].tn; qn = l.a[m].qn;
            if(!strcmp(qI->n.a[qn].name, tI->n.a[tn].name))
            {
                m = k;
                continue;
            }

            for (i = m, b.n = 0; i < k; i++)
            {
                if(opt->p_only && !l.a[i].is_p) continue;
                kv_push(uint64_t, b, l.a[i].ts<<1);
                kv_push(uint64_t, b, l.a[i].te<<1|1);
            }
            ks_introsort_uint64_t(b.n, b.a);


            res.n = 0;
            for (i = 0, dp = 0, start = 0; i < b.n; ++i) 
            {
                old_dp = dp;
                ///if a[j] is qe
                if (b.a[i]&1) 
                {
                    --dp;
                }
                else
                {
                    ++dp;
                } 


                if (old_dp < min_dp && dp >= min_dp) ///old_dp < dp, b.a[j] is qs
                { 
                    ///case 2, a[j] is qs
                    start = b.a[i]>>1;
                } 
                else if (old_dp >= min_dp && dp < min_dp) ///old_dp > min_dp, b.a[j] is qe
                {
                    if(OVL((uint64_t)start, b.a[i]>>1, opt->chrS, opt->chrE) == 0) continue;

                    kv_pushp(ma_sub_t, res, &p);
                    p->s = MAX((uint64_t)start, opt->chrS); p->e = MIN(b.a[i]>>1, opt->chrE);
                    p->dp = old_dp;
                }
            }

            p = NULL;
            if(qI->n.a[qn].name[1] == tI->n.a[tn].name[1])
            {
                kv_pushp(ma_sub_t, contigs, &p);
                p->s = (uint64_t)-1; p->e = 0; p->dp = 0;
            }

            for (i = c_tl = 0; i < res.n; ++i)
            {
                c_tl += res.a[i].e - res.a[i].s;

                if(p)
                {
                    p->s = MIN(p->s, res.a[i].s);
                    p->e = MAX(p->e, res.a[i].e);
                    p->dp += (res.a[i].e - res.a[i].s);
                }
                // fprintf(stderr, "i-%lu, [%lu, %lu), dp-%lu\n", i, res.a[i].s, res.a[i].e, res.a[i].dp);
            } 

            if(p)
            {
                p->s <<= 32; p->s += qn;
                p->e <<= 32; p->e += tn;
            }

            if(qI->n.a[qn].name[1] == '1')
            {
                contig_cov.a[tn].hL[0] += c_tl;
            }

            if(qI->n.a[qn].name[1] == '2')
            {
                contig_cov.a[tn].hL[1] += c_tl;
            }

            contig_cov.a[tn].len = l.a[m].tl;

            m = k;
        }
    }

    for (k = 0; k < contig_cov.n; k++)
    {
        if(tI->n.a[k].name[1] == '1')
        {
            contig_cov.a[k].w = contig_cov.a[k].hL[1] - contig_cov.a[k].hL[0];
        }

        if(tI->n.a[k].name[1] == '2')
        {
            contig_cov.a[k].w = contig_cov.a[k].hL[0] - contig_cov.a[k].hL[1];
        }
        contig_cov.a[k].id = k;
    }
    
    // radix_sort_hapLen(contig_cov.a, contig_cov.a + contig_cov.n);
    qsort(contig_cov.a, contig_cov.n, sizeof(hapLen_t), cmp_hapLen_t);

    for (k = 0; k < contig_cov.n; k++)
    {
        fprintf(stderr, "%s: hap1_len(%ld), hap2_len(%ld), w(%ld), len(%ld)\n", 
                        tI->n.a[contig_cov.a[k].id].name, contig_cov.a[k].hL[0], contig_cov.a[k].hL[1], contig_cov.a[k].w,
                        contig_cov.a[k].len);
    }

    radix_sort_ma(contigs.a, contigs.a+contigs.n);
    for (k = 0; k < contigs.n; k++)
    {
        fprintf(stderr, "**********len: %lu, qn: %s, tn: %s, s: %lu, e: %lu\n", 
        contigs.a[k].dp, qI->n.a[(uint32_t)contigs.a[k].s].name, tI->n.a[(uint32_t)contigs.a[k].e].name,
        contigs.a[k].s>>32, contigs.a[k].e>>32);
    }

    des_cov_main:
    q_idx_t_des(&qI);
    q_idx_t_des(&tI);
    kv_destroy(l);
    kv_destroy(b);
    kv_destroy(res);
    kv_destroy(contigs);
    kv_destroy(contig_cov);
}

uint64_t get_cover_length(interval_t *a, uint64_t n, int64_t min_dp, uint8_t p_only, uint64_t m_s, kv_uint64_t_warp *b)
{
    uint64_t i, c_tl;
    int64_t dp, old_dp, start = 0;
    b->n = 0; 
    for (i = 0; i < n; i++)
    {
        if(p_only && !a[i].is_p) continue;
        if(a[i].tf & m_s)
        {
            kv_push(uint64_t, *b, a[i].qs<<1);
            kv_push(uint64_t, *b, a[i].qe<<1|1);
        }
    }
        
    ks_introsort_uint64_t(b->n, b->a);
    c_tl = 0;
    for (i = 0, dp = 0, start = 0; i < b->n; ++i) 
    {
        old_dp = dp;
        ///if a[j] is qe
        if (b->a[i]&1) 
        {
            --dp;
        }
        else
        {
            ++dp;
        } 


        if (old_dp < min_dp && dp >= min_dp) ///old_dp < dp, b.a[j] is qs
        { 
            ///case 2, a[j] is qs
            start = b->a[i]>>1;
        } 
        else if (old_dp >= min_dp && dp < min_dp) ///old_dp > min_dp, b.a[j] is qe
        {
            c_tl += (b->a[i]>>1) - start;
        }
    }

    return c_tl;
}


uint32_t adjust_contigs(const cov_opt_t *opt, const char* fn, q_idx_t *qI, ref_in_t *res)
{
    uint64_t i, k, m, tf = 1;
    q_idx_t *tI = q_idx_t_init();
    kv_uint64_t_warp b; kv_init(b);
    ist_t l; memset(&l, 0, sizeof(l));
    ist_t ub; memset(&ub, 0, sizeof(ub));
    interval_t *p = NULL, *s = NULL;
    ref_in *r = NULL;
    if(!load_paf(opt, fn, qI, tI, &l, tf, 0, 0, 0))
    {
        q_idx_t_des(&tI);
        kv_destroy(l);
        kv_destroy(ub);
        return 0;
    }

    radix_sort_ist_qn(l.a, l.a+l.n);
    for (k = 1, m = 0; k <= l.n; ++k)
    {
        if (k == l.n || l.a[k].qn != l.a[m].qn) //same qn
        {
            radix_sort_ist_tn(l.a + m, l.a + k);  
            m = k;
        }
    }

    for (k = 1, m = 0; k <= l.n; ++k)
    {
        if (k == l.n || l.a[k].qn != l.a[m].qn || l.a[k].tn != l.a[m].tn) //same qn and tn
        {
            kv_pushp(interval_t, ub, &p);
            *p = l.a[m];

            for (i = m, p->err = 0; i < k; i++)
            {
                if(opt->p_only && !l.a[i].is_p) continue;
                p->err += l.a[i].err;
            }
            p->mlen = get_cover_length(l.a+m, k-m, 1, opt->p_only, tf, &b);
            p->err += p->ql - p->mlen;
            p->qn = l.a[m].qn; p->tn = l.a[m].tn;
            m = k;
        }
    }



    res->n = 0;
    for (k = 1, m = 0; k <= ub.n; ++k)
    {
        if (k == ub.n || ub.a[k].qn != ub.a[m].qn) //same qn, determine which chromesome qn from
        {
            // radix_sort_ist_score(ub.a + m, ub.a + k);
            radix_sort_ist_mlen(ub.a + m, ub.a + k);
            p = s = NULL;
            p = &(ub.a[k-1]);
            s = k - m < 2? NULL:&(ub.a[k-2]);
            if(s && s->mlen >= opt->mis_assembLen)
            {
                fprintf(stderr, "**********[%s] len: %lu, primary[err: %ld, base: %lu]", 
                        qI->n.a[p->qn].name, p->ql, p->err, p->mlen);
                if(s) fprintf(stderr, ", secondary[err: %ld, base: %lu]", s->err, s->mlen);
                fprintf(stderr,"\n");
            }
            // if(strcmp(tI->n.a[p->tn].name, "chr21") == 0)
            // {
            //     fprintf(stderr, "######[%s-to-%s] len: %lu, primary[base: %lu]", 
            //             qI->n.a[p->qn].name, tI->n.a[p->tn].name, p->ql, p->mlen);
            //     if(s) fprintf(stderr, ", secondary[base: %lu]", s->mlen);
            //     fprintf(stderr,"\n");
            // }
            
            kv_pushp(ref_in, *res, &r);
            r->s = 0; r->e = p->ql; r->l = p->ql; r->id = p->qn;
            r->chr = nidx_strdup(tI->n.a[p->tn].name);
            r->chr_id = p->tn;
            r->chr_len = p->tl;
            r->pl = p->mlen;
            m = k;
        }
    }
    

    q_idx_t_des(&tI);
    kv_destroy(l);
    kv_destroy(ub);
    kv_destroy(b);
    return 1;
}


void trio_cov_main(const cov_opt_t *opt)
{
    q_idx_t *qI = q_idx_t_init();
    q_idx_t *tI = q_idx_t_init();
    ist_t l; memset(&l, 0, sizeof(l));
    ref_in_t ref; kv_init(ref);
    ref_cov_t rcov; kv_init(rcov);
    ref_cov *rp = NULL;
    uint64_t i, k, m, m_s[2], hf, c_tl[2], rate[2], rate_a;
    kv_uint64_t_warp b; kv_init(b);
    kvec_t(ma_sub_t) res; kv_init(res);
    kvec_t(ma_sub_t) contigs; kv_init(contigs);
    int64_t s[2], min_dp = opt->covX;
    uint64_t totalBase = 0, totalAnchorBase = 0;

    m_s[0] = 1; m_s[1] = 2;
    if(!adjust_contigs(opt, opt->in, qI, &ref)) goto des;
    if(!load_paf(opt, opt->t_in_1, qI, tI, &l, m_s[0], 0, 0, 0)) goto des;
    if(!load_paf(opt, opt->t_in_2, qI, tI, &l, m_s[1], 0, 0, 0)) goto des;

    radix_sort_ist_qn(l.a, l.a+l.n);
    // radix_sort_ist_tn(l.a, l.a+l.n);
    for (k = 1, m = 0; k <= l.n; ++k)
    {
        if (k == l.n || l.a[k].qn != l.a[m].qn) //same qn
        {           
            b.n = 0; 
            for (i = m, s[0] = s[1] = 0; i < k; i++)
            {
                if(opt->p_only && !l.a[i].is_p) continue;
                if(l.a[i].tf & m_s[0]) s[0] += l.a[i].err;
                if(l.a[i].tf & m_s[1]) s[1] += l.a[i].err;
            }
            c_tl[0] = get_cover_length(l.a+m, k-m, min_dp, opt->p_only, m_s[0], &b);
            c_tl[1] = get_cover_length(l.a+m, k-m, min_dp, opt->p_only, m_s[1], &b);

            s[0] += l.a[m].ql - c_tl[0];
            s[1] += l.a[m].ql - c_tl[1];

            s[0] = (l.a[m].ql - s[0])*opt->match_score - s[0] * opt->mis_penalty;
            s[1] = (l.a[m].ql - s[1])*opt->match_score - s[1] * opt->mis_penalty;

            hf = 0;
            if(s[0] >= s[1]) hf |= m_s[0];
            if(s[1] >= s[0]) hf |= m_s[1];

            char* ref_chr = NULL;
            uint64_t ref_chr_id = (uint64_t)-1, ref_chr_len = (uint64_t)-1;
            for (i = 0; i < ref.n; i++)
            {
                if(ref.a[i].id == l.a[m].qn && ref.a[i].pl >= ref.a[i].l*0.75)
                {
                    ref_chr = ref.a[i].chr;
                    ref_chr_id = ref.a[i].chr_id;
                    ref_chr_len = ref.a[i].chr_len;
                }
            }

            kv_pushp(ref_cov, rcov, &rp);
            rp->tlen = ref_chr_len;
            rp->chr = ref_chr;
            rp->chr_id = ref_chr_id;
            rp->cov_a = rp->cov[0] = rp->cov[1]=0;
            if((hf & m_s[0]) && (hf & m_s[1]))
            {
                rp->cov_a += l.a[m].ql;
            }
            else if((hf & m_s[0]))
            {
                rp->cov[0] += l.a[m].ql;
            }
            else if((hf & m_s[1]))
            {
                rp->cov[1] += l.a[m].ql;
            }
            if(ref_chr && strcmp(ref_chr, "chr7") == 0)
            {
                fprintf(stderr, "[%s] len: %lu: s[0]: %ld, c_tl[0]: %lu, s[1]: %ld, c_tl[1]: %lu, %u%u\n", 
                                            qI->n.a[l.a[m].qn].name, l.a[m].ql, s[0], c_tl[0], s[1], c_tl[1], 
                                            !!(hf & m_s[0]), !!(hf & m_s[1]));
            }

            // fprintf(stderr, "[%s] len: %lu, s[0]: %ld, cl[0]: %lu, s[1]: %ld, cl[1]: %lu, ref_chr: %s\n", 
            //                     qI->n.a[l.a[m].qn].name, l.a[m].ql, s[0], c_tl[0], s[1], c_tl[1], ref_chr?ref_chr:"NA");
            m = k;
        }
    }

    // exit(1);

    radix_sort_ref_cov_chr_id(rcov.a, rcov.a+rcov.n);
    // qsort(rcov.a, rcov.n, sizeof(ref_cov), cmp_ref_cov);
    
    for (k = 1, m = 0; k <= rcov.n; ++k)
    {
        if (k == rcov.n || rcov.a[k].chr_id != rcov.a[m].chr_id) //same chr_id
        {      
            uint64_t cov[2], cov_a, occ = 0;
            cov[0] = cov[1] = cov_a = 0;
            for (i = m; i < k; i++)
            {
                cov[0] += rcov.a[i].cov[0];
                cov[1] += rcov.a[i].cov[1];
                cov_a += rcov.a[i].cov_a;
            }
            char* ref_chr = NULL;
            ref_chr = rcov.a[m].chr;

            occ += (cov[0] > 0 ) + (cov[1] > 0) + (cov_a > 0);
            rate[0] = (cov[0]*100)/rcov.a[m].tlen;
            rate[1] = (cov[1]*100)/rcov.a[m].tlen;
            rate_a = (cov_a*100)/rcov.a[m].tlen;

            // if(occ != 1 || (cov[0] + cov[1] + cov_a) <= rcov.a[m].tlen*0.7)
            // if((cov[0] + cov[1] + cov_a) > rcov.a[m].tlen)
            if (MIN(cov[0], cov[1]) >= MAX(cov[0], cov[1])*0.075)
            {
                fprintf(stderr, "[%s-%lu], # hap1 bases: %lu(%lu%%), # hap2 bases: %lu(%lu%%), # amb bases: %lu(%lu%%)\n", 
                ref_chr?ref_chr:"NA", rcov.a[m].tlen, cov[0], rate[0], cov[1], rate[1], cov_a, rate_a);
            }
            
            totalBase += cov[0] + cov[1] + cov_a;
            if(ref_chr) totalAnchorBase += cov[0] + cov[1] + cov_a;
            m = k;
        }
    }
    fprintf(stderr, "totalBase: %lu, totalAnchorBase: %lu\n", totalBase, totalAnchorBase);



    des:
    q_idx_t_des(&qI);
    q_idx_t_des(&tI);
    kv_destroy(l);
    kv_destroy(b);
    kv_destroy(res);
    kv_destroy(contigs);
}

kv_int8_t_warp *get_trio_flag(const char* fn)
{
    kv_int8_t_warp *a; CALLOC(a, 1);
    rc_t kt; memset(&kt, 0, sizeof(kt));
    char *pch = NULL;
    uint64_t rid;
    int8_t t;
    if(!rc_t_init(&kt, fn))
    {
        rc_t_des(&kt);
        fprintf(stderr, "[M::%s::] ==> Cnnot open paf %s\n", __func__, fn);
        return NULL;
    }
    while (getLine(&kt))
    {
        if(kt.n == 0) continue;
        if(kt.a[0] != 'A') continue;
        t = 0; rid = (uint64_t)-1;
        pch = strtok(kt.a, "\t");
        while (pch != NULL)
        {
            if(strlen(pch) >= 6)
            {
                if(!memcmp(pch, "id:i:", 5))
                {
                    rid = atoll(pch+5);
                }
                else if(!memcmp(pch, "HG:A:", 5))
                {
                    t = pch[5] == 'p'? 1 : t;
                    t = pch[5] == 'm'? -1 : t;
                }
            }

            pch = strtok (NULL, "\t");
        }
        if(rid == (uint64_t)-1) continue;
        if(rid >= a->n)
        {
            kv_resize(int8_t, *a, rid+1);
            memset(a->a + a->n, 0, rid + 1 - a->n);
            a->n = rid+1;
        }
        a->a[rid] = t;
    }
    rc_t_des(&kt);
    return a;
}

ref_cov_t *bin_gfa(const cov_opt_t *opt, const char* cn, const char* rn, q_idx_t *qI)
{
    kv_int8_t_warp *flag = get_trio_flag(rn);
    if(!flag) exit(0);
    rc_t kt; memset(&kt, 0, sizeof(kt));
    char *pch = NULL, *qch = NULL;
    uint64_t rid, i, qid, qlen, k, m;
    if(!rc_t_init(&kt, cn))
    {
        rc_t_des(&kt);
        fprintf(stderr, "[M::%s::] ==> Cnnot open paf %s\n", __func__, cn);
        return NULL;
    }

    kvec_t(uint64_t) R; kv_init(R);
    while (getLine(&kt))
    {
        if(kt.n == 0) continue;
        if(kt.a[0] != 'A' && kt.a[0] != 'S') continue;
        rid = qlen = (uint64_t)-1; i = 0; qch = NULL;
        pch = strtok(kt.a, "\t");
        while (pch != NULL)
        {
            if(i == 1) qch = pch;
            if(strlen(pch) >= 6)
            {
                if(!memcmp(pch, "id:i:", 5))
                {
                    rid = atoll(pch+5);
                }
                else if(!memcmp(pch, "LN:i:", 5))
                {
                    qlen = atoll(pch+5);
                }
            }

            pch = strtok (NULL, "\t");
            i++;
        }
        if(qch && qlen != (uint64_t)-1)///S-line
        {
            qid = nidx_append(qI, qch, 0);
            qI->l.a[qid] = qlen;
            continue;
        }
        if(rid == (uint64_t)-1 || !qch) continue;

        qid = nidx_append(qI, qch, 0);
        qid <<= 2;
        if(rid < flag->n) 
        {
            if(flag->a[rid] > 0) qid += 1;
            if(flag->a[rid] < 0) qid += 2;
        }
        kv_push(uint64_t, R, qid);
    }

    // fprintf(stderr, "# contigs: %u\n", (uint32_t)qI->n.n);

    ref_cov_t *rcov; CALLOC(rcov, 1);
    ref_cov *rp = NULL;
    ks_introsort_uint64_t(R.n, R.a);
    for (k = 1, m = 0; k <= R.n; ++k)
    {
        if (k == R.n || (R.a[k]>>2) != (R.a[m]>>2)) //same qid
        {           
            uint64_t s[3] = {0, 0, 0};
            for (i = m; i < k; i++)
            {
                s[R.a[i]&3]++;
            }
            kv_pushp(ref_cov, *rcov, &rp);
            rp->tlen = qI->l.a[R.a[m]>>2];
            rp->chr = NULL;
            rp->chr_id = R.a[m]>>2;
            rp->cov_a = s[0];
            rp->cov[0] = s[1];
            rp->cov[1] = s[2];
            // fprintf(stderr, "%s:%lu, p-%lu, m-%lu, a-%lu\n", qI->n.a[rp->chr_id].name, rp->tlen, rp->cov[0], rp->cov[1], rp->cov_a);

            m = k;
        }
    }



    kv_destroy(R);
    kv_destroy(*flag); 
    free(flag);
    rc_t_des(&kt);
    return rcov;
}

void trio_cov_main_gfa(const cov_opt_t *opt)
{
    q_idx_t *qI = q_idx_t_init();
    ref_in_t ref; kv_init(ref);
    ref_cov_t *rcov = NULL, chr_info; kv_init(chr_info);
    ref_cov *rp = NULL;
    uint64_t i, k, m, rate[2], rate_a;
    uint64_t totalBase = 0, totalAnchorBase = 0;
    int status[2] = {1, -1};
    if(!adjust_contigs(opt, opt->in, qI, &ref)) goto des;
    rcov = bin_gfa(opt, opt->q_gfa, opt->r_gfa, qI);
    if(!rcov) goto des;

    for (k = 0; k < rcov->n; k++) ///scan all contigs
    {
        char* ref_chr = NULL;
        uint64_t ref_chr_id = (uint64_t)-1, ref_chr_len = (uint64_t)-1;
        for (i = 0; i < ref.n; i++)
        {
            if(ref.a[i].id == rcov->a[k].chr_id && ref.a[i].pl >= ref.a[i].l*0.75)
            {
                ref_chr = ref.a[i].chr;
                ref_chr_id = ref.a[i].chr_id;
                ref_chr_len = ref.a[i].chr_len;
            }
        }


        kv_pushp(ref_cov, chr_info, &rp);
        *rp = rcov->a[k];
        rp->tlen = ref_chr_len;
        rp->chr = ref_chr;
        rp->chr_id = ref_chr_id;
        rp->qid = rcov->a[k].chr_id;
        // if(ref_chr && strcmp(ref_chr, "chr1") == 0)
        // {
        //     fprintf(stderr, "[%s] len: %lu: s[0]-%lu, s[1]-%lu\n", qI->n.a[rcov->a[k].chr_id].name, 
        //     rcov->a[k].tlen, rp->cov[0], rp->cov[1]);
        // }
    }
    

    // exit(1);

    radix_sort_ref_cov_chr_id(chr_info.a, chr_info.a+chr_info.n);
    // qsort(rcov.a, rcov.n, sizeof(ref_cov), cmp_ref_cov);
    
    for (k = 1, m = 0; k <= chr_info.n; ++k)
    {
        if (k == chr_info.n || chr_info.a[k].chr_id != chr_info.a[m].chr_id) //same chr_id
        {      
            uint64_t cov[2], cov_a, occ = 0;
            cov[0] = cov[1] = cov_a = 0;
            for (i = m; i < k; i++)
            {
                cov[0] += chr_info.a[i].cov[0];
                cov[1] += chr_info.a[i].cov[1];
                cov_a += chr_info.a[i].cov_a;
            }
            char* ref_chr = NULL;
            ref_chr = chr_info.a[m].chr;

            occ += (cov[0] > 0 ) + (cov[1] > 0) + (cov_a > 0);
            rate[0] = (cov[0]*100)/(cov[0] + cov[1] + cov_a);
            rate[1] = (cov[1]*100)/(cov[0] + cov[1] + cov_a);
            rate_a = (cov_a*100)/(cov[0] + cov[1] + cov_a);

            // if(occ != 1 || (cov[0] + cov[1] + cov_a) <= rcov.a[m].tlen*0.7)
            // if((cov[0] + cov[1] + cov_a) > rcov.a[m].tlen)
            if(MIN(cov[0], cov[1]) >= MAX(cov[0], cov[1])*0.05 && ref_chr)
            {
                radix_sort_ref_cov_chr_cov(chr_info.a+m, chr_info.a+k);
                fprintf(stderr, "\n[%s-%lu], # hap1 reads: %lu(%lu%%), # hap2 reads: %lu(%lu%%), # amb reads: %lu(%lu%%)\n", 
                ref_chr?ref_chr:"NA", chr_info.a[m].tlen, cov[0], rate[0], cov[1], rate[1], cov_a, rate_a);
                for (i = m; i < k; i++)
                {
                    if(chr_info.a[i].cov[0] == 0 && chr_info.a[i].cov[1] == 0) continue;
                    fprintf(stderr, "[%s-status:%d] len: %lu: s[0]-%lu, s[1]-%lu\n", 
                    qI->n.a[chr_info.a[i].qid].name, status[chr_info.a[i].cov[1] > chr_info.a[i].cov[0]], 
                    qI->l.a[chr_info.a[i].qid], chr_info.a[i].cov[0], chr_info.a[i].cov[1]);
                }
            }
            
            totalBase += cov[0] + cov[1] + cov_a;
            if(ref_chr) totalAnchorBase += cov[0] + cov[1] + cov_a;
            m = k;
        }
    }
    fprintf(stderr, "# Reads: %lu, # AnchorReads: %lu\n", totalBase, totalAnchorBase);

    des:
    q_idx_t_des(&qI);
    kv_destroy(chr_info);
}

typedef struct {
    int64_t *p, *f;
    int64_t m;
} chain_dp;

//ai is the suffix of the aj
int64_t cal_chain_sc(interval_t *ai, interval_t *aj, int64_t fj)
{
    int64_t sc = INT64_MIN, ql, tl;
    if(ai->qs >= aj->qs && ai->qe >= aj->qe && 
            ai->ts >= aj->ts && ai->te >= aj->te) {
        ql = ai->qe-aj->qe; tl = ai->te-aj->te;
        if(ql <= tl) {
            sc = ql - (((double)ql)/((double)(ai->qe-ai->qs)))*ai->err;
        } else {
            sc = tl - (((double)tl)/((double)(ai->te-ai->ts)))*ai->err;
        }
        // sc = MIN(ai->qe-aj->qe, ai->te-aj->te);
        return sc + fj;
    }
    return sc;
}

void gen_aln_chain(interval_t *a, int64_t a_n, int64_t aoff, chain_dp *dp, ist_t *res)
{
    if(a_n <= 0) return;
    // fprintf(stderr, "a_n::%ld\n", a_n);///load all overlaps
    int64_t k, l, i, j, k0, l0, m;
    for (l = 0, k = 1; k <= a_n; k++) {
        if(k == a_n || a[k].rev != a[l].rev) {
            radix_sort_ist_qs(a + l, a + k);
            for (l0 = l, k0 = l + 1; k0 <= k; k0++) {
                if(k0 == k || a[k0].qs != a[l0].qs) {
                    radix_sort_ist_qe(a + l0, a + k0);
                    l0 = k0;
                }
            }
            l = k;
        }
    }
    if(a_n > dp->m) {
        REALLOC(dp->p, a_n); REALLOC(dp->f, a_n);
    }
    int64_t *p = dp->p, *f = dp->f, max_f, st, max_j, sc, msc, msc_i;
    msc = msc_i = -1;
    for (i = 0, st = 0; i < a_n; ++i) {
        max_f = MIN((a[i].qe-a[i].qs), (a[i].te-a[i].ts)); max_f -= a[i].err;
        max_j = -1;
        if(a[i].rev != a[st].rev) st = i;
        for (j = i-1; j >= st; --j) {
            sc = cal_chain_sc(&(a[i]), &(a[j]), f[j]);
            if(sc == INT64_MIN) continue;
            if (sc > max_f) {
                max_f = sc, max_j = j;
            }
        }
        f[i] = max_f; p[i] = max_j;
        if(f[i] > msc) {
            msc = f[i]; msc_i = i;
        }
    }
    if(msc_i < 0) return;
    interval_t z; z = a[msc_i];
    
    i = msc_i; m = 0;
    while (i >= 0) {
        if(a[i].qs < z.qs) z.qs = a[i].qs;
        if(a[i].qe > z.qe) z.qe = a[i].qe;
        if(a[i].ts < z.ts) z.ts = a[i].ts;
        if(a[i].te > z.te) z.te = a[i].te;
        f[m++] = i; i = p[i];
    }
    if(m <= 0) return;

    for (i = m-1, k = 0; i >= 0; i--) a[k++] = a[f[i]];
    assert(k == m);
    if(z.rev) {
        int64_t ts, te;
        ts = z.tl - z.te; 
        te = z.tl - z.ts;
        z.ts = ts; z.te = te; 
        for (i = 0; i < m; i++) {
            ts = a[i].tl - a[i].te; 
            te = a[i].tl - a[i].ts;
            a[i].ts = ts; a[i].te = te; 
        }
    }
    z.mq = aoff; z.mlen = m;
    kv_push(interval_t, *res, z);
}

void srt_ist_t_chain(interval_t *a, int64_t a_n, chain_dp *dp, ist_t *res)
{
    int64_t k, l, i, ts, te; res->n = 0;
    for (k = 0; k < a_n; k++) a[k].qn = ((a[k].qn<<32)|(a[k].tn<<1)|((uint64_t)a[k].rev));
    // for (k = 0; k < a_n; k++) a[k].qn = ((a[k].qn<<32)|(a[k].tn));
    radix_sort_ist_qn(a, a+a_n);
    for (l = 0, k = 1; k <= a_n; k++) {
        if(k == a_n || (a[k].qn>>1) != (a[l].qn>>1)) {
            for (i = l; i < k; i++) {
                a[i].qn >>= 32;
                if(a[i].rev) {
                    ts = a[i].tl - a[i].te; 
                    te = a[i].tl - a[i].ts;
                    a[i].ts = ts; a[i].te = te; 
                }
            }
            gen_aln_chain(a + l, k - l, l, dp, res);
            l = k;
        }
    }
    
}

void t2t_main(const cov_opt_t *opt, char *paf)
{
    q_idx_t *qI = q_idx_t_init();///contig idx
    q_idx_t *tI = q_idx_t_init();///chr idx
    ist_t l; memset(&l, 0, sizeof(l));///all overlaps
    ist_t res; memset(&res, 0, sizeof(res));//chains
    chain_dp dp; memset(&dp, 0, sizeof(dp));

    if(!load_paf(opt, paf, qI, tI, &l, 0, 0, 0, 0)) {
        fprintf(stderr, "Wrong file name: %s\n", paf);///load all overlaps
        exit(1);
    }
    fprintf(stderr, "File name::%s, l.n::%u\n", paf, (uint32_t)l.n);///load all overlaps
    srt_ist_t_chain(l.a, l.n, &dp, &res);
    uint64_t k;
    if(opt->t2tlen == INT64_MAX) {
        for (k = 0; k < res.n; k++) {
            if(((res.a[k].qe-res.a[k].qs) >= (res.a[k].ql*0.95)) && 
                    ((res.a[k].te-res.a[k].ts) >= (res.a[k].tl*0.95))) 
            {
                fprintf(stderr, "%s\t%lu\t%lu\t%lu\t%c\t%s\t%lu\t%lu\t%lu\n", 
                qI->n.a[res.a[k].qn].name, res.a[k].ql, res.a[k].qs, res.a[k].qe,  "+-"[res.a[k].rev],
                tI->n.a[res.a[k].tn].name, res.a[k].tl, res.a[k].ts, res.a[k].te);      
            }
        }
    } else {
        for (k = 0; k < res.n; k++) {
            if((((int64_t)(res.a[k].qe-res.a[k].qs)) >= opt->t2tlen) && 
                                    (((int64_t)(res.a[k].te-res.a[k].ts)) >= opt->t2tlen)) {
                fprintf(stderr, "%s\t%lu\t%lu\t%lu\t%c\t%s\t%lu\t%lu\t%lu\n", 
                qI->n.a[res.a[k].qn].name, res.a[k].ql, res.a[k].qs, res.a[k].qe,  "+-"[res.a[k].rev],
                tI->n.a[res.a[k].tn].name, res.a[k].tl, res.a[k].ts, res.a[k].te);  
                // uint64_t i, a_n; interval_t *a;     
                // a = l.a + res.a[k].mq; a_n = res.a[k].mlen;
                // for (i = 0; i < a_n; i++) {
                //     fprintf(stderr, ">%s\t%lu\t%lu\t%lu\t%c\t%s\t%lu\t%lu\t%lu\n", 
                //     qI->n.a[a[i].qn].name, a[i].ql, a[i].qs, a[i].qe,  "+-"[a[i].rev],
                //     tI->n.a[a[i].tn].name, a[i].tl, a[i].ts, a[i].te);      
                // }
            }
        }
    }
    

    q_idx_t_des(&qI);
    q_idx_t_des(&tI);
    kv_destroy(l);
    kv_destroy(res);
    free(dp.f); free(dp.p); 
}


uint32_t load_bed(const cov_opt_t *opt, const char* fn, q_idx_t *qI, vbed_t *qbed, vec_assemb_t *qseq)
{
    rc_t kt; memset(&kt, 0, sizeof(kt));
    uint64_t i, k, qs, qe, type, qid, os, oe, pe;
    char *pch = NULL, *qn = NULL;
    bed_t *t = NULL;
    if(!rc_t_init(&kt, fn)) {
        rc_t_des(&kt);
        fprintf(stderr, "[M::%s::] ==> Cannot open bed %s\n", __func__, fn);
        return 0;
    }

    while (getLine(&kt)) {
        i = 0; qs = qe = (uint64_t)-1; type = 0;
        if((kt.n >= 5) && (!memcmp(kt.a, "track", 5))) continue;
        pch = strtok(kt.a, "\t");
        while (pch != NULL) {
            if(i == 0) qn = pch;
            if(i == 1) qs = atoll(pch);
            if(i == 2) qe = atoll(pch);
            if((i == 3) && (strlen(pch) == 3)) {
                if(!memcmp(pch, "Hap", 3)) type = 0;
                else if(!memcmp(pch, "Err", 3)) type = 1;
                else if(!memcmp(pch, "Unk", 3)) type = 2;
                else if(!memcmp(pch, "Col", 3)) type = 3;
                else if(!memcmp(pch, "Dup", 3)) type = 4;
            }
            pch = strtok (NULL, "\t");
            i++;
        }
        // if(qe - qs < min_l) continue;
        if(qe <= qs) continue;
        qid = nidx_append(qI, qn, 0);
        if (type && qseq && qseq->a[qid].Ns.n) {
            for (k = 0, pe = qs; k < qseq->a[qid].Ns.n; k++) {
                if(qe <= qseq->a[qid].Ns.a[k].s) break;
                if(qs >= qseq->a[qid].Ns.a[k].e) continue;
                os = MAX(qs, qseq->a[qid].Ns.a[k].s);
                oe = MIN(qe, qseq->a[qid].Ns.a[k].e);
                assert(oe > os);
                if(os > pe) {
                    kv_pushp(bed_t, *qbed, &t);
                    t->qn = qid;
                    t->qs = pe; t->qe = os; t->type = type;
                }
                kv_pushp(bed_t, *qbed, &t);
                t->qn = qid;
                t->qs = os; t->qe = oe; t->type = 5;///N
                pe = oe;
            }
            if(pe < qe) {
                kv_pushp(bed_t, *qbed, &t);
                t->qn = qid;
                t->qs = pe; t->qe = qe; t->type = type;
            }            
        } else {
            kv_pushp(bed_t, *qbed, &t);
            t->qn = qid;
            t->qs = qs; t->qe = qe; t->type = type;
        }
    }

    rc_t_des(&kt); 
    fprintf(stderr, "[M::%s::] ==> Loaded %s\n", __func__, fn);
    return 1;
}

void load_assemble(const char* fn, q_idx_t *qI, vec_assemb_t *res)
{
    gzFile fp; kseq_t *ks = NULL; 
    char *nid; uint32_t cid, s, e, k; int32_t ret; ma_sub_t *x; assemb_t *p;
    if ((fp = gzopen(fn, "r")) == 0) {
        fprintf(stderr, "[M::%s::] ==> Cannot open %s\n", __func__, fn);
        return;
    }
    ks = kseq_init(fp);
    while ((ret = kseq_read(ks)) >= 0) {
        nid = ks->name.s;
        cid = nidx_append(qI, nid, 0);
        kv_pushp(assemb_t, *res, &p);
        memset(p, 0, sizeof((*p)));
        p->id = cid; 
        p->seq.n = p->seq.m = ks->seq.l;
        MALLOC(p->seq.a, ks->seq.l); 
        memcpy(p->seq.a, ks->seq.s, ks->seq.l*(sizeof(*(ks->seq.s))));
        
        for (k = 0, s = e = (uint32_t)-1; k < ks->seq.l; k++) {
            if(ks->seq.s[k] == 'N') {
                if(s == (uint32_t)-1) {
                    s = k; e = k + 1;
                } else {
                    e = k + 1;
                }                
            } else {
                if(s != (uint32_t)-1) {
                    kv_pushp(ma_sub_t, p->Ns, &x);
                    x->dp = cid; x->s = s; x->e = e;
                }
                s = e = (uint32_t)-1;
            }
        }

        if(s != (uint32_t)-1 && e > s) {
            kv_pushp(ma_sub_t, p->Ns, &x);
            x->dp = cid; x->s = s; x->e = e;
        }
    }
    kseq_destroy(ks);
    gzclose(fp);
    fprintf(stderr, "[M::%s::] ==> Loaded %s\n", __func__, fn);
}

void debug_paf_cigar(ist_t *l, vec_assemb_t *qseq, vec_assemb_t *tseq, q_idx_t *qI, q_idx_t *tI)
{
    uint64_t k, z, cn; uint32_t *ca; int64_t offset, ql, tl;
    for (k = 0; k < l->n; k++) {
        fprintf(stderr, "+k::%lu\n", k);
        cn = l->a[k].clen;
        ca = l->cg.a + l->a[k].cidx;
        for (z = offset = 0; z < cn; z++) {
            if((ca[z]&3) == IF) offset -= (int64_t)(ca[z]>>2);///q is longer
            if((ca[z]&3) == DF) offset += (int64_t)(ca[z]>>2);///t is longer
            // fprintf(stderr, "len::%ld, c::%u\n", (int64_t)(ca[z]>>2), (ca[z]&3));
        }
        ql = l->a[k].qe - l->a[k].qs;
        tl = l->a[k].te - l->a[k].ts;
        if((ql + offset) != tl) {
            fprintf(stderr, "%s\t%lu\t%lu\t%lu\t%c\t%s\t%lu\t%lu\t%lu\n", 
            qI->n.a[l->a[k].qn].name, l->a[k].ql, l->a[k].qs, l->a[k].qe, "+-"[l->a[k].rev], 
            tI->n.a[l->a[k].tn].name, l->a[k].tl, l->a[k].ts, l->a[k].te);  
            fprintf(stderr, "ql::%ld, tl::%ld, offset::%ld\n", ql, tl, offset);
            exit(1);
        }
        fprintf(stderr, "-k::%lu\n", k);
    }
}

void idx_vbed_t(vbed_t *in)
{
    uint64_t k, l;
    for (k = in->idx.n = 0; k < in->n; k++) {
        if(in->a[k].qn > in->idx.n) in->idx.n = in->a[k].qn;
    }
    in->idx.n++;
    kv_resize(uint64_t, in->idx, in->idx.n);
    // fprintf(stderr, "in->idx.n::%lu, in->idx.m::%lu\n", (uint64_t)in->idx.n, in->idx.m);
    memset(in->idx.a, 0, sizeof((*(in->idx.a)))*in->idx.n);
    
    radix_sort_bed_qn(in->a, in->a+in->n);
    for (k = 1, l = 0; k <= in->n; ++k) {
        if (k == in->n || in->a[k].qn != in->a[l].qn) {
            radix_sort_bed_qs(in->a+l, in->a+k);
            in->idx.a[in->a[l].qn] = (l<<32)|(k-l);
            l = k;
        }
    }
}

uint64_t cal_map_errors(int64_t qs, int64_t qe, interval_t *mm, ist_t *bmap, vbed_t *tbed, q_idx_t *qI, q_idx_t *tI, uint64_t *rts, uint64_t *rte)
{
    uint32_t *ca; int64_t cn, k, cqs, cqe, cts, cte, ts, te, s, e;
    ca = bmap->cg.a + mm->cidx; cn = mm->clen;
    cqs = mm->qs; cts = mm->ts; cqe = cte = ts = te = -1;
    for (k = 0; k < cn; k++) {
        if((ca[k]&3) == IF) {///q is longer
            cqe = cqs + (int64_t)(ca[k]>>2);
            cte = cts;
        } else if((ca[k]&3) == DF) {///t is longer
            cqe = cqs;
            cte = cts + (int64_t)(ca[k]>>2);
        } else {
            cqe = cqs + (int64_t)(ca[k]>>2);
            cte = cts + (int64_t)(ca[k]>>2);
        }
        if(qs >= cqs && qs < cqe) {
            if((ca[k]&3) == IF) {///q is longer
                ts = cts;
            } else if((ca[k]&3) == DF) {///t is longer
                ts = cts;
            } else {
                ts = cts + qs - cqs;
            }
        }
        if(qe >= cqs && qe < cqe) {
            if((ca[k]&3) == IF) {///q is longer
                te = cte;
            } else if((ca[k]&3) == DF) {///t is longer
                te = cte;
            } else {
                te = cts + qe - cqs;
            }
        }
        cqs = cqe; cts = cte;
        if((ts != -1) && (te != -1)) break;
    }
    if(qs == (int64_t)cqe) ts = cte;
    if(qe == (int64_t)cqe) te = cte;
    // if(!(qe == (int64_t)mm->qe)) {
    //     fprintf(stderr, "%s\t%lu\t%lu\t%lu\t%c\t%s\t%lu\t%lu\t%lu\n", 
    //         qI->n.a[mm->qn].name, mm->ql, mm->qs, mm->qe, "+-"[mm->rev], 
    //         tI->n.a[mm->tn].name, mm->tl, mm->ts, mm->te);  
    //     fprintf(stderr, "q::[%ld, %ld), t::[%ld, %ld)\n", qs, qe, ts, te);  
    // }
    // assert(qe == (int64_t)mm->qe);
    // assert(te == (int64_t)mm->te);
    if(!((ts != -1) && (te != -1))) {
        fprintf(stderr, "%s\t%lu\t%lu\t%lu\t%c\t%s\t%lu\t%lu\t%lu\n", 
            qI->n.a[mm->qn].name, mm->ql, mm->qs, mm->qe, "+-"[mm->rev], 
            tI->n.a[mm->tn].name, mm->tl, mm->ts, mm->te);  
        fprintf(stderr, "q::[%ld, %ld), t::[%ld, %ld)\n", qs, qe, ts, te);  
    }
    
    assert((ts != -1) && (te != -1));
    if(!(mm->rev)) {
        s = ts; e = te;
    } else {
        s = (int64_t)(mm->te) - (te-mm->ts); 
        e = (int64_t)(mm->te) - (ts-mm->ts);
    }
    assert(s <= e);
    if(rts) *rts = s;
    if(rte) *rte = e;
    if(s == e) return 0;
    bed_t *a = tbed->a + (tbed->idx.a[mm->tn]>>32); 
    int64_t a_n = (uint32_t)tbed->idx.a[mm->tn], os, oe, tot = 0;
    for (k = 0; k < a_n; k++) {
        assert(a[k].qn ==mm->tn);
        if(a[k].qs >= (uint64_t)e) break;
        if(a[k].qe <= (uint64_t)s) continue;
        if(a[k].type == 1 || a[k].type == 2 || a[k].type == 3 || a[k].type == 4) {
            os = MAX(a[k].qs, (uint64_t)s); oe = MIN(a[k].qe, (uint64_t)e);
            assert(oe > os);
            tot += oe - os;
        }
    }
    return tot;
}

void print_mask_bed(vbed_t *qbed, vbed_t *tbed, vec_assemb_t *qseq, q_idx_t *qI, q_idx_t *tI, ist_t *bmap, uint64_t *bmap_idx, uint64_t min_l, double cutrate, char *out)
{
    const char *hap = "Hap";
    const char *err = "Err";
    const char *unk = "Unk";
    const char *col = "Col";
    const char *dup = "Dup";
    const char *gap = "Gap";
    const char *raw = "Raw";
    const char *msk = "Msk";
    const char *s0, *s1;

    char *gfa_name; MALLOC(gfa_name, strlen(out)+100);
    sprintf(gfa_name, "%s.mask.bed", out);
    FILE* fn = fopen(gfa_name, "w");
    fprintf(stderr, "[M::%s::] ==> Output masked results to %s\n", __func__, gfa_name);

    char *dbg_name; MALLOC(dbg_name, strlen(out)+100);
    sprintf(dbg_name, "%s.dbg.bed", out);
    FILE* dn = fopen(dbg_name, "w");
    fprintf(stderr, "[M::%s::] ==> Output debug results to %s\n", __func__, dbg_name);

    ///fprintf(stderr, "%s\n", out);
    // fprintf(stderr, "%s\n", gfa_name);
    fprintf(fn, "track name=\"%s\" visibility=1 itemRgb=\"On\"\n", out);
    uint64_t k, z, on, os, oe, mmerr, cerr, raw_tot = 0, mask_tot = 0, rts, rte, is_dbg; interval_t *oa;
    uint64_t *qlen; CALLOC(qlen, qseq->n);
    for (k = 0; k < qseq->n; k++) qlen[qseq->a[k].id] = qseq->a[k].seq.n;
    
    for (k = 0; k < qbed->n; k++) {
        if(qlen[qbed->a[k].qn] < min_l) continue;
        s0 = hap; s1 = raw; is_dbg = 0;
        if(qbed->a[k].type == 1 || qbed->a[k].type == 2 || qbed->a[k].type == 3 || qbed->a[k].type == 4) {
            raw_tot += qbed->a[k].qe - qbed->a[k].qs;
            oa = bmap->a + (bmap_idx[qbed->a[k].qn]>>32); 
            on = (uint32_t)bmap_idx[qbed->a[k].qn];
            mmerr = (uint64_t)-1;
            for (z = 0; z < on; z++) {
                assert(oa[z].qn ==qbed->a[k].qn);
                if(oa[z].qs >= qbed->a[k].qe) break;
                if(oa[z].qe <= qbed->a[k].qs) continue;
                os = MAX(oa[z].qs, qbed->a[k].qs);
                oe = MIN(oa[z].qe, qbed->a[k].qe);
                assert(oe > os);
                if((oe - os) >= ((qbed->a[k].qe - qbed->a[k].qs)*0.5)) {
                    cerr = cal_map_errors(os, oe, &(oa[z]), bmap, tbed, qI, tI, NULL, NULL);
                    if(cerr < mmerr) mmerr = cerr;
                }
            }
            if((mmerr == ((uint64_t)-1)) || ((qbed->a[k].qe-qbed->a[k].qs) >= (mmerr*cutrate))) {///still an error region
                if(qbed->a[k].type == 1) s0 = err;
                if(qbed->a[k].type == 2) s0 = unk;
                if(qbed->a[k].type == 3) s0 = col;
                if(qbed->a[k].type == 4) s0 = dup;
                mask_tot += qbed->a[k].qe - qbed->a[k].qs;
                is_dbg = 1;
            } else {///has been masked
                s1 = msk;
            }
        } else if(qbed->a[k].type == 5) {
            s1 = gap;
        }
        fprintf(fn, "%s\t%lu\t%lu\t%s\t%s\t%lu\n", 
            qI->n.a[qbed->a[k].qn].name, qbed->a[k].qs, qbed->a[k].qe, s0, s1, qbed->a[k].qe-qbed->a[k].qs); 
        if(is_dbg) {
            fprintf(dn, "%s\t%lu\t%lu\t%s\t%s\t%lu\n", 
                qI->n.a[qbed->a[k].qn].name, qbed->a[k].qs, qbed->a[k].qe, s0, s1, qbed->a[k].qe-qbed->a[k].qs); 
            oa = bmap->a + (bmap_idx[qbed->a[k].qn]>>32); 
            on = (uint32_t)bmap_idx[qbed->a[k].qn];
            for (z = 0; z < on; z++) {
                assert(oa[z].qn ==qbed->a[k].qn);
                if(oa[z].qs >= qbed->a[k].qe) break;
                if(oa[z].qe <= qbed->a[k].qs) continue;
                os = MAX(oa[z].qs, qbed->a[k].qs);
                oe = MIN(oa[z].qe, qbed->a[k].qe);
                assert(oe > os);
                if((oe - os) >= ((qbed->a[k].qe - qbed->a[k].qs)*0.5)) {
                    cerr = cal_map_errors(os, oe, &(oa[z]), bmap, tbed, qI, tI, &rts, &rte);
                    fprintf(dn, "#%s\t%lu\t%lu\t%lu\tDbg\n", 
                        tI->n.a[oa[z].tn].name, rts, rte, cerr); 
                    // if(cerr < mmerr) mmerr = cerr;
                }
            }
        } 
    }

    fprintf(stderr, "# raw errors::%lu\n", raw_tot);
    fprintf(stderr, "# mask errors::%lu\n", mask_tot);
    free(gfa_name);
    fclose(fn);
    free(dbg_name);
    fclose(dn);
    free(qlen);
}

void flagger_filter_main(const cov_opt_t *opt, char *fin_bed, char *fin_assemb, char *fin_con, 
char *fmask_bed, char *fmask_assemb, char *fmask_paf, int64_t flen, double fcut)
{
    q_idx_t *qI = q_idx_t_init();///fin idx
    q_idx_t *tI = q_idx_t_init();///fmask idx
    ist_t l; memset(&l, 0, sizeof(l));///all overlaps
    uint64_t k, m, *lidx = NULL; 
    vec_assemb_t qseq, tseq; 
    memset(&qseq, 0, sizeof(qseq)); 
    memset(&tseq, 0, sizeof(tseq));
    vbed_t qbed, tbed; 
    memset(&qbed, 0, sizeof(qbed)); 
    memset(&tbed, 0, sizeof(tbed)); 

    if(fin_assemb) load_assemble(fin_assemb, qI, &qseq);
    if(fin_bed) load_bed(opt, fin_bed, qI, &qbed, &qseq);

    if(fmask_assemb) load_assemble(fmask_assemb, tI, &tseq);
    if(fmask_bed) load_bed(opt, fmask_bed, tI, &tbed, &tseq);

    if(fmask_paf) load_paf(opt, fmask_paf, qI, tI, &l, 0, 0, 0, 1);

    idx_vbed_t(&qbed); idx_vbed_t(&tbed);

    
    CALLOC(lidx, qI->n.n);
    radix_sort_ist_qn(l.a, l.a+l.n);
    for (k = 1, m = 0; k <= l.n; ++k) {
        if (k == l.n || l.a[k].qn != l.a[m].qn) {
            radix_sort_ist_qs(l.a+m, l.a+k);
            lidx[l.a[m].qn] = (m<<32)|(k-m);
            m = k;
        }
    }

    uint64_t oulen = strlen(fin_bed);
    for (k = m = 0; k < oulen; k++) if(fin_bed[k] == '/') m = k+1;
    print_mask_bed(&qbed, &tbed, &qseq, qI, tI, &l, lidx, flen, fcut, fin_bed+m);
    // debug_paf_cigar(&l, &qseq, &tseq, qI, tI);

    q_idx_t_des(&qI);
    q_idx_t_des(&tI);
    kv_destroy(l);
}