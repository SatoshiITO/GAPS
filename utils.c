#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>
  
#include <strings.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <zlib.h>
#include <assert.h>
#include <signal.h>
#include <inttypes.h>
#include <unistd.h>

// Suppress deprecation message for cigar_tab, which we initialise
#include "htslib/hts_defs.h"
#undef HTS_DEPRECATED
#define HTS_DEPRECATED(message)
  
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "cram/cram.h"
#include "hts_internal.h"
#include "sam_internal.h"
#include "htslib/hfile.h"
#include "htslib/hts_endian.h"
#include "htslib/hts_expr.h"
#include "header.h"
  
#include "htslib/khash.h"
KHASH_DECLARE(s2i, kh_cstr_t, int64_t)
KHASH_SET_INIT_INT(tag)
  
#ifndef EFTYPE
#define EFTYPE ENOEXEC
#endif
#ifndef EOVERFLOW
#define EOVERFLOW ERANGE
#endif


static
void GAPS_print_cigar( bam1_t *b )
{
     char buff[1024], op;
     int i, len;

     uint32_t *cigar = bam_get_cigar(b);
     
     //fprintf(stderr,"CIGAR: ");
     for ( i = 0; i < b->core.n_cigar; i++) {
	  len = bam_cigar_oplen( cigar[i] );
	  op  = bam_cigar_opchr( cigar[i] );
	  //len = bam_cigar_oplen( 308 ); //19S
	  //op  = bam_cigar_opchr( 308 ); //19S
	  fprintf(stderr,"cigar = %d\n",(int)cigar[i]);
	  fprintf(stderr,"bam_cigar_op(c) = %d\n",bam_cigar_op(cigar[i]));
	  fprintf(stderr,"len, op = %d %c\n",len, op);
	  snprintf( buff, 1024, "%d%c", len, op );
	  //fprintf(stderr,"%s",buff);
     }
     fprintf(stderr,"\n");
     
     //printf("n_cigar = %d\n", b->core.n_cigar);
     

}

static
void GAPS_print_quality( bam1_t *b )
{
     // format query quality
     uint8_t *tmp_q   = bam_get_qual(b);                 // query quality
     char    *qqual;
     if (tmp_q[0] == 0xff) {
	  qqual=(char *)malloc(2);
	  qqual[0]='*';
	  qqual[1]='\0';
     } else {
	  qqual=(char *)malloc(b->core.l_qseq + 1);
	  for (int i = 0; i < b->core.l_qseq; ++i) qqual[i]=tmp_q[i]+33;
	  qqual[b->core.l_qseq]='\0';
     }
     fprintf(stderr, "(print_qual) Quality: %s\n", qqual);
}

/* These function was copied and modified from bam_read1 in sam.c */
static
void bam_cigar2rqlens(int n_cigar, const uint32_t *cigar,
                              hts_pos_t *rlen, hts_pos_t *qlen)
 {
     int k;
     *rlen = *qlen = 0;
     for (k = 0; k < n_cigar; ++k) {
         int type = bam_cigar_type(bam_cigar_op(cigar[k]));
         int len = bam_cigar_oplen(cigar[k]);
         if (type & 1) *qlen += len;
         if (type & 2) *rlen += len;
     }
 }     
int store_read(uint8_t *cdata, bam1_t *b)
{
     bam1_core_t *c = &b->core;
     int32_t block_len, ret, i;
     uint32_t x[8], new_l_data;
     int nBytes;
     
     b->l_data = 0;

     //if ((ret = bgzf_read(fp, &block_len, 4)) != 4) {
     //    if (ret == 0) return -1; // normal end-of-file
     //    else return -2; // truncated
     //}
     //if (fp->is_be) ed_swap_4p(&block_len);
     //if (block_len < 32) return -4;  // block_len includes core data
     memcpy( &block_len, cdata, 4 );
     nBytes = 4;
     
     //if (bgzf_read(fp, x, 32) != 32) return -3;
     //if (fp->is_be) {
     //    for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
     //}
     memcpy( x, cdata +nBytes, 32 );
     nBytes += 32;
     
     c->tid        = x[0];
     c->pos        = (int32_t)x[1];
     c->bin        = x[2]>>16;
     c->qual       = x[2]>>8&0xff;
     c->l_qname    = x[2]&0xff;
     c->l_extranul = (c->l_qname%4 != 0)? (4 - c->l_qname%4) : 0;
     c->flag       = x[3]>>16;
     c->n_cigar    = x[3]&0xffff;
     c->l_qseq     = x[4];
     c->mtid       = x[5];
     c->mpos       = (int32_t)x[6];
     c->isize      = (int32_t)x[7];

     fprintf(stderr, "POS     : %ld\n", c->pos);     // Clear!!
     fprintf(stderr, "MAPQ    : %d\n",  c->qual);    // Clear!!
     fprintf(stderr, "FLAG    : %d\n",  c->flag);    // Clear!!
     fprintf(stderr, "TID     : %d\n",  c->tid);     // Clear!!
     fprintf(stderr, "l_qname : %d\n",  c->l_qname); // Clear!!
     fprintf(stderr, "MPOS    : %ld\n", c->mpos);    // Clear!!
     fprintf(stderr, "ISIZE: %ld\n", c->isize); 

     // Check length of b->data array. I had better modify here!
     new_l_data = block_len - 32 + c->l_extranul;
     if (new_l_data > INT_MAX || c->l_qseq < 0 || c->l_qname < 1) return -4;
     if ( ((uint64_t) c->n_cigar << 2) + c->l_qname + c->l_extranul
         + ( ((uint64_t) c->l_qseq + 1) >> 1 ) + c->l_qseq > (uint64_t) new_l_data )
         return -4;
     //if (realloc_bam_data(b, new_l_data) < 0) return -4;
     b->data = realloc(b->data, new_l_data); // Because (sam_)realloc_bam_data is not public.
     if ( !b->data ) return -4; 
     b->l_data = new_l_data;
     

     //if (bgzf_read(fp, b->data, c->l_qname) != c->l_qname) return -4;
     memcpy( b->data, cdata +nBytes, b->l_data );
     nBytes += b->l_data;

     
     if (b->data[c->l_qname - 1] != '\0') { // Try to fix missing NUL termination
     //    if (fixup_missing_qname_nul(b) < 0) return -4;
	  fprintf(stderr,"Missing NUL termination is detected.\n");
     }
     
     //GAPS_print_cigar(b); CIGAR Clear!
     for (i = 0; i < c->l_extranul; ++i) b->data[c->l_qname+i] = '\0';
     c->l_qname += c->l_extranul;
     //GAPS_print_cigar(b); CIGAR was broken!!!
     
     //if (b->l_data < c->l_qname ||
     //    bgzf_read(fp, b->data + c->l_qname, b->l_data - c->l_qname) != b->l_data - c->l_qname)
     //    return -4;
     //if (fp->is_be) swap_data(c, b->l_data, b->data, 0);
     memcpy( b->data + c->l_qname, cdata +nBytes, b->l_data - c->l_qname);
     nBytes += b->l_data - c->l_qname;
     if (b->l_data < c->l_qname ) return -4;

     //if (bam_tag2cigar(b, 0, 0) < 0) return -4; // return 0 if CIGAR is untouched; 1 if CIGAR is updated with CG

     fprintf(stderr,"Read Name: %s\n",bam_get_qname(b));
     //fprintf(stderr,"Read seq.: %s\n",bam_get_seq(b));
     
     //GAPS_print_cigar(b);
     //GAPS_print_quality(b);
     
     if (c->n_cigar > 0) { // recompute "bin" and check CIGAR-qlen consistency
         hts_pos_t rlen, qlen;
         bam_cigar2rqlens(c->n_cigar, bam_get_cigar(b), &rlen, &qlen);
         if ((b->core.flag & BAM_FUNMAP) || rlen == 0) rlen = 1;
         b->core.bin = hts_reg2bin(b->core.pos, b->core.pos + rlen, 14, 5);
         // Sanity check for broken CIGAR alignments
         if (c->l_qseq > 0 && !(c->flag & BAM_FUNMAP) && qlen != c->l_qseq) {
             hts_log_error("CIGAR and query sequence lengths differ for %s",
                     bam_get_qname(b));
             return -4;
         }
     }
#if DEBUG           
     return 4 + block_len;
#endif
 }
