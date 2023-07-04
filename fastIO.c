#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "hfile_internal.h"
#include "htslib/thread_pool.h"
#include "htslib/hts_endian.h"

#define  READ_DATA_SIZE 0x8000000 //128M

extern int store_read(uint8_t *cdata, bam1_t *b);


static inline int unpackInt16(const uint8_t *buffer)
{
     return buffer[0] | buffer[1] << 8;
}



static int bgzf_uncompress(uint8_t *dst, size_t *dlen,
                            const uint8_t *src, size_t slen,
                            uint32_t expected_crc) {
     z_stream zs = {
         .zalloc = NULL,
         .zfree = NULL,
         .msg = NULL,
         .next_in = (Bytef*)src,
         .avail_in = slen,
         .next_out = (Bytef*)dst,
         .avail_out = *dlen
     };
  
     int ret = inflateInit2(&zs, -15);
     if (ret != Z_OK) {
	  //hts_log_error("Call to inflateInit2 failed: %s", bgzf_zerr(ret, &zs));
         return -1;
     }
     if ((ret = inflate(&zs, Z_FINISH)) != Z_STREAM_END) {
	  //hts_log_error("Inflate operation failed: %s", bgzf_zerr(ret, ret == Z_DATA_ERROR ? &zs : NULL));
         if ((ret = inflateEnd(&zs)) != Z_OK) {
	      //hts_log_warning("Call to inflateEnd failed: %s", bgzf_zerr(ret, NULL));
         }
         return -1;
     }
     if ((ret = inflateEnd(&zs)) != Z_OK) {
	  //hts_log_error("Call to inflateEnd failed: %s", bgzf_zerr(ret, NULL));
         return -1;
     }
     *dlen = *dlen - zs.avail_out;
     fprintf(stderr,"(bgzf_uncomp) dlen = %ld\n", *dlen);
     
     uint32_t crc = crc32(crc32(0L, NULL, 0L), (unsigned char *)dst, *dlen);
     if (crc != expected_crc) {
         hts_log_error("CRC32 checksum mismatch");
	 fprintf(stderr,"crc, expected_crc = %d %d\n", crc, expected_crc);
         //return -2;
     }
  
     return 0;
 }




int main(int argc, char *argv[]) {
    
    if (argc <= 1 || argc >= 3) {
	 fprintf(stderr, "Error: Invalid arguments are specified.\n");
	 fprintf(stderr, "Usage: mkdup in.bam\n");
	 return 0;
    }
    
    /* ***** Startup messages ***** */
    fprintf(stdout, "GPU accelerated Mark duplicates start.\n");
    fprintf(stdout, "Input: %s\n", argv[1]);

   
    /* *** Open BAM *** */
    char      *f   = argv[1];
    htsFile   *in  = hts_open(f,"rb"); //Clear
    //hfile_set_blksize(hFILE *fp, size_t bufsiz);
    hfile_set_blksize(in->fp.hfile, READ_DATA_SIZE);
    //sam_hdr_t *h   = sam_hdr_read(in);
    bam_hdr_t *h   = sam_hdr_read(in);
    htsFile   *out = hts_open("temp.bam","wb");
    
    if ( sam_hdr_write(out, h) < 0 ){
	fprintf(stderr, "Error: Fail to write header to temp.sam file.\n");
	abort();
    }


    // read bam
    bam1_t *b= bam_init1();
#if DEBUG    
    //kstring_t aux={0, 0, NULL};
    //while (sam_read1(in, h, b) >= 0) {          // end file = -1
    //while (bam_read1( in->fp.bgzf, b ) >= 0) {    // Clear

    //	 if ( sam_write1(out, h, b) < 0 ){
    //      fprintf(stderr, "Error: Fail to write temp.sam file.\n");
    //      abort;
    //   }
    //}
#endif

    void *buff, *data;
    //int   in_size=0, rest_size=0;
    int   i,ret;
    
    buff = (void *)malloc( READ_DATA_SIZE           ); // 128MB
    data = (void *)malloc( READ_DATA_SIZE + 0x200000); // 300MB
    
    if ( !buff && !data ){
	 fprintf(stderr, "Error: Fail to allocate buff memory.\n");
          abort();
     }
    for ( i = 0; i < READ_DATA_SIZE; i++ ) ((uint8_t*)buff)[i] = 0;
    
#define BLOCK_HEADER_LENGTH 18
#define BLOCK_FOOTER_LENGTH 8
    uint8_t header[BLOCK_HEADER_LENGTH];
    
    while(1){
	 //ret = bgzf_read( in->fp.bgzf, buff, READ_DATA_SIZE );
	 //if ( ret <= 0 ) break;
	 //ret = bgzf_write( out->fp.bgzf, buff, ret );

	 ret = bgzf_raw_read( in->fp.bgzf, buff, READ_DATA_SIZE );
	 if ( ret <= 0 ) break;
	 uint8_t *cdata  = (uint8_t *)buff;
	 uint8_t *uncomp = (uint8_t *)data;
	 uint16_t BSIZE;
	 uint32_t CRC;
	 
	 memcpy( header, (uint8_t *)buff, BLOCK_HEADER_LENGTH );
	 fprintf(stderr, "ID1   = %d\n", header[0] );
	 fprintf(stderr, "ID2   = %d\n", header[1] );
	 fprintf(stderr, "CM    = %d\n", header[2] );
	 fprintf(stderr, "FLG   = %d\n", header[3] );
	 fprintf(stderr, "OS    = %d\n", header[9] );
	 //fprintf(stderr, "XLEN  = %d\n", (uint16_t)header[10] );
	 fprintf(stderr, "XLEN  = %d\n", unpackInt16( (uint8_t*)&header[10] ));
	 fprintf(stderr, "SI1   = %d\n", header[12] );
	 fprintf(stderr, "SI2   = %d\n", header[13] );
	 //fprintf(stderr, "SLEN  = %d\n", (uint16_t)header[14] );
	 fprintf(stderr, "SLEN  = %d\n", unpackInt16( (uint8_t*)&header[14] ));
	 //fprintf(stderr, "h[16],h[17] = %d %d\n", header[16], header[17] );
	 BSIZE = unpackInt16( (uint8_t*)&header[16] ) +1;
	 fprintf(stderr, "BSIZE = %d\n", BSIZE );

#ifdef HAVE_LIBDEFLATE
	 fprintf(stderr,"LIBDEFLATE is defined!\n");
#else
	 fprintf(stderr,"LIBDEFLATE is undefined!\n");
#endif
	 //CRC = (uint32_t *)(cdata +BSIZE -8);
	 //CRC = (uint32_t)(buff +BSIZE -8);
	 ////memcpy( &CRC, cdata +BSIZE -8, 4 );
	 //CRC = le_to_u32(cdata + (uint8_t)BSIZE -8);
	 CRC = le_to_u32(cdata + BSIZE -8);    //Clear!
	 //CRC = *((uint32_u *) cdata +BSIZE -8);
	 fprintf(stderr, "CRC   = %d\n", CRC );

	 // Calc. crc for check
	 //CRC = crc32( crc32(0L,NULL,0L), (unsigned char*)uncomp, 0);
	 //fprintf(stderr, "CRC(calc) = %d\n", CRC );

	 size_t dlen = BGZF_MAX_BLOCK_SIZE;
	 fprintf(stderr,"dlen = %d\n", (int)dlen);
	 
	 bgzf_uncompress( uncomp, &dlen, (Bytef*)cdata +18, BSIZE -18, CRC );
	 //bgzf_uncompress( uncomp, &dlen, (Bytef*)cdata +18, BSIZE -18, *(cdata+BSIZE-8) );
	 store_read( uncomp, b);
	 exit(0);

	 //fprintf(stderr,"CDATA = %s", uncomp);
	 

	 
	 ret = bgzf_raw_write( out->fp.bgzf, buff, ret );
    }

	 
    // close bam
    //sam_hdr_destroy(h);
    bam_hdr_destroy(h);
    sam_close(in);
    sam_close(out);
    bam_destroy1(b);
   
    // close threads
    //hts_tpool_destroy(p.pool);
       
    return 0;
}
