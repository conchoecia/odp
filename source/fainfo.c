// fainfo.c — print: <prefix>\t<name>\t<length>
// Uses kseq.h from seqtk and zlib. Supports .fa/.fa.gz or stdin ("-").

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <unistd.h>   // for fileno()

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

static void usage(const char *prog) {
    fprintf(stderr, "Usage: %s <prefix> <fasta[.gz] | ->\n", prog);
}

int main(int argc, char *argv[]) {
    if (argc != 3) { usage(argv[0]); return 1; }

    const char *prefix = argv[1];
    const char *path   = argv[2];

    gzFile fp = NULL;
    if (strcmp(path, "-") == 0) {
        fp = gzdopen(fileno(stdin), "rb");
    } else {
        fp = gzopen(path, "rb"); // transparent for gzip or plain
    }
    if (fp == NULL) {
        fprintf(stderr, "Error: cannot open input '%s'\n", path);
        return 2;
    }

    kseq_t *seq = kseq_init(fp);
    if (!seq) {
        fprintf(stderr, "Error: kseq_init failed\n");
        gzclose(fp);
        return 3;
    }

    while (kseq_read(seq) >= 0) {
        printf("%s\t%s\t%zu\n", prefix, seq->name.s, (size_t)seq->seq.l);
    }

    kseq_destroy(seq);
    gzclose(fp);
    return 0;
}
