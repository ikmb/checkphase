/*SS
 *    Copyright (C) 2018-2025 by Lars Wienbrandt,
 *    Institute of Clinical Molecular Biology, Kiel University
 *
 *    This file is part of Checkphase.
 *
 *    Checkphase is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    Checkphase is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with Checkphase. If not, see <https://www.gnu.org/licenses/>.
 */

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

//#include <omp.h>

#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>

#include "Args.h"
#include "utils.h"

using namespace std;

inline size_t getNumVariantsFromIndex(const string &vcffilename) {
    size_t numvars = 0;

    htsFile* file = hts_open(vcffilename.c_str(), "r");
    if (!file) {
        cerr << "ERROR: Could not read file " << vcffilename << endl;
        exit(EXIT_FAILURE);
    }
    bcf_hdr_t *hdr = bcf_hdr_read(file);
    if (!hdr) {
        cerr << "ERROR: Could not read header from " << vcffilename << endl;
        exit(EXIT_FAILURE);
    }

    tbx_t *tbx = NULL;
    hts_idx_t *idx = NULL;

    if (hts_get_format(file)->format == vcf) {
        tbx = tbx_index_load(vcffilename.c_str());
        if (!tbx) {
            cerr << "ERROR: Could not read reference index file." << endl;
            exit(EXIT_FAILURE);
        }
    }
    else if (hts_get_format(file)->format == bcf)
    {
        idx = bcf_index_load(vcffilename.c_str());
        if (!idx) {
            cerr << "ERROR: Could not read reference index file." << endl;
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        cerr << "ERROR: Could not detect the reference file type as VCF or BCF." << endl;
        exit(EXIT_FAILURE);
    }

    int nseq;
    // get the number and names of sequences stored in the target file
    const char** seq = tbx ? tbx_seqnames(tbx, &nseq) : bcf_index_seqnames(idx, hdr, &nseq);
    // we only pick the first sequence that fits the convention, i.e. a number or literals optionally preceding with "chr" (but no special chars such as "_", ".", ...) and nrecords > 0
    for (int i = 0; i < nseq; i++) {
        uint64_t nrecords, unmapped;
        string seqstring(seq[i]);
//        if (seqstring.substr(0,3).compare("chr") == 0) // seqname starts with "chr"
//            seqstring = seqstring.substr(3);
        if (seqstring.find_first_not_of("0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ") == string::npos) { // only number or literals
            hts_idx_get_stat(tbx ? tbx->idx : idx, i, &nrecords, &unmapped);
//            cerr << "seq: " << i << " seqname: " << seq[i] << " nrecords: " << nrecords << " unmapped: " << unmapped << endl;
            if (nrecords) {
                numvars = nrecords;
                break;
            }
        }
    }
    free(seq); // allocated by HTSlib

    hts_close(file);
    bcf_hdr_destroy(hdr);
    if (tbx)
        tbx_destroy(tbx);
    if (idx)
        hts_idx_destroy(idx);

    return numvars;
}

inline void updateStatus(const char* statfile, float r, float q) {
    if (!statfile)
        return;
    ofstream stat(statfile, ofstream::trunc);
    if (!stat.fail())
        stat << "R: " << r << "\nQ: " << q << endl;
}

double calc_r2_hard(size_t n, size_t sumx, size_t sumx2, size_t sumy, size_t sumy2, size_t sumxy) {
    int64_t r_num = n * sumxy - sumx * sumy;
    int64_t r2_denom_ref = n * sumx2 - sumx * sumx;
    int64_t r2_denom_q = n * sumy2 - sumy * sumy;
    double r2 = (r_num / (double)r2_denom_ref) * (r_num / (double)r2_denom_q);
    return r2;
}

double calc_r2_soft(size_t n, size_t sumx, size_t sumx2, double sumy, double sumy2, double sumxy) {
    double r_num = n * sumxy - sumx * sumy;
    double r2_denom_ref = n * sumx2 - sumx * sumx;
    double r2_denom_q = n * sumy2 - sumy * sumy;
    double r2 = (r_num / r2_denom_ref) * (r_num / r2_denom_q);
    return r2;
}

int main(int argc, char *argv[]) {

    Args args = Args::parseArgs(argc, argv);

    const string& reffile = args.vcfRef;
    const string& queryfile = args.vcfQuery;
    const string& sharedfile = args.vcfShared;
    const char* statfile = args.statfile.empty() ? NULL : args.statfile.c_str();
    bool dump = args.dump;

    if (statfile)
        cout << "Statfile: " << statfile << endl;
    if (dump)
        cout << "--dump enabled. Will dump phase error positions to stderr." << endl;
    cout << endl;

    updateStatus(statfile, 0, 0);

    size_t Mrefidx = getNumVariantsFromIndex(reffile);
    size_t Mqidx = getNumVariantsFromIndex(queryfile);
    size_t Msharedidx = sharedfile.empty() ? 0 : getNumVariantsFromIndex(sharedfile);
    bool useshared = Msharedidx ? true : false;

    bcf_srs_t *sr = bcf_sr_init();

    bcf_sr_set_opt(sr, BCF_SR_PAIR_LOGIC, BCF_SR_PAIR_BOTH_REF);
    bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);

    if (!bcf_sr_add_reader(sr, reffile.c_str())) {
        cerr << "ERROR: Could not open reference for reading: " << bcf_sr_strerror(sr->errnum) << endl;
        exit(EXIT_FAILURE);
    }
    if (!bcf_sr_add_reader(sr, queryfile.c_str())) {
        cerr << "ERROR: Could not open query for reading: " << bcf_sr_strerror(sr->errnum) << endl;
        exit(EXIT_FAILURE);
    }
    if (useshared) {
        if (!bcf_sr_add_reader(sr, sharedfile.c_str())) {
            cerr << "ERROR: Could not open shared file for reading: " << bcf_sr_strerror(sr->errnum) << endl;
            exit(EXIT_FAILURE);
        }
    }

    bcf_hdr_t *ref_hdr = bcf_sr_get_header(sr, 0);
    bcf_hdr_t *q_hdr = bcf_sr_get_header(sr, 1);
    bcf_hdr_t *s_hdr = bcf_sr_get_header(sr, 2);

    size_t Nref = bcf_hdr_nsamples(ref_hdr);
    size_t Nquery = bcf_hdr_nsamples(q_hdr);

    // Read sample IDs
    vector<string> refIDs;
    vector<string> qIDs;
    refIDs.reserve(Nref);
    qIDs.reserve(Nquery);
    for (size_t i = 0; i < Nref; i++)
        refIDs.push_back(ref_hdr->samples[i]);
    for (size_t i = 0; i < Nquery; i++)
        qIDs.push_back(q_hdr->samples[i]);

    cout << "Input files:" << endl;
    cout << "  Reference:            " << reffile << endl;
    cout << "  Query:                " << queryfile << endl;
    if (useshared)
        cout << "  Shared variants file: " << sharedfile << endl;
    cout << "  Reference variants (from index): Mref = " << Mrefidx << endl;
    cout << "  Reference samples:               Nref = " << Nref << endl;
    cout << "  Query variants (from index):   Mquery = " << Mqidx << endl;
    cout << "  Query samples:                 Nquery = " << Nquery << endl;
    if (useshared)
        cout << "  Shared file variants (from index): Ms = " << Msharedidx << endl;
    cout << endl;

    // extract shared samples (shared samples have to be in the same order! all query samples have to be in the reference!)
    vector<size_t> map2ref(Nquery); // index mapping for shared samples from query to reference
    {
        size_t q = 0;
        for (size_t r = 0; r < Nref && q < Nquery; r++) {
            if (qIDs[q].compare(refIDs[r]) == 0) { // equal IDS -> add index to map and continue
                map2ref[q] = r;
                q++;
            }
        }
        if (q < Nquery) {
            cerr << "ERROR: Not all queries could be find in the reference." << endl;
            exit(EXIT_FAILURE);
        }
    }

    cout << "Passed samples check." << endl;

    cout << "Checking shared variants for phase and genotype errors..." << flush;

    size_t Mref = 0; // reference variants
    size_t Mq = 0;   // query variants
    size_t Ms = 0;   // variants from shared file
    size_t Mshared = 0; // shared variants
    size_t Mmaf01 = 0; // shared variants with RefPanelMAF >= 0.1
    size_t Mmaf001 = 0; // shared variants with RefPanelMAF >= 0.01
    size_t Mmaf0001 = 0; // shared variants with RefPanelMAF >= 0.001
    size_t Mmaf00001 = 0; // shared variants with RefPanelMAF >= 0.0001
    size_t Mtyped = 0; // shared variants with TYPED tag

    int mref_gt = 0; void *ref_gt = NULL; // will be allocated once in bcf_get_genotypes() and then reused for each marker (need void* because of htslib)
    int mtgt_gt = 0; void *tgt_gt = NULL; // will be allocated once in bcf_get_genotypes() and then reused for each marker (need void* because of htslib)
    int ndosbuf = 0; void *dosbuf = NULL;

    size_t MrefError = 0;
    size_t MqError = 0;
    size_t MRefAltSwap = 0;
    size_t MStrandFlip = 0;
    size_t MRefAltSwapAndStrandFlip = 0;
    size_t MAlleleDiff = 0;

    // for each tested sample the number of missing and unphased sites
    vector<size_t> MrefMissing (Nquery, 0);
    vector<size_t> MqMissing   (Nquery, 0);
    vector<size_t> MrefUnphasedHet(Nquery, 0);
    vector<size_t> MqUnphasedHet  (Nquery, 0);

    vector<bool> initialized(Nquery, false);
    vector<bool> initialized_typed(Nquery, false);
    vector<bool> switched(Nquery, false);
    vector<bool> switched_typed(Nquery, false);
    vector<size_t> switchErrors(Nquery, 0);
    vector<size_t> switchErrors_typed(Nquery, 0);
    vector<bool> matPatSwitches(Nquery, false);
    vector<bool> matPatSwitches_typed(Nquery, false);
    vector<size_t> gtErrors(Nquery, 0);
    vector<size_t> gtErrors_maf01(Nquery, 0);
    vector<size_t> gtErrors_maf001(Nquery, 0);
    vector<size_t> gtErrors_maf0001(Nquery, 0);
    vector<size_t> gtErrors_maf00001(Nquery, 0);
    vector<double> gtErrorsSoft(Nquery, 0.0);
    vector<double> gtErrorsSoft_maf01(Nquery, 0.0);
    vector<double> gtErrorsSoft_maf001(Nquery, 0.0);
    vector<double> gtErrorsSoft_maf0001(Nquery, 0.0);
    vector<double> gtErrorsSoft_maf00001(Nquery, 0.0);
    // correlation r2
    vector<size_t> gtSumRef(Nquery, 0.0);
    vector<size_t> gt2SumRef(Nquery, 0.0);
    vector<size_t> gtSumQ(Nquery, 0.0);
    vector<size_t> gt2SumQ(Nquery, 0.0);
    vector<size_t> gtSumRefQ(Nquery, 0.0);
    vector<vector<double>> gtDosSumQ;      // bins for a couple of variants, then for each query
    vector<vector<double>> gtDos2SumQ;     // bins for a couple of variants, then for each query
    vector<vector<double>> gtDosSumRefQ;   // bins for a couple of variants, then for each query
    vector<size_t> gtSumRef_maf01(Nquery, 0.0);
    vector<size_t> gt2SumRef_maf01(Nquery, 0.0);
    vector<size_t> gtSumQ_maf01(Nquery, 0.0);
    vector<size_t> gt2SumQ_maf01(Nquery, 0.0);
    vector<size_t> gtSumRefQ_maf01(Nquery, 0.0);
    vector<vector<double>> gtDosSumQ_maf01;
    vector<vector<double>> gtDos2SumQ_maf01;
    vector<vector<double>> gtDosSumRefQ_maf01;
    vector<size_t> gtSumRef_maf001(Nquery, 0.0);
    vector<size_t> gt2SumRef_maf001(Nquery, 0.0);
    vector<size_t> gtSumQ_maf001(Nquery, 0.0);
    vector<size_t> gt2SumQ_maf001(Nquery, 0.0);
    vector<size_t> gtSumRefQ_maf001(Nquery, 0.0);
    vector<vector<double>> gtDosSumQ_maf001;
    vector<vector<double>> gtDos2SumQ_maf001;
    vector<vector<double>> gtDosSumRefQ_maf001;
    vector<size_t> gtSumRef_maf0001(Nquery, 0.0);
    vector<size_t> gt2SumRef_maf0001(Nquery, 0.0);
    vector<size_t> gtSumQ_maf0001(Nquery, 0.0);
    vector<size_t> gt2SumQ_maf0001(Nquery, 0.0);
    vector<size_t> gtSumRefQ_maf0001(Nquery, 0.0);
    vector<vector<double>> gtDosSumQ_maf0001;
    vector<vector<double>> gtDos2SumQ_maf0001;
    vector<vector<double>> gtDosSumRefQ_maf0001;
    vector<size_t> gtSumRef_maf00001(Nquery, 0.0);
    vector<size_t> gt2SumRef_maf00001(Nquery, 0.0);
    vector<size_t> gtSumQ_maf00001(Nquery, 0.0);
    vector<size_t> gt2SumQ_maf00001(Nquery, 0.0);
    vector<size_t> gtSumRefQ_maf00001(Nquery, 0.0);
    vector<vector<double>> gtDosSumQ_maf00001;
    vector<vector<double>> gtDos2SumQ_maf00001;
    vector<vector<double>> gtDosSumRefQ_maf00001;
    vector<double> r2Sum;              // sum of per variant r2, bins for each couple of queries
    vector<double> r2Sum_maf01;        // sum of per variant r2 for variants of MAF>=0.1, bins for each couple of queries
    vector<double> r2Sum_maf001;       // sum of per variant r2 for variants of MAF>=0.01, bins for each couple of queries
    vector<double> r2Sum_maf0001;      // sum of per variant r2 for variants of MAF>=0.001, bins for each couple of queries
    vector<double> r2Sum_maf00001;     // sum of per variant r2 for variants of MAF>=0.0001, bins for each couple of queries
    vector<double> r2SoftSum;          // sum of per variant r2, bins for each couple of queries
    vector<double> r2SoftSum_maf01;    // sum of per variant r2 for variants of MAF>=0.1, bins for each couple of queries
    vector<double> r2SoftSum_maf001;   // sum of per variant r2 for variants of MAF>=0.01, bins for each couple of queries
    vector<double> r2SoftSum_maf0001;  // sum of per variant r2 for variants of MAF>=0.001, bins for each couple of queries
    vector<double> r2SoftSum_maf00001; // sum of per variant r2 for variants of MAF>=0.0001, bins for each couple of queries

    // used only when --dump is set
    vector<vector<size_t>> errPos(Nquery);

    // for counting haploid samples
    bool hapsset_ref = false;
    bool hapsset_tgt = false;
    size_t Nrefhap = 0;
    size_t Nqhap = 0;
    bool havedosages = false;
    bool havetyped = false;

    // for progress
    size_t linesperpercentref = divideRounded(Mrefidx, (size_t)100);
    size_t linesperpercentq = divideRounded(Mqidx, (size_t)100);
    if (!linesperpercentref) linesperpercentref = 1;
    if (!linesperpercentq) linesperpercentq = 1;
    updateStatus(statfile, 0, 0);

    // need to allocate memory for the return value for getting AF and TYPED info from htslib
    int af_size = sizeof(float);
    float *af_ptr = (float*) malloc(af_size);
    int typed_size = sizeof(int);
    int *typed_ptr = (int*) malloc(typed_size);

    // we sum up dosages in bins to compensate for numeric instability
    size_t currbin = 0;

    while (bcf_sr_next_line(sr)) { // read data SNP-wise in positional sorted order from query, reference and optional shared file

        bcf1_t *ref = bcf_sr_get_line(sr, 0); // read one line of reference, if available at current position (otherwise NULL)
        bcf1_t *tgt = bcf_sr_get_line(sr, 1); // read one line of query, if available at current position (otherwise NULL)
        bcf1_t *shd = NULL;
        if (useshared)
            shd = bcf_sr_get_line(sr, 2); // read one line of shared file, if available at current position (otherwise NULL)

        if (ref)
            Mref++;
        if (tgt)
            Mq++;
        if (shd)
            Ms++;

        if ((ref && Mref % linesperpercentref == 0) || (tgt && Mq % linesperpercentq == 0)) {
            if (ref && Mref % linesperpercentref == 0) {
                size_t p = (size_t)(100*Mref/(float)Mrefidx);
                if (p%5 == 0)
                    cout << p << "%.." << flush;
            }
            updateStatus(statfile, Mref/(float)Mrefidx, Mq/(float)Mqidx);
        }

        if (!ref || !tgt || (useshared && !shd)) { // query or reference only variant, or not in shared file
            continue;
        }

        // shared variant
        if (Mshared % 1024 == 0) {
            currbin = Mshared/1024;
            gtDosSumQ.push_back(vector<double>(Nquery, 0.0));
            gtDos2SumQ.push_back(vector<double>(Nquery, 0.0));
            gtDosSumRefQ.push_back(vector<double>(Nquery, 0.0));
            gtDosSumQ_maf01.push_back(vector<double>(Nquery, 0.0));
            gtDos2SumQ_maf01.push_back(vector<double>(Nquery, 0.0));
            gtDosSumRefQ_maf01.push_back(vector<double>(Nquery, 0.0));
            gtDosSumQ_maf001.push_back(vector<double>(Nquery, 0.0));
            gtDos2SumQ_maf001.push_back(vector<double>(Nquery, 0.0));
            gtDosSumRefQ_maf001.push_back(vector<double>(Nquery, 0.0));
            gtDosSumQ_maf0001.push_back(vector<double>(Nquery, 0.0));
            gtDos2SumQ_maf0001.push_back(vector<double>(Nquery, 0.0));
            gtDosSumRefQ_maf0001.push_back(vector<double>(Nquery, 0.0));
            gtDosSumQ_maf00001.push_back(vector<double>(Nquery, 0.0));
            gtDos2SumQ_maf00001.push_back(vector<double>(Nquery, 0.0));
            gtDosSumRefQ_maf00001.push_back(vector<double>(Nquery, 0.0));
            r2Sum.push_back(0.0);
            r2Sum_maf01.push_back(0.0);
            r2Sum_maf001.push_back(0.0);
            r2Sum_maf0001.push_back(0.0);
            r2Sum_maf00001.push_back(0.0);
            r2SoftSum.push_back(0.0);
            r2SoftSum_maf01.push_back(0.0);
            r2SoftSum_maf001.push_back(0.0);
            r2SoftSum_maf0001.push_back(0.0);
            r2SoftSum_maf00001.push_back(0.0);
        }
        Mshared++;

        // exclude monomorphic or multi-allelic markers
        if (ref->n_allele != 2 || tgt->n_allele != 2) {
            if (ref->n_allele != 2)
                MrefError++;
            if (tgt->n_allele != 2)
                MqError++;
            continue;
        }

        // decode genotypes/haplotypes
        size_t nref_gt = bcf_get_genotypes(ref_hdr, ref, &ref_gt, &mref_gt); // calls bcf_unpack() within
        size_t ntgt_gt = bcf_get_genotypes(q_hdr, tgt, &tgt_gt, &mtgt_gt); // calls bcf_unpack() within

        // imputation R2
        float maf = -1.0; // default if no MAF is present
        if (bcf_get_info_float(q_hdr, tgt, "RefPanelAF", (void*)&af_ptr, &af_size) >= 0                      // successfully read RefPanelAF tag from query
            || (useshared && bcf_get_info_float(s_hdr, shd, "RefPanelAF", (void*)&af_ptr, &af_size) >= 0)) { // or successfully read RefPanelAF tag from shared file
            if (*af_ptr <= 0.5)
                maf = *af_ptr; // if there's more than one RefPanelAF entry, take the first.
            else
                maf = 1.0 - *af_ptr;
        }

        // genotyped or not?
        bool typed = false;
        if (bcf_get_info_flag(q_hdr, tgt, "TYPED", (void*)&typed_ptr, &typed_size) > 0) { // successfully read TYPED tag and thus, it was set
//            if (!useshared || bcf_get_info_flag(s_hdr, shd, "TYPED", (void*)&typed_ptr, &typed_size) > 0) { // also successfully read TYPED tag from shared variants file if required
                typed = true;
                havetyped = true;
//            }
        }

        // get dosages, if present
        int nret = bcf_get_format_values(q_hdr,tgt,"ADS",(void**)&dosbuf,&ndosbuf,BCF_HT_REAL); // ADS type
        if (nret <= 0) nret = bcf_get_format_values(q_hdr,tgt,"HDS",(void**)&dosbuf,&ndosbuf,BCF_HT_REAL); // HDS type
        if (nret > 0) {
            // check number of dosages
            if ((size_t)nret != 2*Nquery && (size_t)nret != Nquery) {
               cerr << "ERROR: called dosage number is not as expected! nret = " << nret << " (expected: " << Nquery << " or " << 2*Nquery << ")" << endl;
               exit(EXIT_FAILURE);
            }
            havedosages = true;
        }

        // check alleles
        bool refaltswap = false;
        if (strcmp(tgt->d.allele[0], ref->d.allele[0]) == 0 && strcmp(tgt->d.allele[1], ref->d.allele[1]) == 0) { // all good
            // nothing to do here
        } else if (strcmp(tgt->d.allele[0], ref->d.allele[1]) == 0 && strcmp(tgt->d.allele[1], ref->d.allele[0]) == 0) { // switched alleles
            refaltswap = true;
            MRefAltSwap++;
        } else if (reverseComplement(tgt->d.allele[0]).compare(ref->d.allele[0]) == 0 && reverseComplement(tgt->d.allele[1]).compare(ref->d.allele[1]) == 0) { // strand flip + switched alleles
            MStrandFlip++;
        } else if (reverseComplement(tgt->d.allele[0]).compare(ref->d.allele[1]) == 0 && reverseComplement(tgt->d.allele[1]).compare(ref->d.allele[0]) == 0) { // strand flip + switched alleles
            refaltswap = true;
            MRefAltSwapAndStrandFlip++;
        } else { // different alleles -> next variant
            MAlleleDiff++;
            continue;
        }

        // check alleles also with shared vars file
        if (useshared) {
            if (strcmp(tgt->d.allele[0], shd->d.allele[0]) == 0 && strcmp(tgt->d.allele[1], shd->d.allele[1]) == 0) { // all good
            } else if (strcmp(tgt->d.allele[0], shd->d.allele[1]) == 0 && strcmp(tgt->d.allele[1], shd->d.allele[0]) == 0) { // switched alleles
            } else if (reverseComplement(tgt->d.allele[0]).compare(shd->d.allele[0]) == 0 && reverseComplement(tgt->d.allele[1]).compare(shd->d.allele[1]) == 0) { // strand flip + switched alleles
            } else if (reverseComplement(tgt->d.allele[0]).compare(shd->d.allele[1]) == 0 && reverseComplement(tgt->d.allele[1]).compare(shd->d.allele[0]) == 0) { // strand flip + switched alleles
            } else { // different alleles -> next variant
                MAlleleDiff++;
                continue;
            }
        }

        // check number of called genotypes
        if ((nref_gt != 2*Nref && nref_gt != Nref) || (ntgt_gt != 2*Nquery && ntgt_gt != Nquery)) {
            cerr << "ERROR: called genotype number is not as expected! nref_gt = " << nref_gt << " (expected: " << Nref << " or " << 2*Nref << ") ntgt_gt = " << ntgt_gt << " (expected: " << Nquery << " or " << 2*Nquery << ")" << endl;
            exit(EXIT_FAILURE);
        }

        // check if variant is haploid
        bool haploid_ref = nref_gt == Nref;
        bool haploid_tgt = ntgt_gt == Nquery;
        if  (haploid_ref) {
            if (!hapsset_ref) {
                hapsset_ref = true;
                Nrefhap = Nref;
            } // else could perhaps throw an error if numbers don't match??
        }
        if  (haploid_tgt) {
            if (!hapsset_tgt) {
                hapsset_tgt = true;
                Nqhap = Nquery;
            } // else could perhaps throw an error if numbers don't match??
        }

        // count in MAF categories
        if (maf >= 0.1)
            Mmaf01++;
        if (maf >= 0.01)
            Mmaf001++;
        if (maf >= 0.001)
            Mmaf0001++;
        if (maf >= 0.0001)
            Mmaf00001++;

        // count checked typed variants
        if (typed)
            Mtyped++;

        // for variant-wise correlation (PCC)
        size_t gtsumref = 0;
        size_t gt2sumref = 0;
        size_t gtsumq = 0;
        size_t gt2sumq = 0;
        size_t gtsumrefq = 0;
        double gtdossumq = 0.0;
        double gtdos2sumq = 0.0;
        double gtdossumrefq = 0.0;

        // check samples
        for (size_t q = 0; q < Nquery; q++) {

            // haploid variants will be encoded as homozygous diploid
            int32_t refmat_i = haploid_ref ? ((int32_t*)ref_gt)[map2ref[q]] : ((int32_t*)ref_gt)[2*map2ref[q]];
            int32_t refpat_i = haploid_ref ? refmat_i : ((int32_t*)ref_gt)[2*map2ref[q]+1];
            int32_t qmat_i = haploid_tgt ? ((int32_t*)tgt_gt)[q] : ((int32_t*)tgt_gt)[2*q];
            int32_t qpat_i = haploid_tgt ? qmat_i : ((int32_t*)tgt_gt)[2*q+1];
            float qdosmat = 0.0, qdospat = 0.0;
            if (havedosages) {
                qdosmat = haploid_tgt ? ((float*)dosbuf)[q] : ((float*)dosbuf)[2*q];
                qdospat = haploid_tgt ? qdosmat : ((float*)dosbuf)[2*q+1];
            }

            // check if missing
            if (bcf_gt_is_missing(refmat_i) || bcf_gt_is_missing(refpat_i) || bcf_gt_is_missing(qmat_i) || bcf_gt_is_missing(qpat_i)) {
                if (bcf_gt_is_missing(refmat_i) || bcf_gt_is_missing(refpat_i))
                    MrefMissing[q]++;
                if (bcf_gt_is_missing(qmat_i) || bcf_gt_is_missing(qpat_i))
                    MqMissing[q]++;
                continue;
            }

            // check if unphased or haploid -> only gt check
            bool phased = true;
            bool diploid = true;
            if (refpat_i == bcf_int32_vector_end || qpat_i == bcf_int32_vector_end || !bcf_gt_is_phased(refpat_i) || !bcf_gt_is_phased(qpat_i)) {
                phased = false;

                if (refpat_i == bcf_int32_vector_end) {
                    refpat_i = refmat_i; // treat as homozygous diploid
                    if (!hapsset_ref)
                        Nrefhap++;
                    diploid = false;
                }
                if (qpat_i == bcf_int32_vector_end) {
                    qpat_i = qmat_i; // treat as homozygous diploid
                    qdospat = qdosmat;
                    if (!hapsset_tgt)
                        Nqhap++;
                    diploid = false;
                }
            }

            // decode haplotypes
            bool refmat = bcf_gt_allele(refmat_i);
            bool refpat = bcf_gt_allele(refpat_i);
            bool qmat = bcf_gt_allele(qmat_i);
            bool qpat = bcf_gt_allele(qpat_i);

            // check if dosage is a number, otherwise set dosage corresponding to allele
            if (isnan(qdosmat)) qdosmat = qmat ? 1.0 : 0.0;
            if (isnan(qdospat)) qdospat = qpat ? 1.0 : 0.0;

            // apply ref/alt swap
            if (refaltswap) {
                qmat = !qmat;
                qpat = !qpat;
                qdosmat = 1.0 - qdosmat;
                qdospat = 1.0 - qdospat;
            }

            // gt check
            int refgt = (refmat ? 1 : 0) + (diploid ? (refpat ? 1 : 0) : 0);
            int qgt = (qmat ? 1 : 0) + (diploid ? (qpat ? 1 : 0) : 0);
            // for complete and sample-wise correlation
            gtSumRef[q] += refgt;
            gt2SumRef[q] += refgt*refgt;
            gtSumQ[q] += qgt;
            gt2SumQ[q] += qgt*qgt;
            gtSumRefQ[q] += refgt * qgt;
            if (maf >= 0.1) {
                gtSumRef_maf01[q] += refgt;
                gt2SumRef_maf01[q] += refgt*refgt;
                gtSumQ_maf01[q] += qgt;
                gt2SumQ_maf01[q] += qgt*qgt;
                gtSumRefQ_maf01[q] += refgt * qgt;
            }
            if (maf >= 0.01) {
                gtSumRef_maf001[q] += refgt;
                gt2SumRef_maf001[q] += refgt*refgt;
                gtSumQ_maf001[q] += qgt;
                gt2SumQ_maf001[q] += qgt*qgt;
                gtSumRefQ_maf001[q] += refgt * qgt;
            }
            if (maf >= 0.001) {
                gtSumRef_maf0001[q] += refgt;
                gt2SumRef_maf0001[q] += refgt*refgt;
                gtSumQ_maf0001[q] += qgt;
                gt2SumQ_maf0001[q] += qgt*qgt;
                gtSumRefQ_maf0001[q] += refgt * qgt;
            }
            if (maf >= 0.0001) {
                gtSumRef_maf00001[q] += refgt;
                gt2SumRef_maf00001[q] += refgt*refgt;
                gtSumQ_maf00001[q] += qgt;
                gt2SumQ_maf00001[q] += qgt*qgt;
                gtSumRefQ_maf00001[q] += refgt * qgt;
            }
            // for variant-wise correlation
            gtsumref += refgt;
            gt2sumref += refgt*refgt;
            gtsumq += qgt;
            gt2sumq += qgt*qgt;
            gtsumrefq += refgt * qgt;

            // from here haploids are treated as homozygous diploid
            if (!diploid) {
                refgt += (refpat ? 1 : 0);
                qgt += (qpat ? 1 : 0);
            }
            if (refgt == 1 && !phased)
                MrefUnphasedHet[q]++;
            if (qgt == 1 && !phased)
                MqUnphasedHet[q]++;
            // soft check
            if (havedosages) {
                double err;
                switch (refgt) {
                case 0: // homozygous wild
                    //this is wrong! gtErrorsSoft[q] += diploid ? qdosmat + qdospat : qdosmat;
                    if (!diploid)
                        err = qdosmat;
                    else {
                        double gp = (1.0-qdosmat)*(1.0-qdospat); // gp0
                        err = 1.0 - gp;
                    }
                    break;
                case 2: // homozygous variant
                    //this is wrong! gtErrorsSoft[q] += diploid ? 2.0 - qdosmat - qdospat : 1.0 - qdosmat;
                    if (!diploid)
                        err = 1.0 - qdosmat;
                    else {
                        double gp = qdosmat*qdospat; // gp2
                        err = 1.0 - gp;
                    }
                    break;
                default: // heterozygous (not applicable for haploid)
                    //this is wrong!
                    //float err1 = 1.0 - qdosmat + qdospat; // if mat=1 and pat=0
                    //float err2 = 1.0 - qdospat + qdosmat; // if mat=0 and pat=1
                    //gtErrorsSoft[q] += min(err1, err2); // we count the minimum deviation as we compare only the genotype and do not consider a phase switch here
                    err = (1.0-qdosmat)*(1.0-qdospat) + qdosmat*qdospat; // 1-gp1 = 1-(1-gp0-gp2) = gp0+gp2
                }
                gtErrorsSoft[q] += err;
                if (maf >= 0.1)
                    gtErrorsSoft_maf01[q] += err;
                if (maf >= 0.01)
                    gtErrorsSoft_maf001[q] += err;
                if (maf >= 0.001)
                    gtErrorsSoft_maf0001[q] += err;
                if (maf >= 0.0001)
                    gtErrorsSoft_maf00001[q] += err;
                // genotype dosage = ads0 + ads1
                double gtdos = qdosmat + (diploid ? qdospat : 0);
                gtDosSumQ[currbin][q] += gtdos;
                gtDos2SumQ[currbin][q] += gtdos*gtdos;
                gtDosSumRefQ[currbin][q] += (refgt * gtdos) / (diploid ? 1 : 2); // need to reduce homozygous diploid representation back to haploid in the case of a haploid sample
                if (maf >= 0.1) {
                    gtDosSumQ_maf01[currbin][q] += gtdos;
                    gtDos2SumQ_maf01[currbin][q] += gtdos*gtdos;
                    gtDosSumRefQ_maf01[currbin][q] += (refgt * gtdos) / (diploid ? 1 : 2); // need to reduce homozygous diploid representation back to haploid in the case of a haploid sample
                }
                if (maf >= 0.01) {
                    gtDosSumQ_maf001[currbin][q] += gtdos;
                    gtDos2SumQ_maf001[currbin][q] += gtdos*gtdos;
                    gtDosSumRefQ_maf001[currbin][q] += (refgt * gtdos) / (diploid ? 1 : 2); // need to reduce homozygous diploid representation back to haploid in the case of a haploid sample
                }
                if (maf >= 0.001) {
                    gtDosSumQ_maf0001[currbin][q] += gtdos;
                    gtDos2SumQ_maf0001[currbin][q] += gtdos*gtdos;
                    gtDosSumRefQ_maf0001[currbin][q] += (refgt * gtdos) / (diploid ? 1 : 2); // need to reduce homozygous diploid representation back to haploid in the case of a haploid sample
                }
                if (maf >= 0.0001) {
                    gtDosSumQ_maf00001[currbin][q] += gtdos;
                    gtDos2SumQ_maf00001[currbin][q] += gtdos*gtdos;
                    gtDosSumRefQ_maf00001[currbin][q] += (refgt * gtdos) / (diploid ? 1 : 2); // need to reduce homozygous diploid representation back to haploid in the case of a haploid sample
                }
                gtdossumq += gtdos;
                gtdos2sumq += gtdos*gtdos;
                gtdossumrefq += (refgt * gtdos) / (diploid ? 1 : 2);

            }
            // hard check
            if (refgt != qgt) { // genotype error -> continue with next sample
                gtErrors[q]++;
                if (maf >= 0.1)
                    gtErrors_maf01[q]++;
                if (maf >= 0.01)
                    gtErrors_maf001[q]++;
                if (maf >= 0.001)
                    gtErrors_maf0001[q]++;
                if (maf >= 0.0001)
                    gtErrors_maf00001[q]++;
                continue;
            }

            // phase check
            if (refgt == 1 && phased) { // only proceed on heterozygous site and if site is phased (gt errors are already excluded above)
                if (initialized[q]) { // this is not the first het for this sample -> wrong phase is a switch error!
                    //bool swerr = (refmat != qmat) xor switched[q];
                    bool swerr = switched[q] ? (refmat == qmat) : (refmat != qmat);
                    if (swerr) {
                        switched[q] = !switched[q];
                        switchErrors[q]++;
                        if (dump)
                            errPos[q].push_back(Mq-1); // error position is 0-based
                    }
                } else { // first het site -> simply set the switched-flag according to the current phases in ref and query
                    initialized[q] = true;
                    if (refmat != qmat) {
                        switched[q] = true;
                        matPatSwitches[q] = true;
                        if (dump)
                            errPos[q].push_back(0); // indicates the difference at the first het, for debugging
                    }
                }
                if (typed) {
                    if (initialized_typed[q]) { // this is not the first het for this sample -> wrong phase is a switch error!
                        //bool swerr = (refmat != qmat) xor switched[q];
                        bool swerr = switched_typed[q] ? (refmat == qmat) : (refmat != qmat);
                        if (swerr) {
                            switched_typed[q] = !switched_typed[q];
                            switchErrors_typed[q]++;
                        }
                    } else { // first het site -> simply set the switched-flag according to the current phases in ref and query
                        initialized_typed[q] = true;
                        if (refmat != qmat) {
                            switched_typed[q] = true;
                            matPatSwitches_typed[q] = true;
                        }
                    }
                }
            }
        } // end for every sample

        // variant-wise correlation r2
        double r2 = calc_r2_hard(Nquery, gtsumref, gt2sumref, gtsumq, gt2sumq, gtsumrefq);
        if (!isnan(r2)) { // add only if it is a number, otherwise it's treated as zero
            r2Sum[currbin] += r2;
            if (maf >= 0.1)
                r2Sum_maf01[currbin] += r2;
            if (maf >= 0.01)
                r2Sum_maf001[currbin] += r2;
            if (maf >= 0.001)
                r2Sum_maf0001[currbin] += r2;
            if (maf >= 0.0001)
                r2Sum_maf00001[currbin] += r2;
        }
        double r2soft = 0.0;
        if (havedosages) {
            r2soft = calc_r2_soft(Nquery, gtsumref, gt2sumref, gtdossumq, gtdos2sumq, gtdossumrefq);
            if (!isnan(r2soft)) { // add only if it is a number, otherwise it's treated as zero
                r2SoftSum[currbin] += r2soft;
                if (maf >= 0.1)
                    r2SoftSum_maf01[currbin] += r2soft;
                if (maf >= 0.01)
                    r2SoftSum_maf001[currbin] += r2soft;
                if (maf >= 0.001)
                    r2SoftSum_maf0001[currbin] += r2soft;
                if (maf >= 0.0001)
                    r2SoftSum_maf00001[currbin] += r2soft;
            }
        }

        if (!hapsset_ref) {
            if (Nrefhap)
                hapsset_ref = true;
        } // else could perhaps throw an error if numbers don't match??
        if (!hapsset_tgt) {
            if (Nqhap)
                hapsset_tgt = true;
        } // else could perhaps throw an error if numbers don't match??

    } // end while read line
    cout << " done." << endl;
    updateStatus(statfile, 1, 1);
    free(af_ptr);
    free(typed_ptr);

    size_t Mexclude = MrefError + MqError + MAlleleDiff;
    size_t Mcheck = Mshared - Mexclude;

    cout << "\nSummary:\n" << endl;

    cout << "  Haploid ref samples:      " << Nrefhap << endl;
    cout << "  Haploid query samples:    " << Nqhap << endl;
    cout << endl;

    cout << "  Reference variants:       " << Mref << endl;
    cout << "  Query variants:           " << Mq << endl;
    cout << "  Shared variants:          " << Mshared << endl;
    cout << "  Checked variants:         " << Mcheck << endl;
    cout << "  Checked MAF>=0.1:         " << Mmaf01 << endl;
    cout << "  Checked MAF>=0.01:        " << Mmaf001 << endl;
    cout << "  Checked MAF>=0.001:       " << Mmaf0001 << endl;
    cout << "  Checked MAF>=0.0001:      " << Mmaf00001 << endl;
    cout << "  Checked typed variants:   " << Mtyped << endl;
    cout << endl;

    cout << "  Excluded shared:          " << Mexclude << endl;
    cout << "    Not biallelic in ref:   " << MrefError << endl;
    cout << "    Not biallelic in query: " << MqError << endl;
    cout << "    Alleles do not match:   " << MAlleleDiff << endl;
    cout << endl;

    cout << "  Ref/Alt swaps:            " << MRefAltSwap << endl;
    cout << "  Strand flips:             " << MStrandFlip << endl;
    cout << "  Ref/Alt + Strand flip:    " << MRefAltSwapAndStrandFlip << endl;
    cout << endl;

    size_t totalRefMissing = 0;
    for (size_t v : MrefMissing)
        totalRefMissing += v;

    size_t totalQMissing = 0;
    for (size_t v : MqMissing)
        totalQMissing += v;

    size_t totalRefUnphasedHet = 0;
    for (size_t v : MrefUnphasedHet)
        totalRefUnphasedHet += v;

    size_t totalQUnphasedHet = 0;
    for (size_t v : MqUnphasedHet)
        totalQUnphasedHet += v;

    cout << "  Missing or unphased sites:" << endl;
    cout << "    Total missing in ref:        " << totalRefMissing << endl;
    cout << "    Total missing in query:      " << totalQMissing << endl;
    cout << "    Total unphased het in ref:   " << totalRefUnphasedHet << endl;
    cout << "    Total unphased het in query: " << totalQUnphasedHet << endl;
    cout << endl;

    cout << fixed; // omit scientific notation
    cout << setprecision(8);

    // for correlation (PCC)
    size_t gtsumref = 0;
    size_t gt2sumref = 0;
    size_t gtsumq = 0;
    size_t gt2sumq = 0;
    size_t gtsumrefq = 0;
    double gtdossumq = 0.0;
    double gtdos2sumq = 0.0;
    double gtdossumrefq = 0.0;
    size_t gtsumref_maf01 = 0;
    size_t gt2sumref_maf01 = 0;
    size_t gtsumq_maf01 = 0;
    size_t gt2sumq_maf01 = 0;
    size_t gtsumrefq_maf01 = 0;
    double gtdossumq_maf01 = 0.0;
    double gtdos2sumq_maf01 = 0.0;
    double gtdossumrefq_maf01 = 0.0;
    size_t gtsumref_maf001 = 0;
    size_t gt2sumref_maf001 = 0;
    size_t gtsumq_maf001 = 0;
    size_t gt2sumq_maf001 = 0;
    size_t gtsumrefq_maf001 = 0;
    double gtdossumq_maf001 = 0.0;
    double gtdos2sumq_maf001 = 0.0;
    double gtdossumrefq_maf001 = 0.0;
    size_t gtsumref_maf0001 = 0;
    size_t gt2sumref_maf0001 = 0;
    size_t gtsumq_maf0001 = 0;
    size_t gt2sumq_maf0001 = 0;
    size_t gtsumrefq_maf0001 = 0;
    double gtdossumq_maf0001 = 0.0;
    double gtdos2sumq_maf0001 = 0.0;
    double gtdossumrefq_maf0001 = 0.0;
    size_t gtsumref_maf00001 = 0;
    size_t gt2sumref_maf00001 = 0;
    size_t gtsumq_maf00001 = 0;
    size_t gt2sumq_maf00001 = 0;
    size_t gtsumrefq_maf00001 = 0;
    double gtdossumq_maf00001 = 0.0;
    double gtdos2sumq_maf00001 = 0.0;
    double gtdossumrefq_maf00001 = 0.0;
    vector<double> gtdossumqv(Nquery, 0.0);
    vector<double> gtdos2sumqv(Nquery, 0.0);
    vector<double> gtdossumrefqv(Nquery, 0.0);
    vector<double> gtdossumqv_maf01(Nquery, 0.0);
    vector<double> gtdos2sumqv_maf01(Nquery, 0.0);
    vector<double> gtdossumrefqv_maf01(Nquery, 0.0);
    vector<double> gtdossumqv_maf001(Nquery, 0.0);
    vector<double> gtdos2sumqv_maf001(Nquery, 0.0);
    vector<double> gtdossumrefqv_maf001(Nquery, 0.0);
    vector<double> gtdossumqv_maf0001(Nquery, 0.0);
    vector<double> gtdos2sumqv_maf0001(Nquery, 0.0);
    vector<double> gtdossumrefqv_maf0001(Nquery, 0.0);
    vector<double> gtdossumqv_maf00001(Nquery, 0.0);
    vector<double> gtdos2sumqv_maf00001(Nquery, 0.0);
    vector<double> gtdossumrefqv_maf00001(Nquery, 0.0);
    for (size_t q = 0; q < Nquery; q++) {
       gtsumref  += gtSumRef[q];
       gt2sumref += gt2SumRef[q];
       gtsumq    += gtSumQ[q];
       gt2sumq   += gt2SumQ[q];
       gtsumrefq += gtSumRefQ[q];
       gtsumref_maf01  += gtSumRef_maf01[q];
       gt2sumref_maf01 += gt2SumRef_maf01[q];
       gtsumq_maf01    += gtSumQ_maf01[q];
       gt2sumq_maf01   += gt2SumQ_maf01[q];
       gtsumrefq_maf01 += gtSumRefQ_maf01[q];
       gtsumref_maf001  += gtSumRef_maf001[q];
       gt2sumref_maf001 += gt2SumRef_maf001[q];
       gtsumq_maf001    += gtSumQ_maf001[q];
       gt2sumq_maf001   += gt2SumQ_maf001[q];
       gtsumrefq_maf001 += gtSumRefQ_maf001[q];
       gtsumref_maf0001  += gtSumRef_maf0001[q];
       gt2sumref_maf0001 += gt2SumRef_maf0001[q];
       gtsumq_maf0001    += gtSumQ_maf0001[q];
       gt2sumq_maf0001   += gt2SumQ_maf0001[q];
       gtsumrefq_maf0001 += gtSumRefQ_maf0001[q];
       gtsumref_maf00001  += gtSumRef_maf00001[q];
       gt2sumref_maf00001 += gt2SumRef_maf00001[q];
       gtsumq_maf00001    += gtSumQ_maf00001[q];
       gt2sumq_maf00001   += gt2SumQ_maf00001[q];
       gtsumrefq_maf00001 += gtSumRefQ_maf00001[q];
    }
    if (havedosages) {
        for (size_t bin = 0; bin < gtDosSumQ.size(); bin++) {
            for (size_t q = 0; q < Nquery; q++) {
               gtdossumqv[q]    += gtDosSumQ[bin][q];
               gtdos2sumqv[q]   += gtDos2SumQ[bin][q];
               gtdossumrefqv[q] += gtDosSumRefQ[bin][q];
               gtdossumqv_maf01[q]    += gtDosSumQ_maf01[bin][q];
               gtdos2sumqv_maf01[q]   += gtDos2SumQ_maf01[bin][q];
               gtdossumrefqv_maf01[q] += gtDosSumRefQ_maf01[bin][q];
               gtdossumqv_maf001[q]    += gtDosSumQ_maf001[bin][q];
               gtdos2sumqv_maf001[q]   += gtDos2SumQ_maf001[bin][q];
               gtdossumrefqv_maf001[q] += gtDosSumRefQ_maf001[bin][q];
               gtdossumqv_maf0001[q]    += gtDosSumQ_maf0001[bin][q];
               gtdos2sumqv_maf0001[q]   += gtDos2SumQ_maf0001[bin][q];
               gtdossumrefqv_maf0001[q] += gtDosSumRefQ_maf0001[bin][q];
               gtdossumqv_maf00001[q]    += gtDosSumQ_maf00001[bin][q];
               gtdos2sumqv_maf00001[q]   += gtDos2SumQ_maf00001[bin][q];
               gtdossumrefqv_maf00001[q] += gtDosSumRefQ_maf00001[bin][q];
            }
        }
        for (size_t q = 0; q < Nquery; q++) {
            gtdossumq    += gtdossumqv[q];
            gtdos2sumq   += gtdos2sumqv[q];
            gtdossumrefq += gtdossumrefqv[q];
            gtdossumq_maf01    += gtdossumqv_maf01[q];
            gtdos2sumq_maf01   += gtdos2sumqv_maf01[q];
            gtdossumrefq_maf01 += gtdossumrefqv_maf01[q];
            gtdossumq_maf001    += gtdossumqv_maf001[q];
            gtdos2sumq_maf001   += gtdos2sumqv_maf001[q];
            gtdossumrefq_maf001 += gtdossumrefqv_maf001[q];
            gtdossumq_maf0001    += gtdossumqv_maf0001[q];
            gtdos2sumq_maf0001   += gtdos2sumqv_maf0001[q];
            gtdossumrefq_maf0001 += gtdossumrefqv_maf0001[q];
            gtdossumq_maf00001    += gtdossumqv_maf00001[q];
            gtdos2sumq_maf00001   += gtdos2sumqv_maf00001[q];
            gtdossumrefq_maf00001 += gtdossumrefqv_maf00001[q];
        }
    }

    {
        size_t totalGtErrors = 0;
        size_t totalGtErrors_maf01 = 0;
        size_t totalGtErrors_maf001 = 0;
        size_t totalGtErrors_maf0001 = 0;
        size_t totalGtErrors_maf00001 = 0;
        size_t gtErrorMin = 0xffffffffffffffffull;
        size_t gtErrorMax = 0;
        // calc total errors and identify min and max
        for (size_t err : gtErrors) {
            totalGtErrors += err;
            if (gtErrorMin > err) {
                gtErrorMin = err;
            }
            if (gtErrorMax < err) {
                gtErrorMax = err;
            }
        }
        for (size_t err : gtErrors_maf01)
            totalGtErrors_maf01 += err;
        for (size_t err : gtErrors_maf001)
            totalGtErrors_maf001 += err;
        for (size_t err : gtErrors_maf0001)
            totalGtErrors_maf0001 += err;
        for (size_t err : gtErrors_maf00001)
            totalGtErrors_maf00001 += err;
        double avgterr = totalGtErrors / (double) Nquery;
        double mingterrrate = gtErrorMin / (double)Mcheck;
        double maxgterrrate = gtErrorMax / (double)Mcheck;
        double avgterrrate  = totalGtErrors / (double)(Mcheck * Nquery);
        double avgterrrate_maf01 = totalGtErrors_maf01 / (double)(Mmaf01 * Nquery);
        double avgterrrate_maf001 = totalGtErrors_maf001 / (double)(Mmaf001 * Nquery);
        double avgterrrate_maf0001 = totalGtErrors_maf0001 / (double)(Mmaf0001 * Nquery);
        double avgterrrate_maf00001 = totalGtErrors_maf00001 / (double)(Mmaf00001 * Nquery);
        // calc standard deviation and variance for error rates
        double totgtdev2 = 0.0;
        for (size_t err : gtErrors) {
            double diff = err/(double)Mcheck - avgterrrate;
            totgtdev2 += diff*diff;
        }
        double gervar = totgtdev2 / (double)Nquery;
        double gerdev = sqrt(gervar);

        // correlation r2 (complete)
        double r2 = calc_r2_hard(Mcheck*Nquery, gtsumref, gt2sumref, gtsumq, gt2sumq, gtsumrefq);
        double r2_maf01 = calc_r2_hard(Mmaf01*Nquery, gtsumref_maf01, gt2sumref_maf01, gtsumq_maf01, gt2sumq_maf01, gtsumrefq_maf01);
        double r2_maf001 = calc_r2_hard(Mmaf001*Nquery, gtsumref_maf001, gt2sumref_maf001, gtsumq_maf001, gt2sumq_maf001, gtsumrefq_maf001);
        double r2_maf0001 = calc_r2_hard(Mmaf0001*Nquery, gtsumref_maf0001, gt2sumref_maf0001, gtsumq_maf0001, gt2sumq_maf0001, gtsumrefq_maf0001);
        double r2_maf00001 = calc_r2_hard(Mmaf00001*Nquery, gtsumref_maf00001, gt2sumref_maf00001, gtsumq_maf00001, gt2sumq_maf00001, gtsumrefq_maf00001);

        // correlation r2 (sample-wise)
        double r2_sav = 0.0;
        double r2_sav_maf01 = 0.0;
        double r2_sav_maf001 = 0.0;
        double r2_sav_maf0001 = 0.0;
        double r2_sav_maf00001 = 0.0;
        double r2_smin = 999.0;
        double r2_smax = 0.0;
        // store values for each sample (not required yet, but can be used later)
        vector<double> r2v(Nquery);
        vector<double> r2v_maf01(Nquery);
        vector<double> r2v_maf001(Nquery);
        vector<double> r2v_maf0001(Nquery);
        vector<double> r2v_maf00001(Nquery);
        for (size_t q = 0; q < Nquery; q++) {
            r2v[q] = calc_r2_hard(Mcheck, gtSumRef[q], gt2SumRef[q], gtSumQ[q], gt2SumQ[q], gtSumRefQ[q]);
            r2v_maf01[q] = calc_r2_hard(Mmaf01, gtSumRef_maf01[q], gt2SumRef_maf01[q], gtSumQ_maf01[q], gt2SumQ_maf01[q], gtSumRefQ_maf01[q]);
            r2v_maf001[q] = calc_r2_hard(Mmaf001, gtSumRef_maf001[q], gt2SumRef_maf001[q], gtSumQ_maf001[q], gt2SumQ_maf001[q], gtSumRefQ_maf001[q]);
            r2v_maf0001[q] = calc_r2_hard(Mmaf0001, gtSumRef_maf0001[q], gt2SumRef_maf0001[q], gtSumQ_maf0001[q], gt2SumQ_maf0001[q], gtSumRefQ_maf0001[q]);
            r2v_maf00001[q] = calc_r2_hard(Mmaf00001, gtSumRef_maf00001[q], gt2SumRef_maf00001[q], gtSumQ_maf00001[q], gt2SumQ_maf00001[q], gtSumRefQ_maf00001[q]);
            r2_sav += r2v[q]; // for the average
            r2_sav_maf01    += r2v_maf01[q]; // for the average
            r2_sav_maf001   += r2v_maf001[q]; // for the average
            r2_sav_maf0001  += r2v_maf0001[q]; // for the average
            r2_sav_maf00001 += r2v_maf00001[q]; // for the average
            if (r2v[q] < r2_smin)
                r2_smin = r2v[q];
            if (r2v[q] > r2_smax)
                r2_smax = r2v[q];
        }
        r2_sav /= Nquery;
        r2_sav_maf01 /= Nquery;
        r2_sav_maf001 /= Nquery;
        r2_sav_maf0001 /= Nquery;
        r2_sav_maf00001 /= Nquery;

        // correlation r2 (variant-wise)
        double r2_vav = 0.0;
        double r2_vav_maf01 = 0.0;
        double r2_vav_maf001 = 0.0;
        double r2_vav_maf0001 = 0.0;
        double r2_vav_maf00001 = 0.0;
        for (size_t bin = 0; bin < r2Sum.size(); bin++) {
            r2_vav += r2Sum[bin];
            r2_vav_maf01 += r2Sum_maf01[bin];
            r2_vav_maf001 += r2Sum_maf001[bin];
            r2_vav_maf0001 += r2Sum_maf0001[bin];
            r2_vav_maf00001 += r2Sum_maf00001[bin];
        }
        r2_vav /= Mcheck;
        r2_vav_maf01 /= Mmaf01;
        r2_vav_maf001 /= Mmaf001;
        r2_vav_maf0001 /= Mmaf0001;
        r2_vav_maf00001 /= Mmaf00001;

        cout << "  Genotype errors (hard):" << endl;
        cout << "    Total genotype errors:                    " << totalGtErrors << endl;
        cout << "    Total genotype errors (MAF>=0.1):         " << totalGtErrors_maf01 << endl;
        cout << "    Total genotype errors (MAF>=0.01):        " << totalGtErrors_maf001 << endl;
        cout << "    Total genotype errors (MAF>=0.001):       " << totalGtErrors_maf0001 << endl;
        cout << "    Total genotype errors (MAF>=0.0001):      " << totalGtErrors_maf00001 << endl;
        cout << "    Minimum genotype errors:                  " << gtErrorMin << endl;
        cout << "    Maximum genotype errors:                  " << gtErrorMax << endl;
        cout << "    Average gt err per sample:                " << avgterr << endl;
        cout << "    Minimum gt error rate:                    " << mingterrrate << endl;
        cout << "    Maximum gt error rate:                    " << maxgterrrate << endl;
        cout << "    Average gt error rate:                    " << avgterrrate << endl;
        cout << "    Average gt error rate (MAF>=0.1):         " << avgterrrate_maf01 << endl;
        cout << "    Average gt error rate (MAF>=0.01):        " << avgterrrate_maf001 << endl;
        cout << "    Average gt error rate (MAF>=0.001):       " << avgterrrate_maf0001 << endl;
        cout << "    Average gt error rate (MAF>=0.0001):      " << avgterrrate_maf00001 << endl;
        cout << "    Standard GER deviation:                   " << gerdev << endl;
        cout << "    correlation r2 (complete):                " << r2 << endl;
        cout << "    correlation r2 (complete, MAF>=0.1):      " << r2_maf01 << endl;
        cout << "    correlation r2 (complete, MAF>=0.01):     " << r2_maf001 << endl;
        cout << "    correlation r2 (complete, MAF>=0.001):    " << r2_maf0001 << endl;
        cout << "    correlation r2 (complete, MAF>=0.0001):   " << r2_maf00001 << endl;
        cout << "    correlation r2 (sample av.):              " << r2_sav << endl;
        cout << "    correlation r2 (sample av., MAF>=0.1):    " << r2_sav_maf01 << endl;
        cout << "    correlation r2 (sample av., MAF>=0.01):   " << r2_sav_maf001 << endl;
        cout << "    correlation r2 (sample av., MAF>=0.001):  " << r2_sav_maf0001 << endl;
        cout << "    correlation r2 (sample av., MAF>=0.0001): " << r2_sav_maf00001 << endl;
        cout << "    min correlation r2 (sample):              " << r2_smin << endl;
        cout << "    max correlation r2 (sample):              " << r2_smax << endl;
        cout << "    correlation r2 (var av.):                 " << r2_vav << endl;
        cout << "    correlation r2 (var av., MAF>=0.1):       " << r2_vav_maf01 << endl;
        cout << "    correlation r2 (var av., MAF>=0.01):      " << r2_vav_maf001 << endl;
        cout << "    correlation r2 (var av., MAF>=0.001):     " << r2_vav_maf0001 << endl;
        cout << "    correlation r2 (var av., MAF>=0.0001):    " << r2_vav_maf00001 << endl;
        cout << endl;
    }

    if (havedosages) {
        double totalGtErrorsSoft = 0;
        double totalGtErrorsSoft_maf01 = 0;
        double totalGtErrorsSoft_maf001 = 0;
        double totalGtErrorsSoft_maf0001 = 0;
        double totalGtErrorsSoft_maf00001 = 0;
        double gtErrorMinSoft = 0xffffffffffffffffull;
        double gtErrorMaxSoft = 0;
        // calc total errors and identify min and max
        for (double err : gtErrorsSoft) {
            totalGtErrorsSoft += err;
            if (gtErrorMinSoft > err) {
                gtErrorMinSoft = err;
            }
            if (gtErrorMaxSoft < err) {
                gtErrorMaxSoft = err;
            }
        }
        for (double err : gtErrorsSoft_maf01)
            totalGtErrorsSoft_maf01 += err;
        for (double err : gtErrorsSoft_maf001)
            totalGtErrorsSoft_maf001 += err;
        for (double err : gtErrorsSoft_maf0001)
            totalGtErrorsSoft_maf0001 += err;
        for (double err : gtErrorsSoft_maf00001)
            totalGtErrorsSoft_maf00001 += err;
        double avgterr = totalGtErrorsSoft / (double) Nquery;
        double mingterrrate = gtErrorMinSoft / (double)Mcheck;
        double maxgterrrate = gtErrorMaxSoft / (double)Mcheck;
        double avgterrrate  = totalGtErrorsSoft / (double)(Mcheck * Nquery);
        double avgterrrate_maf01 = totalGtErrorsSoft_maf01 / (double)(Mmaf01 * Nquery);
        double avgterrrate_maf001 = totalGtErrorsSoft_maf001 / (double)(Mmaf001 * Nquery);
        double avgterrrate_maf0001 = totalGtErrorsSoft_maf0001 / (double)(Mmaf0001 * Nquery);
        double avgterrrate_maf00001 = totalGtErrorsSoft_maf00001 / (double)(Mmaf00001 * Nquery);
        // calc standard deviation and variance for error rates
        double totgtdev2 = 0.0;
        for (double err : gtErrorsSoft) {
            double diff = err/(double)Mcheck - avgterrrate;
            totgtdev2 += diff*diff;
        }
        double gervar = totgtdev2 / (double)Nquery;
        double gerdev = sqrt(gervar);

        // correlation r2 (complete)
        double r2 = calc_r2_soft(Mcheck*Nquery, gtsumref, gt2sumref, gtdossumq, gtdos2sumq, gtdossumrefq);
        double r2_maf01 = calc_r2_soft(Mmaf01*Nquery, gtsumref_maf01, gt2sumref_maf01, gtdossumq_maf01, gtdos2sumq_maf01, gtdossumrefq_maf01);
        double r2_maf001 = calc_r2_soft(Mmaf001*Nquery, gtsumref_maf001, gt2sumref_maf001, gtdossumq_maf001, gtdos2sumq_maf001, gtdossumrefq_maf001);
        double r2_maf0001 = calc_r2_soft(Mmaf0001*Nquery, gtsumref_maf0001, gt2sumref_maf0001, gtdossumq_maf0001, gtdos2sumq_maf0001, gtdossumrefq_maf0001);
        double r2_maf00001 = calc_r2_soft(Mmaf00001*Nquery, gtsumref_maf00001, gt2sumref_maf00001, gtdossumq_maf00001, gtdos2sumq_maf00001, gtdossumrefq_maf00001);

        // correlation r2 (sample-wise)
        double r2_sav = 0.0;
        double r2_sav_maf01 = 0.0;
        double r2_sav_maf001 = 0.0;
        double r2_sav_maf0001 = 0.0;
        double r2_sav_maf00001 = 0.0;
        double r2_smin = 999.0;
        double r2_smax = 0.0;
        // store values for each sample (not required yet, but can be used later)
        vector<double> r2v(Nquery);
        vector<double> r2v_maf01(Nquery);
        vector<double> r2v_maf001(Nquery);
        vector<double> r2v_maf0001(Nquery);
        vector<double> r2v_maf00001(Nquery);
        for (size_t q = 0; q < Nquery; q++) {
            r2v[q] = calc_r2_soft(Mcheck, gtSumRef[q], gt2SumRef[q], gtdossumqv[q], gtdos2sumqv[q], gtdossumrefqv[q]);
            r2v_maf01[q] = calc_r2_soft(Mmaf01, gtSumRef_maf01[q], gt2SumRef_maf01[q], gtdossumqv_maf01[q], gtdos2sumqv_maf01[q], gtdossumrefqv_maf01[q]);
            r2v_maf001[q] = calc_r2_soft(Mmaf001, gtSumRef_maf001[q], gt2SumRef_maf001[q], gtdossumqv_maf001[q], gtdos2sumqv_maf001[q], gtdossumrefqv_maf001[q]);
            r2v_maf0001[q] = calc_r2_soft(Mmaf0001, gtSumRef_maf0001[q], gt2SumRef_maf0001[q], gtdossumqv_maf0001[q], gtdos2sumqv_maf0001[q], gtdossumrefqv_maf0001[q]);
            r2v_maf00001[q] = calc_r2_soft(Mmaf00001, gtSumRef_maf00001[q], gt2SumRef_maf00001[q], gtdossumqv_maf00001[q], gtdos2sumqv_maf00001[q], gtdossumrefqv_maf00001[q]);
            r2_sav += r2v[q]; // for the average
            r2_sav_maf01    += r2v_maf01[q]; // for the average
            r2_sav_maf001   += r2v_maf001[q]; // for the average
            r2_sav_maf0001  += r2v_maf0001[q]; // for the average
            r2_sav_maf00001 += r2v_maf00001[q]; // for the average
            if (r2v[q] < r2_smin)
                r2_smin = r2v[q];
            if (r2v[q] > r2_smax)
                r2_smax = r2v[q];
        }
        r2_sav /= Nquery;
        r2_sav_maf01 /= Nquery;
        r2_sav_maf001 /= Nquery;
        r2_sav_maf0001 /= Nquery;
        r2_sav_maf00001 /= Nquery;

        // correlation r2 (variant-wise)
        double r2_vav = 0.0;
        double r2_vav_maf01 = 0.0;
        double r2_vav_maf001 = 0.0;
        double r2_vav_maf0001 = 0.0;
        double r2_vav_maf00001 = 0.0;
        for (size_t bin = 0; bin < r2SoftSum.size(); bin++) {
            r2_vav += r2SoftSum[bin];
            r2_vav_maf01 += r2SoftSum_maf01[bin];
            r2_vav_maf001 += r2SoftSum_maf001[bin];
            r2_vav_maf0001 += r2SoftSum_maf0001[bin];
            r2_vav_maf00001 += r2SoftSum_maf00001[bin];
        }
        r2_vav /= Mcheck;
        r2_vav_maf01 /= Mmaf01;
        r2_vav_maf001 /= Mmaf001;
        r2_vav_maf0001 /= Mmaf0001;
        r2_vav_maf00001 /= Mmaf00001;

        cout << "  Genotype errors (soft):" << endl;
        cout << "    Total genotype errors (soft):                    " << totalGtErrorsSoft << endl;
        cout << "    Total genotype errors (MAF>=0.1) (soft):         " << totalGtErrorsSoft_maf01 << endl;
        cout << "    Total genotype errors (MAF>=0.01) (soft):        " << totalGtErrorsSoft_maf001 << endl;
        cout << "    Total genotype errors (MAF>=0.001) (soft):       " << totalGtErrorsSoft_maf0001 << endl;
        cout << "    Total genotype errors (MAF>=0.0001) (soft):      " << totalGtErrorsSoft_maf00001 << endl;
        cout << "    Minimum genotype errors (soft):                  " << gtErrorMinSoft << endl;
        cout << "    Maximum genotype errors (soft):                  " << gtErrorMaxSoft << endl;
        cout << "    Average gt err per sample (soft):                " << avgterr << endl;
        cout << "    Minimum gt error rate (soft):                    " << mingterrrate << endl;
        cout << "    Maximum gt error rate (soft):                    " << maxgterrrate << endl;
        cout << "    Average gt error rate (soft):                    " << avgterrrate << endl;
        cout << "    Average gt error rate (MAF>=0.1) (soft):         " << avgterrrate_maf01 << endl;
        cout << "    Average gt error rate (MAF>=0.01) (soft):        " << avgterrrate_maf001 << endl;
        cout << "    Average gt error rate (MAF>=0.001) (soft):       " << avgterrrate_maf0001 << endl;
        cout << "    Average gt error rate (MAF>=0.0001) (soft):      " << avgterrrate_maf00001 << endl;
        cout << "    Standard GER deviation (soft):                   " << gerdev << endl;
        cout << "    correlation r2 (complete) (soft):                " << r2 << endl;
        cout << "    correlation r2 (complete, MAF>=0.1) (soft):      " << r2_maf01 << endl;
        cout << "    correlation r2 (complete, MAF>=0.01) (soft):     " << r2_maf001 << endl;
        cout << "    correlation r2 (complete, MAF>=0.001) (soft):    " << r2_maf0001 << endl;
        cout << "    correlation r2 (complete, MAF>=0.0001) (soft):   " << r2_maf00001 << endl;
        cout << "    correlation r2 (sample av.) (soft):              " << r2_sav << endl;
        cout << "    correlation r2 (sample av., MAF>=0.1) (soft):    " << r2_sav_maf01 << endl;
        cout << "    correlation r2 (sample av., MAF>=0.01) (soft):   " << r2_sav_maf001 << endl;
        cout << "    correlation r2 (sample av., MAF>=0.001) (soft):  " << r2_sav_maf0001 << endl;
        cout << "    correlation r2 (sample av., MAF>=0.0001) (soft): " << r2_sav_maf00001 << endl;
        cout << "    min correlation r2 (sample) (soft):              " << r2_smin << endl;
        cout << "    max correlation r2 (sample) (soft):              " << r2_smax << endl;
        cout << "    correlation r2 (var av.) (soft):                 " << r2_vav << endl;
        cout << "    correlation r2 (var av., MAF>=0.1) (soft):       " << r2_vav_maf01 << endl;
        cout << "    correlation r2 (var av., MAF>=0.01) (soft):      " << r2_vav_maf001 << endl;
        cout << "    correlation r2 (var av., MAF>=0.001) (soft):     " << r2_vav_maf0001 << endl;
        cout << "    correlation r2 (var av., MAF>=0.0001) (soft):    " << r2_vav_maf00001 << endl;
        cout << endl;
    }

    {
        size_t totalSwErrors = 0;
        size_t swErrorMin = 0xffffffffffffffffull;
        size_t swErrorMax = 0;
        for (size_t err : switchErrors) {
            totalSwErrors += err;
            if (swErrorMin > err) {
                swErrorMin = err;
            }
            if (swErrorMax < err) {
                swErrorMax = err;
            }
        }
        double avswerr = totalSwErrors/(double) Nquery;
        double minswerrrate = swErrorMin / (double)Mcheck;
        double maxswerrrate = swErrorMax / (double)Mcheck;
        double avswerrrate  = totalSwErrors / (double)(Mcheck * Nquery);
        // calc standard deviation and variance for error rates
        double totswdev2 = 0.0;
        for (size_t err : switchErrors) {
            double diff = err/(double)Mcheck - avswerrrate;
            totswdev2 += diff*diff;
        }
        double servar = totswdev2 / (double)Nquery;
        double serdev = sqrt(servar);

        // tab-delimited list of queryID, mat/pat switched?, # swerrs, comma-separated list of sw error positions
        if (dump)
            cout << "Switch error positions (0-based) to cerr..." << endl;
        size_t errfree = 0;
        size_t matpatswitches = 0;
        for (size_t q = 0; q < Nquery; q++) {
            // switch error free targets
            if (switchErrors[q] == 0)
                errfree++;
            bool matpatsw = matPatSwitches[q];
            if (matpatsw) { // mat/pat switch
                matpatswitches++;
            }
            if (dump) {
                const auto &ep = errPos[q];
                cerr << q << "\t" << (matpatsw ? 1 : 0) << "\t" << (matpatsw ? (ep.size()-1) : ep.size()) << "\t";
                auto epit = ep.begin();
                if (matpatsw)
                    epit++; // jump over the zero encoding the matpat switch
                for (; epit != ep.end(); epit++)
                    cerr << *epit << ",";
                cerr << endl;
            }
        }
        if (dump)
            cout << " done.\n" << endl;

        cout << "  Switch errors:" << endl;
        cout << "    Total switch errors:       " << totalSwErrors << endl;
        cout << "    Minimum switch errors:     " << swErrorMin << endl;
        cout << "    Maximum switch errors:     " << swErrorMax << endl;
        cout << "    Average sw err per sample: " << avswerr << endl;
        cout << "    Minimum sw error rate:     " << minswerrrate << endl;
        cout << "    Maximum sw error rate:     " << maxswerrrate << endl;
        cout << "    Average sw error rate:     " << avswerrrate << endl;
        cout << "    Standard SER deviation:    " << serdev << endl;
        cout << "    Switch error free targets: " << errfree << endl;
        cout << "    Mat/Pat switches:          " << matpatswitches << endl;
        cout << endl;
    }

    if (havetyped) {
        size_t totalSwErrors = 0;
        size_t swErrorMin = 0xffffffffffffffffull;
        size_t swErrorMax = 0;
        for (size_t err : switchErrors_typed) {
            totalSwErrors += err;
            if (swErrorMin > err) {
                swErrorMin = err;
            }
            if (swErrorMax < err) {
                swErrorMax = err;
            }
        }
        double avswerr = totalSwErrors/(double) Nquery;
        double minswerrrate = swErrorMin / (double)Mtyped;
        double maxswerrrate = swErrorMax / (double)Mtyped;
        double avswerrrate  = totalSwErrors / (double)(Mtyped * Nquery);
        // calc standard deviation and variance for error rates
        double totswdev2 = 0.0;
        for (size_t err : switchErrors_typed) {
            double diff = err/(double)Mtyped - avswerrrate;
            totswdev2 += diff*diff;
        }
        double servar = totswdev2 / (double)Nquery;
        double serdev = sqrt(servar);

        size_t errfree = 0;
        size_t matpatswitches = 0;
        for (size_t q = 0; q < Nquery; q++) {
            // switch error free targets
            if (switchErrors_typed[q] == 0)
                errfree++;
            bool matpatsw = matPatSwitches_typed[q];
            if (matpatsw) { // mat/pat switch
                matpatswitches++;
            }
        }

        cout << "  Switch errors (typed):" << endl;
        cout << "    Total switch errors (typed):       " << totalSwErrors << endl;
        cout << "    Minimum switch errors (typed):     " << swErrorMin << endl;
        cout << "    Maximum switch errors (typed):     " << swErrorMax << endl;
        cout << "    Average sw err per sample (typed): " << avswerr << endl;
        cout << "    Minimum sw error rate (typed):     " << minswerrrate << endl;
        cout << "    Maximum sw error rate (typed):     " << maxswerrrate << endl;
        cout << "    Average sw error rate (typed):     " << avswerrrate << endl;
        cout << "    Standard SER deviation (typed):    " << serdev << endl;
        cout << "    Switch error free targets (typed): " << errfree << endl;
        cout << "    Mat/Pat switches (typed):          " << matpatswitches << endl;
        cout << endl;
    } // END if (typed)

    // dump sample-wise information to file
    string samplefile(queryfile);
    samplefile += ".checkphase.samples";
    ofstream ofs(samplefile);
    ofs << "Checked variants:\t" << Mcheck << endl;
    if (havetyped)
        ofs << "Checked typed variants:\t" << Mtyped << endl;
    // header
    ofs << "Idx\tSwerr\tGterr_hard\tGterr_soft" << endl;
    for (size_t q = 0; q < Nquery; q++) {
        ofs << q << "\t" << switchErrors[q] << "\t" << gtErrors[q] << "\t" << gtErrorsSoft[q] << endl;
    }
    ofs.close();

	return EXIT_SUCCESS;
}
