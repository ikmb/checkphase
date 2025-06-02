/*
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

// convert NULL terminated C string to a reverse complemented std::string
// only charactes ATCG are supported, others are kept as they are
inline string reverseComplement(const char *s) {
    string rc(s);
    auto rcit = rc.begin();
    for (int i = rc.size()-1; i >= 0; i--) {
        char c = s[i];
        switch (c) {
        case 'A':
            c = 'T';
            break;
        case 'T':
            c = 'A';
            break;
        case 'C':
            c = 'G';
            break;
        case 'G':
            c = 'C';
            break;
        default:
            break;
        }
        *rcit = c;
        rcit++;
    }
    return rc;
}

inline void updateStatus(const char* statfile, float r, float q) {
    if (!statfile)
        return;
    ofstream stat(statfile, ofstream::trunc);
    if (!stat.fail())
        stat << "R: " << r << "\nQ: " << q << endl;
}

inline void printUsageAndDie(const char* argv0) {
    cerr << "ERROR: Please provide two input VCFs/BCFs.\nUsage: " << argv0 << " <reference> <query> [--stat <statfile>] [--dump]" << endl;
    exit(EXIT_FAILURE);
}

int main(int argc, char *argv[]) {

    if (argc < 3 || argc > 6) {
        printUsageAndDie(argv[0]);
    }

    char* reffile = argv[1];
    char* queryfile = argv[2];
    char* statfile = NULL;
    bool dump = false;

    if (argc == 5 && strcmp(argv[3], "--stat") == 0 && strcmp(argv[4], "--dump") != 0) // only --stat option
        statfile = argv[4];
    else if (argc == 4 && strcmp(argv[3], "--dump") == 0) // only ---dump option
        dump = true;
    else if (argc == 6) { // --dump and --stat option
        if (strcmp(argv[3], "--dump") == 0 && strcmp(argv[4], "--stat") == 0) {
            dump = true;
            statfile = argv[5];
        } else if (strcmp(argv[3], "--stat") == 0 && strcmp(argv[5], "--dump") == 0) {
            dump = true;
            statfile = argv[4];
        } else
            printUsageAndDie(argv[0]);
    } else if (argc > 3)
        printUsageAndDie(argv[0]);

    if (statfile)
        cout << "Statfile: " << statfile << endl;
    if (dump)
        cout << "--dump enabled. Will dump phase error positions to stderr." << endl;
    cout << endl;

    updateStatus(statfile, 0, 0);

    size_t Mrefidx = getNumVariantsFromIndex(reffile);
    size_t Mqidx = getNumVariantsFromIndex(queryfile);

    bcf_srs_t *sr = bcf_sr_init();

    bcf_sr_set_opt(sr, BCF_SR_PAIR_LOGIC, BCF_SR_PAIR_BOTH_REF);
    bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);

    if (!bcf_sr_add_reader(sr, reffile)) {
        cerr << "ERROR: Could not open reference for reading: " << bcf_sr_strerror(sr->errnum) << endl;
        exit(EXIT_FAILURE);
    }
    if (!bcf_sr_add_reader(sr, queryfile)) {
        cerr << "ERROR: Could not open query for reading: " << bcf_sr_strerror(sr->errnum) << endl;
        exit(EXIT_FAILURE);
    }

    bcf_hdr_t *ref_hdr = bcf_sr_get_header(sr, 0);
    bcf_hdr_t *q_hdr = bcf_sr_get_header(sr, 1);

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
    cout << "  Reference: " << reffile << endl;
    cout << "  Query:     " << queryfile << endl;
    cout << "  Reference variants (from index): Mref = " << Mrefidx << endl;
    cout << "  Reference samples:               Nref = " << Nref << endl;
    cout << "  Query variants (from index):   Mquery = " << Mqidx << endl;
    cout << "  Query samples:                 Nquery = " << Nquery << endl;
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

    size_t Mshared = 0;    // shared variants
    size_t Mref = 0; // reference variants
    size_t Mq = 0;   // query variants

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
    vector<bool> switched(Nquery, false);
    vector<size_t> switchErrors(Nquery, 0);
    vector<size_t> gtErrors(Nquery, 0);
    vector<double> gtErrorsSoft(Nquery, 0.0);
    vector<vector<size_t>> errPos(Nquery);

    // for counting haploid samples
    bool hapsset_ref = false;
    bool hapsset_tgt = false;
    size_t Nrefhap = 0;
    size_t Nqhap = 0;
    bool havedosages = false;

    // for progress
    size_t linesperpercentref = Mrefidx/100;
    size_t linesperpercentq = Mqidx/100;
    if (!linesperpercentref) linesperpercentref = 1;
    if (!linesperpercentq) linesperpercentq = 1;
    updateStatus(statfile, 0, 0);

    while (bcf_sr_next_line(sr)) { // read data SNP-wise in positional sorted order from target and reference

        bcf1_t *ref = bcf_sr_get_line(sr, 0); // read one line of reference, if available at current position (otherwise NULL)
        bcf1_t *tgt = bcf_sr_get_line(sr, 1); // read one line of target, if available at current position (otherwise NULL)

        if (ref)
            Mref++;
        if (tgt)
            Mq++;

        if (Mref % linesperpercentref == 0 || Mq % linesperpercentq == 0)
            updateStatus(statfile, Mref/(float)Mrefidx, Mq/(float)Mqidx);

        if (!ref || !tgt) { // query or reference only variant
            continue;
        }

        // shared variant
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

        // get dosages, if present
        int nret = bcf_get_format_values(q_hdr,tgt,"ADS",(void**)&dosbuf,&ndosbuf,BCF_HT_REAL); // ADS type
        if (!nret) nret = bcf_get_format_values(q_hdr,tgt,"HDS",(void**)&dosbuf,&ndosbuf,BCF_HT_REAL); // HDS type
        if (nret) {
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
            int refgt = (refmat ? 1 : 0) + (refpat ? 1 : 0);
            int qgt = (qmat ? 1 : 0) + (qpat ? 1 : 0);
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
            }
            // hard check
            if (refgt != qgt) { // genotype error -> continue with next sample
                gtErrors[q]++;
                continue;
            }

            // phase check
            if (refgt == 1 && phased) { // only proceed on heterozygous site and if site is phased
                if (initialized[q]) { // this is not the first het for this sample -> wrong phase is a switch error!
                    //bool swerr = (refmat != qmat) xor switched[q];
                    bool swerr = switched[q] ? (refmat == qmat) : (refmat != qmat);
                    if (swerr) {
                        switched[q] = !switched[q];
                        switchErrors[q]++;
                        errPos[q].push_back(Mq-1); // error position is 0-based
                    }
                } else { // first het site -> simply set the switched-flag according to the current phases in ref and query
                    initialized[q] = true;
                    if (refmat != qmat) {
                        switched[q] = true;
                        errPos[q].push_back(0); // indicates the difference at the first het, for debugging
                    }
                }
            }
        } // end for every sample

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
    cout << endl;

    cout << "  Excluded shared:         " << Mexclude << endl;
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

    {
        size_t totalGtErrors = 0;
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
        double avgterr = round(totalGtErrors/(double) Nquery);
        double mingterrrate = gtErrorMin / (double)Mcheck;
        double maxgterrrate = gtErrorMax / (double)Mcheck;
        double avgterrrate  = totalGtErrors / (double)(Mcheck * Nquery);
        // calc standard deviation and variance
        double totgtdev = 0.0;
        double totgtdev2 = 0.0;
        for (size_t err : gtErrors) {
            double diff = ((double)err) - avgterr;
            totgtdev += abs(diff);
            totgtdev2 += diff*diff;
        }
        double gerdev = totgtdev / (double)(Mcheck * Nquery);
        double gervar = sqrt(totgtdev2) / (double)(Mcheck * Nquery);

        cout << "  Genotype errors (hard):" << endl;
        cout << "    Total genotype errors:     " << totalGtErrors << endl;
        cout << "    Minimum genotype errors:   " << gtErrorMin << endl;
        cout << "    Maximum genotype errors:   " << gtErrorMax << endl;
        cout << "    Average gt err per sample: " << avgterr << endl;
        cout << "    Minimum gt error rate:     " << mingterrrate << endl;
        cout << "    Maximum gt error rate:     " << maxgterrrate << endl;
        cout << "    Average gt error rate:     " << avgterrrate << endl;
        cout << "    Standard GER deviation:    " << gerdev << endl;
        cout << "    GER variance:              " << gervar << endl;
        cout << endl;
    }

    if (havedosages) {
        double totalGtErrorsSoft = 0;
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
        double avgterr = round(totalGtErrorsSoft/(double) Nquery);
        double mingterrrate = gtErrorMinSoft / (double)Mcheck;
        double maxgterrrate = gtErrorMaxSoft / (double)Mcheck;
        double avgterrrate  = totalGtErrorsSoft / (double)(Mcheck * Nquery);
        // calc standard deviation and variance
        double totgtdev = 0.0;
        double totgtdev2 = 0.0;
        for (double err : gtErrorsSoft) {
            double diff = ((double)err) - avgterr;
            totgtdev += abs(diff);
            totgtdev2 += diff*diff;
        }
        double gerdev = totgtdev / (double)(Mcheck * Nquery);
        double gervar = sqrt(totgtdev2) / (double)(Mcheck * Nquery);

        cout << "  Genotype errors (soft):" << endl;
        cout << "    Total genotype errors (soft):     " << totalGtErrorsSoft << endl;
        cout << "    Minimum genotype errors (soft):   " << gtErrorMinSoft << endl;
        cout << "    Maximum genotype errors (soft):   " << gtErrorMaxSoft << endl;
        cout << "    Average gt err per sample (soft): " << avgterr << endl;
        cout << "    Minimum gt error rate (soft):     " << mingterrrate << endl;
        cout << "    Maximum gt error rate (soft):     " << maxgterrrate << endl;
        cout << "    Average gt error rate (soft):     " << avgterrrate << endl;
        cout << "    Standard GER deviation (soft):    " << gerdev << endl;
        cout << "    GER variance (soft):              " << gervar << endl;
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
        double avswerr = round(totalSwErrors/(double) Nquery);
        double minswerrrate = swErrorMin / (double)Mcheck;
        double maxswerrrate = swErrorMax / (double)Mcheck;
        double avswerrrate  = totalSwErrors / (double)(Mcheck * Nquery);
        // calc standard deviation and variance
        double totswdev = 0.0;
        double totswdev2 = 0.0;
        for (size_t err : switchErrors) {
            double diff = ((double)err) - avswerr;
            totswdev += abs(diff);
            totswdev2 += diff*diff;
        }
        double serdev = totswdev / (double)(Mcheck * Nquery);
        double servar = sqrt(totswdev2) / (double)(Mcheck * Nquery);

        // tab-delimited list of queryID, mat/pat switched?, # swerrs, comma-separated list of sw error positions
        if (dump)
            cout << "Switch error positions (0-based) to cerr..." << endl;
        size_t errfree = 0;
        size_t matpatswitches = 0;
        for (size_t q = 0; q < Nquery; q++) {
            const auto &ep = errPos[q];
            if (ep.empty() || (ep.size() == 1 && ep[0] == 0)) {
                errfree++;
            }
            bool matpatsw = false;
            if (!ep.empty() && ep[0] == 0) {
                matpatsw = true;
                matpatswitches++;
            }
            if (dump) {
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
        cout << "    SER variance:              " << servar << endl;
        cout << "    Switch error free targets: " << errfree << endl;
        cout << "    Mat/Pat switches:          " << matpatswitches << endl;
        cout << endl;
    }

    // dump sample-wise information to file
    string samplefile(queryfile);
    samplefile += ".checkphase.samples";
    ofstream ofs(samplefile);
    ofs << "Checked variants:\t" << Mcheck << endl;
    // header
    ofs << "Idx\tSwerr\tGterr_hard\tGterr_soft" << endl;
    for (size_t q = 0; q < Nquery; q++) {
        ofs << q << "\t" << switchErrors[q] << "\t" << gtErrors[q] << "\t" << gtErrorsSoft[q] << endl;
    }
    ofs.close();

	return EXIT_SUCCESS;
}
