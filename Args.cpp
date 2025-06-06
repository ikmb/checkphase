/*
 *    Copyright (C) 2025 by Lars Wienbrandt,
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

#include <iostream>
#include <iomanip>
#include <cmath>
#include <thread>

#include <boost/program_options.hpp>

#include "Args.h"
#include <version.h>

namespace bpo = boost::program_options;

using namespace std;
using namespace bpo;

/**
 * @brief Constructs, parses and verifies all command-line arguments
 *
 * This function constructs the Args object and parses all given command-line
 * arguments. If unknown arguments and/or missing obligatory arguments are
 * detected, this function does not return and instead prints an appropriate
 * help message and calls exit(), executing all exit handlers registered
 * up to the parseArgs call.
 *
 * @param argc argument count
 * @param argv argument vector
 * @return an rvalue reference of a new Args object containing all defined or defaulted CLI arguments
 */
/*static*/ Args Args::parseArgs(int argc, char *argv[]) {
    Args args {argc, argv};

    if (argc <= 1 || args.count("help")) {
        args.printHelp(argv[0], cout);
        exit(EXIT_SUCCESS);
    }

    if (args.count("version")) {
        args.printVersion();
        exit(EXIT_SUCCESS);
    }

    // set variables
    args.parseVars();

    if (!args.count("query")) {
        cerr << "ERROR: No query file specified. Try --help." << endl;
        exit(EXIT_FAILURE);
    }

    if (!args.count("ref")) {
        cerr << "ERROR: No reference file specified. Try --help." << endl;
        exit(EXIT_FAILURE);
    }

    return args;
}

ostream &operator<<(ostream &out, const Args &args) {
    variables_map::const_iterator it;

    long name_width = 0;
    long value_width = 0;
    long orig_width = out.width();

    // collect field widths
    for(it = args.vars.begin(); it != args.vars.end(); ++it) {
        long this_width = static_cast<long>(it->first.length());

        if(this_width > name_width) {
            name_width = this_width;
        }

        this_width = 0;

        if(it->second.value().type() == typeid(string)) {
            this_width = static_cast<long>(it->second.as<string>().length());
        }

        if(it->second.value().type() == typeid(int)) {
            this_width = static_cast<long>(log10(it->second.as<int>()));
        }

        if(this_width > value_width) {
            value_width = this_width;
        }
    }

    // dump option values
    out.setf(ios::adjustfield, ios::left);

    for(it = args.vars.begin(); it != args.vars.end(); ++it) {
        out.width(name_width + 2); // more space for the colon
        out << (it->first + ":") << endl;

        out.width(value_width);

        if(it->second.value().type() == typeid(string)) {
            out << it->second.as<string>() << endl;
        } else if(it->second.value().type() == typeid(int)) {
            out << it->second.as<int>() << endl;
        } else {
            out << "(unknown)" << endl;
        }
    }

    resetiosflags(ios::adjustfield);
    out.width(orig_width);

    return out;
}

Args::Args(int argc, char *argv[]) :
    opts_regular("Program options"),
	opts_hidden("Hidden options (only visible in debug mode)")
	{

    opts_regular.add_options()
    ("help,h", "Produces this help message and exits.")
    ("version,v", "Prints version information and exits.")
    ("query,q", value<string>(&vcfQuery), "Indexed compressed VCF/BCF file for query genotypes that should be checked against a reference.")
    ("ref,r", value<string>(&vcfRef), "Indexed compressed VCF/BCF file for the reference that the query is checked against. If RefPanelAF tag exists, MAF categories are enabled.")
    ("shared,s", value<string>(&vcfShared), "Tabix-indexed compressed VCF/BCF file containing variants that should exclusively be taken for the comparison. Does not need to contain genotypes (are ignored anyway). If RefPanelAF tag exists, MAF categories are enabled. (optional)")
    ("stat", value<string>(&statfile), "File for status output. (optional)")
    ;

    opts_hidden.add_options()
    ("dump", "Dumps switch error positions to stderr.") // TODO remove?
    ("debug", "Produce lots of debug output.")
    ;

    parse(argc, argv);
}

void Args::parse(int argc, char *argv[]) {
    bpo::options_description all_options;

    // combine all options
    all_options.add(opts_regular);
    all_options.add(opts_hidden);

    // do the actual parsing
    store(command_line_parser(argc, argv).options(all_options).run(), vars);
    notify(vars);

}

void Args::parseVars() {

    if (vars.count("debug"))
        debug = true;

    if (vars.count("dump"))
        dump = true;

    if (!vcfQuery.compare(vcfRef)) {
        cerr << "ERROR: Query and reference cannot be the same file." << endl;
        exit(EXIT_FAILURE);
    }

    if (!vcfShared.empty() && (!vcfShared.compare(vcfRef) || !vcfShared.compare(vcfQuery))) {
        cerr << "ERROR: Shared variants file cannot be equal to query or reference." << endl;
        exit(EXIT_FAILURE);
    }

    if (!statfile.empty()) {
        bool ok = true;
        size_t pos = statfile.rfind(".vcf.gz");
        if (pos != string::npos && pos == statfile.size()-7)
            ok = false;
        pos = statfile.rfind(".vcf");
        if (pos != string::npos && pos == statfile.size()-4)
            ok = false;
        pos = statfile.rfind(".bcf");
        if (pos != string::npos && pos == statfile.size()-4)
            ok = false;
        if(!ok) {
            cerr << "ERROR: Status file must not end with .vcf.gz/.vcf/.bcf" << endl;
            exit(EXIT_FAILURE);
        }
    }

}

bool Args::isDefined(const string &optname) const {
    bool found = false;
    found = !!this->opts_regular.find_nothrow(optname, false); // return null (-> false) if option has not been found
    found |= !!this->opts_hidden.find_nothrow(optname, false);
    return found;
}

void Args::printHelp(const string &progname, ostream &out) const {
    out << "Usage: " << progname << " [options]" << endl << endl;
    out << opts_regular << endl;
    if (debug) {
        out << opts_hidden << endl;
    }
    out << endl;

    printVersion();
}

/* static */
void Args::printVersion() {
    cout << "This is version " << prog_version << ", compiled on " << prog_timestamp << endl;
    cout << "Send bugs to " << prog_bugaddress << endl;
}
