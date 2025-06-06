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

#ifndef ARGS_H_
#define ARGS_H_

#include <string>
#include <vector>

#include <boost/program_options.hpp>

namespace bpo = boost::program_options;

using namespace std;

/**
 * Class for storing and retrieving command-line arguments.
 */
class Args {
public:

	static Args parseArgs(int argc, char *argv[]);

    /**
     * Returns the argument with the given (long) name. The template argument
     * specifies the return type you want the parameter casted to.
     * @param name argument name
     * @return the variable value
     */
    template<typename T> T get(const std::string &name) const {
        auto where = vars.find(name);
        if(where == std::end(vars)) {
            if(!isDefined(name))
                throw std::invalid_argument("Option undefined: " + name + " (This is a bug)");
            else
                throw std::out_of_range("Option has not been specified and does not have a default value associated: " + name);
        }
        return where->second.as<T>();
    }

    /**
     * Counts the number of argument occurences. Mainly useful for boolean switches
     * and/or counting flags, such as verbosity or debug mode.
     * @param name argument name
     * @return argument value
     */
    unsigned int count(const std::string &name) const {
        if(!isDefined(name))
            throw std::invalid_argument("Option undefined: " + name + " (This is a bug)");
        return vars.count(name);
    }

    bool operator()(const std::string &name) const {
        return count(name) > 0;
    }

    /**
     * Prints a help message.
     * @param progname the program name, usually argv[0]
     * @param out output stream, e.g. cout
     */
    void printHelp(const std::string &progname, std::ostream &out) const;
    static void printVersion();

    Args(Args&& other) = default;

    string vcfRef;
    string vcfQuery;
    string vcfShared;
    string statfile;

    bool debug = false;
    bool dump = false;

protected:
    /** Constructs the arguments list and adds all defined options */
    Args();
    Args(Args const &);
    void operator=(Args const &);

    void parse(int argc, char *argv[]);
    void parseVars();
    bool isDefined(const std::string &optname) const;

    bpo::options_description opts_regular;        /**< regular options, shown on help output */
    bpo::options_description opts_hidden;         /**< hidden options */

    bpo::variables_map vars;    /**< this is where the option values go */

    /** dump all options to the given stream */
    friend std::ostream &operator<<(std::ostream &out, const Args &args);

private:
    /** parses the main() options into the variable map */
    Args(int argc, char *argv[]);

};

#endif /* ARGS_H_ */
