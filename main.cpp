/**
 *
 */
#include <stdexcept>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include "boost/iostreams/stream.hpp"
#include "boost/iostreams/device/null.hpp"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "fcsp.hpp"
#include "ctab.h"

using namespace std;
using boost::lexical_cast;
namespace bio = boost::iostreams;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

bio::stream<bio::null_sink> nullOstream((bio::null_sink()));

void readReplacements(ifstream& inp, vector<Replacement>& repls)
{
	auto sdfs = readSdf(inp);
	for_each(sdfs.begin(), sdfs.end(), [&repls](SDF& sdf){
		int dc = lexical_cast<int>(sdf.props["DC"][0]);
		int couple = lexical_cast<int>(sdf.props["COUPLING"][0]);
		repls.emplace_back(toGraph(sdf.mol), dc, couple);
		/*ofstream out(lexical_cast<string>(repls.size())+".dot");
		dumpGraph(repls.back().piece, out);*/
	});
}

int main(int argc, const char* argv[])
{
	bool plot = false, linear = false, cyclic = false, replace = false;
	string descriptors;
	vector<string> inputs;
	po::options_description visible("Options");
	visible.add_options()
		("cyclic,c", po::bool_switch(), "Dump cyclic FCSP descriptors.")
		("help,h", po::bool_switch(), "Show help message.")
		("plot,p", po::bool_switch(), "Dump Graph-viz dot of molecule.")
		("linear,l", po::bool_switch(), "Dump linear FCSP descriptors.")
		("replacement,r", po::bool_switch(), "Dump replacements FCSP descriptors.")
		("descriptors,d", po::value<string>(), "Directory with descriptor lists.")
		;
	po::options_description implicit("Implicit options");
	implicit.add_options()
		("inputs", po::value<vector<string>>())
		;
	po::options_description whole;
	whole.add(visible).add(implicit);
	po::positional_options_description posOpt;
	posOpt.add("inputs", -1);
	po::variables_map vars;
	try{
		po::store(po::command_line_parser(argc, argv)
			.options(whole).positional(posOpt).run(), vars);
		if (!vars.count("inputs") /*|| vars.count("help")*/)
		{
			cout << visible << endl;
			return 1;
		}
		inputs = vars["inputs"].as<vector<string>>();
		plot = vars.count("plot") ? vars["plot"].as<bool>() : false;
		linear = vars.count("linear") ? vars["linear"].as<bool>() : false;
		replace = vars.count("replacement") ? vars["replacement"].as<bool>() : false;
		descriptors = vars.count("descriptors") ? vars["descriptors"].as<string>() : ".";
	}
	catch(const po::invalid_command_line_syntax &e) {
        switch (e.kind()) {
        case po::invalid_syntax::missing_parameter:
            cout << "Missing argument for option '" << e.tokens() << "'.\n";
            break;
        default:
            cout << "Syntax error: " << (int)e.kind() << "\n";
            break;
        };
        return 1;
    }
    catch (const po::unknown_option &e) {
        cout << "Unknown option '" << e.get_option_name() << "'\n";
        return 1;
    }
    
    fs::path base(descriptors);
    ifstream descr1((base / "descr1.csv").c_str());
    if(!descr1){
    	cerr << "Failed to open descr1.csv, check your -d option" << endl;
    	return 1;
    }
    ifstream descr2((base / "descr2.sdf").c_str());
    if(!descr2){
		cerr << "Failed to open descr2.sdf, check -d option" << endl;
		return 1;
	}
	ifstream repl((base / "replacement.sdf").c_str());
	try{
		vector<LevelOne> order1;
		vector<LevelTwo> order2;
		vector<Replacement> replacements;
		read1stOrder(descr1, order1);
		read2ndOrder(descr2, order2);
		readReplacements(repl, replacements);
		FCSP fcsp(move(order1), move(order2), move(replacements));
		for(auto& arg : inputs)
		{
			ifstream f(arg);
			if(!f){
				cerr << "ERROR: cannot open '" << arg << "'\n";
				continue;
			}
			try{
				fcsp.load(f);
				if (plot){
					ofstream dot(arg+".dot");
					fcsp.dumpGraph(dot);
				}
				if (linear)
					fcsp.linear(cout);
				if (replace)
					fcsp.replacement(cout);
			}
			catch(std::exception &e)
			{
				cerr << e.what() << endl;
			}
		}
	}
	catch(std::exception& e)
	{
		cerr << e.what() << endl;
	}
	return 0;
}
