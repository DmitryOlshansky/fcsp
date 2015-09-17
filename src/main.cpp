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

void processFile(FCSP& fcsp, const string& path, bool plot)
{
    cerr << "Reading " << path << endl;
    ifstream f(path);
    if(!f){
        cerr << "ERROR: cannot open '" << path << "'\n";
        return;
    }
    try{
        fcsp.load(f);
        if (plot){
            ofstream dot(path+".dot");
            fcsp.dumpGraph(dot);
        }
        fcsp.process(cout);
    }
    catch(std::exception &e)
    {
        cerr << e.what() << endl;
    }
}

int main(int argc, const char* argv[])
{
	bool plot = false;
	string descriptors;
	vector<string> inputs;
	po::options_description visible("Options");
	visible.add_options()
		("help,h", po::bool_switch(), "Show help message.")
		("plot,p", po::bool_switch(), "Dump Graph-viz dot of molecule.")
		("descriptors,d", po::value<string>(), "Directory with descriptor database files.")
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
		if (vars.count("inputs"))
		{
            inputs = vars["inputs"].as<vector<string>>();	
		}
		plot = vars.count("plot") ? vars["plot"].as<bool>() : false;
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
        if(inputs.empty()){
            fcsp.load(cin);
            fcsp.process(cout);
        }
		else
            for(auto& arg : inputs) processFile(fcsp, arg, plot);
	}
	catch(std::exception& e)
	{
		cerr << e.what() << endl;
	}
	return 0;
}