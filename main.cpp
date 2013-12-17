/**
 *
 */
#include <stdexcept>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "boost/iostreams/stream.hpp"
#include "boost/iostreams/device/null.hpp"
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "fcsp.hpp"

using namespace std;
namespace bio = boost::iostreams;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

bio::stream<bio::null_sink> nullOstream((bio::null_sink()));

int main(int argc, const char* argv[])
{
	bool plot = false, linear = false, help = false;
	string descriptors=".";
	vector<string> inputs;
	po::options_description visible("Options");
	visible.add_options()
		("help,h", po::bool_switch(), "Show help message.")
		("plot,p", po::bool_switch(), "Dump Graph-viz dot of molecule.")
		("linear,l", po::bool_switch(), "Dump linear FCSP descriptors.")
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
		inputs = vars["inputs"].as<vector<string>>();
		help = vars["help"].as<bool>();
		plot = vars["plot"].as<bool>();
		linear = vars["linear"].as<bool>();
		descriptors = vars["descriptors"].as<string>();
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
    if(!inputs.size() || help)
    {
    	cout << visible << endl;
		return 1;
    }
    fs::path base(descriptors);
    ifstream descr1((base / "descr1.csv").c_str());
    ifstream descr2((base / "descr2.sdf").c_str());
	try{
		auto order1 = read1stOrder(descr1);
		auto order2 = read2ndOrder(descr2);
		FCSP fcsp(move(order1), move(order2));
		for(auto& arg : inputs)
		{
			ifstream f(arg);
			if(!f){
				cerr << "ERROR: cannot open '" << arg << "'\n";
				continue;
			}
			try{
				fcsp.load(f);
				if(plot){
					ofstream dot(arg+".dot");
					fcsp.dumpGraph(dot);
				}
				if(linear)
					fcsp.linear(cout);
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
