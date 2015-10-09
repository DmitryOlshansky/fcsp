/**
 *
 */
#include <stdexcept>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "fcsp.hpp"
#include "ctab.hpp"
#include "log.hpp"

using namespace std;

#ifdef __linux__
#include <unistd.h>

// array of default search paths for descriptor database
// on linux also try ../share/fcss-2a/descr
vector<string> descrPaths()
{
	char buf[2048];
	int len = readlink("/proc/self/exe", buf, sizeof(buf));
	auto defaults = vector<string>({"."});
	if(len)
	{
		auto exe = string(buf, len);
		auto p = exe.find_last_of("/");
		auto share = exe.substr(0, p) + "/../share/fcss-2a/descr";
		defaults.insert(defaults.begin(), share);
	}
	return defaults;
}

#else

// just deafult to current directory
vector<string> descrPaths()
{
	return vector<string>({"."});
}

#endif

namespace po = boost::program_options;
namespace fs = boost::filesystem;

void processFile(FCSP& fcsp, const string& path, bool plot)
{
    LOG(INFO) << "Reading " << path << endline;
    fs::path fpath(path);
    ifstream f(path);
    if(!f){
        LOG(ERROR) << "ERROR: cannot open '" << path << "'\n";
        return;
    }
    try{
        fcsp.load(f);
        if (plot){
            ofstream dot(path+".dot");
            fcsp.dumpGraph(dot);
        }
        fcsp.process(cout, fpath.filename().string());
    }
    catch(std::exception &e)
    {
        LOG(ERROR) << e.what() << endline;
    }
}

FCSPFMT toFCSPFMT(string fmt)
{
	if(fmt == "json") return FCSPFMT::JSON;
	if(fmt == "csv") return FCSPFMT::CSV;
	if(fmt == "txt") return FCSPFMT::TXT;
	throw logic_error("No such format "+fmt);
}

FCSPOptions configure(vector<string> paths, bool long41, FCSPFMT fmt)
{
	auto default_ex = logic_error("DB not found in any of search paths");
	string found = "";
	exception& ex = default_ex;
	vector<LevelOne> order1;
	vector<LevelTwo> order2;
	vector<Replacement> replacements;
  for(auto p : paths)
	{
		LOG(DEBUG) << "Trying to load DCs from " << p << endline;
  	try{
  		auto descrBase = fs::path(p);
  		ifstream descr1((descrBase / "descr1.csv").c_str());
		  if(!descr1)
		  	throw logic_error("Failed to open descr1.csv DB, check your -d option");
		  ifstream descr2((descrBase / "descr2.sdf").c_str());
		  if(!descr2)
				throw logic_error("Failed to open descr2.sdf DB, check -d option");
			ifstream repl((descrBase / "replacement.sdf").c_str());
			if(!repl)
				throw logic_error("Failed to open replacement.sdf DB, check -d option");
			order1 = read1stOrder(descr1);
			order2 = read2ndOrder(descr2);
			replacements = readReplacements(repl);
			found = p;
			break;
		}
		catch(exception &te){
			ex = te;
		}
	}
	if(found.empty())
		throw ex;
	else
		LOG(INFO) << "Loaded DCs from "<< found << endline;
	return {order1, order2, replacements, long41, fmt};
}

int main(int argc, const char* argv[])
{
	bool plot = false;
	bool long41 = true;
	string descriptors;
	FCSPFMT fmt = FCSPFMT::JSON;
	vector<string> inputs;
	po::options_description visible("Options");
	visible.add_options()
		("help,h", po::bool_switch(), "Show help message.")
		("long41", po::bool_switch(), "Enable old FCSS-2 treatment of DC #41 - adding +1 to the chain")
		("plot,p", po::bool_switch(), "Dump Graph-viz dot of molecule.")
		("verbosity,v", po::value<int>(), "Set verbosity level (0-6).")
		("format", po::value<string>(), "Choose output format: json - JSON with bindings (default), csv - CSV with just codes, txt - simple line-based format")
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
		if (vars.count("format"))
		{
			fmt = toFCSPFMT(vars["format"].as<string>());
		}
		if (vars.count("verbosity"))
		{
			logLevel = vars["verbosity"].as<int>();
		}
		plot = vars.count("plot") ? vars["plot"].as<bool>() : false;
		descriptors = vars.count("descriptors") ? vars["descriptors"].as<string>() : ".";
	}
	catch(const po::invalid_command_line_syntax &e){
    switch (e.kind()) {
      case po::invalid_syntax::missing_parameter:
          cerr << "Missing argument for option '" << e.tokens() << "'.\n";
          break;
      default:
          cerr << "Syntax error: " << (int)e.kind() << "\n";
          break;
      };
      return 1;
  }
  catch (const po::unknown_option &e) {
      cerr << "Unknown option '" << e.get_option_name() << "'\n";
      return 1;
  }
	try{
		auto paths = descrPaths();
		if(descriptors != ".")
			paths.insert(paths.begin(), descriptors);
		FCSP fcsp(configure(paths, long41, fmt));
    if(inputs.empty())
    {
        fcsp.load(cin);
        fcsp.process(cout);
    }
		else
		{
      for(auto& arg : inputs)
    		processFile(fcsp, arg, plot);
		}
	}
	catch(std::exception& e)
	{
		cerr << e.what() << endline;
	}
	return 0;
}
