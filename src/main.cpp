/**
 *
 */
#include <stdexcept>
#include <boost/filesystem.hpp>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <thread>
#include "cxxopts.hpp"
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

// just default to current directory
vector<string> descrPaths()
{
    return vector<string>({"."});
}

#endif

namespace fs = boost::filesystem;

void processFile(FCSP& fcsp, const string& path, ostream& out)
{
    LOG(INFO) << "Reading " << path << endline;
    fs::path fpath(path);
    ifstream f(path);
    if(!f){
        LOG(ERROR) << "ERROR: cannot open '" << path << "'\n";
        return;
    }
    try {
        fcsp.load(f);
        fcsp.process(out, fpath.filename().string());
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

int main(int argc, char* argv[])
{
    bool long41 = true;
    int threads = 0;
    string descriptors;
    FCSPFMT fmt = FCSPFMT::JSON;
    vector<string> inputs;
    cxxopts::Options options(argv[0], " - example command line options");
    options.add_options()
    ("h,help", "Print help")
    ("input", "List of MOL files to encode", cxxopts::value<vector<string>>())
    ("t,threads", "Number of threads to use", cxxopts::value<int>(), "0")
    ("v,verbosity", "Level of verbosity", cxxopts::value<int>(), "0")
    ("f,format", "Output format: txt, csv, json", cxxopts::value<string>(), "json");

    try {
        options.parse_positional("input");
        options.parse(argc, argv);
        if (options.count("help"))
        {
            cout << options.help({""}) << endl;
            return 0;
        }
        if (options.count("input"))
        {
            inputs = options["input"].as<vector<string>>();    
        }
        if (options.count("format"))
        {
            fmt = toFCSPFMT(options["format"].as<string>());
        }
        if (options.count("threads"))
        {
            threads = options["threads"].as<int>();
        }
        if (options.count("verbosity"))
        {
            logLevel = options["verbosity"].as<int>();
        }
        descriptors = options.count("descriptors") ? options["descriptors"].as<string>() : ".";
    }
    catch(cxxopts::OptionException &e) {
        cerr << "Argument parsing error: ";
        cerr << e.what() << endline;
        return 1;
    }
    catch(exception & e) {
        cerr << "Error during parameter parsing: ";
        cerr << e.what() << endline;
        return 1;
    }
    try {
        auto paths = descrPaths();
        if(descriptors != ".")
            paths.insert(paths.begin(), descriptors);
        auto conf = configure(paths, long41, fmt);
        
        if(inputs.empty()) {
            FCSP fcsp(conf);
            fcsp.load(cin);
            fcsp.process(cout);
        }
        else {
            size_t n = threads <= 0 ? thread::hardware_concurrency() : threads;
            size_t batch = (inputs.size() + n - 1) / n;
            size_t batches = (inputs.size()+ batch - 1)/ batch;
            LOG(INFO) << "CPUs: " << n << " batch-size: " << batch << endline;
            vector<thread> threads(batches);
            vector<stringstream> streams(batches);
            for (size_t i = 0; i < batches; i ++) {
                size_t start = i * batch;
                size_t stop = start + batch > inputs.size() ? inputs.size() : start + batch;
                auto t = thread([&inputs, &streams, &conf, i, start, stop](){
                    FCSP fcsp(conf);
                    for_each(inputs.begin()+start, inputs.begin()+stop, [&streams, &fcsp, i](string& inp){
                        processFile(fcsp, inp, streams[i]);
                    });
                });
                threads[i] = move(t);
            }
            for (size_t i = 0; i < batches; i++) {
                threads[i].join();
                cout << streams[i].str();
            }
        }
    }
    catch(std::exception& e) {
        cerr << e.what() << endline;
    }
    return 0;
}
