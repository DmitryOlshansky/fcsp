import glob, os, sys

prefix = "/usr/local"

AddOption('--release', dest='release', action='store_true', help='release build')

release = GetOption('release')
boost = ARGUMENTS.get('boost')
env = Environment()
#env["CXX"] = "powerpc64le-linux-gnu-g++"
if env['CC'] == 'cl':
    if not boost:
        print "Needs path to boost on Windows e.g. scons boost=<path-to-boost>"
        sys.exit(1)
    boost_lib = boost+"\\stage\\lib"
    copts = "/EHsc /MD /I%s /link %s" % (boost, boost+"\\stage\\lib")
    if release:
        copts += " /O2"
    libs = ["libboost_program_options-*", "libboost_system-*", "libboost_filesystem-*"]
    libs = [glob.glob(boost_lib+'\\'+lib)[0] for lib in libs]
elif 'g++' in env['CXX'] or 'clang++' in env['CXX']:
    copts = "-std=gnu++11 -static -Wno-sign-compare -Wall "
    if release:
        copts += "-O2 "
    else:
        copts += "-g "
    libs = ["boost_system", "boost_filesystem"]

env.Append(CCFLAGS=copts, LINKFLAGS="-static -L../boost_1_64_0/stage/lib/")
src = Glob("src/*.cpp")
prog = env.Program('fcss-2a', src, LIBS=libs);
env.Alias("install", env.Install(os.path.join(prefix, "bin"), prog))
env.Alias("install", env.Install(os.path.join(prefix, "bin"), 'fcss-comp'))
env.Alias("install", env.Install(os.path.join(prefix, "share/fcss-2a/descr"), ['descr1.csv', 'descr2.sdf', 'replacement.sdf']))
