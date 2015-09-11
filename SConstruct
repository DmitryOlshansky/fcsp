import os
prefix = "/usr/local"

env = Environment()
if env['CXX'] == 'cl':
    env.Append(CCFLAGS="/EHsc")
elif env['CXX'] == 'g++':
	env.Append(CCFLAGS="-std=gnu++11 -static -O2")
libs = ["boost_program_options", "boost_system", "boost_filesystem"]
src = Glob("src/*.cpp")
prog = env.Program('fcsp', src, LIBS=libs);
env.Alias("install", env.Install(os.path.join(prefix, "bin"), prog))