import os
prefix = "/usr/local"

AddOption('--release', dest='release', action='store_true', help='release build')
release = GetOption('release')
env = Environment()
if env['CXX'] == 'cl':
    copts = "/EHsc"
elif env['CXX'] == 'g++' or env['CXX'] == 'clang++':
	copts = "-std=gnu++11 -static -Wall "
	if release:
		copts += "-O2 "
	else:
		copts += "-g "

env.Append(CCFLAGS=copts)
libs = ["boost_program_options", "boost_system", "boost_filesystem"]
src = Glob("src/*.cpp")
prog = env.Program('fcss-2a', src, LIBS=libs);
env.Alias("install", env.Install(os.path.join(prefix, "bin"), prog))
env.Alias("install", env.Install(os.path.join(prefix, "bin"), 'fcss-comp'))
