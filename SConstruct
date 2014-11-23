env = Environment()
if env['CXX'] == 'cl':
    env.Append(CCFLAGS="/EHsc")
elif env['CXX'] == 'g++':
	env.Append(CCFLAGS="-std=gnu++11")
libs = ["boost_program_options", "boost_system", "boost_filesystem"]
src = Glob("src/*.cpp")
env.Program('fcsp', src, LIBS=libs);