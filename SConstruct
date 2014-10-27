env = Environment()
if env['CXX'] == 'cl':
    env.Append(CCFLAGS="/EHsc")
elif env['CXX'] == 'g++':
	env.Append(CCFLAGS="-std=gnu++11")
libs = ["boost_program_options", "boost_system", "boost_filesystem"]
src = Glob("*.cpp")
env.Program(b'fcsp', src, LIBS=libs);