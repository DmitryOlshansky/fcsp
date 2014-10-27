env = Environment()
if env['CC'] == 'cl':
    env.Append(CCFLAGS="/EHsc")
src = Glob("*.cpp")
env.Program(b'fcsp', src);