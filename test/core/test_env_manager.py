import os
from enzy_htp.core import EnvironmentManager


def test_check_environment():
	os.environ['___VAR1___'] = 'VALUE1'
	os.environ['___VAR2___'] = 'VALUE2'
	os.environ['___VAR3___'] = 'VALUE3'
	em = EnvironmentManager()
	em.add_env_var('___VAR1___')
	em.add_env_var('___VAR2___')
	em.add_env_var('___VAR3___')
	em.add_executable('ls')
	em.add_executable('less')
	em.add_executable('cat')
	em.check_environment()
	assert not em.is_missing()
	
	em = EnvironmentManager(
		env_vars=['___VAR1___','___VAR2___','___VAR3___'],
		executables=['ls','less','cat']
	)
	em.check_environment()
	assert not em.is_missing()
	os.unsetenv('___VAR1___')
	os.unsetenv('___VAR2___')
	os.unsetenv('___VAR3___')


def test_missing_env_var():
	em = EnvironmentManager()
	em.add_env_var('___VAR4___')
	em.check_environment()
	assert em.is_missing()
	em.reset()
	assert not em.is_missing()

def test_missing_executable():
	em = EnvironmentManager()
	em.add_executable('does-not-exist')
	em.check_environment()
	assert em.is_missing()
	em.reset()
	assert not em.is_missing()

