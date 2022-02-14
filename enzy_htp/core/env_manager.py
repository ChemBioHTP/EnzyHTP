""""""
import os
import shutil
import logging
from .exception import MissingEnvironmentElement
class EnvironmentManager:
	def __init__( self, **kwargs ):
		self.env_vars_ = kwargs.get('env_vars', [] )
		self.executables_ = kwargs.get('executables', [] )
		self.missing_env_vars_ = []
		self.missing_executables_ = []

	def add_executable( self, exe_name ):
		self.executables_.append( exe_name )

	def add_env_var(self, env_var):
		self.env_vars_.append( env_var )
	
	def check_env_vars(self):
		for env_var in self.env_vars_:
			print(os.getenv(env_var))
			if os.getenv(env_var) is None:
				self.missing_env_vars_.append( env_var )

	def check_executables(self):
		def exe_exists( exe_name ):
			full_path = os.path.expandvars( exe_name )
			return shutil.which(exe_name) is not None
		
		for exe in self.executables_:
			if not exe_exists( exe ):
				self.missing_executables_.append( exe )

	def display_missing(self):
		if not self.is_missing():
			logging.info('Environment has all required elements!')
			return	
		logging.warning('Environment is missing some required elements...')

		if len(self.missing_executables_):
			logging.warning('\tMissing excecutables:')
			for me in self.missing_executables_:
				logging.warning(f'\t\t{me}')
		
		if len(self.missing_env_vars_):
			logging.warning('\tMissing environment variables:')
			for mev in self.missing_env_vars_:
				logging.warning(f'\t\t{mev}')


	def check_environment(self):
		logging.info('Checking environment for required elements...')
		self.check_env_vars()
		self.check_executables()
		self.display_missing()
		logging.info('Environment check completed!')

	def reset(self):
		self.executables_ = []
		self.missing_executables_ = []
		self.env_vars_ = []
		self.missing_env_vars_ = []

	def is_missing(self):
		return len(self.missing_executables_) or len(self.missing_env_vars_)


