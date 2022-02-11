import os
import shutil


def safe_rm( fname : str ) -> None:
	if os.path.exists( fname ):
		os.remove( fname )

def safe_rmdir( dirname : str ) -> None:
	if os.path.isdir( dirname ):	
		shutil.rmtree( dirname )
