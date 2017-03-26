
import os
import errno


def make_directory(path):
	"""Make a directory reliably"""
	try:
		os.makedirs(path)
	except OSError as exc:
		if exc.errno == errno.EEXIST and os.path.isdir(path):
			pass
		else:
			raise


def safe_open(filename, read_write):
	"""Open a file and create subdirectories as necessary."""
	make_directory(os.path.dirname(filename))
	return open(filename, read_write)
