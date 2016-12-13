import os
import errno


def make_directory(path):
	try:
		os.makedirs(path)
	except OSError as exc:
		if exc.errno == errno.EEXIST and os.path.isdir(path):
			pass
		else:
			raise


def safe_open(filename, read_write):
	make_directory(os.path.dirname(filename))
	return open(filename, read_write)
