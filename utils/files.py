"""File module."""
from os.path import dirname, join
MAIN_DIRECTORY = dirname(join(dirname(__file__), '..', '..'))


def get_full_path(*path):
    """Return the full path."""
    return join(MAIN_DIRECTORY, *path)
