from pathlib import Path

MAIN_DIRECTORY = Path(__file__).parent.parent.parent.parent


def get_full_path(*path):
    return Path(MAIN_DIRECTORY, *path)
