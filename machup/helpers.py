"""Collection of helper functions and classes."""

import os


def check_valid_filepath(filepath):
    """Check if the file path is valid.

    Parameters
    ----------
    filepath : str
        The path of the file in question.

    Raises
    ------
    IOError
        If filepath is not valid.

    """
    if not os.path.isfile(filepath):
        print('Error: Connot find file "{0}". Make sure'.format(filepath))
        print(' the path is correct and the file is accessible.')
        raise IOError(filepath)
