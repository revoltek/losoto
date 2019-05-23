import os
import sys

from .. import logging
import numpy as np
import pytest

TEST_FOLDER = os.path.abspath('./test_output')
os.makedirs(TEST_FOLDER,exist_ok=True)

def clean_test_output():
    logging.debug("Removing {}".format(TEST_FOLDER))
    os.unlink(TEST_FOLDER)


@pytest.fixture
def ex_h5parm():
    return h5parm(os.path.join(TEST_FOLDER,'test_h5parm.h5'))
