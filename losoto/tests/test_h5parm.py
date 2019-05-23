from .common_setup import *

from ..h5parm import h5parm

def test_h5parm():
    H = h5parm(os.path.join(TEST_FOLDER,'test_h5parm.h5'))
    assert isinstance(H, h5parm)
    assert isinstance(str(H), str)
    H.makeSolset('test_solset')
    assert 'test_solset' in H.getSolsets(0
    H.close()
    #TODO: assert closed somehow

