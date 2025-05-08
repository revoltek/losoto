from losoto.h5parm import h5parm
from losoto.operations.reweight import run

def test_reweight():
    """
    """
    with h5parm('test_reweight.h5', readonly=False) as h5:
        ss=h5.getSolset('sol000')
        st=ss.getSoltab('amplitude000')
        run(st, mode='window')


if __name__ == "__main__":
    test_reweight()