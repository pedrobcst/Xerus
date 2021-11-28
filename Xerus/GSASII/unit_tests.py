########### SVN repository information ###################
# $Date: 2017-11-30 13:38:49 +0900 (木, 30 11月 2017) $
# $Author: toby $
# $Revision: 3166 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/unit_tests.py $
# $Id: unit_tests.py 3166 2017-11-30 04:38:49Z toby $
########### SVN repository information ###################
'''
*unit_tests: Self-test Module*
------------------------------

A script that can be run to test a series of self-tests in GSAS-II. At present,
only modules ``GSASIIspc`` and ``GSASIIlattice`` have self-tests. 

'''

import GSASIIspc
import GSASIIlattice
def test_GSASIIspc():
    '''Test registered self-tests in ``GSASIIspc``.
    Takes no input and returns nothing. Throws an Exception if a test fails.
    '''
    #GSASIIspc.selftestquiet = False
    for test in GSASIIspc.selftestlist:
        test()
def test_GSASIIlattice():
    '''Test registered self-tests in ``GSASIIlattice``.
    Takes no input and returns nothing. Throws an Exception if a test fails.
    '''
    #GSASIIlattice.selftestquiet = False
    for test in GSASIIlattice.selftestlist:
        test()

if __name__ == '__main__':
    test_GSASIIspc()
    test_GSASIIlattice()
    print("OK")
