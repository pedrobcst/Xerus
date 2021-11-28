#!/usr/bin/env python
# Installs GSAS-II from network using subversion and creates platform-specific shortcuts.
# works for Mac & Linux & Windows
from __future__ import division, print_function
import os, stat, sys, platform, subprocess, datetime

version = "$Id: bootstrap.py 4645 2020-11-04 00:19:08Z toby $"
g2home = 'https://subversion.xray.aps.anl.gov/pyGSAS/'
path2GSAS2 = os.path.dirname(os.path.abspath(os.path.expanduser(__file__)))

skipInstallChecks = False  # if set to True, continue installation even if current Python lacks packages and
                           # skips platform-dependent installation steps
skipDownloadSteps = False  # if False, svn is used to download GSAS-II files and binaries 
skipProxy = False          # if False then the script creates a prompt asking for proxy info
showWXerror = False        # if True, errors are shown in wxPython windows (used on Windows only)
help = False               # if True, a help message is shown
allBinaries = False        # if True, causes all binaries to be installed
binaryVersion=None         # if specified, gives a specific numpy version to use (such as 1.18)
                           # for selecting a set of binaries to use
numpyVersion=None

for a in sys.argv[1:]:
    if 'noinstall' in a.lower():
        skipInstallChecks = True
        if sys.platform.startswith('win'): showWXerror = True
    elif 'nonet' in a.lower():
        skipDownloadSteps = True
        skipProxy = True
    elif 'noproxy' in a.lower():
        skipProxy = True
    elif 'allbin' in a.lower() or 'server' in a.lower():
        allBinaries = True
    elif 'binary' in a.lower():
        numpyVersion = a.split('=')[1]
        skipInstallChecks = True
        if sys.platform.startswith('win'): showWXerror = True
        skipProxy = True
    else:
        help = True
    if 'help' in a.lower():
        help = True

if help:
    print('''
  bootstrap.py options:

    --noinstall   skip post-install, such as creating run shortcuts, turns on 
                  wxpython error display in Windows
    --noproxy     do not ask for proxy information
    --server      load all binary versions 
    --allbin      load all binary versions (same as --server)
    --help        this message
    --nonet       skip steps requiring internet
    --binary=ver  specifies a specific numpy binary directory to use (such as 
                  --binary=1.18); turns on --noinstall and --noproxy
''')
    sys.exit()

now = str(datetime.datetime.now())
print('Running bootstrap from {} at {}\n\tId: {}'.format(path2GSAS2,now,version))
print ("Python:     %s"%sys.version.split()[0])
try:
    import numpy as np
    print ("numpy:      %s"%np.__version__)
except:
    pass
fp = open(os.path.join(path2GSAS2,'bootstrap.log'),'a')
fp.write('Running bootstrap from {} at {}\n\tId: {}\n'.format(path2GSAS2,now,version))
fp.close()
        
################################################################################
################################################################################
def BailOut(msg):
    '''Exit with an error message. Use a GUI to show the error when
    showWXerror == True (on windows when --noinstall is specified)
    '''
    print(msg)
    if showWXerror:
        import wx
        app = wx.App()
        app.MainLoop()
        dlg = wx.MessageDialog(None,msg,'GSAS-II installation error', 
                wx.OK | wx.ICON_ERROR | wx.STAY_ON_TOP)
        dlg.Raise()
        dlg.ShowModal()
        dlg.Destroy()
    else:
        print("Recreate error this using command:")
        print("     {} {} {}".format(
            os.path.abspath(sys.executable),
            os.path.abspath(os.path.expanduser(__file__)),
            ' '.join(sys.argv[1:]),
        ),file=sys.stderr)
        print("\nBOOTSTRAP WARNING: ",file=sys.stderr)
        for line in msg.split('\n'):
            print(line,file=sys.stderr)
    sys.exit()
        
def GetConfigValue(*args): return True
# routines copied from GSASIIpath.py
proxycmds = []
'Used to hold proxy information for subversion, set if needed in whichsvn'
svnLocCache = None
'Cached location of svn to avoid multiple searches for it'

def MakeByte2str(arg):
    '''Convert output from subprocess pipes (bytes) to str (unicode) in Python 3.
    In Python 2: Leaves output alone (already str). 
    Leaves stuff of other types alone (including unicode in Py2)
    Works recursively for string-like stuff in nested loops and tuples.

    typical use::

        out = MakeByte2str(out)

    or::

        out,err = MakeByte2str(s.communicate())
    
    '''
    if isinstance(arg,str): return arg
    if isinstance(arg,bytes):
        try:
            return arg.decode()
        except:
            print('Decode error')
            return arg
    if isinstance(arg,list):
        return [MakeByte2str(i) for i in arg]
    if isinstance(arg,tuple):
        return tuple([MakeByte2str(i) for i in arg])
    return arg

def getsvnProxy():
    '''Loads a proxy for subversion from the file created by bootstrap.py
    '''
    global proxycmds
    proxycmds = []
    proxyinfo = os.path.join(os.path.expanduser('~/.G2local/'),"proxyinfo.txt")
    if not os.path.exists(proxyinfo):
        proxyinfo = os.path.join(path2GSAS2,"proxyinfo.txt")
    if not os.path.exists(proxyinfo):
        return '','',''
    fp = open(proxyinfo,'r')
    host = fp.readline().strip()
    # allow file to begin with comments
    while host.startswith('#'):
        host = fp.readline().strip()
    port = fp.readline().strip()
    etc = []
    line = fp.readline()
    while line:
        etc.append(line.strip())
        line = fp.readline()
    fp.close()
    setsvnProxy(host,port,etc)
    return host,port,etc

def setsvnProxy(host,port,etc=[]):
    '''Sets the svn commands needed to use a proxy
    '''
    global proxycmds
    proxycmds = []
    host = host.strip()
    port = port.strip()
    if host: 
        proxycmds.append('--config-option')
        proxycmds.append('servers:global:http-proxy-host='+host)
        if port:
            proxycmds.append('--config-option')
            proxycmds.append('servers:global:http-proxy-port='+port)
    for item in etc:
        proxycmds.append(item)
        
def whichsvn():
    '''Returns a path to the subversion exe file, if any is found.
    Searches the current path after adding likely places where GSAS-II
    might install svn. 

    :returns: None if svn is not found or an absolute path to the subversion
      executable file.
    '''
    # use a previosuly cached svn location
    global svnLocCache
    if svnLocCache: return svnLocCache
    # prepare to find svn
    is_exe = lambda fpath: os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    svnprog = 'svn'
    if sys.platform.startswith('win'): svnprog += '.exe'
    host,port,etc = getsvnProxy()
    if GetConfigValue('debug') and host:
        print('DBG_Using proxy host {} port {}'.format(host,port))
    # add likely places to find subversion when installed with GSAS-II
    pathlist = os.environ["PATH"].split(os.pathsep)
    pathlist.insert(0,os.path.split(sys.executable)[0])
    pathlist.insert(1,path2GSAS2)
    for rpt in ('..','bin'),('..','Library','bin'),('svn','bin'),('svn',),('.'):
        pt = os.path.normpath(os.path.join(path2GSAS2,*rpt))
        if os.path.exists(pt):
            pathlist.insert(0,pt)    
    # search path for svn or svn.exe
    for path in pathlist:
        exe_file = os.path.join(path, svnprog)
        if is_exe(exe_file):
            try:
                p = subprocess.Popen([exe_file,'help'],stdout=subprocess.PIPE)
                res = p.stdout.read()
                p.communicate()
                svnLocCache = os.path.abspath(exe_file)
                return svnLocCache
            except:
                pass        
    svnLocCache = None

def svnVersion(svn=None):
    '''Get the version number of the current subversion executable

    :returns: a string with a version number such as "1.6.6" or None if
      subversion is not found.

    '''
    if not svn: svn = whichsvn()
    if not svn: return

    cmd = [svn,'--version','--quiet']
    s = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = MakeByte2str(s.communicate())
    if err:
        print ('subversion error!\nout=%s'%out)
        print ('err=%s'%err)
        s = '\nsvn command:  '
        for i in cmd: s += i + ' '
        print(s)
        return None
    return out.strip()

def svnVersionNumber(svn=None):
    '''Get the version number of the current subversion executable

    :returns: a fractional version number such as 1.6 or None if
      subversion is not found.

    '''
    ver = svnVersion(svn)
    if not ver: return 
    M,m = ver.split('.')[:2]
    return int(M)+int(m)/10.

def showsvncmd(cmd):
    s = '\nsvn command:  '
    for i in cmd: s += i + ' '
    print(s)

def svnChecksumPatch(svn,fpath,verstr):
    '''This performs a fix when svn cannot finish an update because of
    a Checksum mismatch error. This seems to be happening on OS X for 
    unclear reasons. 
    '''
    print('\nAttempting patch for svn Checksum mismatch error\n')
    cmd = [svn,'cleanup',fpath]
    showsvncmd(cmd)        
    s = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = MakeByte2str(s.communicate())
    #if err: print('error=',err)
    cmd = ['svn','update','--set-depth','empty',
               os.path.join(fpath,'bindist')]
    showsvncmd(cmd)        
    s = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = MakeByte2str(s.communicate())
    #if err: print('error=',err)
    try:
        import GSASIIpath
        print('import of GSASIIpath completed')
    except Exception as err:
        msg = 'Failed with import of GSASIIpath. This is unexpected.'
        msg += '\nGSAS-II will not run without correcting this. Contact toby@anl.gov'
        BailOut(msg)
    cmd = ['svn','switch',g2home+'/trunk/bindist',
               os.path.join(fpath,'bindist'),
               '--non-interactive', '--trust-server-cert', '--accept',
               'theirs-conflict', '--force', '-rHEAD', '--ignore-ancestry']
    showsvncmd(cmd)        
    s = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = MakeByte2str(s.communicate())
    GSASIIpath.DownloadG2Binaries(g2home,verbose=True)
    cmd = ['svn','update','--set-depth','infinity',
               os.path.join(fpath,'bindist')]
    showsvncmd(cmd)        
    s = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = MakeByte2str(s.communicate())
    #if err: print('error=',err)
    cmd = [svn,'update',fpath,verstr,
                       '--non-interactive',
                       '--accept','theirs-conflict','--force']
    if svnVersionNumber() >= 1.6:
        cmd += ['--trust-server-cert']
    if proxycmds: cmd += proxycmds
    #print(cmd)
    showsvncmd(cmd)        
    s = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = MakeByte2str(s.communicate())
    #if err: print('error=',err)
    return err

################################################################################
################################################################################
print(70*'*')
#testing for incorrect locale code'
try:
    import locale
    locale.getdefaultlocale()
except ValueError:
    print('Your location is not set properly. This causes problems for matplotlib')
    print('  (see https://github.com/matplotlib/matplotlib/issues/5420.)')
    print('Will try to bypass problem by setting LC_ALL to en_US.UTF-8 (US English)')
    os.environ['LC_ALL'] = 'en_US.UTF-8'
    locale.getdefaultlocale()
print('Preloading matplotlib to build fonts...')
try:
    import matplotlib
except:
    pass
print('Checking python packages...')
missing = []
for pkg in ['numpy','scipy','matplotlib','wx','OpenGL',]:
    try:
        exec('import '+pkg)
    except:
        missing.append(pkg)

if missing and not skipInstallChecks:
    msg = """Sorry, this version of Python cannot be used
for GSAS-II. It is missing the following package(s):
\t"""
    for pkg in missing: msg += " "+pkg
    msg += "\nPlease install these package(s) and try running bootstrap.py again."
    #print("Showing first error: ")
    #for pkg in ['numpy','scipy','matplotlib','wx','OpenGL',]:
    #    exec('import '+pkg)
    BailOut(msg)

if not skipDownloadSteps:
    host = None
    port = '80'
    print('\nChecking for subversion...')
    svn = whichsvn() # resets host & port if proxyinfo.txt is found
    if not svn:
        msg ="Sorry, subversion (svn) could not be found on your system."
        msg += "\nPlease install this or place in path and rerun this."
        BailOut(msg)
    else:
        print(' found svn image: '+svn)

#if install_with_easyinstall:  	 	 
#    print('\nInstalling PyOpenGL. Lots of warnings will follow... ')
#    install_with_easyinstall('PyOpenGl')  	 	 
#    print('done.')
    
print('Ready to bootstrap GSAS-II from repository\n\t'+g2home+'\nto '+path2GSAS2)
proxycmds = []
host,port,etc = getsvnProxy()
if sys.version_info[0] == 2:
    getinput = raw_input
else:
    getinput = input

# get proxy setting from environment variable
key = None
for i in os.environ.keys():
    if 'https_proxy' == i.lower():
        key = i
        break
else:
    for i in os.environ.keys():
        if 'http_proxy' == i.lower():
            key = i
            break
val = ''
if key:
    val = os.environ[key].strip()
if val:
    if val[-1] == '/':
        val = val[:-1]
    if len(val.split(':')) > 2:
        host = ':'.join(val.split(':')[:-1])
        port = val.split(':')[-1]
    else:
        host = ':'.join(val.split(':')[:-1])
        port = val.split(':')[-1]

# get proxy from user, if terminal available
try:
    if skipProxy:
        host = ""
    elif host:
        print('\n'+75*'*')
        ans = getinput("Enter the proxy address (type none to remove) ["+host+"]: ").strip()
        if ans.lower() == "none": host = ""
    else:
        ans = getinput("Enter your proxy address [none needed]: ").strip()
        if ans: host = ans
    if host:
        ans = getinput("Enter the proxy port ["+port+"]: ").strip()
        if ans == "": ans=port
        port = ans
        print('If your site needs additional svn commands (such as \n\t',
                  '--config-option servers:global:http-proxy-username=*account*','\n\t',
                  '--config-option servers:global:http-proxy-password=*password*',
                  '\nenter them now:')
        if etc:
            prevetc = ' '.join(etc)
            print('\nDefault for next input is "{}"'.format(prevetc))
            prompt = "Enter additional svn options (if any) [use previous]: "
        else:
            prompt = "Enter additional svn options (if any) [none]: "
        ans = 'start'
        etcstr = ''
        while ans:
            ans = getinput(prompt).strip()
            prompt = "more svn options (if any): "
            if etcstr: etcstr += ' '
            etcstr += ans
        if etcstr.strip():
           etc = etcstr.split()
except EOFError:
    host = ""
    port = ""
    etc = []
setsvnProxy(host,port,etc)
# delete old proxy files
localproxy = os.path.join(os.path.expanduser('~/.G2local/'),"proxyinfo.txt")
for proxyinfo in localproxy,os.path.join(path2GSAS2,"proxyinfo.txt"):
    if os.path.exists(proxyinfo):
        try:
            os.remove(proxyinfo)
            print('Deleted file {}'.format(proxyinfo))
        except:
            pass
if host:
    try:
        fp = open(proxyinfo,'w')
    except:
        fp = open(localproxy,'w')
        proxyinfo = localproxy
    try:
        fp.write(host.strip()+'\n')
        fp.write(port.strip()+'\n')
        for i in etc:
            if i.strip():
                fp.write(i.strip()+'\n')
        fp.close()
        msg = 'Proxy info written: {} port {} etc {}\n'.format(host,port,etc)
        print(msg)
        fp = open(os.path.join(path2GSAS2,'bootstrap.log'),'a')
        fp.write(msg)
        fp.close()
    except Exception as err:
        print('Error writing file {}:\n{}'.format(proxyinfo,err))
        
if not skipDownloadSteps:
    # patch: switch GSAS-II location if linked to XOR server (relocated May/June 2017)
    cmd = [svn, 'info']
    p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    res,err = p.communicate()
    if '.xor.' in str(res):
        print('Switching previous install with .xor. download location to\n\thttps://subversion.xray.aps.anl.gov/pyGSAS')
        cmd = [svn, 'switch','--relocate','https://subversion.xor.aps.anl.gov/pyGSAS',
               'https://subversion.xray.aps.anl.gov/pyGSAS']
        if proxycmds: cmd += proxycmds
        p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        res,err = p.communicate()
        if err:
            print('Please report this error to toby@anl.gov:')
            print(err)
            print(res)
    # patch: switch GSAS-II location if switched to 2frame version (removed August 2017)
    if '2frame' in str(res):
        print('Switching previous 2frame install to trunk\n\thttps://subversion.xray.aps.anl.gov/pyGSAS')
        cmd = [svn, 'switch',g2home + '/trunk',path2GSAS2,
               '--non-interactive','--trust-server-cert',
               '--accept','theirs-conflict','--force','--ignore-ancestry']
        if proxycmds: cmd += proxycmds
        p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        res,err = p.communicate()
        if err:
            print('Please report this error to toby@anl.gov:')
            print(err)
            print(res)

    print('\n'+75*'*')
    print('Now preparing to install GSAS-II')
    tryagain = True
    err = False
    firstPass = 0
    while(tryagain):
        tryagain = False
        if err:
            print('Retrying after a cleanup...')
            cmd = [svn, 'cleanup', path2GSAS2]
            s = subprocess.Popen(cmd,stderr=subprocess.PIPE)
            out,err = MakeByte2str(s.communicate())
            if err:
                print('subversion returned an error:')
                print(out)
                print(err)
        cmd = [svn, 'co', g2home+ 'trunk/', path2GSAS2, '--non-interactive', '--trust-server-cert']
        if proxycmds: cmd += proxycmds
        msg = 'svn load command: '
        for item in cmd: msg += " "+item
        print(msg)
        s = subprocess.Popen(cmd,stderr=subprocess.PIPE)
        print('\nsubversion output:')
        out,err = MakeByte2str(s.communicate())
        if 'Checksum' in err:  # deal with Checksum problem
            err = svnChecksumPatch(svn,path2GSAS2,'-rHEAD')
            if err:
                print('error from svnChecksumPatch\n\t',err)
        elif err:
            print('subversion returned an error:')
            print(out)
            print(err)
            if firstPass == 0: tryagain = True
        firstPass += 1
    if err:
        print('Retrying with a command for older svn version...')
        cmd = [svn, 'co', g2home+ 'trunk/', path2GSAS2]
        if proxycmds: cmd += proxycmds
        msg = ""
        for item in cmd: msg += " " + item
        print(msg)
        s = subprocess.Popen(cmd,stderr=subprocess.PIPE)
        out,err = MakeByte2str(s.communicate())
        if err:
            msg = 'subversion returned an error:\n'
            msg += err
            if os.path.exists(os.path.join(path2GSAS2,"makeBat.py")):
                msg += '\nGSAS-II appears to be installed but failed to be updated. A likely reason'
                msg += '\nis a network access problem. If your web browser works, but this update did'
                msg += '\nnot, the most common reason is you need to use a network proxy. Please'
                msg += '\ncheck with a network administrator or use http://www.whatismyproxy.com/'
                msg += '\n\nIf installing from the gsas2full dist and your computer is not connected'
                msg += '\nto the internet, this error message can be ignored.\n'
            else:
                # this will happen only with initial installs where all files
                # are to be downloaded (not gsas2full or updates)
                msg += '\n\n  *** GSAS-II failed to be installed. A likely reason is a network access'
                msg += '\n  *** problem, most commonly because you need to use a network proxy. Please'
                msg += '\n  *** check with a network administrator or use http://www.whatismyproxy.com/\n'
            BailOut(msg)
    print('\n'+75*'*')

# subsequent commands require GSASIIpath which better be here now, import it
try:
    import GSASIIpath
    print('import of GSASIIpath completed')
except Exception as err:
    msg = 'Failed with import of GSASIIpath. This is unexpected.'
    msg += '\nGSAS-II will not run without correcting this. Contact toby@anl.gov'
    BailOut(msg)

if skipDownloadSteps:
    pass
elif allBinaries:
    print('Loading all binaries with command...')
    if not GSASIIpath.svnSwitchDir('AllBinaries','',g2home+ 'Binaries/',None,True):
        msg = 'Binary load failed. Subversion problem? Please seek help'
        BailOut(msg)
elif numpyVersion:
    binaryVersion = GSASIIpath.GetBinaryPrefix()+'_n'+numpyVersion
    if not GSASIIpath.svnSwitchDir('bindist','',g2home+ 'Binaries/'+binaryVersion,None,True):
        msg = 'Binary load failed with '+binaryVersion+'. Subversion problem? Please seek help'
        BailOut(msg)
else:
    GSASIIpath.DownloadG2Binaries(g2home)
        
#===========================================================================
# test if the compiled files load correctly
#===========================================================================
GSASIItested = False
if not skipInstallChecks:
    script = """  
# commands that test each module can at least be loaded & run something in pyspg
try:
    import GSASIIpath
    GSASIIpath.SetBinaryPath(loadBinary=False)
    import pyspg
    import histogram2d
    import polymask
    import pypowder
    import pytexture
    pyspg.sgforpy('P -1')
    print('==OK==')
except Exception as err:
    print(err)
"""
    p = subprocess.Popen([sys.executable,'-c',script],stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                         cwd=path2GSAS2)
    res,err = MakeByte2str(p.communicate())
    if '==OK==' not in str(res) or p.returncode != 0:
        #print('\n'+75*'=')
        msg = 'Failed when testing the GSAS-II compiled files. GSAS-II will not run'
        msg += ' without correcting this.\n\nError message:\n'
        if res: 
            msg += res
            msg += '\n'
        if err:
            msg += err
        #print('\nAttempting to open a web page on compiling GSAS-II...')
        msg += '\n\nPlease see web page\nhttps://subversion.xray.aps.anl.gov/trac/pyGSAS/wiki/CompileGSASII if you wish to compile for yourself (usually not needed for windows and Mac, but sometimes required for Linux.)'
        BailOut(msg)
        #import webbrowser
        #webbrowser.open_new('https://subversion.xray.aps.anl.gov/trac/pyGSAS/wiki/CompileGSASII')
        #print(75*'=')
    #    if '86' in platform.machine() and (sys.platform.startswith('linux')
    #                                        or sys.platform == "darwin"
    #                                        or sys.platform.startswith('win')):
    #        print('Platform '+sys.platform+' with processor type '+platform.machine()+' is supported')
    #    else:
    #        print('Platform '+sys.platform+' with processor type '+platform.machine()+' not is supported')
    else:
        print('Successfully tested compiled routines')
        GSASIItested = True
#===========================================================================
# import all .py files so that .pyc files get created
if not skipInstallChecks:
    print('Byte-compiling all .py files...')
    import compileall
    compileall.compile_dir(path2GSAS2,quiet=True)
    print('done')
#===========================================================================
# do platform-dependent stuff
#===========================================================================
if sys.version_info[0] > 2:
    def execfile(file):
        with open(file) as source_file:
            exec(source_file.read())

if skipInstallChecks:
    pass
#===========================================================================
# on Windows, make a batch file with Python and GSAS-II location hard-coded
elif sys.platform.startswith('win') and os.path.exists(
    os.path.join(path2GSAS2,"makeBat.py")):
    execfile(os.path.join(path2GSAS2,"makeBat.py"))
#===========================================================================
# on a Mac, make an applescript 
elif sys.platform.startswith('darwin') and os.path.exists(
    os.path.join(path2GSAS2,"makeMacApp.py")):
    sys.argv = [os.path.join(path2GSAS2,"makeMacApp.py")]
    print(u'running '+sys.argv[0])
    execfile(sys.argv[0])
#===========================================================================
# On linux, make desktop icon
elif sys.platform.startswith('linux') and os.path.exists(
    os.path.join(path2GSAS2,"makeLinux.py")):
    sys.argv = [os.path.join(path2GSAS2,"makeLinux.py")]
    print(u'running '+sys.argv[0])
    execfile(sys.argv[0])

