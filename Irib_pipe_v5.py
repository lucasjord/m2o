
from AIPS import AIPS
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage
import os
import sys
import subprocess

# This pipeline does: 1: FITLD, INDXR, 2: listr, prtan, plots CL1, makes a VPLOT, 3: POSSMs all srcs, Gnuplots all POSSMs, finds the peak chan for each source
# Run it using:
#    ParselTongue Irib_pipe.py x
# where x is the tmask to start from (defaults to 1 if not sepcified)



############ FITS file

# Get the fits file from the current directory
fits = subprocess.check_output("ls *IDI*", shell=True).strip()
fits = fits.decode()
print("Fits =", fits)

# Or set the fits file explicitly
#fits = 'irib5.IDI'


############# Antab file
cat = "cat *antabfs > ./antab"
os.system(cat)

############# Flag file
cat = "cat *uvflgfs > ./uvflg"
os.system(cat)

# The AIPS ID 
AIPS.userno = 15091




# Get the Tmask either from command or from this script
if len(sys.argv) > 1:
	tmask =int(sys.argv[1])
else:
	tmask = 1

refant = 1

#################### This is where you set the restfrequency of the maser
tfreq = 6668519200
freq_h = float(int(tfreq*1e-6))*1e6
freq_l = tfreq - freq_h


# Open a file called "essentials" where some useful info gets printed
essentials = open("essentials.txt","w")


##### Print some info at the start
print("  Doing from Tmask:", tmask)
print("  AIPS ID:", AIPS.userno)
print("  Refant:", refant)
print("  FreqRef:", tfreq, freq_h, freq_l)

##### Define the different data types we will use
uvdata = AIPSUVData('MULTI','UVDATA',1,1)
splatdata1 = AIPSUVData('MULTI','SPLAT',1,1)
splatdata2 = AIPSUVData('MULTI','SPLAT',1,2)
splatdata11 = AIPSUVData('MULTI','SPLAT',1,3)
splatdata22 = AIPSUVData('MULTI','SPLAT',1,4)



############ Start loading and sorting

if tmask <= 1:
    if uvdata.exists():
        print("Zapping old data")
        uvdata.zap()
    fitld = AIPSTask('FITLD')
    datain = 'PWD:' + fits
    fitld.datain = datain
    #fitld.outdata = uvdata
    fitld.clint = 0.166
    print("Fits is:", datain)
    fitld.doconcat = 1
    fitld.ncount = 2
    fitld.go()
    #########Index data
    indxr = AIPSTask('indxr')
    indxr.indata = uvdata
    indxr.cparm[3] = 0.166
    indxr.go()
    ########

##### Define list of targets in the full project source list (Currently just three but there will be 30 or so when the source list is finalised by Artis. Call this 'all_masers'
##### Then define a list called 'masers' which are the sources that are in both the 'all_masers' list and also in the uvdata in this particlar fits file.
##### Then define a list of all the sources that are in the uvdata but are not masers. These will be the calibrators, call this 'conts'.
all_masers = ['G85.411', 'CEPA', 'G24.33+0.14']
masers = set(uvdata.sources) & set(all_masers)
conts = list(set(uvdata.sources).difference(masers))

# Or if you want to explicitly set the maser and cont list:
#masers = ['G85.411', 'CEPA']
#conts = ['3C345', '3C454.3']

essentials.write("Maser sources are " + str(masers) + "\n")
essentials.write("Calibrator sources are " + str(conts) + "\n")


########## Have the Pipeline determine the number of channels, and calculate the 10% band edge
N_edgeCh = int((uvdata.header['naxis'][2])/10)
print (' Number of chans, and the 10% of it', (uvdata.header['naxis'][2]), N_edgeCh)

############################## Listr and Prtan
if tmask <= 2:
    essentials.write("Reached tmask 2 \n")
    if os.path.exists('listr.txt'):
         print("  Removing old listr.txt")
         os.remove('listr.txt')
            
    listr = AIPSTask('LISTR')
    listr.indata = uvdata
    listr.optype = 'SCAN'
    listr.outprint = './listr.txt'
    listr.go()
    #######
    if os.path.exists('prtan.txt'):
            print("  Removing old prtan.txt")
            os.remove('prtan.txt')
            
    prtan = AIPSTask('PRTAN')
    prtan.indata = uvdata
    prtan.outprint = './prtan.txt'
    prtan.go()
    #######

    while uvdata.table_highver('AIPS PL') > 0:
            print("  Removing old PL tables")
            uvdata.zap_table('AIPS PL', -1)

    if os.path.exists('VPLOT.eps'):
         print("  Removing old VPLOT.eps")
         os.remove('VPLOT.eps')



    #######
    snplt = AIPSTask('SNPLT')
    snplt.indata = uvdata
    snplt.inext = 'CL'
    snplt.inver = 1
    snplt.dotv = -1
    snplt.nplots = 8
    snplt.optype = 'AMP'
    snplt.go()
    #######
    if os.path.exists('CL1.eps'):
         print("removing old CL1.eps")
         os.remove('CL1.eps')

    lwpla = AIPSTask('lwpla')
    lwpla.indata = uvdata
    lwpla.plver = 1
    lwpla.inver = uvdata.table_highver('AIPS PL')
    lwpla.outfile = 'PWD:CL1.eps'
    lwpla.go()
    while uvdata.table_highver('AIPS PL') > 0:
            uvdata.zap_table('AIPS PL', 0)


    ###### Plotting CL1 VPLOT
    vplot = AIPSTask('vplot')
    vplot.indata = uvdata
    vplot.solint = 0.166
    vplot.bparm[2] = -1
    #vplot.avgchan = 1
    #vplot.avgif = 1
    vplot.dotv = -1
    vplot.refant = refant
    #vplot.nplots = 12
    # Also want to make bchan echan more automated
    vplot.bchan = 1
    vplot.echan = (uvdata.header['naxis'][2])
    vplot.crowded = 3
    vplot.do3col = 1
    vplot.go()
    print("bchan echan is", vplot.bchan, vplot.echan)
    #######

    if os.path.exists('VPLOT_CL1.eps'):
         print("removing old VPLOT_CL1.eps")
         os.remove('VPLOT_CL1.eps')

    lwpla = AIPSTask('lwpla')
    lwpla.indata = uvdata
    lwpla.plver = 1
    lwpla.inver = uvdata.table_highver('AIPS PL')
    lwpla.outfile = 'PWD:VPLOT_CL1.eps'
    lwpla.go()
    while uvdata.table_highver('AIPS PL') > 0:
            uvdata.zap_table('AIPS PL', 0)


################ Do ACCOR, ANTAB, APCAL
if tmask <= 3:
    essentials.write("Reached tmask 3 \n")
    #Start by exdesting any old SN tables
    while uvdata.table_highver('AIPS SN') > 0:
            uvdata.zap_table('AIPS SN', 0)
    #######
    accor = AIPSTask('accor')
    accor.indata = uvdata
    accor.solint = 1
    accor.go()
    #######
    clcal = AIPSTask('CLCAL')
    clcal.indata = uvdata
    clcal.snver = 1
    clcal.gainver = 1
    clcal.gainuse = 2
    clcal.interpol = 'self'
    clcal.go()
    ####### ANTAB
    #Start by exdesting any old GC tables
    while uvdata.table_highver('AIPS TY') > 0:
            uvdata.zap_table('AIPS TY', 0)
    while uvdata.table_highver('AIPS GC') > 0:
            uvdata.zap_table('AIPS GC', 0)
    antab = AIPSTask('ANTAB')
    antab.indata = uvdata
    antab.calin = 'PWD:antab'
    antab.go()
    #######
    apcal = AIPSTask('APCAL')
    apcal.indata = uvdata
    apcal.dofit[1] = -1
    print(apcal.dofit)
    apcal.go()
    #######
    clcal = AIPSTask('CLCAL')
    clcal.indata = uvdata
    clcal.snver = 2
    clcal.gainver = 2
    clcal.gainuse = 3
    clcal.interpol = 'self'
    clcal.go()
    #######
    #######  Plot the CL table
    while uvdata.table_highver('AIPS PL') > 0:
        uvdata.zap_table('AIPS PL', 0)
    snplt = AIPSTask('SNPLT')
    snplt.indata = uvdata
    snplt.inext = 'CL'
    snplt.inver = 3
    snplt.dotv = -1
    snplt.nplots = 8
    snplt.optype = 'AMP'
    snplt.go()
    #######
    if os.path.exists('CL3.eps'):
         print("removing old CL3.eps")
         os.remove('CL3.eps')

    lwpla = AIPSTask('lwpla')
    lwpla.indata = uvdata
    lwpla.plver = 1
    lwpla.inver = uvdata.table_highver('AIPS PL')
    lwpla.outfile = 'PWD:CL3.eps'
    lwpla.go()
    
    #######
    while uvdata.table_highver('AIPS PL') > 0:
            uvdata.zap_table('AIPS PL', 0)
    vplot = AIPSTask('vplot')
    vplot.indata = uvdata
    vplot.solint = 0.166
    vplot.bparm[2] = -1
    #vplot.avgchan = 1
    #vplot.avgif = 1
    vplot.dotv = -1
    vplot.refant = refant
    #vplot.nplots = 12
    # Also want to make bchan echan more automated
    vplot.bchan = 1
    vplot.echan = (uvdata.header['naxis'][2])
    vplot.crowded = 3
    vplot.do3col = 1
    vplot.docal = 1
    vplot.gainuse = 3
    vplot.stokes = 'RR'
    vplot.go()
    vplot.stokes = 'LL'
    vplot.go()
    print("bchan echan is", vplot.bchan, vplot.echan)
    #######
    if os.path.exists('VPLOT_CL3.eps'):
         print("removing old VPLOT_CL3.eps")
         os.remove('VPLOT_CL3.eps')

    lwpla = AIPSTask('lwpla')
    lwpla.indata = uvdata
    lwpla.plver = 1
    lwpla.inver = 2
    lwpla.outfile = 'PWD:VPLOT_CL3.eps'
    lwpla.go()


################ Do SETJY BPASS and POSSSM 
if tmask <= 4:
    essentials.write("Reached tmask 4 \n")
    ############ Flag using the flag file provided
    uvflg = AIPSTask('uvflg')
    uvflg.indata = uvdata
    uvflg.intext = 'PWD:uvflg'
    uvflg.opcode = 'flag'
    uvflg.go()
    ######### Flag the band edges (10%)
    uvflg = AIPSTask('uvflg')
    uvflg.indata = uvdata
    uvflg.opcode = 'flag'
    uvflg.bchan = 1
    uvflg.echan = N_edgeCh
    uvflg.go()
    uvflg.bchan = (uvdata.header['naxis'][2])-N_edgeCh
    uvflg.echan = (uvdata.header['naxis'][2])
    uvflg.go()
    #########
    setjy = AIPSTask('setjy')
    setjy.indata = uvdata
    #setjy.sources[1:] = ''
    setjy.restfreq[1:] = [freq_h,freq_l]
    setjy.veltyp = 'LSR'
    setjy.optype = 'VCAL'
    setjy.veldef = 'RADIO'
    setjy.go()
    ########## POSSM all IFs and Stokes for all srcs, no cal, no flag, no BP
    while uvdata.table_highver('AIPS PL') > 0:
            uvdata.zap_table('AIPS PL', 0)
    for source in uvdata.sources:
        possm = AIPSTask('possm')
        possm.indata = uvdata
        possm.solint = 666
        possm.sources[1:] = [source]
        possm.stokes = 'HALF'
        possm.flagver = -1 
        possm.nplots = 2
        possm.dotv = -1
        #possm.doband = 1
        possm.docal = -1
        possm.aparm[7] = 0
        possm.aparm[9] = 1
        possm.aparm[1] = -1
        # Fixed scale for phase
        possm.aparm[2] = 1
        possm.aparm[5] = -180
        possm.aparm[6] = 180
        possm.go()

        if os.path.exists('POSSM_CROSS_CL1.eps'):
             print("removing old POSSM_CROSS_CL1.eps")
             os.remove('POSSM_CROSS_CL1.eps')

        lwpla = AIPSTask('lwpla')
        lwpla.indata = uvdata
        lwpla.plver = 1
        lwpla.inver = 12
        lwpla.outfile = 'PWD:POSSM_CROSS_CL1.eps'
        lwpla.go()
    ########## POSSM all IFs and Stokes for all srcs, no cal, no flag, with BP
    while uvdata.table_highver('AIPS BP') > 0:
            uvdata.zap_table('AIPS BP', 0)
    bpass = AIPSTask('BPASS')
    bpass.indata = uvdata
    bpass.bpassprm[1] = 1
    bpass.docal = 1
    bpass.gainuse = 3
    bpass.bpassprm[1] = 1
    #bpass.bpassprm[10] = 1
    #bpass.bpassprm[5] = 1
    bpass.calsour[1:] = conts
    bpass.go()
    while uvdata.table_highver('AIPS PL') > 0:
            uvdata.zap_table('AIPS PL', 0)
    for source in uvdata.sources:
        possm = AIPSTask('possm')
        possm.indata = uvdata
        possm.solint = 666
        possm.sources[1:] = [source]
        possm.stokes = 'HALF'
        possm.flagver = -1 
        possm.nplots = 4
        possm.dotv = -1
        possm.doband = -1
        possm.docal = -1
        possm.aparm[7] = 0
        possm.aparm[9] = 1
        possm.aparm[1] = -1
        possm.aparm[8] = 1
        # Fixed scale for phase
        possm.aparm[2] = 1
        possm.aparm[5] = -180
        possm.aparm[6] = 180
        possm.go()

        if os.path.exists('POSSM_AUTO_BP.eps'):
             print("removing old POSSM_AUTO_BP.eps")
             os.remove('POSSM_AUTO_BP.eps')

        lwpla = AIPSTask('lwpla')
        lwpla.indata = uvdata
        lwpla.plver = 1
        lwpla.inver = 12
        lwpla.outfile = 'PWD:POSSM_AUTO_BP.eps'
        lwpla.go()
    while uvdata.table_highver('AIPS PL') > 0:
            uvdata.zap_table('AIPS PL', 0)
    for source in uvdata.sources:
        possm = AIPSTask('possm')
        possm.indata = uvdata
        possm.solint = 666
        possm.sources[1:] = [source]
        possm.stokes = 'HALF'
        possm.flagver = 1 
        possm.nplots = 4
        possm.dotv = -1
        possm.doband = -1
        possm.docal = 1
        possm.gainuse = 3
        possm.aparm[7] = 2
        possm.aparm[9] = 1
        possm.aparm[8] = 1
        # Fixed scale for phase
        possm.aparm[2] = 1
        possm.aparm[5] = -180
        possm.aparm[6] = 180
        possm.go()

        if os.path.exists('POSSM_AUTO_CL3.eps'):
             print("removing old POSSM_AUTO_CL3.eps")
             os.remove('POSSM_AUTO_CL3.eps')

        lwpla = AIPSTask('lwpla')
        lwpla.indata = uvdata
        lwpla.plver = 1
        lwpla.inver = 12
        lwpla.outfile = 'PWD:POSSM_AUTO_CL3.eps'
        lwpla.go()
    while uvdata.table_highver('AIPS PL') > 0:
            uvdata.zap_table('AIPS PL', 0)
    for source in uvdata.sources:
        possm = AIPSTask('possm')
        possm.indata = uvdata
        possm.solint = 666
        possm.sources[1:] = [source]
        possm.stokes = 'HALF'
        possm.aparm[1] = -1
        possm.flagver = 1 
        possm.nplots = 2
        possm.dotv = -1
        possm.doband = 1
        possm.docal = 1
        possm.gainuse = 3
        possm.aparm[7] = 2
        possm.aparm[9] = 1
        possm.aparm[1] = -1
        # Fixed scale for phase
        possm.aparm[2] = 1
        possm.aparm[5] = -180
        possm.aparm[6] = 180
        possm.go()

        if os.path.exists('POSSM_CROSS_CL3_BP.eps'):
             print("removing old POSSM_CROSS_CL3_BP.eps")
             os.remove('POSSM_CROSS_CL3_BP.eps')

        lwpla = AIPSTask('lwpla')
        lwpla.indata = uvdata
        lwpla.plver = 1
        lwpla.inver = 12
        lwpla.outfile = 'PWD:POSSM_CROSS_CL3_BP.eps'
        lwpla.go()

################ Do FRING (ALIGN IFs SO THEY CAN BE AVGD) zero Rates
if tmask <= 5:
        essentials.write("Reached tmask 5 \n")
        print("Tmaks 5")
        while uvdata.table_highver('AIPS SN') > 2:
                uvdata.zap_table('AIPS SN', 0)
        fring = AIPSTask('FRING')
        fring.indata = uvdata
        fring.calsour[1:] = conts
        fring.aparm = AIPSList([ 2,-1,0,0,0,0,3,0,0])
        fring.dparm = AIPSList([ 1,200,300,1,0,0,0,1])
        fring.doband = -1
        fring.refant = 1
        fring.docal = 1
        fring.gainuse = 3
        fring.solint = 10
        #Makes SN3
        fring.go()
        #######
        while uvdata.table_highver('AIPS CL') > 3:
                uvdata.zap_table('AIPS CL', 0)
        clcal = AIPSTask('CLCAL')
        clcal.indata = uvdata
        clcal.snver = 3
        clcal.gainver = 3
        clcal.gainuse = 4
        #clcal.sources = AIPSList(['3C345'])
        clcal.interpol = 'SIMP'
        clcal.bparm = AIPSList([0,0,0,1])
        clcal.dobtween = 1
        clcal.go()
        #######  Plot the CL table
        while uvdata.table_highver('AIPS PL') > 0:
             uvdata.zap_table('AIPS PL', 0)
        snplt = AIPSTask('SNPLT')
        snplt.indata = uvdata
        snplt.inext = 'CL'
        snplt.inver = 4
        snplt.dotv = -1
        snplt.nplots = 8
        snplt.optype = 'DELA'
        snplt.go()
        #######
        if os.path.exists('CL4.eps'):
             print("removing old CL4.eps")
             os.remove('CL4.eps')

        lwpla = AIPSTask('lwpla')
        lwpla.indata = uvdata
        lwpla.plver = 1
        lwpla.inver = uvdata.table_highver('AIPS PL')
        lwpla.outfile = 'PWD:CL4.eps'
        lwpla.go()
        
        ###### Plotting CL4 VPLOT
        while uvdata.table_highver('AIPS PL') > 0:
            uvdata.zap_table('AIPS PL', 0)
        
        vplot = AIPSTask('vplot')
        vplot.indata = uvdata
        vplot.solint = 0.166
        vplot.bparm[2] = -1
        #vplot.avgchan = 1
        #vplot.avgif = 1
        vplot.dotv = -1
        vplot.refant = refant
        #vplot.nplots = 12
        # Also want to make bchan echan more automated
        vplot.bchan = 1
        vplot.echan = (uvdata.header['naxis'][2])
        vplot.docal = 1
        vplot.gainuse = 4
        vplot.crowded = 3
        vplot.do3col = 1
        vplot.stokes = 'RR'
        vplot.go()
        vplot.stokes = 'LL'
        vplot.go()
        print("bchan echan is", vplot.bchan, vplot.echan)
        #######

        if os.path.exists('VPLOT_CL4.eps'):
             print("removing old VPLOT_CL4.eps")
             os.remove('VPLOT_CL4.eps')

        lwpla = AIPSTask('lwpla')
        lwpla.indata = uvdata
        lwpla.plver = 1
        lwpla.inver = uvdata.table_highver('AIPS PL')
        lwpla.outfile = 'PWD:VPLOT_CL4.eps'
        lwpla.go()
        while uvdata.table_highver('AIPS PL') > 0:
                uvdata.zap_table('AIPS PL', 0)

################ Do FRING (GET SINGLE SOL FOR ALL IFS)
if tmask <= 6:
        essentials.write("Reached tmask 6 \n")
        print("Tmaks 6")
        while uvdata.table_highver('AIPS SN') > 3:
                uvdata.zap_table('AIPS SN', 0)
        fring = AIPSTask('FRING')
        fring.indata = uvdata
        fring.calsour[1:] = conts
        fring.aparm = AIPSList([ 2,-1,0,0,1,0,3,0,0])
        fring.dparm = AIPSList([ 1,200,300,1,0,0,0,0])
        fring.doband = -1
        fring.refant = 1
        fring.docal = 1
        fring.gainuse = 4
        fring.solint = 1
        #fring.solsub = 4
        fring.go()
        #######
        while uvdata.table_highver('AIPS CL') > 4:
                uvdata.zap_table('AIPS CL', 0)
        clcal = AIPSTask('CLCAL')
        clcal.indata = uvdata
        clcal.snver = 4
        clcal.gainver = 4
        clcal.gainuse = 5
        #clcal.sources = AIPSList(['3C345'])
        clcal.interpol = 'AMBG'
        clcal.bparm = AIPSList([0,0,0,1])
        clcal.dobtween = 1
        clcal.go()
        
        ###### Plotting CL5 VPLOT
        vplot = AIPSTask('vplot')
        vplot.indata = uvdata
        vplot.solint = 0.166
        vplot.bparm[2] = -1
        #vplot.avgchan = 1
        #vplot.avgif = 1
        vplot.dotv = -1
        vplot.refant = refant
        #vplot.nplots = 12
        # Also want to make bchan echan more automated
        vplot.bchan = 1
        vplot.echan = (uvdata.header['naxis'][2])
        vplot.docal = 1
        vplot.gainuse = 5
        vplot.crowded = 3
        vplot.do3col = 1
        vplot.stokes = 'RR'
        vplot.go()
        vplot.stokes = 'LL'
        vplot.go()
        print("bchan echan is", vplot.bchan, vplot.echan)
        #######

        if os.path.exists('VPLOT_CL5.eps'):
             print("removing old VPLOT_CL5.eps")
             os.remove('VPLOT_CL5.eps')

        lwpla = AIPSTask('lwpla')
        lwpla.indata = uvdata
        lwpla.plver = 1
        lwpla.inver = uvdata.table_highver('AIPS PL')
        lwpla.outfile = 'PWD:VPLOT_CL5.eps'
        lwpla.go()
        while uvdata.table_highver('AIPS PL') > 0:
                uvdata.zap_table('AIPS PL', 0)

        
        
###################  FRING the maser for rate sols

if tmask <= 7:
    essentials.write("Reached tmask 7 \n")
    print("Tmaks 7")
    while uvdata.table_highver('AIPS CL') > 5:
        uvdata.zap_table('AIPS CL', 0)
    
    ########## Superarray
    while uvdata.table_highver('AIPS PL') > 0:
            uvdata.zap_table('AIPS PL', 0)
    for source in masers:
         possm = AIPSTask('possm')
         possm.indata = uvdata
         possm.solint = 666
         possm.sources = AIPSList([source])
         possm.nplots = 0
         possm.dotv = -1
         possm.doband = 1
         possm.docal = 1
         possm.gainuse = 5
         possm.aparm[7] = 2
         possm.aparm[1] = -1
         possm.outtext = './' + source + '.TXT'
         sa_spec       = './' + source + '.TXT'
         if os.path.exists(sa_spec):
              print("removing old plot", sa_spec)
              os.remove(sa_spec)
         possm.go()
         gnuplot = "gnuplot -e \"set term postscript eps; set out '" + source + "_SA.eps'; plot '" + sa_spec + "' u 5:6 w lines\""
         os.system(gnuplot)
         if os.path.exists('Spectra'):
             mv = "mv " + sa_spec + " Spectra"
             os.system (mv)
             mv = "mv " + source + "_SA.eps Spectra"
             os.system (mv)
         else:
             os.system('mkdir Spectra')
             mv = "mv " + sa_spec + " Spectra"
             os.system (mv)
             mv = "mv " + source + "_SA.eps Spectra"
             os.system (mv)

         #  for source in ['IRAS18056']:
         print("   we look for maschan")
         cmd = "sed -E -n '/LL|RR/p' Spectra/" + source + ".TXT | sort -k6 -n | tail -n 1"
         maspeak = subprocess.check_output(cmd, shell=True)
         print(cmd)
         print("Maspeak =")
         print(maspeak)
         mas = maspeak.split()
         maschan = int(mas[0])
         print("   Source and maser ch:", source, " ", maschan)
         #Write in the maser peak channels into Essentials.txt
         essentials.write("Maser peak channels are \n")
         essentials.write(source + " " + str(maschan))
         essentials.write( '\n')

         ####### Zoomed possm
         possm = AIPSTask('possm')
         possm.indata = uvdata
         possm.solint = 666
         possm.sources = AIPSList([source])
         possm.nplots = 0
         possm.dotv = -1
         possm.doband = -1
         possm.docal = 1
         possm.gainuse = 5
         possm.aparm[7] = 2
         possm.aparm[1] = -1
         possm.stokes = 'I'
         if maschan > 50:
             possm.bchan = maschan - 50
         else:
             possm.bchan = 1
         if maschan + 50 > uvdata.header['naxis'][2]:
             possm.echan = uvdata.header['naxis'][2]
         else:
             possm.echan = maschan + 50
         possm.bif = int(mas[1])
         possm.eif = int(mas[1])
         possm.outtext = './' + source + '.zoom.TXT'
         sa_spec       = './' + source + '.zoom.TXT'
         if os.path.exists(sa_spec):
             print("removing old plot", sa_spec)
             os.remove(sa_spec)
         possm.go()
         
         # Gnuplot the TXT spectra
         gnuplot = "gnuplot -e \"set term postscript eps; set out '" + source + "_SA_zoom.eps'; set xlabel 'Velocity [km/s]'; set ylabel 'Flux density [Jy]'; set yrange [-5:]; plot '" + sa_spec + "' u 5:6 w lines title 'IF " + str(possm.bif) + "' \""
         os.system(gnuplot)
         if os.path.exists('Spectra'):
             mv = "mv " + source + "_SA_zoom.eps Spectra"
             mv = "mv " + source + ".zoom.TXT Spectra"
             os.system (mv)
         else:
             os.system('mkdir Spectra')
             mv = "mv " + source + "_SA_zoom.eps Spectra"
             mv = "mv " + source + ".zoom.TXT Spectra"
             os.system (mv)
         ############# Maser FRING
         # Need to figure out how to deal with SN numbers here. Maybe make, SNCOP, use, delete
         while uvdata.table_highver('AIPS SN') > 4:
             uvdata.zap_table('AIPS SN', 0)
         fring = AIPSTask('FRING')
         fring.indata = uvdata
         fring.calsour = AIPSList([source])
         fring.aparm = AIPSList([ 2,0,3,0,0,0,1,0,0])
         fring.dparm = AIPSList([ 0,-1,200,0,0,0,0,0])
         fring.doband = -1
         fring.refant = 1
         fring.docal = 1
         fring.gainuse = 5
         fring.bif = int(mas[1])
         fring.eif = int(mas[1])
         fring.bchan = maschan
         fring.echan = maschan
         fring.solint = 0.333
         fring.go()
         ####### makes SN5
         
         ### Copy SN5 masIF sols to all IFs
         while uvdata.table_highver('AIPS SN') > 5:
                 uvdata.zap_table('AIPS SN', 0)
         clcop = AIPSTask('CLCOP')
         clcop.indata = uvdata
         clcop.inext = 'SN'
         clcop.invers = 5
         clcop.optype = 'AVER'
         clcop.fparm = AIPSList([ int(mas[1])])
         clcop.vparm = AIPSList([ 1,2,4,5,6,7,8])
         clcop.go()
         ### Makes SN6
         
         
         ### Plot SN5, kill PL files first
         while uvdata.table_highver('AIPS PL') > 0:
             uvdata.zap_table('AIPS PL', 0)
         
         if os.path.exists('SN5_' + source + '.eps'):
              print("removing old SN5.eps")
              os.remove('SN5_' + source + '.eps')

         #######
         snplt = AIPSTask('SNPLT')
         snplt.indata = uvdata
         snplt.inext = 'SN'
         snplt.inver = 6
         snplt.dotv = -1
         snplt.nplots = 4
         snplt.do3col = 1
         snplt.opcode = 'ALIF'
         snplt.optype = 'PHAS'
         snplt.go()
         #######

         lwpla = AIPSTask('lwpla')
         lwpla.indata = uvdata
         lwpla.plver = 1
         lwpla.inver = uvdata.table_highver('AIPS PL')
         lwpla.outfile = 'PWD:SN5_' + source + '.eps'
         lwpla.go()
         #### Finish plotting SN5 for each source
         
         
         clcal = AIPSTask('CLCAL')
         clcal.indata = uvdata
         clcal.snver = 6
         clcal.gainver = 5
         clcal.gainuse = 6
         clcal.sources = AIPSList([source])
         clcal.calsour = AIPSList([source])
         clcal.interpol = 'AMBG'
         clcal.go()
         

    ### Plot CL6, kill PL files first
    while uvdata.table_highver('AIPS PL') > 0:
        uvdata.zap_table('AIPS PL', 0)
    
    if os.path.exists('CL6.eps'):
         print("removing old CL6.eps")
         os.remove('CL6.eps')
    
    ### Of there are no masers there will be no CL6 at this point
    # at which case lets duplicate CL4 to CL6
    while uvdata.table_highver('AIPS CL') < 6:
        tacop = AIPSTask('TACOP')
        tacop.indata = uvdata
        tacop.outdata = uvdata
        tacop.inext = 'CL'
        tacop.invers = 5
        tacop.outver = 6
        tacop.go()

    #######
    snplt = AIPSTask('SNPLT')
    snplt.indata = uvdata
    snplt.inext = 'CL'
    snplt.inver = 6
    snplt.do3col = 1
    snplt.opcode = 'ALIF'
    snplt.dotv = -1
    snplt.nplots = 4
    snplt.optype = 'PHAS'
    snplt.go()
    #######

    lwpla = AIPSTask('lwpla')
    lwpla.indata = uvdata
    lwpla.plver = 1
    lwpla.inver = uvdata.table_highver('AIPS PL')
    lwpla.outfile = 'PWD:CL6.eps'
    lwpla.go()
    #### Finish plotting SN5 for each source
    
    while uvdata.table_highver('AIPS PL') > 0:
            uvdata.zap_table('AIPS PL', 0)
    vplot = AIPSTask('vplot')
    vplot.indata = uvdata
    vplot.solint = 0.166
    vplot.bparm[2] = -1
    #vplot.avgchan = 1
    #vplot.avgif = 1
    vplot.dotv = -1
    vplot.refant = refant
    #vplot.nplots = 12
    # Also want to make bchan echan more automated
    vplot.bchan = 1
    vplot.echan = (uvdata.header['naxis'][2])
    vplot.crowded = 3
    vplot.do3col = 1
    vplot.docal = 1
    vplot.gainuse = 6
    vplot.stokes = 'RR'
    vplot.go()
    vplot.stokes = 'LL'
    vplot.go()
    print("bchan echan is", vplot.bchan, vplot.echan)
    #######
    if os.path.exists('VPLOT_CL6.eps'):
         print("removing old VPLOT_CL6.eps")
         os.remove('VPLOT_CL6.eps')

    lwpla = AIPSTask('lwpla')
    lwpla.indata = uvdata
    lwpla.plver = 1
    lwpla.inver = uvdata.table_highver('AIPS PL')
    lwpla.outfile = 'PWD:VPLOT_CL6.eps'
    lwpla.go()
    
    #### Plot CL6 POSSM
    while uvdata.table_highver('AIPS PL') > 0:
            uvdata.zap_table('AIPS PL', 0)
    for source in uvdata.sources:
        print('doing source' + source)
        possm = AIPSTask('possm')
        possm.indata = uvdata
        possm.solint = 666
        possm.sources = AIPSList([source])
        possm.stokes = 'HALF'
        possm.aparm[1] = -1
        possm.flagver = 1
        possm.nplots = 2
        possm.dotv = -1
        possm.doband = 1
        possm.docal = 1
        possm.gainuse = 6
        possm.aparm[7] = 2
        possm.aparm[9] = 1
        possm.aparm[1] = -1
        # Fixed scale for phase
        possm.aparm[2] = 1
        possm.aparm[5] = -180
        possm.aparm[6] = 180
        possm.go()

        if os.path.exists('POSSM_CROSS_CL6.eps'):
             print("removing old POSSM_CROSS_CL6.eps")
             os.remove('POSSM_CROSS_CL6.eps')

        lwpla = AIPSTask('lwpla')
        lwpla.indata = uvdata
        lwpla.plver = 1
        lwpla.inver = uvdata.table_highver('AIPS PL')
        lwpla.outfile = 'PWD:POSSM_CROSS_CL6.eps'
        lwpla.go()


if tmask <= 8:
    ########## FRING all sources: long solint, IFs indivi
        essentials.write("Reached tmask 8 \n")
        while uvdata.table_highver('AIPS SN') > 6:
            uvdata.zap_table('AIPS SN', 0)
        fring = AIPSTask('FRING')
        fring.indata = uvdata
        #fring.calsour = AIPSList([ uvdata.sources])
        fring.aparm = AIPSList([ 2,-1,0,0,0,0,3,0,0])
        fring.dparm = AIPSList([ 1,200,300,1,0,0,0,1])
        fring.doband = -1
        fring.refant = 1
        fring.docal = 1
        fring.gainuse = 6
        #fring.bif = int(mas[1])
        #fring.eif = int(mas[1])
        #fring.bchan = maschan
        #fring.echan = maschan
        fring.solint = 10
        fring.go()
        ####### makes SN7
        while uvdata.table_highver('AIPS CL') > 6:
                uvdata.zap_table('AIPS CL', 0)

if tmask <= 9:
        essentials.write("Reached tmask 9 \n")
        clcal = AIPSTask('CLCAL')
        clcal.indata = uvdata
        clcal.snver = 7
        clcal.gainver = 6
        clcal.gainuse = 7
        #clcal.sources =
        #clcal.calsour = AIPSList([source])
        clcal.interpol = 'SELF'
        clcal.go()
        
        while uvdata.table_highver('AIPS PL') > 0:
                uvdata.zap_table('AIPS PL', 0)
        vplot = AIPSTask('vplot')
        vplot.indata = uvdata
        vplot.solint = 0.166
        vplot.bparm[2] = -1
        #vplot.avgchan = 1
        #vplot.avgif = 1
        vplot.dotv = -1
        vplot.refant = refant
        #vplot.nplots = 12
        # Also want to make bchan echan more automated
        vplot.bchan = 1
        vplot.echan = (uvdata.header['naxis'][2])
        vplot.crowded = 3
        vplot.do3col = 1
        vplot.docal = 1
        vplot.gainuse = 7
        vplot.stokes = 'RR'
        vplot.go()
        vplot.stokes = 'LL'
        vplot.go()
        print("bchan echan is", vplot.bchan, vplot.echan)
        #######
        if os.path.exists('VPLOT_CL7.eps'):
             print("removing old VPLOT_CL7.eps")
             os.remove('VPLOT_CL7.eps')

        lwpla = AIPSTask('lwpla')
        lwpla.indata = uvdata
        lwpla.plver = 1
        lwpla.inver = uvdata.table_highver('AIPS PL')
        lwpla.outfile = 'PWD:VPLOT_CL7.eps'
        lwpla.go()
        
        #### Plot CL6 POSSM
        while uvdata.table_highver('AIPS PL') > 0:
                uvdata.zap_table('AIPS PL', 0)
        for source in uvdata.sources:
            possm = AIPSTask('possm')
            possm.indata = uvdata
            possm.solint = 666
            possm.sources = AIPSList([source])
            possm.stokes = 'HALF'
            possm.aparm[1] = -1
            possm.flagver = 1
            possm.nplots = 2
            possm.dotv = -1
            possm.doband = 1
            possm.docal = 1
            possm.gainuse = 7
            possm.aparm[7] = 2
            possm.aparm[9] = 1
            possm.aparm[1] = -1
            # Fixed scale for phase
            possm.aparm[2] = 1
            possm.aparm[5] = -180
            possm.aparm[6] = 180
            possm.go()

            if os.path.exists('POSSM_CROSS_CL7.eps'):
                 print("removing old POSSM_CROSS_CL7.eps")
                 os.remove('POSSM_CROSS_CL7.eps')

            lwpla = AIPSTask('lwpla')
            lwpla.indata = uvdata
            lwpla.plver = 1
            lwpla.inver = uvdata.table_highver('AIPS PL')
            lwpla.outfile = 'PWD:POSSM_CROSS_CL7.eps'
            lwpla.go()

if tmask <= 10:
    essentials.write("Reached tmask 10 \n")
    for source in uvdata.sources:
        if os.path.exists('UVFIT_' + source + '.txt'):
            print("  Removing old UVFIT results text file")
            os.remove('UVFIT_' + source + '.txt')
        uvfit = AIPSTask('uvfit')
        uvfit.indata = uvdata
        uvfit.source = AIPSList([source])
        uvfit.docal = 1
        uvfit.gainuse = 7
        uvfit.stokes = 'I'
        uvfit.prtlev = 1
        #uvfit.dopos = -1
        uvfit.solmode = ''
        uvfit.fitout = 'PWD:UVFIT_' + source + '.txt'
        uvfit.go()



#### Concatenate all the CL VPLOTS into one
if os.path.exists('VPLOTS_ALL.pdf'):
    print("removing old VPLOTS_ALL.pdf")
    os.remove('VPLOTS_ALL.pdf')

gs = "gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=VPLOTS_ALL.pdf VPLOT*.eps"

os.system(gs)
###########################################

#### Concatenate all the CL POSSMS into one
if os.path.exists('POSSM_CROSS_ALL.pdf'):
    print("removing old POSSM_CROSS_ALL.pdf")
    os.remove('POSSM_CROSS_ALL.pdf')

gs = "gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=POSSM_CROSS_ALL.pdf POSSM_CROSS*.eps"

os.system(gs)
###########################################



############ extract just the flux values from the results for essentials
essentials.write("Final Flux and error values are \n")
essentials.close()

awk = "awk '$3 ~ /\*/ {print FILENAME, $5,$6} FNR==2 && $3 !~ /\*/ {print FILENAME, $7,$8}' ./UVFIT_* cat >> essentials.txt"
os.system(awk)

