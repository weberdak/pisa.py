#!/usr/bin/env python

VERSION='1.0'
REVISION='Nov 20 2019'


import pandas as pd
import numpy as np
import math
import itertools
import argparse
import time
import random
import multiprocessing
import concurrent.futures


def initData():
    df = pd.DataFrame({'resname': [], 'rho': [], 'shift_exp': [], 
                      'coupling_exp': [], 'shift_calc': [], 'coupling_calc': [], 
                      'score': [], 'shift_raw': [], 'coupling_raw': []})
    return df


def initFromSequence(df,seq,nres,seq_start,rho_start):
    '''Initialize residue names into dataframe by input sequence.
    
    Parameters
    ----------
    df: Pandas dataframe
        An existing data set.
    seq: str
        Amino acid sequence.
    nres: int
        Number of residues. Ignored if seq used.
    seq_start: int
        Residue number of first amino acid.
    rho_start: int
        Residue of rho0.
        
    Returns
    -------
    updated pandas dataframe
    '''
    if seq:
        nres=len(seq)
        seq = list(seq)
    indices = [ i+seq_start-rho_start for i in range(nres) ]
    resnums = [ i+seq_start for i in  range(nres) ]
    if seq:
        resnames = [ s+str(i) for s,i in zip(seq,resnums) ]
    else:
        resnames = [ str(i) for i in resnums ]
    dict_temp = {'resname': pd.Series(resnames, index=indices)}
    return pd.concat([df,pd.DataFrame(dict_temp)])


def initFromFile(df,infile,rho_start):
    '''Initialize residue names and experimental data from file.
    
    Parameters
    ---------
    df: Pandas dataframe
        An existing data set.
    infile: TSV file
        Tab-delimited file of format: assignment, shifts (ppm),  couplings (kHz).
        Residues without a residue number will be ignored.
    rho_start: int
        Residue of rho0.
        
    Returns
    -------
    df: Pandas dataframe
        Updated dataframe.
    '''
    tnames = np.genfromtxt(infile,usecols=(0),dtype=str)
    tshifts = np.genfromtxt(infile,usecols=(1),dtype=float)
    tcouplings = np.genfromtxt(infile,usecols=(2),dtype=float)
    # Read names, shifts and couplings into dataframe. Handle errrors if only 1 assignment.
    try:
        for name,shift,coupling in zip(tnames,tshifts,tcouplings):
            s, i = resSplit(name)
            if i:
                index = i-rho_start
                df.at[index,'resname'] = s+str(i)
                df.at[index,'shift_exp'] = shift
                df.at[index,'coupling_exp'] = np.absolute(coupling)
    except TypeError:
        s, i = resSplit(str(tnames))
        if i:
            index = i-rho_start
            df.at[index,'resname'] = s+str(i)
            df.at[index,'shift_exp'] = tshifts
            df.at[index,'coupling_exp'] = np.absolute(tcouplings)
    return df


def filterData(df,resnums,rho_start,exclude=False,include=False):
    '''Remove rows belonging to specified residue numbers.

    Parameters
    ----------
    df: Pandas dataframe
       An existing data set.
    resnums: list of ints
       List of residue numbers to strip or keep.
    rho_start: int
       This residue number is indexed to zero in PISA calculations.
    exclude: bool
       Remove these residue numbers from dataframe
    include: bool
       Keep only these residues numbers in the dataframe
    
    Returns
    -------
    df: dataframe
       Dataframe without unwanted residues.
    '''
    if exclude:
        remove = [ i-rho_start for i in resnums ]
    if include:
        keep = [ i-rho_start for i in resnums ]
        print(keep)
        remove = []
        for i in df.index.values:
            if i not in keep:
                remove.append(i)
    return df.drop(remove)


def calcShiftCoupling(tau,rho,order,flip,pas,beta,dc,matrixA):
    '''Compute the 15N shift and 15N-1H dipolar coupling for residue.
    
    Parameters
    ----------
    tau: float
        Helix tilt angle (radians).
    rho: float
        Azimuthal angle of this residue (radians).
    order: float
        Order parameter (0 to 1).
    flip: float
        Flip angle of membrane normal to magnetic field (radians).
    pas: list of three floats
        15N principal axis system (ppm).
    beta: float
        Second polar angle relating the 15N PAS to NH dipolar vector.
    dc: float
        Full 15N-1H dipolar coupling.
    matrixA: numpy array
        Matrix relating peptide structure to helical axis.

    Returns
    -------
    N, NH: list of floats
        15N chemical shift (ppm) and 15N-1H dipolar coupling (kHz)
    '''
    iso = np.mean(pas)
    X = np.array([[ np.cos(rho)*np.sin(tau)], 
         [ np.sin(rho)*np.sin(tau)], 
         [ np.cos(tau)]])
    x,y,z = np.matmul(matrixA,X); # Generate PAF coordinates for rho angle
    n = pas[0]*x**2 + pas[1]*y**2 + pas[2]*z**2 
    n = (n - iso) * order * 0.5 * (3*(math.cos(flip))**2-1) + iso
    nh = dc * 0.5 * (3*(math.sin(beta)*x + math.cos(beta)*z)**2 - 1)
    nh = nh * order * 0.5 * (3 * (math.cos(flip))**2-1)
    return n[0],nh[0]


def calcMatrixA(phi,psi,beta,aCaCN,aCNCa,aNCaC,aCaNH,bCaC,bCN,bNCa):
    '''Derive Matrix A relating helical axis frame (HAF) to
    15N principal axis frame (PAF).
    
    Parameters
    ----------
    phi: float
        C(i-1),N(i),CA(i),C(i) dihedral angle (radians)
    psi: float
        N(i),CA(i),C(i),N(i+1) dihedral angle (radians)
    beta: float
        As per calcShiftCoupling().
    aCaCN: float
        CA,C,N bond angle (radians)
    aCNCa: float
        C,N,CA bond angle (radians)
    aNCaC: float
        N,CA,C bond angle (radians)
    aCaNH: float
        CA,N,NH bond angle (radians)
    bCaC: float
        CA,C bond length (angstroms)
    bCN: float
        C,N bond length (angstroms)
    bNCa: float
        N,Ca bond length (angstroms)
    
    Returns
    -------
    matrixA: numpy array
    '''
    omega = 180 * 0.017453
    eAlpha = (180 * 0.017453) - aCaCN
    eBeta = (180 * 0.017453) - aCNCa
    eGamma = (180 * 0.017453) - aNCaC
    v = np.array([ bCN*np.cos(eBeta) + bCaC*np.cos(eAlpha-eBeta) + bNCa,
                   bCN*-np.sin(eBeta) + bCaC*np.sin(eAlpha-eBeta),
                   0 ])
    R1phi = np.array([[ 1, 0, 0],
            [ 0, np.cos(phi), -np.sin(phi) ],
            [ 0, np.sin(phi), np.cos(phi) ]])
    R3gamma = np.array([[ np.cos(eGamma), -np.sin(eGamma), 0],
            [ np.sin(eGamma), np.cos(eGamma), 0 ],
            [ 0, 0, 1 ]])
    R1psi = np.array([[ 1, 0, 0],
            [ 0, np.cos(psi), -np.sin(psi) ],
            [ 0, np.sin(psi), np.cos(psi) ]])
    R3alpha = np.array([[ np.cos(eAlpha), -np.sin(eAlpha), 0],
            [ np.sin(eAlpha), np.cos(eAlpha), 0 ],
            [ 0, 0, 1 ]])
    R1omega = np.array([[ 1, 0, 0],
            [ 0, np.cos(omega), -np.sin(omega) ],
            [ 0, np.sin(omega), np.cos(omega) ]])
    R3beta = np.array([[ np.cos(eBeta), -np.sin(eBeta), 0],
               [ np.sin(eBeta), np.cos(eBeta), 0 ],
               [ 0, 0, 1 ]])
    R3 = np.array([[ np.cos(beta + aCaNH), -np.sin(beta + aCaNH), 0 ],
          [ np.sin(beta + aCaNH),  np.cos(beta + aCaNH), 0 ],
          [ 0, 0, 1]])
    C = R1phi @ R3gamma @ R1psi @ R3alpha @ R1omega @ R3beta
    a1 = C[2][1] - C[1][2]
    a2 = C[0][2] - C[2][0]
    a3 = C[1][0] - C[0][1]
    a = np.array([a1 / np.sqrt(a1**2+a2**2+a3**2),
                  a2 / np.sqrt(a1**2+a2**2+a3**2),
                  a3 / np.sqrt(a1**2+a2**2+a3**2)])
    p = 0.5 * ((1/np.tan(50*0.017453)) * np.cross(a,v) + v - np.dot(a,v) * a )
    r = np.array([-p[0]/np.sqrt(p[0]**2+p[1]**2+p[2]**2),
                  -p[1]/np.sqrt(p[0]**2+p[1]**2+p[2]**2),
                  -p[2]/np.sqrt(p[0]**2+p[1]**2+p[2]**2)])
    HAF = np.array([r,np.cross(a,r),a]).T
    return np.array([[0,1,0],[0,0,1],[1,0,0]]) @ R3 @ HAF


def calcWheel(df,tau,rho0,period,order,flip,pas,beta,dc,matrixA):
    '''Fill dataframe with calculated rhos, 15N chemical shifts and 15N-1H dipolar couplings.
    
    Parameters
    ----------
    df: Pandas dataframe
        An initialized dataframe for the sequence.
    rho0: float
        Azimuthal angle of reference residue (radians).
    period: float
        Number of residues per turn of helix.
    tau, order, flip, pas, beta, dc, matrixA
        See calcShiftCoupling()
        
    Returns
    -------
    df: Pandas dataframe
        An updated dataframe.
    '''
    for index in df.index.values:
        rho = (rho0-((360/period)*index) % 360)
        if rho < 0:
            rho = rho + 360
        n,nh = calcShiftCoupling(tau,rho*0.017453,order,flip,pas,beta,dc,matrixA)
        df.at[index,'rho'] = rho
        df.at[index,'shift_calc'] = n
        df.at[index,'coupling_calc'] = np.absolute(nh)
    return df


def calcWave(df,tau,rho0,period,order,flip,pas,beta,dc,rho_start,matrixA,margin=5.0,increment=0.1):
    '''Compute continuous (interpolated) wheel of PISA data points.

    Parameters
    ----------
    df, rho0, period, order, flip, pas, beta, dc, matrixA: See calcWheel()
    rho_start: int
       Residue number that rho is referenced to zero.
    margin: float
       Compute values this number of residues beyond data.
    increment: float
       Increment rho values by this for each point on wheel.

    Returns
    -------
    resnums: list of floats
       Interpolated list of residue numbers.
    rhos: list of floats
       Interpolated list of rho indices.
    shifts: list of floats
       Calculated chemical shifts for rhos.
    couplings: list of floats
       Calculated dipolar couplings for rhos.
    '''
    low_index = min(df.index.values) - margin
    high_index = max(df.index.values) + margin
    resnums = []
    rhos = []
    shifts = []
    couplings = []
    for i in np.arange(low_index,high_index+increment,increment):
        rho = (rho0-((360/period)*i) % 360)
        if rho < 0:
            rho = rho + 360
        n,nh = calcShiftCoupling(tau,rho*0.017453,order,flip,pas,beta,dc,matrixA)
        resnums.append(i+rho_start)
        rhos.append(rho)
        shifts.append(n)
        couplings.append(np.absolute(nh))
    return resnums,rhos,shifts,couplings
        

def calcScores(df,scalar):
    '''Fill dataframe with RMSD fit scores.
    
    Parameters
    ----------
    df: Pandas dataframe
        An initialized dataframe for sequence.
    scalar: float
        Factor to scale up dipolar couplings to match 15N shift ppm scale.
    
    Returns
    -------
    df: Pandas dataframe
        An updated dataframe.
    score: float
        Total fit score (sum of all residues).
    '''
    score = 0
    for index in df.index.values:
        if pd.isnull(df.at[index,'shift_exp']):
            n_dev = 0
        else:
            n_dev = np.sqrt((df.at[index,'shift_exp'] - df.at[index,'shift_calc'])**2)
        if pd.isnull(df.at[index,'coupling_exp']):
            nh_dev = 0
        else:
            nh_dev = scalar * np.sqrt((df.at[index,'coupling_exp'] - df.at[index,'coupling_calc'])**2)
        sub_score = n_dev + nh_dev
        df.at[index,'score'] = sub_score
        score += sub_score
    return df,score


def randomizeShiftsCouplings(df,lw_shifts,lw_couplings):
    '''Add random values to experimental shifts and couplings in dataframe.
    
    Parameters
    ----------
    df: Pandas dataframe
        An initialized dataframe for sequence.
    lw_shifts: float
        Average linewidths of chemical shifts in dataframe (in PPM).
    lw_couplings: float
        Average linewidths of dipolar couplings in dataframe (in kHz)
    
    Returns
    -------
    df: Pandas dataframe
        A dataframe with modified shifts and couplings.
    '''
    for index in df.index.values:
        if not pd.isnull(df.at[index,'shift_exp']):
            df.at[index,'shift_exp'] = df.at[index,'shift_exp'] + random.gauss(0, lw_shifts)
        if not pd.isnull(df.at[index,'coupling_exp']):
            df.at[index,'coupling_exp'] = df.at[index,'coupling_exp'] + random.gauss(0, lw_couplings)
    return df


def explore_comb(df,period,flip,pas,beta,dc,matrixA,scalar,comb,verbose=True):
    '''Explore combinations of tau, rho0 and order for lowest fit score.
    
    Parameters
    ----------
    df: Pandas dataframe
        An initialized dataframe for sequence.
    comb: list of list of floats
        Combinations of [tau,rho0,order] to be tested.
    verbose: Bool
        Report progress.
    flip, pas, beta, dc, matrixA, scalar
        See calcShiftCoupling()
    
    Returns
    -------
    results: list of list of floats
        [tau, rho, order, score] for all combinations tested.
    '''
    ncomb = len(comb)
    c = 0
    results = []
    for t,r,o in comb:
        sim = calcWheel(df,t*0.017453,r,period,o,flip,pas,beta,dc,matrixA)
        df_t, score = calcScores(sim,scalar)
        results.append((t,r,o,score))
        if verbose:
            c+=1
            print('# Progress {0:.2f}%'.format(100*c/ncomb), end='\r')
    if verbose:
        print('')
    return results


def quickfit_comb(minShift,maxShift,minCouple,maxCouple,rho0,period,flip,pas,beta,dc,matrixA,scalar,comb,verbose=True):
    '''Explore combinations of tau and oder for best fit to data. Ignores assignments.
    
    Parameters
    ----------
    minShift: float
        Minimum experimental 15N shift (ppm).
    maxShift: float
        Maximum experimental 15N shift (ppm).
    minCoupling: float
        Minimum experimental 15N-1H dipolar coupling (kHz).
    maxCoupling: float
        Maximum experimental 15N-1H dipolar coupling (kHz).
    rho0, period, flip, pas, beta, dc, matrixA, scalar
        See calcShiftCoupling()
    
    Returns
    -------
    results: list of list of floats
        [tau, order, score] for all combinations tested.
    '''
    ncomb = len(comb)
    df = initData()
    df = initFromSequence(df,'',18,1,1)
    c = 0
    results = []
    for t,o in comb:
        df_temp = calcWheel(df,t*0.017453,rho0,period,o,flip,pas,beta,dc,matrixA)
        maxShift_calc = max(df_temp['shift_calc'])
        minShift_calc = min(df_temp['shift_calc'])
        maxCouple_calc = max(df_temp['coupling_calc'])
        minCouple_calc = min(df_temp['coupling_calc'])
        if minShift:
            n_dev = np.sqrt((maxShift-maxShift_calc)**2 + (minShift-minShift_calc)**2)
        else:
            n_dev = 0
        if minCouple:
            nh_dev = scalar*np.sqrt((maxCouple-maxCouple_calc)**2 + (minCouple-minCouple_calc)**2)
        else:
            nh_dev = 0
        score = n_dev + nh_dev
        results.append((t,o,score))
        if verbose:
            c+=1
            print('# Progress {0:.2f}%'.format(100*c/ncomb), end='\r')
    if verbose:
        print('')
    return results



### THREADING

def chunks(l, n):
    """Yield successive n-sized chunks from l. From stackoverflow 312443"""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def explore_comb_multiproc(df,period,flip,pas,beta,dc,matrixA,scalar,comb,procs=2):
    """Parallel version of explore_comb.

    Parameters
    ----------
    df, period, flip, pas, beta, dc, matrixA, scalar, comb
       See explore_comb()
    processors: int
       Number of CPUs to use
    """
    ncomb = len(comb)
    chunk_size = int(ncomb/(procs))
    ncomb_chunks = chunks(comb,chunk_size)
    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=procs) as executor:
        future_to_result = { executor.submit(explore_comb,df,period,flip,pas,beta,dc,matrixA,scalar,l,verbose=False): l for l in ncomb_chunks }
        for future in concurrent.futures.as_completed(future_to_result):
            results = results + future.result()
            complete = len(results)
            print('# Progress {0:.2f}%'.format(100*complete/ncomb), end='\r')
    return results


def quickfit_comb_multiproc(minShift,maxShift,minCouple,maxCouple,rho0,period,flip,pas,beta,dc,matrixA,scalar,comb,procs=2):
    """Parallel version of quickfit_comb.

    Parameters
    ----------
    df, period, flip, pas, beta, dc, matrixA, scalar, comb
       See explore_comb()
    processors: int
       Number of CPUs to use
    """
    ncomb = len(comb)
    chunk_size = int(ncomb/(procs))
    ncomb_chunks = chunks(comb,chunk_size)
    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=procs) as executor:
        future_to_result = { executor.submit(quickfit_comb,minShift,maxShift,minCouple,maxCouple,rho0,period,flip,pas,beta,dc,matrixA,scalar,l,verbose=False): l for l in ncomb_chunks }
        for future in concurrent.futures.as_completed(future_to_result):
            results = results + future.result()
            complete = len(results)
            print('# Progress {0:.2f}%'.format(100*complete/ncomb), end='\r')
    return results

### END THREADING


def calcRaw(df,order,flip,pas):
    '''Fill dataframe with experimental shifts and couplings that reverse effects
    of disorder and non-zero flip angles.
    
    Parameters
    ----------
    df: Pandas dataframe
        An initialized dataframe for sequence.
    order, flip, pas
        See calcShiftCoupling()
        
    Returns
    -------
    df: Pandas dataframe
        An updated dataframe.
    '''
    iso = np.mean(pas)
    for index in df.index.values:
        df.at[index,'shift_raw'] = (df.at[index,'shift_exp']-iso)/order/0.5/(3*(math.cos(flip))**2-1)+iso
        df.at[index,'coupling_raw'] = np.absolute((df.at[index,'coupling_exp']/order/0.5/(3*(math.cos(flip))**2-1)))
    return df


def log(text,f):
    '''Very simple logger. Writes to file object and terminal.

    Parameters
    ----------
    text: str
       Text to write to file.
    f: file object
       File to write lof text to.
    '''
    f.write(text+'\n')
    print(text)

    
def resSplit(text):
    '''Find integers in text and concat. For processing residue names.

    Parameters
    ---------
    text: str
       String of text and numbers.
    
    Returns
    -------
    strs, ints: str, int
        Strings and integers concatenated together.
    '''
    strs=''
    ints=''
    for c in text:
        if is_int(c):
            ints+=c
        else:
            strs+=c
    return strs,int(ints)


def is_int(s):
    '''Return True if string is an integer.
    '''
    try: 
        int(s)
        return True
    except ValueError:
        return False


def is_float(s):
    '''Return True is string is float.
    '''
    try: 
        float(s)
        return True
    except ValueError:
        return False


def parse_args():
    parser = argparse.ArgumentParser(description='A python tool for simulating and fitting PISA wheels to oriented solid-state NMR data.',
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        '-t','--tau', type=float,
        help='Helical tilt angle in degrees (20.0).',
        default=20.0
    )
    parser.add_argument(
        '-r','--rho0', type=float,
        help='Azimuthal angle of reference residue (0.0).',
        default=0.0
    )
    parser.add_argument(
        '-o','--order', type=float,
        help='Order parameter (1.0).',
        default=1.0
    )
    parser.add_argument(
        '--phi', type=float,
        help='Phi dihedral angle for an ideal helix in degrees (-63.0).',
        default=-63.0
    )
    parser.add_argument(
        '--psi', type=float,
        help='Psi dihedral angle for an ideal helix in degrees (-42.0).',
        default=-42.0
    )
    parser.add_argument(
        '--beta', type=float,
        help='Beta rotation of the 15N PAS with respect to the amide plane in degrees (17.0).',
        default=17.0
    )
    parser.add_argument(
        '--dc', type=float,
        help='Maximum 15N-1H dipolar coupling in kHz (10.735).',
        default=10.735
    )
    parser.add_argument(
        '--pas', type=float, nargs='+',
        help='N15 principal axis system components in ppm (57.3 81.2 228.1).',
        default=(57.3,81.2,228.1)
    )
    parser.add_argument(
        '-s','--seq', type=str,
        help='Single-letter amino acid sequence of the protein (none).',
        default=''
    )
    parser.add_argument(
        '--seq_start', type=int,
        help='Starting residue number in the sequence (1).',
        default=1
    )
    parser.add_argument(
        '--rho_start', type=int,
        help='Residue number to reference as zero in all PISA calculations. I.e., residue of rho0 angle (1).',
        default=1
    )
    parser.add_argument(
        '-n','--nres', type=int,
        help='Number of residues to simulate (18 or determined from --seq)',
        default=18
    )
    parser.add_argument(
        '--flip', type=float,
        help='Angle of membrane normal to magnetic field in degrees (0.0).',
        default=0.0
    )
    parser.add_argument(
        '--period', type=float,
        help='Number of residues per helical rotation (3.6).',
        default=3.6
    )
    parser.add_argument(
        '--negative_dc',
        help='Output all dipolar couplings as negative values in output files.',
        action='store_true'
    )
    parser.add_argument(
        '--aCaCN', type=float,
        help='Angle of Ca-C-N in degrees (117.5).',
        default=117.5
    )
    parser.add_argument(
        '--aCNCa', type=float,
        help='Angle of C-N-Ca in degrees (124.0).',
        default=124.0
    )
    parser.add_argument(
        '--aNCaC', type=float,
        help='Angle of N-Ca-C in degrees (107.4).',
        default=107.4
    )
    parser.add_argument(
        '--aCaNH', type=float,
        help='Angle of Ca-N-H in degrees (116.0).',
        default=116.0
    )
    parser.add_argument(
        '--bCaC', type=float,
        help='Length of Ca-C bond in angstroms (1.52).',
        default=1.52
    )
    parser.add_argument(
        '--bCN', type=float,
        help='Length of C-N bond in angstroms (1.35).',
        default=1.35
    )
    parser.add_argument(
        '--bNCa', type=float,
        help='Length of N-Ca bond in angstroms (1.45).',
        default=1.45
    )
    parser.add_argument(
        '-f','--peaks', type=str,
        help='Input peak list file.',
        default=''
    )
    parser.add_argument(
        '--quickfit', nargs='+',
        help='Min. CS, max. CS, min. DC, max. DC. Explore combinations of tau angles and order parameters to best fit these bounds.',
        default=[]
    )
    parser.add_argument(
        '--explore',
        help='Specify to explore combinations of tau, rho0 and order to best fit input peak list.',
        action='store_true'
    )
    parser.add_argument(
        '--fit_tau', type=float, nargs='+',
        help='Min., max. and step to fit the tilt angle (0.0 90.0 1.0).',
        default=(0.0,90.0,1.0)
    )
    parser.add_argument(
        '--fit_rho0', type=float, nargs='+',
        help='Min., max. and step to fit the azimuthal angle (0.0 360.0 4.0).',
        default=(0.0,360.0,4.0)
    )
    parser.add_argument(
        '--fit_order', type=float, nargs='+',
        help='Min., max. and step to fit the order parameter (0.85 0.85 0.1).',
        default=(0.85, 0.85, 0.1)
    )
    parser.add_argument(
        '--scalar', type=float,
        help='Value to scale up dipolar couplings to match chemical shift dispersion (10.0)',
        default=10.0
    )
    parser.add_argument(
        '--fit_only', type=int, nargs='+',
        help='Only fit wheel to these residues (none).',
        default=[]
    )
    parser.add_argument(
        '--fit_exclude', type=int, nargs='+',
        help='Exclude these residues from fitting (none).',
        default=[]
    )
    parser.add_argument(
        '--show_only', type=int, nargs='+',
        help='Display only these residues in the output files (none).',
        default=[]
    )
    parser.add_argument(
        '--out_log', type=str,
        help='File name/path of output log file (pisa_log.dat).',
        default='pisa_log.dat'
    )
    parser.add_argument(
        '--out_wave', type=str,
        help='File name/path of output wave file (pisa_wave.dat).',
        default='pisa_wave.dat'
    )
    parser.add_argument(
        '--out_fit', type=str,
        help='File name/path of output fit file (pisa_fit.dat).',
        default='pisa_fit.dat'
    )
    parser.add_argument(
        '--errors', type=float, nargs='+',
        help='Peroform error analysis for explore function. Specify three numbers: number of replicates (int),  avg. linewidths in CS (ppm) and DC (kHz) dimensions.',
        default=[]
    )
    parser.add_argument(
        '--procs', type=int,
        help='Number of CPUs to use in explore and quickfit functions (32).',
        default=32
    )
    args = parser.parse_args()
    return args


def main():
    
    # Handle simulation arguments
    args = parse_args()
    tau = args.tau
    rho0 = args.rho0
    order = args.order
    flip = args.flip
    dc = args.dc
    pas = args.pas
    negative_dc = args.negative_dc

    # Handle fitting arguments
    procs = args.procs
    fit_tau = args.fit_tau
    fit_rho0 = args.fit_rho0
    fit_order = args.fit_order
    scalar = args.scalar
    explore = args.explore
    quickfit = args.quickfit
    fit_only = args.fit_only
    fit_exclude = args.fit_exclude
    show_only = args.show_only
    # TO DO: Fix this part to handle mis-entered values better.
    if len(args.errors) == 3:
        errors = True
        num_replicates = int(args.errors[0])
        lw_shifts = args.errors[1]
        lw_couplings = args.errors[2]
    
    # Handle sequence/data input arguments
    seq = args.seq
    nres = args.nres
    seq_start = args.seq_start
    rho_start = args.rho_start
    peaks = args.peaks

    # Handle output arguments
    out_log = args.out_log
    out_wave = args.out_wave
    out_fit = args.out_fit
    
    # Handle geometric arguments
    aCaCN = args.aCaCN
    aCNCa = args.aCNCa
    aNCaC = args.aNCaC
    aCaNH = args.aCaNH
    bCaC = args.bCaC
    bCN = args.bCN
    bNCa = args.bNCa
    phi = args.phi
    psi = args.psi
    beta = args.beta
    period = args.period
    
    # Initialize log file and log parameters
    f = open(out_log, 'w')
    log('# Job executed using pisa.py',f)
    log('# Version: {}'.format(VERSION),f)
    log('# Last revision: {}'.format(REVISION),f)
    log('# The following geometric parameters will be used:',f)
    log('# --aCaCN {}'.format(aCaCN),f)
    log('# --aCNCa {}'.format(aCNCa),f)
    log('# --aNCaC {}'.format(aNCaC),f)
    log('# --aCaNH {}'.format(aCaNH),f)
    log('# --bCaC {}'.format(bCaC),f)
    log('# --bCN {}'.format(bCN),f)
    log('# --bNCa {}'.format(bNCa),f)
    log('# --phi {}'.format(phi),f)
    log('# --psi {}'.format(psi),f)
    log('# --beta {}'.format(beta),f)
    log('# --period {}'.format(period),f)
    
    # Convert degrees to radians
    tau = tau * 0.017453
    flip = flip * 0.017453
    aCaCN = aCaCN * 0.017453
    aCNCa = aCNCa * 0.017453
    aNCaC = aNCaC * 0.017453
    aCaNH = aCaNH * 0.017453
    phi = phi * 0.017453
    psi = psi * 0.017453
    beta = beta * 0.017453

    # Initialize sequence/data
    df = initData()
    
    # Import sequence from text
    if seq:
        log('# Importing sequence: {}'.format(seq),f)
        log('# Sequence starts at residue: {}'.format(seq_start),f)
        log('# Rho indices referenced to residue: {}'.format(rho_start),f)
        df = initFromSequence(df,seq,nres,seq_start,rho_start)

    # Import sequence, shifts and couplings from file
    if peaks:
        log('# Importing sequence and data from file: {}'.format(peaks),f)
        log('# Rho indices referenced to residue: {}'.format(rho_start),f)            
        df = initFromFile(df,peaks,rho_start)

    # Generate integer sequence if neither seq nor peaks file specified
    if not seq and not peaks:
        log('# Neither --seq or --peaks specified. Generating {} residues.'.format(nres),f)
        log('# Sequence starts at residue: {}'.format(seq_start),f)
        log('# Rho indices referenced to residue: {}'.format(rho_start),f)
        df = initFromSequence(df,seq,nres,seq_start,rho_start)

    # Generate Matrix A
    matrixA = calcMatrixA(phi,psi,beta,aCaCN,aCNCa,aNCaC,aCaNH,bCaC,bCN,bNCa)
    log('# Computed Matrix A as:',f)
    for row in matrixA:
        log('# {}'.format(row),f)

    # Explore combinations of tau, rho0 and order for best fit to data
    if explore:
        assert (not quickfit), 'Options --explore and --quickfit cannot be specified together!'
        taus = np.arange(fit_tau[0],fit_tau[1]+fit_tau[2],fit_tau[2])
        rhos = np.arange(fit_rho0[0],fit_rho0[1]+fit_rho0[2],fit_rho0[2])
        ords = np.arange(fit_order[0],fit_order[1]+fit_order[2],fit_order[2])
        comb = list(itertools.product(*[taus,rhos,ords]))
        ncomb = len(comb)           
        log('# Exploring combinations of tau, rho0 and order for best fit to peak list.',f)
        log('# --fit_tau {}'.format(fit_tau),f)
        log('# --fit_rho0 {}'.format(fit_rho0),f)
        log('# --fit_order {}'.format(fit_order),f)
        log('# {} values of tau will be tested.'.format(len(taus)),f)
        log('# {} values of rho0 will be tested.'.format(len(rhos)),f)
        log('# {} values of order parameters be tested.'.format(len(ords)),f)
        log('# A total of {} combinations will be tested.'.format(ncomb),f)
        log('# {} CPUs will be used.'.format(procs),f)
        
        if fit_exclude:
            assert (not fit_only), 'Cannot use --fit_exclude and --fit_only flags together!'
            log('# --fit_exclude {}'.format(fit_exclude),f)
            df_fit = filterData(df,fit_exclude,rho_start,exclude=True,include=False)
        elif fit_only:
            assert (not fit_exclude), 'Cannot use --fit_exclude and --fit_only flags together!'
            log('# --fit_only {}'.format(fit_only),f)
            df_fit = filterData(df,fit_only,rho_start,exclude=False,include=True)
        else:
            df_fit = df

        # Report which residues will be used for scoring. Exclude residues missing assignments.
        log_residues = []
        for index in df_fit.index.values:
            check = 0
            if not pd.isnull(df_fit.at[index,'shift_exp']):
                 check+=1
            if not pd.isnull(df_fit.at[index,'coupling_exp']):
                 check+=1
            if check > 0:
                log_residues.append(df_fit.at[index,'resname'])
        log('# Residues to fit: {}'.format(log_residues),f)

        time_start = time.time()
        # Perform error analysis if specified.
        if args.errors:
            log('# Performing error analysis using {} replicates with {} ppm and {} kHz linewidths.'.format(num_replicates,lw_shifts,lw_couplings),f)
            tau_list = []
            rho0_list = []
            ord_list = []
            replicate = 1
            while replicate <= num_replicates:
                df_err = df_fit.copy()
                df_err = randomizeShiftsCouplings(df_err,lw_shifts,lw_couplings)
                if procs > 1:
                    results = explore_comb_multiproc(df_err,period,flip,pas,beta,dc,matrixA,scalar,comb,procs=procs)
                else:
                    results = explore_comb(df_err,period,flip,pas,beta,dc,matrixA,scalar,comb)
                min_score = 100000
                for t,r,o,s in results:
                     if s < min_score:
                         min_score = s
                         min_tau = t
                         min_rho0 = r
                         min_ord = o
                log('# Replicate {0}: tau = {1:.2f}, rho0 = {2:.2f}, ord = {3:.2f}, score = {4:.3f}'.format(replicate,min_tau,min_rho0,min_ord,min_score),f)
                tau_list.append(min_tau)
                rho0_list.append(min_rho0)
                ord_list.append(min_ord)
                replicate += 1
            avg_tau = np.mean(tau_list)
            std_tau = np.std(tau_list)
            avg_rho0 = np.mean(rho0_list)
            std_rho0 = np.std(rho0_list)
            avg_ord = np.mean(ord_list)
            std_ord = np.std(ord_list)
            log('# Average tau: {0:.2f} +/- {1:.2f}'.format(avg_tau,std_tau),f)
            log('# Average rho0: {0:.2f} +/- {1:.2f}'.format(avg_rho0,std_rho0),f)
            log('# Average ord: {0:.2f} +/- {1:.2f}'.format(avg_ord,std_ord),f)
            tau = avg_tau * 0.017453
            rho0 = avg_rho0
            order = avg_ord

        # Or just fit directly to data
        else:
            if procs > 1:
                results = explore_comb_multiproc(df_fit,period,flip,pas,beta,dc,matrixA,scalar,comb,procs=procs)
            else:
                results = explore_comb(df_fit,period,flip,pas,beta,dc,matrixA,scalar,comb)
            log('# Writing score for each combination to: {}'.format(out_fit),f)
            f_fit = open(out_fit, 'w')
            f_fit.write('# tau rho ord score\n')
            min_score = 100000
            for t,r,o,s in results:
                f_fit.write('{0:.3f} {1:.3f} {2:.3f} {3}\n'.format(t,r,o,s))
                if s < min_score:
                    min_score = s
                    min_tau = t
                    min_rho0 = r
                    min_ord = o
            f_fit.close()
            log('# Minimum score: {0:.3f}'.format(min_score),f)
            log('# Fitted tau: {0:.3f}'.format(min_tau),f)
            log('# Fitted rho0: {0:.3f}'.format(min_rho0),f)
            log('# Fitted order: {0:.3f}'.format(min_ord),f)
            tau = min_tau * 0.017453
            rho0 = min_rho0
            order = min_ord

        time_stop = time.time()
        log('# Fitting completed in {0:.2f} minutes ({1:.2f} seconds).'.format((time_stop-time_start)/60,time_stop-time_start),f)

        
    # Explore combinations of tau and order for best fit to minimum and maximum bounds (--quickfit)
    if quickfit:
        
        # Read minimum and maximum points if specified. Convert to positive values.
        log('# Exploring combinations of tau and order for best fit to bounds specified.',f)
        assert (len(quickfit) == 4), 'The quickfit option should contain 4 variables!'
        if is_float(quickfit[0]):
            minShift = np.absolute(float(quickfit[0]))
            log('# Minimum chemical shift: {}.'.format(minShift),f)
        else:
            minShift = False
            log('# Minimum chemical shift not assigned.',f)
        if is_float(quickfit[1]):
            maxShift = np.absolute(float(quickfit[1]))
            log('# Maximum chemical shift: {}.'.format(maxShift),f)
            assert (minShift), 'A minimum chemical shift must also be specified!'
            assert (maxShift > minShift), 'The absolute maximum chemical shift must be greater than the absolute minimum shift.'
        else:
            maxShift = False
            log('# Maximum chemical shift not assigned.',f)
        if is_float(quickfit[2]):
            minCouple = np.absolute(float(quickfit[2]))
            log('# Minimum dipolar coupling: {}.'.format(minCouple),f)
        else:
            minCouple = False
            log('# Minimum  dipolar coupling not assigned.',f)
        if is_float(quickfit[3]):
            maxCouple = np.absolute(float(quickfit[3]))
            log('# Maximum dipolar coupling: {}.'.format(maxCouple),f)
            assert (minCouple), 'A minimum dipolar coupling must also be specified!'
            assert (maxCouple > minCouple), 'The absolute maximum chemical shift must be greater than the absolute minimum shift.'
        else:
            maxCouple = False
            log('# Maximum dipolar coupling not assigned.',f)
        
        taus = np.arange(fit_tau[0],fit_tau[1]+fit_tau[2],fit_tau[2])
        ords = np.arange(fit_order[0],fit_order[1]+fit_order[2],fit_order[2])
        comb = list(itertools.product(*[taus,ords]))
        ncomb = len(comb)
        log('# --fit_tau {}'.format(fit_tau),f)
        log('# --fit_order {}'.format(fit_order),f)
        log('# {} values of tau will be tested.'.format(len(taus)),f)
        log('# {} values of order parameters be tested.'.format(len(ords)),f)
        log('# A total of {} combinations will be tested.'.format(ncomb),f)
        time_start = time.time()
        if procs > 1:
            results = quickfit_comb_multiproc(minShift,maxShift,minCouple,maxCouple,rho0,period,flip,pas,beta,dc,matrixA,scalar,comb,procs=procs)
        else:
            results = quickfit_comb(minShift,maxShift,minCouple,maxCouple,rho0,period,flip,pas,beta,dc,matrixA,scalar,comb)
        time_stop = time.time()
        log('# Fitting completed in {0:.2f} minutes ({1:.2f} seconds).'.format((time_stop-time_start)/60,time_stop-time_start),f)
        log('# Writing score for each combination to: {}'.format(out_fit),f)
        f_fit = open(out_fit, 'w')
        f_fit.write('# tau rho ord score\n')
        min_score = 100000
        for t,o,s in results:
            f_fit.write('{0:.3f} {1:.3f} {2}\n'.format(t,o,s))
            if s < min_score:
                min_score = s
                min_tau = t
                min_ord = o
        f_fit.close()
        log('# Minimum score: {0:.3f}'.format(min_score),f)
        log('# Fitted tau: {0:.3f}'.format(min_tau),f)
        log('# Fitted order: {0:.3f}'.format(min_ord),f)
        tau = min_tau * 0.017453
        order = min_ord

    # Simulate wheel with final parameters
    log('# Simulating PISA wheel with the following parameters:',f)
    log('# --flip {}'.format(flip/0.017453),f)
    log('# --dc {}'.format(dc),f)
    log('# --pas {}'.format(pas),f)
    log('# --tau {0:.2f}'.format(tau/0.017453),f)
    log('# --rho0 {0:.2f}'.format(rho0),f)
    log('# --order {0:.2f}'.format(order),f)
    df = calcWheel(df,tau,rho0,period,order,flip,pas,beta,dc,matrixA)
    df, score = calcScores(df,scalar)
    df = calcRaw(df,order,flip,pas)

    # Specify of negative_dc is swithed on
    if negative_dc:
        log('# All dipolar couplings will be output as negative values.',f)
    
    # Write final results to log file
    if show_only:
        log('# --show_only specified. Outputting data only for residues {}'.format(show_only),f)
        df_out = filterData(df,show_only,rho_start,exclude=False,include=True)
    else:
        df_out = df
        
    log('# resname, index, rho, shift_calc, coupling_calc, shift_exp, coupling_exp, score, shift_raw, coupling_raw', f)
    for index in df_out.index.values:
        if negative_dc:
            coupling_calc = df_out.at[index,'coupling_calc'] * -1
            coupling_exp = df_out.at[index,'coupling_exp'] * -1
        else:
            coupling_calc = df_out.at[index,'coupling_calc']
            coupling_exp = df_out.at[index,'coupling_exp']
        log('{0} {1} {2:.2f} {3:.3f} {4:.3f} {5:.3f} {6:.3f} {7:.3f} {8:.3f} {9:.3f}'.format(
            df_out.at[index,'resname'],
            index,
            df_out.at[index,'rho'],
            df_out.at[index,'shift_calc'],
            coupling_calc,
            df_out.at[index,'shift_exp'],
            coupling_exp,
            df_out.at[index,'score'],
            df_out.at[index,'shift_raw'],
            df_out.at[index,'coupling_raw']
            ),f)
    #log('# Final fitting score: {}'.format(score),f)

    # Write wave file
    log('# Writing interpolated wave data to: {}'.format(out_wave),f)
    wf = open(out_wave, 'w')
    wave = calcWave(df_out,tau,rho0,period,order,flip,pas,beta,dc,rho_start,matrixA,margin=5.0,increment=0.1)
    for i in range(len(wave[0])):
        if negative_dc:
            coupling_wave = wave[3][i] * -1
        else:
            coupling_wave = wave[3][i]
        wf.write('{0:.1f} {1:.1f} {2:.3f} {3:.3f}\n'.format(wave[0][i],wave[1][i],wave[2][i],coupling_wave))
    wf.close()
        
    # Close files
    f.close()
    
if __name__ == '__main__':
    main()
