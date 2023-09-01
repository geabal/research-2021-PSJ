# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 16:12:40 2021

@author: sunjin
"""

import pandas as pd
import re
import numpy as np
import pyteomics.mzxml as mzxml

#txt 파일 파싱

#txt파일에서 C2N PSM을 파싱함. rank 1 C2N만 있는 경우
#specnum도 가져오도록 수정
def getC2N(filePath):
    f = open(filePath)
    
    specnums = []
    filenums = []
    rank = []
    score = []
    QVal = []
    sequence = []
    numlons = []
    protein = []
    c2nmass = []
    

    startidx = []
    idx = 0
    specnum = 0
    filelist = []
    fileidx = 0
    isc2n = 1
    filenum = 0
    prevms2num = 20000
    for i, line in enumerate(f):
        spl = line.split()
    
        if len(spl)==0:
            continue
        elif len(spl[0])>3:
            specnum = int(spl[2])
            mass = float(spl[-3])
            z = int(spl[-2])
            ms2num = int(spl[1])
            if prevms2num - ms2num > 100:
                filenum = filenum + 1
                fileName = spl[0].split('\\')[-1]
                filelist.append([fileName, fileidx])
            if isc2n ==1:
                startidx.append([specnum, idx, mass, z])
                fileidx = fileidx + 1
            isc2n = 0
            prevms2num = ms2num
            continue
    
        if  spl[0][0]=='A':
            isc2n = 1
            idx = idx+1
            specnums.append(specnum)
            filenums.append(filenum)
            c2nmass.append(float(spl[1]))
            rank.append(1)
            QVal.append(float(spl[4]))
            score.append(int(spl[3]))
            sequence.append(spl[5])
            protein.append(spl[6])
            numlons.append(int(spl[7]))
    
    if startidx[-1][1] == idx:
        del startidx[-1]
    
    f.close()
    dict1 = {"filenum":filenums,"SPECNUM":specnums,"RANK":rank,"MASS":c2nmass,"QScore":QVal,"SCORE":score, "SEQUENCE":sequence, "PROTEIN":protein, "NUMIONS":numlons }
    c2ndf = pd.DataFrame(dict1)
    
    return c2ndf, startidx, filelist

#A1, A2..., A10까지 rank1외의 C2N rank도 주어진 경우
def getC2N_rank(filePath):
    f = open(filePath)
    
    rank = []
    score = []
    QVal = []
    sequence = []
    numlons = []
    protein = []
    c2nmass = []

    startidx = []
    idx = 0
    specnum = 0
    filelist = []
    fileidx = 0
    isc2n = 1
    for i, line in enumerate(f):
        spl = line.split()
    
        if len(spl)==0:
            continue
        elif len(spl[0])>3:
            specnum = int(spl[2])
            mass = float(spl[-3])
            z = int(spl[-2])
            if spl[1] =='1':
                fileName = spl[0].split('\\')[-1]
                filelist.append([fileName, fileidx])
            if isc2n ==1:
                startidx.append([specnum, idx, mass, z])
                fileidx = fileidx + 1
            isc2n = 0
            
            continue
    
        if  spl[0][0]=='A':
            isc2n = 1
            idx = idx+1
            c2nmass.append(float(spl[1]))
            rank.append(int(spl[0][1:]))
            QVal.append(float(spl[4]))
            score.append(int(spl[3]))
            sequence.append(spl[5])
            protein.append(spl[6])
            numlons.append(int(spl[7]))
    
    if startidx[-1][1] == idx:
        del startidx[-1]
    
    f.close()
    dict1 = {"MASS":c2nmass,"RANK":rank,"QScore":QVal,"SCORE":score, "SEQUENCE":sequence, "PROTEIN":protein, "NUMIONS":numlons }
    c2ndf = pd.DataFrame(dict1)
    
    return c2ndf, startidx, filelist

#txt파일에서 NOR PSM을 파싱함
def getNOR(filePath):
    f = open(filePath)
    
    specnums = []
    filenums = []    
    rank = []
    score = []
    QVal = []
    sequence = []
    numlons = []
    protein = []
    normass = []

    startidx = []
    idx = 0
    specnum = 0
    filelist = []
    fileidx = 0
    filenum = 0
    isnor = 1
    prevms2num = 20000
    for i, line in enumerate(f):
        spl = line.split()
    
        if len(spl)==0:
            continue
        elif len(spl[0])>3:
            specnum = int(spl[2])
            mass = float(spl[-3])
            z = int(spl[-2])
            ms2num = int(spl[1])
            if prevms2num - ms2num > 100:
                filenum = filenum + 1
                fileName = spl[0].split('\\')[-1]
                filelist.append([fileName, fileidx])
            if isnor ==1:
                startidx.append([specnum, idx, mass, z])
                fileidx = fileidx + 1
            isnor = 0
            prevms2num = ms2num

            continue
    
        if  spl[0][0]!='A' and spl[0]!='T' and spl[0]!='F':
            isnor = 1
            idx = idx+1
            specnums.append(specnum)
            filenums.append(filenum)            
            rank.append(int(spl[0]))
            normass.append(float(spl[1]))
            QVal.append(float(spl[4]))
            score.append(int(spl[3]))
            sequence.append(spl[5])
            protein.append(spl[6])
            numlons.append(int(spl[7]))
            
    if startidx[-1][1] == idx:
        del startidx[-1]
    
    f.close()
    dict1 = {"filenum":filenums,"SPECNUM":specnums,"RANK":rank,"MASS":normass, "QScore":QVal,"SCORE":score, "SEQUENCE":sequence, "PROTEIN":protein, "NUMIONS":numlons }
    nordf = pd.DataFrame(dict1)
    
    return nordf, startidx, filelist

#전체 PSM 파싱
def getAllPSM(filePath):
    f = open(filePath)
    specnums = []
    filenums = []    
    rank = []
    score = []
    QVal = []
    sequence = []
    numlons = []
    protein = []
    isnor = []

    startidx = []
    idx = 0
    specnum = 0
    filelist = []
    fileidx = 0
    filenum = 0
    prevms2num = 20000
    for i, line in enumerate(f):
        spl = line.split()
    
        if len(spl)==0:
            continue
        elif len(spl[0])>3:
            specnum = int(spl[2])
            mass = float(spl[-3])
            z = int(spl[-2])
            startidx.append([specnum, idx, mass, z])
            fileidx = fileidx + 1
            ms2num = int(spl[1])
            if prevms2num-ms2num > 100:
                filenum = filenum + 1
                fileName = spl[0].split('\\')[-1]
                filelist.append([fileName, fileidx])
            prevms2num = ms2num
    
        if len(spl[0])<3:
            idx = idx+1
            if  spl[0][0]!='A' and spl[0]!='T' and spl[0]!='F':
                rank.append(int(spl[0]))
                isnor.append(1)
            else:
                rank.append(spl[0])
                isnor.append(0)
            specnums.append(specnum)
            filenums.append(filenum)
            QVal.append(float(spl[4]))
            score.append(int(spl[3]))
            sequence.append(spl[5])
            protein.append(spl[6])
            numlons.append(int(spl[7]))
            
    if startidx[-1][1] == idx:
        del startidx[-1]
    
    f.close()
    dict1 = {"filenum":filenums,"SPECNUM":specnums,"isNOR":isnor,"RANK":rank,"QScore":QVal,"SCORE":score, "SEQUENCE":sequence, "PROTEIN":protein, "NUMIONS":numlons }
    df = pd.DataFrame(dict1)
    
    return df, startidx, filelist

#percolator 관련 함수
#.pin 파일 읽기
def readPin(filePath):
    f = open(filePath)
    cols = f.readline().split()
    print(cols)
    
    rows=[]

    for i, line in enumerate(f):
        s = line.split()
        row = s[:len(cols)-1]
        row.append(s[len(cols)-1:])
        rows.append(row)

    f.close()
    df = pd.DataFrame(rows, columns=cols)
    
    return df
#.pin 파일 쓰
def toPin(df, filePath):
    f = open(filePath, 'w', encoding='utf-8')
    
    vals = df.values.tolist()
    line = '\t'.join(df.columns.tolist())+'\n'
    f.write(line)
    
    for i in range(len(df)):
        if i == len(df)-1:
            line = '\t'.join(map(str,vals[i][:-1])) +'\t' + '\t'.join(vals[i][-1])
            f.write(line)
        else:
            line = '\t'.join(map(str,vals[i][:-1])) +'\t' + '\t'.join(vals[i][-1])+'\n'
            f.write(line)

    # 파일 객체 닫기
    f.close()
    
    return

def makeProteinlist(df):
    lst = []
    proteins = df['Proteins']
    cols = df.columns.tolist()
    vals = df.values.tolist()
    didx = cols.index('Label')
    
    for j in range(len(proteins)):
        protein = proteins[j]
        plst = protein.split(';')
        nplst = []
        for i in range(len(plst)):
            #target에는 decoy protein이 없고, decoy PSM은 decoy protein만 갖도록 함
            p=plst[i]
            if vals[j][didx]==1:
                if p[:3]=="XXX":
                    continue            
            else: 
                if p[:3]!="XXX":
                    continue
            p = p.split('[')[0]
            nplst.append(p)
        
        lst.append(nplst)
    
    s = pd.Series(lst)     
    return s

def makeFeatures(d,df):
    df['SCORE'] = d['SCORE']
    df['NUMIONS'] = d['NUMIONS']
    df['NTT'] = d['NTT']
    df['pepLen'] = d['pepLen']
    
    
    return df


def makePinput(d):
    df = pd.DataFrame(index=range(0,0), columns = ['SpecId', 'Label', 'ScanNr'])
    df['ScanNr']=d['SPECNUM']
    df['SpecId']=d['SpecId']
    df['Label']=d['isDecoy']
    df['Label'] = df['Label'].replace(1,-1)
    df['Label'] = df['Label'].replace(0,1)
    
    
    df = makeFeatures(d,df)
    
    df['Peptide'] = d['SEQUENCE']
    
    df['Proteins'] = d['PROTEIN']
    df.reset_index(drop=True,inplace=True)
    df['Proteins'] =makeProteinlist(df)
    
    return df

#기타 전처리

def getNTT(df):
    NTTdf = df.copy()
    NTTdf['NTT'] = 0
    cols = NTTdf.columns.tolist()
    vals = NTTdf.values.tolist()
    pcol = cols.index("PROTEIN")
    ncol = cols.index('NTT')
    
    for i in range(len(vals)):
        #vals[i][3]은 PROTEIN column
        spl = vals[i][pcol].split(';')
        ntt = 0
        for j in range(len(spl)):
            tmp = int(spl[j][-3])
            if tmp > ntt:
                ntt = tmp
                if ntt ==2:
                    break
        vals[i][ncol] = ntt
    
    NTTdf = pd.DataFrame(vals, columns=cols)

    return NTTdf

def getIsDecoy(df):
    result = df.copy()
    result['isDecoy'] = 0
    cols = result.columns.tolist()
    vals = result.values.tolist()
    didx = cols.index('isDecoy')
    pidx = cols.index('PROTEIN')
    
    for i in range(len(vals)):
        spls = vals[i][pidx].split(';')
        idflag = 0
        tpeps = []
        dpeps = []
        for j in range(len(spls)):
            if spls[j][0:3]=='XXX':
                idflag = idflag +1
                dpeps.append(j)
            else:
                tpeps.append(j)
        
        if idflag == len(spls):
            vals[i][didx] = 1
        elif idflag > 0:
            tntt = int(spls[tpeps[0]].split('(')[1][0])
            dntt = int(spls[dpeps[0]].split('(')[1][0])
            if tntt < dntt:
                vals[i][didx]=1
            
    result = pd.DataFrame(vals, columns=cols)
    return result

#PTM 수를 세어 column으로 추가
def getNPTM(df):
    tmpdf = df.copy()
    tmpdf['NPTM'] = 0
    
    cols = tmpdf.columns.tolist()
    vals = tmpdf.values.tolist()    
    colnp = cols.index('NPTM')
    colp = cols.index('PTM')
    
    if type(vals[0][colp]) == list:
        for i in range(len(vals)):
            vals[i][colnp] = len(vals[i][colp])
    else:
        for i in range(len(vals)):
            vals[i][colnp] = vals[i][colp].count('(')
            
    result = pd.DataFrame(vals, columns=cols)
    return result

def getBaseForm(df):
    tmpdf = df.copy()
    tmpdf['BaseForm'] = df['SEQUENCE']
    
    cols = tmpdf.columns.tolist()
    vals = tmpdf.values.tolist()
    colbf = cols.index('BaseForm')
    
    for i in range(len(vals)):
        vals[i][colbf] = getstr(vals[i][colbf])
    
    result = pd.DataFrame(vals, columns=cols)
    return result

# unmodified form PSM이 짝으로 존재하는 modified form PSM과 그 쌍인 unmodified PSM을 1로 표기
def getisPair(df):
    tmpdf = df.copy()
    if 'NPTM' not in tmpdf.columns:
        tmpdf = getNPTM(tmpdf)
    if 'BaseForm' not in tmpdf.columns:
        tmpdf = getBaseForm(tmpdf)
        
    tmpdf['isPair'] = 0
    
    cols = tmpdf.columns.tolist()
    vals = tmpdf.values.tolist()
    colbf = cols.index('BaseForm')
    colnp = cols.index('NPTM')
    colip = cols.index('isPair')
 
    vals.sort(key = lambda x:(x[colbf], x[colnp]))
    
    pairFlag = 0
    for i in range(len(vals)-1):
        if vals[i][colnp]==0:
            pairFlag = 1
            
        if vals[i][colbf] == vals[i+1][colbf] and pairFlag == 1:
            vals[i][colip] = 1
            vals[i+1][colip] = 1
        else:
            pairFlag = 0
    
    result = pd.DataFrame(vals, columns=cols)
    return result
#unmodified form PSM이 존재하지 않고, modified form PSM끼리 짝인 것을 1이라 표시
def getisModPair(df):
    tmpdf = df.copy()
    tmpdf['isModPair'] = 0
    
    cols = tmpdf.columns.tolist()
    vals = tmpdf.values.tolist()
    colbf = cols.index('BaseForm')
    colnp = cols.index('NPTM')
    colip = cols.index('isModPair')
 
    vals.sort(key = lambda x:(x[colbf], x[colnp]))
    
    pairFlag = 1
    for i in range(len(vals)-1):
        if vals[i][colnp]==0:
            pairFlag = 0
        elif (i!=0) and (vals[i][colbf] != vals[i-1][colbf]):
            pairFlag = 1
            
        if vals[i][colbf] == vals[i+1][colbf] and pairFlag == 1:
            vals[i][colip] = 1
            vals[i+1][colip] = 1
        
    result = pd.DataFrame(vals, columns=cols)
    return result

#peptide length 컬럼 추가
def getpeplen(df):
    pldf = df.copy()
    pldf['pepLen'] = 0
    cols = pldf.columns.tolist()
    vals = pldf.values.tolist()
    scol = cols.index("SEQUENCE")
    pcol = cols.index('pepLen')
    
    for i in range(len(vals)):
        seq = vals[i][scol]
        #-2:sequence에 N, C-term도 있으므로 그 길이만큼 
        vals[i][pcol] = len(getstr(seq))-2
    
    pldf = pd.DataFrame(vals, columns=cols)

    return pldf

def getdif(seq1, seq2):
    diff = 0
    islendif = 0
    for i in range(len(seq1)):
            if seq1[i]!=seq2[i]:
                if seq1[i] == '_' or seq2[i] =='_':
                    islendif = 1
                else:
                    diff = diff+1
    if islendif ==1:
        diff = diff+1
    return diff

#sequence distance를 구하는 함수 구현
def getseqdist(seq1, seq2):
    dist = 0
    if len(seq1) == len(seq2):
        dist = getdif(seq1, seq2)
    elif len(seq1) > len(seq2):
        lendif = len(seq1)-len(seq2)
        tmp1 = seq2 + "_"*lendif
        tmp2 = "_"*lendif + seq2 
        dist = min(getdif(seq1, tmp1), getdif(seq1, tmp2))
    else:
        lendif = len(seq2)-len(seq1)
        tmp1 = seq1 + "_"*lendif
        tmp2 = "_"*lendif + seq1 
        dist = min(getdif(seq2, tmp1), getdif(seq2, tmp2))
    
    return dist

def getstr(seq):
    result = re.sub(r'[^A-Z]','', seq)
    return result

def getRank1(df,stList):
    result = df.copy()
    cols = result.columns.tolist()
    vals = result.values.tolist()
    newVals = []
    
    for st in stList:
        newVals.append(vals[st[1]])
    
    result = pd.DataFrame(newVals, columns=cols)
    
    return result

#delta score를 계산해 반환함.
#NOR/C2N PSM이 따로 있는 경우 사용. rank 1만 반환
#penalty: seqdiff 기준을 만족하지 못하는 PSM의 delta score에 적용하는 페널티값.
def getDeltaScore_rank1(df, startidx, seqdiff=2, threshold=3, penalty=0):
    df['DeltaScore'] = df['SCORE']
    cols = df.columns.tolist()
    vals = df.values.tolist()
    rank1s = []
    
    colScore = cols.index('SCORE')
    colSeq = cols.index('SEQUENCE')
    
    for i in range(len(startidx)):
        start = startidx[i][1]
        
        if i == len(startidx)-1:
            end = len(df)
        else:
            end = startidx[i+1][1]
        
        tmpdf = vals[start:end]
        
        seq1 = getstr(tmpdf[0][colSeq])
        diff = 0
        #rank1의 score가 음수면 버림
        if tmpdf[0][colScore]<0:
            continue
        for j in range(1, len(tmpdf)):
            seq2 = getstr(tmpdf[j][colSeq])
            diff = getseqdist(seq1, seq2)
            if diff >seqdiff:
                delta = tmpdf[0][colScore] - tmpdf[j][colScore]
                tmpdf[0][-1] = delta
                break
        
        if diff < seqdiff:
            tmpdf[0][-1] = tmpdf[0][-1]-penalty
                
        if tmpdf[0][-1]<threshold:
            continue
        else:
            rank1s.append(tmpdf[0])
    rank1df = pd.DataFrame(rank1s, columns=cols)
    return rank1df

def calfdr(t, d):
    return (d+1)/t

#fdr 계산 기준: calfdr
#col: 정렬 기준이 될 column. ex) score, delta score 등. string으로 주어짐.
#TD:결과로 반환하는 fdr이 target만 반환하는 경우 1, decoy만 반환하는 경우 2, 둘 다 반환하는 경우 0
def FDR(df, col, threshold=0.01, TD=0):
    fdrdf = df.copy()
    fdrdf = fdrdf.sort_values(by = col, ascending = False)
    
    cols= fdrdf.columns.tolist()
    vals = fdrdf.values.tolist()
    didx = cols.index('isDecoy')
    
    ntarget = len(fdrdf[fdrdf['isDecoy']==0])
    ndecoy = len(fdrdf[fdrdf['isDecoy']==1])
    print('T: ', ntarget, 'D: ', ndecoy)
    lastPSM = 0
    
    for i in reversed(range(len(vals))):
        if vals[i][didx] ==0:
            ntarget = ntarget -1
        else:
            ndecoy = ndecoy -1
        fdr = calfdr(ntarget, ndecoy)
        if fdr <= threshold:
            lastPSM = i
            break
    vals = vals[:lastPSM]
    fdrdf = pd.DataFrame(vals, columns=cols)
    if TD ==1:
        fdrdf = fdrdf[fdrdf['isDecoy']==0]
    elif TD==2:
        fdrdf = fdrdf[fdrdf['isDecoy']==1]
    
    return fdrdf

# MS1 관련 코드



#specs: mzxml 객체
#specnum: ms2의 specnum. 이 ms2의 ms1 해당하는 specnum 반환
def findms1idx(specs, specnum):
    tmp = int(specnum)
    msL = specs.get_by_id(specnum)['msLevel']
    while msL==2:
        msL = specs.get_by_id(str(tmp))['msLevel']
        tmp = tmp-1
    return str(tmp)

def findSpec(mzarray,premz, et=10):
    mzs = mzarray.tolist()
    ET = premz*(et/(10**6))
    
    #이진탐색 알고리즘 수랭
    left = 0
    right = len(mzs)-1
    cur = (left+right)//2
    upflag = 0
    lowflag = 0
    while (left!=cur) and (right!=cur):
        upflag = 0
        lowflag = 0
        if mzs[cur]<premz+ET:
            upflag = 1
        if mzs[cur]>premz-ET:
            lowflag = 1
            
        if (lowflag==1) and (upflag==1):
            return cur
        elif (lowflag==1) and (upflag==0):
            right = cur
        elif (lowflag==0) and (upflag==1):
            left = cur
        cur = (left+right)//2
        
    
    if (lowflag==0) and (upflag==1):
        cur = right
    
    if (mzs[cur]>premz-ET) and (mzs[cur]<premz+ET):
        return cur
    else:
        return -1
    

def front5MS1Match(df, st, files):
    result = df.copy()
    result['MS1match_front5'] = 0
    cols = result.columns.tolist()
    vals = result.values.tolist()
    colmass = cols.index('MASS')
    midx = cols.index('MS1match_front5')
    
    mzxmlPaths = ['mzxml/b1906_293T_proteinID_01A_QE3_122212.mzXML','mzxml/b1922_293T_proteinID_02A_QE3_122212.mzXML'
                ,'mzxml/b1923_293T_proteinID_03A_QE3_122212.mzXML','mzxml/b1924_293T_proteinID_04A_QE3_122212.mzXML'
                ,'mzxml/b1925_293T_proteinID_05A_QE3_122212.mzXML','mzxml/b1926_293T_proteinID_06A_QE3_122212.mzXML'
                ,'mzxml/b1927_293T_proteinID_07A_QE3_122212.mzXML','mzxml/b1928_293T_proteinID_08A_QE3_122212.mzXML'
                ,'mzxml/b1929_293T_proteinID_09A_QE3_122212.mzXML','mzxml/b1930_293T_proteinID_10A_QE3_122212.mzXML'
                ,'mzxml/b1931_293T_proteinID_11A_QE3_122212.mzXML','mzxml/b1932_293T_proteinID_12A_QE3_122212.mzXML'
                ,'mzxml/b1937_293T_proteinID_01B_QE3_122212.mzXML','mzxml/b1938_293T_proteinID_02B_QE3_122212.mzXML'
                ,'mzxml/b1939_293T_proteinID_03B_QE3_122212.mzXML','mzxml/b1940_293T_proteinID_04B_QE3_122212.mzXML'
                ,'mzxml/b1941_293T_proteinID_05B_QE3_122212.mzXML','mzxml/b1942_293T_proteinID_06B_QE3_122212.mzXML'
                ,'mzxml/b1943_293T_proteinID_07B_QE3_122212.mzXML','mzxml/b1944_293T_proteinID_08B_QE3_122212.mzXML'
                ,'mzxml/b1945_293T_proteinID_09B_QE3_122212.mzXML','mzxml/b1946_293T_proteinID_10B_QE3_122212.mzXML'
                ,'mzxml/b1947_293T_proteinID_11B_QE3_122212.mzXML','mzxml/b1948_293T_proteinID_12B_QE3_122212.mzXML']
                
    ms1hit = 0
    for i in range(len(files)):
        mzxmlPath = mzxmlPaths[i]
        specs = mzxml.read(mzxmlPath)
        print("MS1 Match Start: ",mzxmlPath)
        
        firstIdx = files[i][1]
        
        j = firstIdx
        ms1s = []
        for spec in specs:
            msL = spec['msLevel']
        
            if msL==1:
                ms1s.append(spec)
                if len(ms1s)==6:
                    del ms1s[0]
            
            if j==len(st):
                break
            specnum = str(st[j][0])
            #MS1과 MS2 매치 시작
            if specnum==spec['num']:
                for ms1 in ms1s:
                    matched = 0
                    for isotope in range(2):
                        exmass = vals[j][colmass]
                        z = st[j][3]
                        #1.008664: neutron mass
                        premz = (exmass + z*1.007276 +isotope*1.008664)/z
                        ms1specidx = findSpec(ms1['m/z array'],premz, 10)
                        if ms1specidx !=-1:
                            matched = matched +1
                    if matched == 2:
                        vals[j][midx] = 1
                        ms1hit = ms1hit+1
                        break
                j= j+1
    
        print("MS1 Match End, hits: ", ms1hit)
        
    result = pd.DataFrame(vals, columns=cols)    
    print(ms1hit)
    return result
    
def rear5MS1Match(df, st, files):
    result = df.copy()
    result['MS1match_rear5'] = 0
    cols = result.columns.tolist()
    vals = result.values.tolist()
    colmass = cols.index('MASS')
    midx = cols.index('MS1match_rear5')
    
    mzxmlPaths = ['mzxml/b1906_293T_proteinID_01A_QE3_122212.mzXML','mzxml/b1922_293T_proteinID_02A_QE3_122212.mzXML'
                ,'mzxml/b1923_293T_proteinID_03A_QE3_122212.mzXML','mzxml/b1924_293T_proteinID_04A_QE3_122212.mzXML'
                ,'mzxml/b1925_293T_proteinID_05A_QE3_122212.mzXML','mzxml/b1926_293T_proteinID_06A_QE3_122212.mzXML'
                ,'mzxml/b1927_293T_proteinID_07A_QE3_122212.mzXML','mzxml/b1928_293T_proteinID_08A_QE3_122212.mzXML'
                ,'mzxml/b1929_293T_proteinID_09A_QE3_122212.mzXML','mzxml/b1930_293T_proteinID_10A_QE3_122212.mzXML'
                ,'mzxml/b1931_293T_proteinID_11A_QE3_122212.mzXML','mzxml/b1932_293T_proteinID_12A_QE3_122212.mzXML'
                ,'mzxml/b1937_293T_proteinID_01B_QE3_122212.mzXML','mzxml/b1938_293T_proteinID_02B_QE3_122212.mzXML'
                ,'mzxml/b1939_293T_proteinID_03B_QE3_122212.mzXML','mzxml/b1940_293T_proteinID_04B_QE3_122212.mzXML'
                ,'mzxml/b1941_293T_proteinID_05B_QE3_122212.mzXML','mzxml/b1942_293T_proteinID_06B_QE3_122212.mzXML'
                ,'mzxml/b1943_293T_proteinID_07B_QE3_122212.mzXML','mzxml/b1944_293T_proteinID_08B_QE3_122212.mzXML'
                ,'mzxml/b1945_293T_proteinID_09B_QE3_122212.mzXML','mzxml/b1946_293T_proteinID_10B_QE3_122212.mzXML'
                ,'mzxml/b1947_293T_proteinID_11B_QE3_122212.mzXML','mzxml/b1948_293T_proteinID_12B_QE3_122212.mzXML']
                
    ms1hit = 0
    for i in range(len(files)):
        mzxmlPath = mzxmlPaths[i]
        specs = mzxml.read(mzxmlPath)
        print("MS1 Match Start: ",mzxmlPath)
        
        firstIdx = files[i][1]
        
        j = firstIdx
        #ms1s: ms1 list
        #matchInfo : [idx, ms2Mass, ms2Charge, [ms1s]]의 list
        ms1s = []
        matchInfo = []
        MIqueue = []
        for spec in specs:
            msL = spec['msLevel']
            
            if (msL == 1) and (len(MIqueue)!=0):
                d = 0
                for t in range(len(MIqueue)) :
                    MIqueue[t][3].append(spec)
                    if len(MIqueue[t][3])==5:
                        matchInfo.append(MIqueue[t])
                        d = d+1
                del MIqueue[:d]
            
            if j==len(st):
                break
            specnum = str(st[j][0])
            #+5MS1 저장
            if specnum==spec['num']:
                exmass = vals[j][colmass]
                z = st[j][3]
                tmp = [j, exmass, z, []]
                MIqueue.append(tmp)
                j= j+1
                
        #MS1, MS2 match 시작
        for mi in matchInfo:
            ms1s = mi[3]
            for ms1 in ms1s:
                matched = 0
                for isotope in range(2):
                    exmass = mi[1]
                    z = mi[2]
                    #1.008664: neutron mass
                    premz = (exmass + z*1.007276 +isotope*1.008664)/z
                    ms1specidx = findSpec(ms1['m/z array'],premz, 10)
                    if ms1specidx !=-1:
                            matched = matched +1
                if matched == 2:
                    vals[mi[0]][midx] = 1
                    ms1hit = ms1hit+1
                    break
    
        print("MS1 Match End, hits: ", ms1hit)
        
    result = pd.DataFrame(vals, columns=cols)    
    print(ms1hit)
    return result

def MS1Match(df, st, files):
    result = df.copy()
    result['MS1match'] = 0
    cols = result.columns.tolist()
    vals = result.values.tolist()
    colmass = cols.index('MASS')
    midx = cols.index('MS1match')
    
    mzxmlPaths = ['mzxml/b1906_293T_proteinID_01A_QE3_122212.mzXML','mzxml/b1922_293T_proteinID_02A_QE3_122212.mzXML'
                ,'mzxml/b1923_293T_proteinID_03A_QE3_122212.mzXML','mzxml/b1924_293T_proteinID_04A_QE3_122212.mzXML'
                ,'mzxml/b1925_293T_proteinID_05A_QE3_122212.mzXML','mzxml/b1926_293T_proteinID_06A_QE3_122212.mzXML'
                ,'mzxml/b1927_293T_proteinID_07A_QE3_122212.mzXML','mzxml/b1928_293T_proteinID_08A_QE3_122212.mzXML'
                ,'mzxml/b1929_293T_proteinID_09A_QE3_122212.mzXML','mzxml/b1930_293T_proteinID_10A_QE3_122212.mzXML'
                ,'mzxml/b1931_293T_proteinID_11A_QE3_122212.mzXML','mzxml/b1932_293T_proteinID_12A_QE3_122212.mzXML'
                ,'mzxml/b1937_293T_proteinID_01B_QE3_122212.mzXML','mzxml/b1938_293T_proteinID_02B_QE3_122212.mzXML'
                ,'mzxml/b1939_293T_proteinID_03B_QE3_122212.mzXML','mzxml/b1940_293T_proteinID_04B_QE3_122212.mzXML'
                ,'mzxml/b1941_293T_proteinID_05B_QE3_122212.mzXML','mzxml/b1942_293T_proteinID_06B_QE3_122212.mzXML'
                ,'mzxml/b1943_293T_proteinID_07B_QE3_122212.mzXML','mzxml/b1944_293T_proteinID_08B_QE3_122212.mzXML'
                ,'mzxml/b1945_293T_proteinID_09B_QE3_122212.mzXML','mzxml/b1946_293T_proteinID_10B_QE3_122212.mzXML'
                ,'mzxml/b1948_293T_proteinID_12B_QE3_122212.mzXML','mzxml/b1948_293T_proteinID_12B_QE3_122212.mzXML']
    ms1hit = 0
    for i in range(len(files)):
        mzxmlPath = mzxmlPaths[i]
        specs = mzxml.read(mzxmlPath)
        print("MS1 Match Start: ",mzxmlPath)
        
        firstIdx = files[i][1]
        
        j = firstIdx
        ms1 = None
        for spec in specs:
            msL = spec['msLevel']
        
            if msL==1:
                ms1 =  spec
            
            if j==len(st):
                break
            specnum = str(st[j][0])
            matched = 0
            #MS1과 MS2 매치 시작
            if specnum==spec['num']:
                for isotope in range(2):
                    exmass = vals[j][colmass]
                    z = st[j][3]
                    #1.008664: neutron mass
                    premz = (exmass + z*1.007276 +isotope*1.008664)/z
                    ms1specidx = findSpec(ms1['m/z array'],premz, 10)
                    if ms1specidx !=-1:
                        matched = matched +1
                if matched == 2:
                    vals[j][midx] = 1
                    ms1hit = ms1hit+1
                j= j+1
    
        print("MS1 Match End, hits: ", ms1hit)
        
    result = pd.DataFrame(vals, columns=cols)    
    print(ms1hit)
    return result

#학습 관련 함수들
#df: fdr을 통과한 PSMs
#mdf: merged df.
#첫 입력으로 들어가는 training data를 만드는 함수.
def makeTraindf_initial(df, mdf):
    df1 = df.copy()
    df1 = df1.sort_values(by = ['filenum','SPECNUM','isNOR'], ascending =True)
    mdf1 = mdf.sort_values(by = ['filenum','SPECNUM','isNOR'], ascending =True)
    
    
    cols= df1.columns.tolist()
    vals = df1.values.tolist()
    mvals = mdf1.values.tolist()
    
    nvals = []
    ncols = ['filenum','SPECNUM','DeltaScore_C2N','DeltaScore_NOR',
             'MS1match_rear5_C2N','MS1match_front5_C2N','MS1match_C2N','MS1match_rear5_NOR','MS1match_front5_NOR','MS1match_NOR',
             'seqDist','pepLen_C2N','pepLen_NOR','isCoelute']
    
    nidx = cols.index('isNOR')
    fidx = cols.index('filenum')
    snidx = cols.index('SPECNUM')
    
    dsidx = cols.index('DeltaScore')
    sdidx = cols.index('seqDist')
    plidx = cols.index('pepLen')
    mridx = cols.index('MS1match_rear5')
    mfidx = cols.index('MS1match_front5')
    midx = cols.index('MS1match')
    icidx = cols.index('isCoelute')
    
    j = 0
    for i in range(len(mvals)):
        nrow=[]
        if (mvals[i][snidx]==vals[j][snidx]):
            if mvals[i][nidx] == 1:
                nrow = [mvals[i][fidx],mvals[i][snidx],mvals[i-1][dsidx],mvals[i][dsidx],
                       mvals[i-1][mridx],mvals[i-1][mfidx],mvals[i-1][midx],mvals[i][mridx],mvals[i][mfidx],mvals[i][midx],
                       mvals[i-1][sdidx],mvals[i-1][plidx],mvals[i][plidx],vals[j][icidx]]
            else:
                nrow = [mvals[i][fidx],mvals[i][snidx],mvals[i][dsidx],mvals[i+1][dsidx],
                       mvals[i][mridx],mvals[i][mfidx],mvals[i][midx],mvals[i+1][mridx],mvals[i+1][mfidx],mvals[i+1][midx],
                       mvals[i][sdidx],mvals[i][plidx],mvals[i+1][plidx],vals[j][icidx]]
            j = j + 1
            if j == len(vals):
                break
        else:
            continue
                
        if len(nvals)!=0:
            if nvals[-1][1] == nrow[1]:
                continue
        nvals.append(nrow)
    
    result = pd.DataFrame(nvals, columns=ncols)
            
    return result

#selectPSM_1과 마찬가지로 기존 기준대로 PSM 선택. st(start list, 각 spectnum 별 시작 인덱스 정보가 있음)가 없을 때 사용
def selectPSM_01(modpdf):
    result = modpdf.copy()
    cols = result.columns.tolist()
    vals = result.values.tolist()
    newVals = []
    
    colScore = cols.index('SCORE')
    colni = cols.index('NUMIONS')
    colisnor = cols.index('isNOR')
    colsn = cols.index('SPECNUM')
    
    presn = -1
    for i in range(len(vals)):
        
        if i+1 != len(vals):
            if (vals[i][colisnor]==1) and (vals[i+1][colsn] == vals[i][colsn]):
                presn = vals[i][colsn]
                continue
        if vals[i][colisnor]==1:
            newVals.append(vals[i])
            presn = vals[i][colsn]
            continue
        
        if (vals[i][colisnor]==0) and (presn == vals[i][colsn]):
            if (vals[i][colScore] >= vals[i-1][colScore]+6) and (vals[i][colni] >= vals[i-1][colni]+1):
                newVals.append(vals[i])
            else:
                newVals.append(vals[i-1])
            presn = vals[i][colsn]
        elif (vals[i][colisnor]==0) and (presn != vals[i][colsn]):
            newVals.append(vals[i])
            presn = vals[i][colsn]
            
    
    result = pd.DataFrame(newVals, columns=cols)
    
    return result

def selectPSM_1(modpdf, modpst):
    result = modpdf.copy()
    cols = result.columns.tolist()
    vals = result.values.tolist()
    newVals = []
    
    colScore = cols.index('SCORE')
    colni = cols.index('NUMIONS')
    colrank = cols.index('RANK')
    
    
    for i in range(len(modpst)):
        stidx = modpst[i][1]
        if i == len(modpst)-1:
            endidx = len(vals)
        else:
            endidx = modpst[i+1][1]
            
        nor = []
        c2n = []
        for j in range(stidx,endidx):
            if (vals[j][colrank] == 1) and (nor==[]):
                nor = vals[j]
            if (vals[j][colrank] == 'A') and (c2n==[]):
                c2n = vals[j]
        #스펙트럼에 어떤 매칭도 없는 경우
        if (nor == []) and (c2n == []):
            continue
        #스펙트럼에 NOR 매치가 없는 경우
        if nor == []:
            newVals.append(c2n)
        #C2N 매치가 없는 경우
        elif c2n == []:
            newVals.append(nor)
        else:
            #기존 기준대로 최종 PSM 선택
            if (c2n[colScore] >= nor[colScore]+6) and (c2n[colni] >= nor[colni]+1):
                newVals.append(c2n)
            else:
                newVals.append(nor)
    
    result = pd.DataFrame(newVals, columns=cols)
    
    return result

#coeluting peptide를 고려해서 PSM 선택.
def selectPSM_2(modpdf):
    result = modpdf.copy()
    cols = result.columns.tolist()
    vals = result.values.tolist()
    newVals = []
    
    colScore = cols.index('SCORE')
    colni = cols.index('NUMIONS')
    colrank = cols.index('RANK')
    colisnor = cols.index('isNOR')
    colsn = cols.index('SPECNUM')
    colsd = cols.index('seqDist')
    colmm = cols.index('MS1match')
    
    nnum = 0
    cnum = 0
    #nor, c2n 둘 다 있을 때 둘 중 하나 선택. select nor or select c2n
    snnum = 0
    scnum = 0
    coelutePep = 0
    presn = -1
    for i in range(len(vals)):
        
        if i+1 != len(vals):
            if (vals[i][colisnor]==1) and (vals[i+1][colsn] == vals[i][colsn]):
                presn = vals[i][colsn]
                continue
        if vals[i][colisnor]==1:
            newVals.append(vals[i])
            presn = vals[i][colsn]
            nnum = nnum +1
            continue
        
        #nor, c2n 매치가 둘 다 있는 경우
        if (vals[i][colisnor]==0) and (presn == vals[i][colsn]):
            #coeluting peptide
            if (vals[i][colmm]==1) and (vals[i][colsd]>2):
                newVals.append(vals[i])
                newVals.append(vals[i-1])
                coelutePep = coelutePep +1
                presn = vals[i][colsn]
                continue
            #둘 중 하나 선택
            if (vals[i][colScore] >= vals[i-1][colScore]+6) and (vals[i][colni] >= vals[i-1][colni]+1):
                newVals.append(vals[i])
                scnum = scnum +1
            else:
                newVals.append(vals[i-1])
                snnum = snnum +1
            presn = vals[i][colsn]
        elif (vals[i][colisnor]==0) and (presn != vals[i][colsn]):
            newVals.append(vals[i])
            presn = vals[i][colsn]
            cnum = cnum +1
            
    print("only NOR: ", nnum,"\nonly C2N: ", cnum,"\nseleicting NOR or C2N: ",snnum, scnum, "\nCoeluting peptides: ", coelutePep)
    
    result = pd.DataFrame(newVals, columns=cols)
    
    return result