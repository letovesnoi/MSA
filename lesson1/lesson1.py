import matplotlib
matplotlib.use('Agg')

__author__ = 'lenk'

from pylab import *

import math


def readDNA(fin):
    DNA = ''
    while True:
        temp = fin.readline()
        if '//' in temp:
            break
        tempDNA = temp[10:-1]
        for fragment in tempDNA.split(' '):
            DNA += fragment
    return DNA

def reverseComplement(s):
    n = len(s)
    i = n - 1
    rcs = ''
    while i != -1:
        if s[i] == 'a':
            rcs += 't'
        if s[i] == 'c':
            rcs += 'g'
        if s[i] == 'g':
            rcs += 'c'
        if s[i] == 't':
            rcs += 'a'
        if s[i] == 'n':
            rcs += 'n'
        i -= 1
    return rcs

def getRegion(temp):
    region= {}
    iBegin = temp.rfind('(')
    iEnd = temp.find(')')
    if iBegin == -1:
        iBegin = temp.rfind(' ')
        iEnd = temp.rfind('\r\n')
    region['intervals'] = temp[iBegin + 1:iEnd].split(',')
    for i in range(len(region['intervals'])):
        region['intervals'][i] = region['intervals'][i].split('..')
    for i in range(len(region['intervals'])):
        region['intervals'][i][0] = int(region['intervals'][i][0])
        if len(region['intervals'][i]) != 1:
            region['intervals'][i][1] = int(region['intervals'][i][1])
    if 'complement' in temp:
        region['complement'] = True
    else:
        region['complement'] = False
    return region

def getRCRegions(n, regions):
    for region in regions:
        if region['complement'] == True:
            for i in range(len(region['intervals'])):
                region['intervals'][i][0] = n - region['intervals'][i][0]
                if len(region['intervals'][i]) != 1:
                    region['intervals'][i][1] = n - region['intervals'][i][1]
                    temp = region['intervals'][i][0]
                    region['intervals'][i][0] = region['intervals'][i][1]
                    region['intervals'][i][1] = temp

            region['intervals'].reverse()
    return regions

def getStrandTranslation(region, DNA):
    countNucl = 0
    translation = ''
    j = 0
    begin = region['intervals'][j][0] - 10
    translation += DNA[begin:begin + 11]
    countNucl += 11
    current = begin + 11
    while countNucl < 21:
        if len(region['intervals'][j]) != 1:
            for nucl in DNA[current:min(region['intervals'][j][1] + 1, current + 21 - countNucl)]:
                if nucl != 'n':
                    translation += nucl
                    countNucl += 1
        else:
            nucl = DNA[region['intervals'][j][0]]
            if nucl != 'n':
                translation += nucl
                countNucl += 1

        #translation += DNA[current:min(region['intervals'][j][1] + 1, current + 21 - countNucl)]
        #countNucl = min(countNucl + region['intervals'][j][1] - current + 1, 21)
        j += 1
        if countNucl != 21:
            current = region['intervals'][j][0]
    return translation

def getTranslation(regions, DNA, rcDNA):
    for i in range(len(regions)):
        if regions[i]['complement'] == False:
            regions[i]['translation'] = getStrandTranslation(regions[i], '0' + DNA)
        if regions[i]['complement'] == True:
            regions[i]['translation'] = getStrandTranslation(regions[i], rcDNA)
    return regions

def divideSets(n, regions):
    set = {}
    set['trainingSet'] = []
    set['testSet'] = []
    for region in regions:
        if region['intervals'][0][0] - 10 < n / 2:
            set['trainingSet'].append(region)
        else:
            set['testSet'].append(region)
    return set

def getPFM(trainingSet):
    PFM = {'a':[], 'c':[], 'g':[], 't':[]}
    for symbol in PFM:
        for i in range(21):
            PFM[symbol].append(0)
    for region in trainingSet:
        for i in range(len(region['translation'])):
            for symbol in PFM:
                if region['translation'][i] == symbol:
                    PFM[symbol][i] += 1
    for symbol in PFM:
        for i in range(len(region['translation'])):
            PFM[symbol][i] = round(PFM[symbol][i] * 1.0 / len(trainingSet), 2)
    return PFM

def getPFMPsevdoCounts(trainingSet):
    PFMPC = {'a':[], 'c':[], 'g':[], 't':[]}
    for symbol in PFMPC:
        for i in range(21):
            PFMPC[symbol].append(1)
    for region in trainingSet:
        for i in range(21):
            for symbol in PFMPC:
                if region['translation'][i] == symbol:
                    PFMPC[symbol][i] += 1
    count = []
    for i in range(21):
        count.append(0)
        for symbol in PFMPC:
            count[-1] += PFMPC[symbol][i]
    for symbol in PFMPC:
        for i in range(21):
            PFMPC[symbol][i] = PFMPC[symbol][i] * 1.0 / count[i]

    '''s = []
    for i in range(21):
         s.append(0)
    for i in range(21):
        for symbol in PFMPC:
            s[i] += PFMPC[symbol][i]
    print 's = ', s'''

    #print 'PFMC = ', PFMPC

    return PFMPC

def getISubstring(i, DNA):
    substring = ''
    j = i
    count = 0
    while count != 21:
        if DNA[j] != 'n':
            substring += DNA[j]
            count += 1
        j += 1
    return substring

def getP(DNA):
    symbols = ['a', 'c', 'g', 't']
    p = {}
    for symbol in symbols:
        p[symbol] = 0.0
    for i in range(len(DNA) - 21 + 1):
        substring = getISubstring(i, DNA)
        for symbol in substring:
            p[symbol] += 1
    for symbol in symbols:
        p[symbol] = p[symbol] / (len(DNA) - 21 + 1) / 21
    return p

def getPSSM(trainingSet, n, DNA):
    PSSM = {'a':[], 'c':[], 'g':[], 't':[]}
    PFMPC = getPFMPsevdoCounts(trainingSet)

    '''p = {}
    for symbol in PFMPC:
        p[symbol] = 0
        for i in range(21):
            p[symbol] += PFMPC[symbol][i] / 21'''

    p = getP(DNA)

    '''s = 0
    for symbol in PFMPC:
        s += p[symbol]
    print 's = ', s'''

    print 'p = ', p, '\r\n'

    for symbol in PSSM:
        for i in range(21):
            PSSM[symbol].append(round(math.log(((PFMPC[symbol][i]) / p[symbol]), math.exp(1)), 1))
    return PSSM

def getScoreSet(PSSM, set):
    scores = []
    for region in set:
        score = 0.0
        for i in range(21):
            score += PSSM[region['translation'][i]][i]
        scores.append(score)
    return scores

def getScoreAllSubstrings(PSSM, DNA, scores):
    for i in range(len(DNA) - 21 + 1):
        score = 0.0
        substring = getISubstring(i, DNA)
        for i in range(21):
            score += PSSM[substring[i]][i]
        scores.append(score)
    return scores

def getMean(y):
    mean = 0.0
    for i in range(len(y)):
        mean += y[i]
    mean /= len(y)
    return mean

def getDispersion(y, mean):
    disp = 0.0
    for i in range(len(y)):
        disp += (y[i] - mean) * (y[i] - mean)
    disp /= len(y)
    return disp

def getGistScores(scores, fileName, boolN, xName, yName):
    minS = scores[0]
    maxS = scores[0]
    for score in scores:
        if score > maxS:
            maxS = score
        elif score < minS:
            minS = score
    x = []
    y = []
    for i in range(int(round(minS * 10, 1)), int(round(maxS * 10, 1)) + 1):
        x.append(i * 1.0 / 10)
        y.append(0)
    for score in scores:
        y[int(round((score - minS) * 10, 1))] += 1

    maxY = y[0]
    for i in range(len(y)):
        if y[i] > maxY:
            maxY = y[i]

    figure()

    if boolN == True:
        mean = getMean(scores)
        sigma = sqrt(getDispersion(scores, mean))
        normDistribution = []
        for i in range(len(x)):
            normDistribution.append(len(scores) * 0.1 / (sigma * sqrt(2 * math.pi)) *
                                    math.exp((-(x[i] - mean) * (x[i] - mean)) *
                                             1.0 / (2 * sigma * sigma)))
        plot(x, normDistribution, color='r')
        text(-20, 60000, r'$\mu=' + str(round(mean,2)) + ',\ \sigma=' + str(round(sigma, 2)) + '$')

    bar(x, y, align='center', width=0.1)
    title('Histogram of score')
    xlabel(xName)
    ylabel(yName)
    savefig(fileName)

    ans = {}
    ans['x'] = x
    ans['y'] = y

    return ans

def getPlotFPFN(all, CDS, fileName):
    x = []
    yFP = []
    yFN = []
    for iThreshold in range(len(all['x'])):
        threshold = all['x'][iThreshold]
        sLessCDS = 0
        sMoreCDS = 0
        sCDS = 0
        sMoreAll = 0
        sAll = 0
        for j in range(len(CDS['x'])):
            if CDS['x'][j] < threshold:
                sLessCDS += CDS['y'][j]
            else:
                sMoreCDS += CDS['y'][j]
            sCDS += CDS['y'][j]
        for j in range(len(all['x'])):
            if all['x'][j] > threshold:
                sMoreAll += all['y'][j]
            sAll += all['y'][j]
        yFP.append((sMoreAll - sMoreCDS) * 1.0 / (sAll - sCDS))
        yFN.append(sLessCDS * 1.0 / sCDS)
        x.append(threshold)
    figure()
    title('red - FP, blue - FN')
    xlabel('threshold')
    ylabel('FP, FN')
    plot(x, yFP)
    plot(x, yFN, color='r')
    savefig(fileName)

def findCDSofSingleStrand(PSSM, threshold, strand):
    CDSofSS = []
    scores = []
    scores = getScoreAllSubstrings(PSSM, strand, scores)
    for i in range(len(scores)):
        if scores[i] > threshold:
            substring = getISubstring(i, strand)
            CDSofSS.append(substring)
    return CDSofSS

def findCDSofSecondPart(PSSM, threshold, DNA):
    secondPart = DNA[:len(DNA) / 2]
    rcSecondPart = reverseComplement(secondPart)
    CDSofSP = []
    CDSofSP.extend(findCDSofSingleStrand(PSSM, threshold, secondPart))
    CDSofSP.extend(findCDSofSingleStrand(PSSM, threshold, rcSecondPart))
    return CDSofSP

def getAnalysisofResults(CDSofSP, testSet, nSP):
    analysis = {}
    analysis['TP'] = 0.0
    analysis['TN'] = 0.0
    analysis['FP'] = 0.0
    analysis['FN'] = 0.0

    testTranslation = []
    for region in testSet:
        testTranslation.append(region['translation'])

    for cds in CDSofSP:
        if cds in testTranslation:
            analysis['TP'] += 1
        else:
            analysis['FP'] += 1
    analysis['FN'] = len(testSet) - analysis['TP']
    analysis['TN'] = nSP - len(CDSofSP) - analysis['FN']

    analysis['sensitivity'] = analysis['TP'] / (analysis['TP'] + analysis['FN'])
    analysis['specificity'] = analysis['TN'] / (analysis['TN'] + analysis['FP'])
    analysis['PPV'] = analysis['TP'] / (analysis['TP'] + analysis['FP'])

    return analysis

def parseCDS(fileName):
    CDS = {}
    regions = []
    with open(fileName, 'r') as fin:
        while True:
            temp = fin.readline()
            if 'CDS ' in temp:
                temp1 = ''
                while '/' not in temp1:
                    temp1 = fin.readline()
                    if '/' not in temp1:
                        temp += temp1
                if '>' in temp or '<' in temp:
                    continue
                regions.append(getRegion(temp))
            if 'ORIGIN' in temp:
                DNA = readDNA(fin)
                break
    #print regions
    CDS['regions'] = regions
    CDS['DNA'] = DNA
    return CDS

def main():
    CDS = parseCDS('NC_004317.gbk')
    DNA = CDS['DNA']
    regions = CDS['regions']

    rcDNA = reverseComplement(DNA)
    regions = getRCRegions(len(DNA), regions)

    regions = getTranslation(regions, DNA, rcDNA)

    with open('CDS.txt', 'w') as fout:
        for region in regions:
            if 'translation' in region:
                fout.write(region['translation'] + '\r\n')

    set = divideSets(len(DNA), regions)
    PFM = getPFM(set['trainingSet'])
    print 'PFM = ', PFM, '\r\n'

    PSSM = getPSSM(set['trainingSet'], len(regions), DNA)
    print 'PSSM = ', PSSM
    print '\r\n'

    scores = getScoreSet(PSSM, set['trainingSet'])
    train = getGistScores(scores, 'gistTrainingScore.png', False, 'score for training set', 'count')

    scores = []
    scores = getScoreAllSubstrings(PSSM, DNA, scores)
    scores = getScoreAllSubstrings(PSSM, rcDNA, scores)
    all = getGistScores(scores,'gistAllSubsScore.png', True, 'score of all substrings', 'count')

    getPlotFPFN(all, train, 'graphFPFN.png')

    CDSofSP = findCDSofSecondPart(PSSM, 1, DNA)
    with open('CDSofSP.txt', 'w') as fout:
        for CDS in CDSofSP:
            fout.write(CDS + '\r\n')

    analysis = getAnalysisofResults(CDSofSP, set['testSet'], (len(DNA) / 2 - 21 + 1))
    print analysis

    '''for region in regions:
        if region['translation'][9:12] != 'atg':
            print region'''

main()