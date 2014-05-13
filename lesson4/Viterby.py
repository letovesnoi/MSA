__author__ = 'lenk'

import sys
import copy
import math

def viterbi(obs, states, start_p, trans_p, emit_p):
    V = [[]]
    countStates = {0: [0, 0], 1: [0, 0]}
    iSegmentsB = {0: [[], []], 1: [[], []]}
    iSegmentsE = {0: [[], []], 1: [[], []]}
    for y in states:
        V[0].append(start_p[y] + emit_p[y][obs[0]])
        countStates[y][y] += 1
        iSegmentsB[y][y].append(0)
    for t in range(1, len(obs)):
        newCountStates = {}
        newISegmentsB = {}
        newISegmentsE = {}
        V.append([])
        for y in states:
            (prob, state) = max((V[t - 1][y0] + trans_p[y0][y] + emit_p[y][obs[t]], y0)
                            for y0 in range(len(states)))
            V[t].append(prob)
            if y != state:
                newISegmentsB[y] = copy.deepcopy(iSegmentsB[state])
                newISegmentsB[y][y].append(t)
                newISegmentsE[y] = copy.deepcopy(iSegmentsE[state])
                newISegmentsE[y][state].append(t - 1)
            else:
                newISegmentsB[y] = iSegmentsB[state]
                newISegmentsE[y] = iSegmentsE[state]
            newCountStates[y] = copy.deepcopy(countStates[state])
            newCountStates[y][y] += 1
        countStates = newCountStates
        if newISegmentsB != {}:
            iSegmentsB = newISegmentsB
        if newISegmentsE != {}:
            iSegmentsE = newISegmentsE
    (prob, state) = max((V[len(obs) - 1][y], y) for y in states)
    iSegmentsE[state][state].append(len(obs) - 1)

    ans = {}
    ans['countStates'] = countStates[state]
    ans['iGC1SegmentsB'] = iSegmentsB[state]
    ans['iGC1SegmentsE'] = iSegmentsE[state]
    ans['countSegments'] = []
    ans['countSegments'].append(len(iSegmentsB[state][0]))
    ans['countSegments'].append(len(iSegmentsB[state][1]))
    #print ans
    return ans

def getTransation_probability(countStates, countSegments):
    transition_probability = [[], []]

    transition_probability[0].append(math.log((countStates[0] + 1) * 1.0 /
                                              (countStates[0] + countSegments[1] + 2)))
    transition_probability[0].append(math.log((countSegments[1] + 1) * 1.0 /
                                              (countStates[0] + countSegments[1] + 2)))
    transition_probability[1].append(math.log((countSegments[0] + 1) * 1.0 /
                                              (countStates[1] + countSegments[0] + 2)))
    transition_probability[1].append(math.log((countStates[1] + 1) * 1.0 /
                                              (countStates[1] + countSegments[0] + 2)))
    return transition_probability

def main():
    states = [0, 1]

    start_probability = [math.log(0.996), math.log(0.004)]

    emission_probability = [{'A': math.log(0.291), 'C': math.log(0.209),
                             'G': math.log(0.209), 'T': math.log(0.291)},
                            {'A': math.log(0.169), 'C': math.log(0.331),
                             'G': math.log(0.331), 'T': math.log(0.169)}]

    observations = ''
    with open(sys.argv[1], 'r') as fin:
        seqName = fin.readline().strip()
        for line in fin:
            observations += line.strip()
    #print len(observations)
    intervals = []
    nameIntervals = []
    with open(sys.argv[2], 'r') as fin:
        line = ''
        while line != '//':
            line = fin.readline().strip()
            if 'gene' in line and 'complement' not in line and '..' in line:
                line = line.replace('gene', '').strip()
                intervals.append(line.split('..'))
                nameIntervals.append(fin.readline().strip())
            elif 'gene' in line and 'complement' in line and '..' in line:
                line = line[line.find('(') + 1:line.find(')')]
                #temp = line.split('..')
                #intervals.append([len(observations) - int(temp[1]), len(observations) - int(temp[0])])
                intervals.append(line.split('..'))
                nameIntervals.append(fin.readline().strip())
    transition_probability = [[[math.log(0.999), math.log(0.001)], [math.log(0.01), math.log(0.99)]]]
    countStates = []
    countSegments = []
    for i in range(1, 11):
        print 'i =', i
        iterationAns = viterbi(observations, states, start_probability, transition_probability[-1],
                                  emission_probability)
        countStates.append(iterationAns['countStates'])
        countSegments.append(iterationAns['countSegments'])
        transition_probability.append(getTransation_probability(countStates[-1], countSegments[-1]))

    with open('output.txt', 'w') as fout:
        fout.write(sys.argv[1] + '\r\n')
        fout.write(seqName + '\r\n')
        for i in range(1, 11):
            #print 'i =', i
            fout.write('i: ' + str(i) + '\r\n')
            fout.write('countStates: ' + str(countStates[i - 1]) + '\r\n')
            fout.write('countSegments: ' + str(countSegments[i - 1]) + '\r\n')
            fout.write('transition_probability: ')
            for state1 in range(len(transition_probability[i])):
                for state2 in range(len(transition_probability[i][state1])):
                    fout.write(str(state1) + '->' + str(state2) + ': %.4f '
                               % math.exp(transition_probability[i][state1][state2]))
            fout.write('\r\n')

        iSegmentsB = iterationAns['iGC1SegmentsB']
        iSegmentsE = iterationAns['iGC1SegmentsE']
        print iSegmentsB[1], iSegmentsE[1]
        fout.write('GC rich segments:\r\n')
        for i in range(len(iSegmentsB[1])):
            fout.write(observations[iSegmentsB[1][i]:iSegmentsE[1][i]])
            fout.write('\r\n')

        i = 0
        #iGene = []
        while i in range(len(iSegmentsB[1])) and i < 5:
            for j in range(len(intervals)):
                interval = intervals[j]
                if (int(interval[0]) >= iSegmentsB[1][i] and int(interval[0]) <= iSegmentsE[1][i]) \
                    or (int(interval[1]) >= iSegmentsB[1][i] and int(interval[1]) <= iSegmentsE[1][i]) \
                    or (int(interval[0]) <= iSegmentsB[1][i] and int(interval[1]) >= iSegmentsE[1][i]):
                    fout.write(nameIntervals[j] + '\r\n')
                    print 'yes'
            i += 1
        #fout.write('coordinates of gene: ' + str(iGene) + '\r\n')

main()


