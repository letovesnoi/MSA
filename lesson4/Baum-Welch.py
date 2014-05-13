__author__ = 'lenk'

import sys
import math

'''def baumwelch(obs, states, start_p, trans_p, emit_p):
    V = [[]]
    for y in states:
        V[0].append(math.log(start_p[y], 2) + math.log(emit_p[y][obs[0]], 2))
    for t in range(1, len(obs)):
        V.append([])
        for y in states:
            (prob, state) = max(((V[t - 1][y0] + math.log(trans_p[y0][y], 2) + math.log(emit_p[y][obs[t]], 2)), y0)
                                for y0 in range(len(states)))
            V[t].append(prob)
    (prob, state) = max((V[len(obs) - 1][y], y) for y in states)
    return prob'''

def getAlphaProbability(obs, states, start_p, trans_p, emit_p):
    alpha = []
    alpha.append([])
    for i in range(len(states)):
        alpha[0].append(2 ** (math.log(start_p[i], 2) + math.log(emit_p[i][obs[0]], 2)))
    for t in range(1, len(obs)):
        alpha.append([])
        for j in range(len(states)):
            #p = alpha[t - 1][0] * trans_p[0][j]
            #q = alpha[t - 1][1] * trans_p[1][j]
            alpha[t].append(math.log(emit_p[j][obs[t]], 2) + alpha[t - 1][0] + math.log(trans_p[0][j], 2) +
                                  math.log(1 + 2 ** (alpha[t - 1][1] + math.log(trans_p[1][j], 2) - alpha[t - 1][0] -
                                                     math.log(trans_p[0][j], 2)), 2))
    return alpha

def getBetaProbability(obs, states, trans_p, emit_p):
    beta = []
    beta.append([])
    for i in range(len(states)):
        beta[0].append(1)
    for t in range(1, len(obs)):
        beta.append([])
        for i in range(len(states)):
            #p = beta[t - 1][0] * trans_p[i][0] * emit_p[0][obs[len(obs) - 1 - (t - 1)]]
            #q = beta[t - 1][1] * trans_p[i][1] * emit_p[1][obs[len(obs) - 1 - (t - 1)]]
            beta[t].append(beta[t - 1][0] + math.log(trans_p[i][0], 2) + math.log(emit_p[0][obs[len(obs) - 1 - (t - 1)]], 2) +
                                 math.log(1 + 2 ** (beta[t - 1][1] + math.log(trans_p[i][1], 2) +
                                                    math.log(emit_p[1][obs[len(obs) - 1 - (t - 1)]], 2) -
                                                beta[t - 1][0] - math.log(trans_p[i][0], 2) -
                                                    math.log(emit_p[0][obs[len(obs) - 1 - (t - 1)]], 2)), 2))
    beta.reverse()
    return beta

def getGammaProbability(obs, states, alpha, beta):
    gamma = []
    for t in range(len(obs)):
        gamma.append([])
        for i in range(len(states)):
            #p = alpha[t][0] * beta[t][0]
            #q = alpha[t][1] * beta[t][1]
            gamma[t].append(alpha[t][i] + beta[t][i] - alpha[t][0] - beta[t][0] - math.log(1 +
                                       2 ** (alpha[t][1] + beta[t][1] - alpha[t][0] - beta[t][0]), 2))
    return gamma

def getPower(states, obs, alpha, beta, gamma, ksi):
    for t in range(len(obs)):
        for i in range(len(states)):
            alpha[t][i] = 2 ** alpha[t][i]
            beta[t][i] = 2 ** beta[t][i]
            gamma[t][i] = 2 ** gamma[t][i]
            if t != len(obs) - 1:
                for j in range(len(states)):
                    ksi[t][i][j] = 2 ** ksi[t][i][j]
    return alpha, beta, gamma, ksi

def getKsiProbability(obs, states, trans_p, emit_p, alpha, beta):
        ksi = []
        for t in range(len(obs) - 1):
            ksi.append([])
            for k in range(len(states)):
                ksi[t].append([])
                for l in range(len(states)):
                    expression0 = alpha[t][0] + math.log(trans_p[0][1], 2) + beta[t + 1][1] + math.log(emit_p[1][obs[t + 1]], 2) - \
                                  (alpha[t][0] + math.log(trans_p[0][0], 2) + beta[t + 1][0] + math.log(emit_p[0][obs[t + 1]], 2))
                    expression2 = alpha[t][1] + math.log(trans_p[1][1], 2) + beta[t + 1][1] + math.log(emit_p[1][obs[t + 1]], 2) - \
                                  (alpha[t][1] + math.log(trans_p[1][0], 2) + beta[t + 1][0] + math.log(emit_p[0][obs[t + 1]], 2))
                    expression3 = alpha[t][0] + math.log(trans_p[0][1], 2) + beta[t + 1][1] + math.log(emit_p[1][obs[t + 1]], 2) - \
                                  (alpha[t][0] + math.log(trans_p[0][0], 2) + beta[t + 1][0] + math.log(emit_p[0][obs[t + 1]], 2))
                    expression1 = alpha[t][1] + math.log(trans_p[1][0], 2) + beta[t + 1][0] + math.log(emit_p[0][obs[t + 1]], 2) + \
                                  math.log(1 + 2 ** expression2, 2) - (alpha[t][0] + math.log(trans_p[0][0], 2) + beta[t + 1][0] +
                                                                       math.log(emit_p[0][obs[t + 1]], 2) + math.log(1 + 2 ** expression3, 2))
                    denominator = alpha[t][0] + math.log(trans_p[0][0], 2) + beta[t + 1][0] + math.log(emit_p[0][obs[t + 1]], 2) + \
                                  math.log(1 + 2 ** expression0, 2) + math.log(1 + 2 ** expression1, 2)
                    ksi[t][k].append(alpha[t][k] + math.log(trans_p[k][l], 2) + beta[t + 1][l] + math.log(emit_p[l][obs[t + 1]], 2) -
                                     denominator)
        return ksi

def getNewProbability(obs, different_obs, states, start_p, trans_p, emit_p):
    ans = {}
    new_start_p = []
    new_trans_p = []
    new_emit_p = []

    alpha = getAlphaProbability(obs, states, start_p, trans_p, emit_p)
    beta = getBetaProbability(obs, states, trans_p, emit_p)
    gamma = getGammaProbability(obs, states, alpha, beta)
    ksi = getKsiProbability(obs, states, trans_p, emit_p, alpha, beta)

    ans['probability'] = alpha[len(obs) - 1][0] + alpha[len(obs) - 1][1]
    alpha, beta, gamma, ksi = getPower(states, obs, alpha, beta, gamma, ksi)

    for i in range(len(states)):
        new_start_p.append(gamma[0][i])
    for i in range(len(states)):
        new_trans_p.append([])
        new_emit_p.append({})
        for j in range(len(states)):
            sum_tn = 0.0
            sum_td = 0.0
            for t in range(len(obs) - 1):
                    sum_tn += ksi[t][i][j]
                    sum_td += gamma[t][i]
            new_trans_p[i].append(2 ** (math.log(sum_tn, 2) -
                                        math.log(sum_td, 2)))
        for j in range(len(different_obs)):
            sum_en = 0.0
            sum_ed = 0.0
            for t in range(len(obs)):
                if obs[t] == different_obs[j]:
                    sum_en += gamma[t][i]
                sum_ed += gamma[t][i]
            new_emit_p[i][different_obs[j]] = 2 ** (math.log(sum_en, 2) - math.log(sum_ed, 2))
    ans['start_probability'] = new_start_p
    ans['emission_probability'] = new_emit_p
    ans['transition_probability'] = new_trans_p
    print 'prob', ans['probability']
    return ans

def main():
    states = [0, 1]
    different_obs = ['A', 'C', 'G', 'T']

    start_probability = [0.996, 0.004]

    emission_probability = [{'A': 0.291, 'C': 0.209, 'G': 0.209, 'T': 0.291},
                            {'A': 0.169, 'C': 0.331, 'G': 0.331, 'T': 0.169}]

    transition_probability = [[0.999, 0.001], [0.01, 0.99]]

    observations = ''
    with open(sys.argv[1], 'r') as fin:
        seqName = fin.readline().strip()
        for line in fin:
            observations += line.strip()

    # !!!!
    count_iterations = 0
    new_probability = 0.0
    prev_probability = 0.0
    while new_probability - prev_probability > 0.1 or count_iterations < 3:
        print 'i =', count_iterations
        count_iterations += 1
        newProbabilities = getNewProbability(observations, different_obs, states, start_probability,
                                             transition_probability, emission_probability)
        prev_probability = new_probability
        new_probability = newProbabilities['probability']
        start_probability = newProbabilities['start_probability']
        transition_probability = newProbabilities['transition_probability']
        emission_probability = newProbabilities['emission_probability']
        print 'start_p=' + str(start_probability), 'trans_p=' + str(transition_probability), \
            'emit_p=' + str(emission_probability), 'p = ' + str(prev_probability) + ', ' + str(new_probability)
    print 'delta_p = ' + str(new_probability - prev_probability)

    with open('output.txt', 'w') as fout:
        fout.write(sys.argv[1] + '\r\n')
        fout.write(seqName + '\r\n')
        fout.write('count_iterations = ' + str(count_iterations) + '\r\n')
        fout.write('final_log-likelihood = ' + str(new_probability) + '\r\n')
        fout.write('transition_probability:\r\n')
        for state1 in range(len(transition_probability)):
            for state2 in range(len(transition_probability[state1])):
                fout.write(str(state1) + '->' + str(state2) + ': %.4f '
                           % transition_probability[state1][state2])
                fout.write('\r\n')
        fout.write('emission_probability:\r\n')
        for state in states:
            for d in different_obs:
                fout.write('e_' + str(d) + str(state) + ' = ' + str(emission_probability[state][d]) + '\r\n')

main()
